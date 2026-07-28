# find the optimal solution determined by the priority ordering `by`:
# a package => version-index dict, or nothing when the requirements are
# not jointly satisfiable
#
# The descent repeatedly pins the highest-priority *front-forced* package --
# one contained in every Pareto-optimal tight solution consistent with the
# pins so far -- at the best version among those solutions, until no unpinned
# package is front-forced. "Tight" means the solution is exactly the
# dependency closure of the requirements, and dominance is coverage-monotone
# version dominance; see docs/design/front-necessity.md for the semantics
# and the proofs backing this implementation.
#
# Implementation shape: candidate solutions are *constructed*, not searched
# for -- the dependency-layered greedy runs under per-query assumptions and
# lands on provably Pareto-maximal solutions whenever its answer is junk-free
# (a dominator would have to agree with every greedy pin before its earliest
# strict improvement, contradicting that pin's optimality). Pins are carried
# as solver assumptions, never clauses, because dominance queries range over
# all tight solutions, not just pin-consistent ones.
function resolve_core(
    sat  :: SAT{P},
    reqs :: SetOrVec{P} = keys(sat.info);
    by   :: Function = identity, # priority ordering
    restore :: Bool = true, # restore the SAT instance's state before returning
) where {P}
    # a solution exists iff the requirements are jointly satisfiable; the
    # check only assumes, so an unsatisfiable resolve leaves no state behind
    is_satisfiable(sat, reqs) || return nothing
    sol = Dict{P,Int}() # scratch model
    pins = Dict{P,Int}() # the solution
    front = Dict{P,Int}[] # certified front members consistent with the pins
    # packages present in every pin-consistent solution -- requirements and
    # dependencies of pinned versions; sound necessity fast paths since
    # forced-everywhere implies forced-on-the-front
    depforced = Set{P}(reqs)

    assume_pins() = sat_assume(sat, pins)

    # the dependency closure of the requirements under the model's versions
    function closure(m::Dict{P,Int})
        seen = Set{P}(reqs)
        queue = P[q for q in reqs]
        while !isempty(queue)
            q = pop!(queue)
            info_q = sat.info[q]
            i = m[q]
            for (j, r) in enumerate(info_q.depends)
                info_q.conflicts[i, j] || continue
                r in seen && continue
                push!(seen, r)
                push!(queue, r)
            end
        end
        return seen
    end
    # prune a model to its tight core (models may carry unjustified packages)
    prune(m::Dict{P,Int}) =
        let c = closure(m)
            Dict{P,Int}(q => i for (q, i) in m if q in c)
        end

    # reverse dependency edges: package => [(dependent, version index)];
    # used to require forced-present packages to be depended upon, which
    # prunes unjustified-junk models from query spaces (tight solutions
    # always satisfy these clauses, so nothing sound is lost)
    rdeps = Dict{P,Vector{Tuple{P,Int}}}()
    for (q, info_q) in sat.info
        for (j, r) in enumerate(info_q.depends)
            edges = get!(() -> Tuple{P,Int}[], rdeps, r)
            for i = 1:length(info_q.versions)
                info_q.conflicts[i, j] && push!(edges, (q, i))
            end
        end
    end
    reqset = Set{P}(q for q in reqs)

    # a present package must be depended upon by some chosen version --
    # requirements are roots and exempt. Restricting all queries to these
    # *supported* models loses no tight solution (every non-root member of a
    # tight solution is depended upon) and reduces junk to dependency cycles,
    # which the pruning and blocking machinery still handles.
    function support_clause(x::P)
        x in reqset && return
        sat_add(sat, x, not=true)
        for (q, i) in get(() -> Tuple{P,Int}[], rdeps, x)
            sat_add(sat, q, i)
        end
        sat_add(sat)
    end

    # exclude every model whose tight core is exactly `s`: a model agreeing
    # with s on supp(s) has core s -- s is dependency-closed and contains the
    # requirements, so the closure can never escape it -- hence one small
    # clause disposes of the entire junk-variant family of a processed core
    function block_core(s::Dict{P,Int})
        for (q, i) in s
            sat_add(sat, q, i, not=true)
        end
        sat_add(sat)
    end

    # exclude the dominance cone of t: every solution it covers at least as
    # well is either t itself or strictly dominated by it
    function block_cone(t::Dict{P,Int})
        for q in keys(sat.info)
            if !haskey(t, q)
                sat_add(sat, q)
            else
                for k = 1:t[q]-1
                    sat_add(sat, q, k)
                end
            end
        end
        sat_add(sat)
    end

    # is s dominated-or-equal to a cached front member? (cheap, outside SAT)
    function cache_dominated(s::Dict{P,Int})
        for f in front
            ok = true
            for (q, i) in s
                if get(f, q, typemax(Int)) > i
                    ok = false
                    break
                end
            end
            ok && return true
        end
        return false
    end

    # the dependency-layered greedy over the models that satisfy the ambient
    # clauses and the per-solve assumptions: optimize each root package to
    # its best feasible rank (given the greedy's own pins, kept as units in
    # a nested scope), then expand along the chosen versions' dependencies,
    # layer by layer, in priority order. Returns the greedy assignment (its
    # layer closure), or nothing when the space is empty -- which is
    # conclusive, since candidates are models. A junk-free result (equal to
    # its own requirement closure) is Pareto-maximal *within the space*:
    # an in-space dominator would have to agree with every greedy pin before
    # its earliest strict improvement, contradicting that pin's optimality.
    function greedy(roots::Set{P}; assume::Function = () -> nothing)
        g = Dict{P,Int}() # greedy pins
        assume()
        is_satisfiable(sat) || return nothing
        extract_solution!(sat, sol)
        with_temp_clauses(sat) do
            opts = roots
            seen = copy(opts)
            while true
                for q in sort!(collect(opts); by)
                    optimize_version!(sat, sol, q; assume)
                    g[q] = sol[q]
                    sat_add(sat, q, g[q])
                    sat_add(sat)
                end
                opts′ = empty(opts)
                for q in opts
                    info_q = sat.info[q]
                    i = g[q]
                    for (j, r) in enumerate(info_q.depends)
                        info_q.conflicts[i, j] || continue
                        r ∉ seen && push!(opts′, r)
                    end
                end
                isempty(opts′) && break
                union!(seen, opts′)
                opts = opts′
            end
        end
        return g
    end

    # does the solution agree exactly with the pins?
    pin_exact(f::Dict{P,Int}) = all(get(f, q, 0) == i for (q, i) in pins)

    # find a strict tight dominator of s, or nothing (conclusively). The
    # search space -- models covering s at versions at least as good, strictly
    # better somewhere -- is set up as clauses, and the greedy constructs a
    # maximal candidate in it; a junk-free result is a genuine dominator (its
    # roots include supp(s), so junk-freeness implies its core covers s), and
    # by the greedy-maximality argument it is itself a front member. Junky
    # results are excluded by core and retried.
    function dominator(s::Dict{P,Int})
        any(i > 1 for i in values(s)) || return nothing
        result = nothing
        with_temp_clauses(sat) do
            for (q, i) in s
                sat_add(sat, q)
                sat_add(sat)
                for k = 1:i
                    sat_add(sat, q, k)
                end
                sat_add(sat)
            end
            for (q, i) in s, k = 1:i-1
                sat_add(sat, q, k)
            end
            sat_add(sat)
            roots = union(Set{P}(q for q in reqs), keys(s))
            while true
                g = greedy(roots)
                g === nothing && break # empty space: no dominator at all
                t = prune(g)
                if length(t) == length(g)
                    result = t
                    break
                end
                # an improvement orphaned part of the space; exclude the
                # offending core (and its junk variants) and retry
                block_core(t)
            end
        end
        return result
    end

    # is package p contained in every front member consistent with the pins?
    # A witness (a pin-exact front member without p) is *constructed*: the
    # greedy runs under the assumptions pins ∧ ¬p; if even its first solve is
    # unsatisfiable there is no candidate at all and p is forced. A junk-free
    # candidate that no tight solution dominates is a front member -- the
    # witness. When a dominator exists it provably contains p and agrees with
    # the pins (a dominator lacking p would sit in the search space and beat
    # the greedy; one improving a pin would hand its own maximal dominator a
    # pin improvement, contradicting Theorem A) -- so it is cached and its
    # entire dominance cone is excluded, and the search continues.
    function forced(p::P)
        p in depforced && return true
        # forced by top-level propagation from the requirement units
        sat_forced_toplevel(sat, p) && return true
        # a certified front member without p settles it for free
        any(f -> !haskey(f, p), front) && return false
        without_p = () -> begin
            assume_pins()
            sat_assume(sat, p, not=true)
        end
        roots = union(Set{P}(q for q in reqs), keys(pins))
        blocks = Dict{P,Int}[] # cores to exclude (junk variants included)
        cones = Dict{P,Int}[] # cone tops (front members containing p)
        while true
            g = with_temp_clauses(sat) do
                for c in blocks
                    block_core(c)
                end
                for t in cones
                    block_cone(t)
                end
                greedy(roots, assume = without_p)
            end
            g === nothing && return true # no candidate at all: forced
            s = prune(g)
            if length(s) != length(g) || !pin_exact(s) || cache_dominated(s)
                # junky, or the core dropped a pin, or a cached front member
                # already dominates it: not a usable candidate
                push!(blocks, s)
                continue
            end
            t = dominator(s)
            if t === nothing
                # undominated: a pin-exact front member without p
                push!(front, s)
                return false
            end
            @assert haskey(t, p)
            @assert pin_exact(t)
            push!(front, t)
            push!(cones, t)
        end
    end

    # pin front-forced p at the best version any pin-consistent front member
    # gives it: bisect to the best pin-consistent model rank as a lower
    # bound, then walk ranks upward; at each rank the greedy runs under
    # pins ∧ p@rank, and a junk-free pin-exact result is certified with no
    # further solves -- a dominator would have p at an exhausted better rank
    # or improve a pin, both impossible, so it would sit in the search space
    # and beat the greedy. Rank advance happens only when no model exists at
    # the rank, which is conclusive.
    function pin!(p::P)
        assume_pins()
        sat_assume(sat, p)
        is_satisfiable(sat) ||
            error("internal error: no model for forced package")
        extract_solution!(sat, sol)
        optimize_version!(sat, sol, p, assume = () -> begin
            assume_pins()
            sat_assume(sat, p)
        end)
        r = sol[p]
        n = length(sat.info[p].versions)
        # a cached front member already at rank r certifies the pin for free:
        # it is junk-free, agrees with the pins, and adding p@r keeps it exact
        cached = findfirst(f -> f[p] == r, front)
        if cached !== nothing
            pins[p] = r
            @goto pinned
        end
        roots = union(Set{P}(q for q in reqs), keys(pins), (p,))
        blocks = Dict{P,Int}[]
        while true
            at_rank = () -> begin
                assume_pins()
                sat_assume(sat, p, r)
            end
            g = with_temp_clauses(sat) do
                for c in blocks
                    block_core(c)
                end
                greedy(roots, assume = at_rank)
            end
            if g === nothing
                # no model at this rank at all: advance
                r += 1
                r ≤ n ||
                    error("internal error: no front version for forced package")
                empty!(blocks)
                continue
            end
            s = prune(g)
            if length(s) == length(g) && pin_exact(s) && get(s, p, 0) == r
                pins[p] = r
                push!(front, s)
                break
            end
            push!(blocks, s)
        end
        @label pinned
        # cached members that disagree with the pin are no longer
        # pin-consistent; the pinned version's dependencies become forced
        filter!(f -> f[p] == pins[p], front)
        info_p = sat.info[p]
        for (j, q) in enumerate(info_p.depends)
            info_p.conflicts[pins[p], j] && push!(depforced, q)
        end
    end

    function descend!()
        # force all requirements: every tight solution covers them, so the
        # dominance queries want them too
        for q in reqs
            sat_add(sat, q)
            sat_add(sat)
        end
        # restrict every query to supported models
        for x in keys(sat.info)
            support_clause(x)
        end
        # seed the certificate cache with a front member: the greedy over the
        # unconstrained space, retried past junky results
        seed_roots = Set{P}(q for q in reqs)
        seed_blocks = Dict{P,Int}[]
        while true
            g = with_temp_clauses(sat) do
                for c in seed_blocks
                    block_core(c)
                end
                greedy(seed_roots)
            end
            g === nothing && error("internal error: no seed solution")
            s = prune(g)
            if length(s) == length(g)
                push!(front, s)
                break
            end
            push!(seed_blocks, s)
        end
        while true
            # candidate packages: those in every certified front member --
            # anything else is certifiably avoidable
            cand = reduce(intersect, (Set{P}(keys(f)) for f in front))
            next = nothing
            for p in sort!(collect(cand); by)
                haskey(pins, p) && continue
                if forced(p)
                    next = p
                    break
                end
            end
            next === nothing && break
            pin!(next)
        end
    end

    # `restore = true` rolls the instance's clauses back (reusable-instance
    # contract); `restore = false` runs at the top level, leaving the
    # instance constrained by the requirements -- measurably faster (picosat
    # keeps its learned clauses), for single-use instances only.
    if restore
        with_temp_clauses(descend!, sat)
    else
        descend!()
    end

    return pins
end

function resolve(
    sat  :: SAT{P,V},
    reqs :: SetOrVec{P} = keys(sat.info);
    by   :: Function = identity, # priority ordering
    restore :: Bool = true, # restore the SAT instance's state before returning
) where {P,V}
    sol = resolve_core(sat, reqs; by, restore)
    sol === nothing && return nothing
    return Dict{P,V}(p => sat.info[p].versions[i] for (p, i) in sol)
end

# convenience entry points

# NOTE: the data-level entry points build the *unreduced* problem. The
# front-necessity semantics is relative to the problem being solved, and a
# reduction is only admissible if every Pareto-front member survives it as a
# model; pkg_info's reachability filter does not guarantee that (it can drop
# front members outright), so it must not run here. A provably
# front-preserving reduction is future work (docs/design/front-necessity.md).
function resolve(
    deps :: DepsProvider{P},
    reqs :: SetOrVec{P} = deps.packages;
    by   :: Function = identity, # package ordering
) where {P}
    info = pkg_info(deps, reqs; filter=false)
    resolve(info, reqs; by)
end

function resolve(
    data :: AbstractDict{P,<:PkgData{P}},
    reqs :: SetOrVec{P} = keys(data);
    by   :: Function = identity, # package ordering
) where {P}
    info = pkg_info(data, reqs; filter=false)
    resolve(info, reqs; by)
end

function resolve(
    info :: Dict{P,PkgInfo{P,V}},
    reqs :: SetOrVec{P} = keys(info);
    by   :: Function = identity, # package ordering
) where {P,V}
    sat = SAT(info)
    # the instance is single-use, so don't bother restoring its state
    try resolve(sat, reqs; by, restore=false)
    finally
        finalize(sat)
    end
end
