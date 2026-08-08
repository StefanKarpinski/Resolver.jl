mutable struct SAT{P,V}
    info :: Dict{P,PkgInfo{P,V}}
    # per package, the version index each class stands for (0: deactivated) —
    # what turns a solution over classes back into one over versions
    reps :: Dict{P,Vector{Int}}
    pico :: Ptr{Cvoid}
    vars :: Dict{P,Int}
    # the literals forbidding the deactivated classes, sorted; internal.
    # they are asserted as unit clauses inside a push frame, so one sat_pop
    # reactivates every class at once — after which any subset of them can be
    # deactivated again by assumption, no clause rewriting involved
    deact :: Vector{Int}
    # how many push frames are open. `sat_pop` retracts whatever was pushed
    # last, so every helper that opens a frame is relying on the frames below
    # being the ones it thinks they are; this is what lets it say so. picosat
    # reports the current context as a literal rather than a level, and recycles
    # those literals, so counting is the way to know
    depth :: Int
end

function Base.show(io::IO, sat::SAT)
    show(io, typeof(sat))
    p = length(sat.vars)
    v = PicoSAT.var_count(sat.pico)
    c = PicoSAT.clause_count(sat.pico)
    print(io,
        "(packages: ", p,
        ", classes: ", v-p,
        ", clauses: ", c, ")")
end

# variable indices:
#   p     vars[p]     package p chosen
#   p@c   vars[p]+c   class c of p chosen
#   p@≤k  lads[p]+k   some class ≤ k of p chosen (prefix ladder,
#                     only for packages with several classes)
#
# A variable is a *class*, not a version: class members are indistinguishable
# to everything in the registry, so no solution can depend on which member
# stands for a class, and the instance is the smaller for it. `sat.reps` names
# the member at the end.
#
# returns the (sorted) names, the two index maps, and the last used variable
function sat_variables(
    info :: Dict{P, <: PkgInfo{P}},
) where {P}
    # sort names for predictability
    names = sort!(collect(keys(info)))
    N = 1
    vars = Dict{P,Int}()
    lads = Dict{P,Int}()
    for p in names
        vars[p] = N
        n_p = nclasses(info[p])
        # 1 variable for package
        # n_p variables for classes
        N += 1 + n_p
        if n_p ≥ 2
            # n_p prefix-ladder variables at lads[p]+1 : lads[p]+n_p
            lads[p] = N - 1
            N += n_p
        end
    end
    return names, vars, lads, N - 1
end

# append the literals for "no class in lo:hi of the package
# with variable base v, ladder base l, and n classes is chosen"
function run_lits!(lits::Vector{Int}, v::Int, l::Int, lo::Int, hi::Int, n::Int)
    if lo == hi
        push!(lits, -(v + lo))
    elseif lo == 1 && hi == n
        push!(lits, -v)
    elseif lo == 1
        push!(lits, -(l + hi))
    else
        push!(lits, -(l + hi), l + lo - 1)
    end
    return lits
end

function SAT(
    univ :: Universe{P,V},
) where {P,V}
    info = univ.info
    names, vars, lads, N = sat_variables(info)

    # instantiate picosat solver
    pico = PicoSAT.init() # TODO: use jl_malloc?
    try # free memory on error
        PicoSAT.adjust(pico, N)

        lits = Int[]
        # default unconstrained variables to false: models then carry
        # fewer spuriously-true packages ("junk"), which makes solves
        # faster and improvement steps land on better versions
        PicoSAT.set_global_default_phase(pico, 0)
        # ... except each package's best class, which defaults to
        # true: when a package must be chosen the solver tries its best
        # class first, so models land at or near the optimum and the
        # descent's improvement loops mostly never need to run
        for v_p in values(vars)
            PicoSAT.set_default_phase_lit(pico, v_p + 1, 1)
        end

        # generate SAT problem
        for p in names
            info_p = info[p]
            n_p = nclasses(info_p)
            v_p = vars[p]

            # package implies some class
            #   p => OR_i p@i
            PicoSAT.add(pico, -v_p)
            for i = 1:n_p
                PicoSAT.add(pico, v_p + i)
            end
            PicoSAT.add(pico, 0)

            # class implies its package
            #   p@i => p
            for i = 1:n_p
                PicoSAT.add(pico, -(v_p + i))
                PicoSAT.add(pico, v_p)
                PicoSAT.add(pico, 0)
            end

            # prefix ladder: L_k holds iff some class ≤ k is chosen.
            # the two upward directions make chosen classes force their
            # ladder suffix true; the completion direction, together with
            # at-most-one below, forces the ladder prefix *below* the
            # chosen class false — which the interval conflict literals
            # rely on when they occur positively
            if n_p ≥ 2
                l_p = lads[p]
                for k = 1:n_p
                    # p@k => L_k
                    PicoSAT.add(pico, -(v_p + k))
                    PicoSAT.add(pico, l_p + k)
                    PicoSAT.add(pico, 0)
                end
                for k = 2:n_p
                    # L_(k-1) => L_k
                    PicoSAT.add(pico, -(l_p + k - 1))
                    PicoSAT.add(pico, l_p + k)
                    PicoSAT.add(pico, 0)
                end
                # completion: L_k => p@k OR L_(k-1)
                PicoSAT.add(pico, -(l_p + 1))
                PicoSAT.add(pico, v_p + 1)
                PicoSAT.add(pico, 0)
                for k = 2:n_p
                    PicoSAT.add(pico, -(l_p + k))
                    PicoSAT.add(pico, v_p + k)
                    PicoSAT.add(pico, l_p + k - 1)
                    PicoSAT.add(pico, 0)
                end
                # classes are mutually exclusive:
                #   p@k => no class < k (linear via the ladder)
                for k = 2:n_p
                    PicoSAT.add(pico, -(v_p + k))
                    PicoSAT.add(pico, -(l_p + k - 1))
                    PicoSAT.add(pico, 0)
                end
            end

            # dependencies, one clause per maximal run of classes
            # sharing the dependency (dep sets change rarely across
            # versions, so runs are long):
            #   (some class in run chosen) => q
            l_p = get(lads, p, 0)
            for (k, q) in enumerate(info_p.depends)
                v_q = vars[q]
                i = 1
                while i ≤ n_p
                    if info_p.conflicts[i, k]
                        lo = i
                        while i < n_p && info_p.conflicts[i + 1, k]
                            i += 1
                        end
                        empty!(lits)
                        run_lits!(lits, v_p, l_p, lo, i, n_p)
                        for x in lits
                            PicoSAT.add(pico, x)
                        end
                        PicoSAT.add(pico, v_q)
                        PicoSAT.add(pico, 0)
                    end
                    i += 1
                end
            end

            # conflicts are added below, rectangle-encoded
        end

        # conflicts, encoded as rectangles: within each interacting pair's
        # conflict bitmap, group p's classes by their conflict pattern
        # against q (identical patterns are the norm, since conflicts come
        # from shared compat entries) and split each pattern into maximal
        # runs of q's classes. each (class group) × (run) rectangle
        # becomes one clause forbidding "some class in the group chosen
        # AND some class in the run chosen". a side contributes its
        # class literal when it is a singleton, its package variable
        # when every class is covered, prefix-ladder literals when it
        # is an interval — ¬L_hi ∨ L_(lo-1), correct because a chosen
        # class in the interval forces L_hi true and (by at-most-one
        # plus ladder completion) L_(lo-1) false — and otherwise a shared
        # auxiliary trigger defined by one implication per class.
        # triggers occur only negatively, so their one-directional
        # definitions suffice for exact model projection.
        psets = Dict{Tuple{Int,Vector{Int}},Int}() # (v_p, classes) => var
        pat = UInt64[]
        runs = Tuple{Int,Int}[]

        # literal for "some class in S of the package with variable
        # base v and n classes is chosen" (non-contiguous S only)
        function set_trigger(v::Int, S::Vector{Int})
            get!(() -> begin
                t = PicoSAT.inc_max_var(pico)
                for i in S
                    PicoSAT.add(pico, -(v + i))
                    PicoSAT.add(pico, t)
                    PicoSAT.add(pico, 0)
                end
                t
            end, psets, (v, S))
        end

        for p in names
            info_p = info[p]
            n_p = nclasses(info_p)
            v_p = vars[p]
            l_p = get(lads, p, 0)
            for (q, b) in info_p.interacts
                p < q || continue # conflicts are symmetrical
                info_q = info[q]
                n_q = nclasses(info_q)
                v_q = vars[q]
                l_q = get(lads, q, 0)
                Y = info_q.conflicts
                c = info_q.interacts[p]
                W = col_words(Y)
                resize!(pat, W)
                # group p's classes by conflict pattern: the pattern of
                # p@i is the contiguous column c+i of q's matrix
                groups = Dict{Vector{UInt64},Vector{Int}}()
                for i = 1:n_p
                    col_copy!(pat, Y, c + i)
                    clear_rows_above!(pat, n_q)
                    all(iszero, pat) && continue
                    push!(get!(() -> Int[], groups, copy(pat)), i)
                end
                for (pattern, Sg) in groups
                    # p-side literals: interval via the ladder, otherwise
                    # a shared trigger
                    empty!(lits)
                    if Sg[end] - Sg[1] + 1 == length(Sg)
                        run_lits!(lits, v_p, l_p, Sg[1], Sg[end], n_p)
                    else
                        push!(lits, -set_trigger(v_p, Sg))
                    end
                    np_lits = length(lits)
                    # maximal runs of the pattern
                    empty!(runs)
                    prev = -2
                    @inbounds for w = 1:W
                        z = pattern[w]
                        while !iszero(z)
                            j = ((w - 1) << 6) + trailing_zeros(z) + 1
                            z &= z - 1
                            if j == prev + 1
                                runs[end] = (runs[end][1], j)
                            else
                                push!(runs, (j, j))
                            end
                            prev = j
                        end
                    end
                    for (lo, hi) in runs
                        resize!(lits, np_lits)
                        run_lits!(lits, v_q, l_q, lo, hi, n_q)
                        for x in lits
                            PicoSAT.add(pico, x)
                        end
                        PicoSAT.add(pico, 0)
                    end
                end
            end
        end
    catch # on error free picosat solver
        PicoSAT.reset(pico)
        rethrow()
    end
    sat = finalizer(finalize, SAT(info, univ.reps, pico, vars, Int[], 0))
    try deactivate_classes!(sat)
    catch
        finalize(sat)
        rethrow()
    end
    return sat
end

# the structural instance for a bare artifact: every class stands for its best
# member and nothing is deactivated, which is what an unconstrained query makes
# of it (see `Universe`)
SAT(info :: Dict{P,PkgInfo{P,V}}) where {P,V} = SAT(Universe(info))

# Forbid the deactivated classes — the ones this query admits no member of.
#
# This is the *whole* of how a user constraint reaches the solver. It never
# forbids a version, because a version is not something the instance can talk
# about: the members of a class are indistinguishable, so a constraint that
# spares one of them has left the class choosable and there is nothing to say.
# What it can do is empty a class, and an empty class is one unit clause.
#
# The clauses go in a push frame of their own, so production solves see them at
# level 0 (no assumptions) while a single `sat_pop` reactivates every class at
# once — after which any subset can be deactivated again by assuming its
# literal, with no clause rewritten. Nothing pops the frame in production; the
# frame exists so that what is built on top of this can, which is
# `with_classes_relaxed`.
function deactivate_classes!(sat::SAT{P,V}) where {P,V}
    pico = sat.pico
    for (p, reps_p) in sat.reps
        v_p = sat.vars[p]
        best = 0
        for (i, r) in enumerate(reps_p)
            if iszero(r)
                push!(sat.deact, forbidden_lit(sat, p, i))
            elseif best == 0
                best = i
            end
        end
        # the structural instance defaults every package's best class to true,
        # so that a package that must be chosen is tried at its best class
        # first. when this query has emptied that class, point the phase at the
        # best class it does admit instead — otherwise the solver's first guess
        # is infeasible for every such package at once
        if best != 1
            PicoSAT.set_default_phase_lit(pico, v_p + 1, 0)
            best == 0 || PicoSAT.set_default_phase_lit(pico, v_p + best, 1)
        end
    end
    isempty(sat.deact) && return sat
    sort!(sat.deact; by = abs) # dict order must not reach the solver
    return push_deactivations!(sat)
end

# Assert `sat.deact` as unit clauses in a push frame of its own: one `sat_pop`
# retracts all of them at once, and asserting them again is what puts the frame
# back (`with_classes_relaxed`). Which classes are forbidden and where each
# package's phase hint points are properties of the query, true whether or not
# the frame is in force, so `deactivate_classes!` settles them once, around this.
function push_deactivations!(sat::SAT)
    isempty(sat.deact) && return sat
    sat_push(sat)
    for l in sat.deact
        PicoSAT.add(sat.pico, l)
        PicoSAT.add(sat.pico, 0)
    end
    return sat
end

# Run `body` with every class of the universe choosable — the deactivating
# clauses retracted — and put the frame back afterwards, whatever `body` does,
# so the instance is as reusable as it was. Returns `body`'s value.
#
# Inside, activation is entirely a matter of assumption: assuming
# `forbidden_lit(sat, p, c)` forbids class `c` of `p` for one solve, and
# assuming all of `sat.deact` states exactly what the frame stated, so any
# subset of the query's deactivations can be imposed one solve at a time with no
# clause added and none rewritten. Restoring re-asserts the same unit clauses,
# because popping a frame discards the clauses asserted in it.
#
# The deactivation frame has to be the innermost one, since `sat_pop` retracts
# whatever was pushed last: this cannot run inside `with_temp_clauses`, though
# `body` may open one. Both halves of that are checked rather than trusted —
# lifting one frame too deep retracts someone else's clauses, and a `body` that
# pushes without popping leaves the instance a frame deeper than it found it,
# and neither says so on its own.
function with_classes_relaxed(body::Function, sat::SAT)
    depth = sat.depth
    @assert depth == (isempty(sat.deact) ? 0 : 1) """
        the deactivation frame is not the innermost one — $depth frames are \
        open where $(isempty(sat.deact) ? 0 : 1) is expected"""
    isempty(sat.deact) || sat_pop(sat) # nothing forbidden: no frame to lift
    try body()
    finally
        push_deactivations!(sat)
        @assert sat.depth == depth balance_error(sat, depth)
    end
end

function finalize(sat::SAT)
    pico = sat.pico
    pico == C_NULL && return
    sat.pico = C_NULL
    PicoSAT.reset(pico)
end

sat_new_variable(sat::SAT) = PicoSAT.inc_max_var(sat.pico)

sat_add_var(sat::SAT, v::Integer) = PicoSAT.add(sat.pico, v)
sat_assume_var(sat::SAT, v::Integer) = PicoSAT.assume(sat.pico, v)

sat_add(sat::SAT) = sat_add_var(sat, 0)
sat_add(sat::SAT{P}, p::P, i::Integer=0; not::Bool=false) where {P} =
    sat_add_var(sat, (-1)^not*(sat.vars[p] + i))

sat_assume(sat::SAT{P}, p::P, i::Integer=0; not::Bool=false) where {P} =
    sat_assume_var(sat, (-1)^not*(sat.vars[p] + i))
sat_assume(sat::SAT{P}, v::SetOrVec{P}) where {P} =
    foreach(p -> sat_assume(sat, p), v)
sat_assume(sat::SAT{P}, d::Dict{P,<:Integer}) where {P} =
    for (p, i) in d
        sat_assume(sat, p, i)
    end

# The two facts about this instance a caller can put to it as a raw literal, in
# the `Int` form `sat_assume_var` and the UnsatCores primitives take:

# "package p is installed"
installed_lit(sat::SAT{P}, p::P) where {P} = sat.vars[p]

# "class c of package p is forbidden" — the same literal `deactivate_classes!`
# asserts as a unit clause for a class this query admits no member of, so
# assuming it says what that clause says, for one solve
forbidden_lit(sat::SAT{P}, p::P, c::Integer) where {P} = -(sat.vars[p] + c)

is_satisfiable(sat::SAT) =
    PicoSAT.sat(sat.pico) == PicoSAT.SATISFIABLE

function is_satisfiable(sat::SAT{P}, reqs::Union{P,SetOrVec{P}}) where {P}
    sat_assume(sat, reqs)
    is_satisfiable(sat)
end

const is_unsatisfiable = !is_satisfiable

# iterate the current solution's package => class-index assignments by
# dereferencing the solver's assignment without re-solving: only valid right
# after a solve that returned satisfiable, with no clauses added since --
# ensuring that is the caller's responsibility
function each_solution_index(f::Function, sat::SAT)
    for (p, v_p) in sat.vars
        PicoSAT.deref(sat.pico, v_p) < 0 && continue
        i = 1
        while true
            if PicoSAT.deref(sat.pico, v_p + i) > 0
                # guaranteed to happen by SAT construction:
                # v_p => v_p + i for some i = 1:n_p
                f(p, i)
                break
            end
            i += 1
        end
    end
end

# extract the solver's current solution; same contract as
# each_solution_index: the last solve must have returned satisfiable,
# with no clauses added since
function extract_solution!(sat::SAT{P}, sol::Dict{P,Int}) where {P}
    empty!(sol)
    each_solution_index(sat) do p, i
        sol[p] = i
    end
    return sol
end

function solution(sat::SAT{P,V}) where {P,V}
    sol = Dict{P,V}()
    is_satisfiable(sat) || return sol
    each_solution_index(sat) do p, i
        sol[p] = sat.info[p].versions[sat.reps[p][i]]
    end
    return sol
end

sat_push(sat::SAT) = (sat.depth += 1; PicoSAT.push(sat.pico))
sat_pop(sat::SAT) = (sat.depth -= 1; PicoSAT.pop(sat.pico))

# What an unbalanced body did, said the same way wherever it is caught. An
# exception thrown from a `finally` while another is in flight does not lose it:
# Julia carries both and reports the original as the cause.
balance_error(sat::SAT, depth::Int) =
    "the body left $(sat.depth - depth) push frame(s) open"

function with_temp_clauses(body::Function, sat::SAT)
    depth = sat.depth
    sat_push(sat)
    try body()
    finally
        sat_pop(sat)
        @assert sat.depth == depth balance_error(sat, depth)
    end
end

# improve package p to its best feasible class: repeatedly demand some
# strictly better class until that becomes unsatisfiable, keeping the last
# good model in sol
function optimize_version!(
    sat :: SAT{P},
    sol :: Dict{P,Int},
    p   :: P,
) where {P}
    # cheap first try: after filtering, the best class is feasible more
    # often than not, and probing it by assumption needs no temp clauses
    # at all — one solve replaces the whole improvement loop
    if sol[p] > 1
        sat_assume(sat, p, 1)
        if is_satisfiable(sat)
            extract_solution!(sat, sol)
            @assert sol[p] == 1
        end
    end
    # sol[p] == 1 is already optimal: the improvement clause below would be
    # empty, forcing a guaranteed-UNSAT solve
    while sol[p] > 1
        improved = with_temp_clauses(sat) do
            # some strictly better class of p
            for i = 1:sol[p]-1
                sat_add(sat, p, i)
            end
            sat_add(sat)
            is_satisfiable(sat) || return false
            extract_solution!(sat, sol)
            return true
        end
        improved || break
    end
end

# fix a given set of versions — as the classes that hold them, which is all
# the instance can express (and all that is needed: a version's class is
# interchangeable with it)

function fix_versions(
    sat  :: SAT{P,V},
    pkgs :: AbstractVector{P},
    vers :: AbstractVector{Union{Nothing,V}},
) where {P,V}
    for (p, v) in zip(pkgs, vers)
        v === nothing && continue
        info_p = sat.info[p]
        i = findfirst(==(v), info_p.versions)
        i === nothing && throw(ArgumentError("package $p: unknown version $v"))
        sat_add(sat, p, info_p.classes[i])
        sat_add(sat)
    end
end
