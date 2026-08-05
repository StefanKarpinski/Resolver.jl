# find the optimal solution determined by the priority ordering `by`:
# a package => class-index dict, or nothing when the requirements are
# not jointly satisfiable
function resolve_core(
    sat  :: SAT{P},
    reqs :: SetOrVec{P} = keys(sat.info);
    by   :: Function = identity, # priority ordering
    restore :: Bool = true, # restore the SAT instance's state before returning
) where {P}
    # a requirement the instance doesn't know has no installable version —
    # the filter drops such packages — so it cannot be satisfied
    all(p -> haskey(sat.info, p), reqs) || return nothing
    # a solution exists iff the requirements are jointly satisfiable; the
    # check only assumes, so an unsatisfiable resolve leaves no state behind
    is_satisfiable(sat, reqs) || return nothing
    sol = Dict{P,Int}() # the solution
    # capture the solution (it covers all the requirements) before the
    # clauses added by the descent invalidate the solver's assignment
    extract_solution!(sat, sol)

    # optimize each package in `opts` to its best feasible version wrt quality,
    # pinning each as it goes (in priority order); return the newly-reachable
    # dependencies of the chosen versions that still need optimizing
    function optimize!(opts::Set{P}, seen::Set{P})
        layer = sort!(collect(opts); by)
        # optimistic joint probe: assume the best class of every
        # package in the layer that the current model doesn't already
        # have at its best. when jointly feasible — the common case —
        # one solve pins the whole layer, and joint feasibility of the
        # best classes witnesses each of the sequential optimizations
        # below, so the layered answer is unchanged. a failed probe
        # tells us nothing and the sequential path proceeds as usual
        # (skipped for a single package, whose own probe covers it)
        todo = [p for p in layer if sol[p] > 1]
        if length(todo) > 1
            for p in todo
                sat_assume(sat, p, 1)
            end
            is_satisfiable(sat) && extract_solution!(sat, sol)
        end
        for p in layer
            optimize_version!(sat, sol, p)
            # fix optimized class
            sat_add(sat, p, sol[p])
            sat_add(sat)
        end
        # next optimization set: new dependencies of the chosen classes
        opts′ = empty(opts)
        for p in opts
            info_p = sat.info[p]
            i = sol[p]
            for (j, q) in enumerate(info_p.depends)
                info_p.conflicts[i, j] || continue
                q ∉ seen && push!(opts′, q)
            end
        end
        return opts′
    end

    function descend!()
        # force all requirements
        for p in reqs
            sat_add(sat, p)
            sat_add(sat)
        end
        # optimize quality in priority order, layer by layer
        opts = Set{P}(reqs)
        seen = copy(opts)
        while true
            @assert opts ⊆ keys(sol)
            opts = optimize!(opts, seen)
            isempty(opts) && break
            union!(seen, opts)
        end
        # prune to the packages the descent reached: the SAT model may set
        # variables of unreachable packages true
        filter!(kv -> first(kv) in seen, sol)
    end

    # `restore = true` rolls the instance's clauses back (reusable-instance
    # contract); `restore = false` runs at the top level, leaving the instance
    # constrained by the requirements and the solution's pins — measurably
    # faster (picosat keeps its learned clauses), for single-use instances only.
    if restore
        with_temp_clauses(descend!, sat)
    else
        descend!()
    end

    return sol
end

"""
    resolve(data, prob::Problem; by = identity) -> Union{Dict{P,V}, Nothing}
    resolve(data, reqs; by = identity, compat = ..., pins = ...)

Resolve the requirements `reqs` against the package universe described by
`data`, returning the optimal solution as a dict mapping each needed package
to its chosen version, or `nothing` when the requirements are not jointly
satisfiable. The solution covers every requirement and is exactly the
dependency closure of the requirements under its own chosen versions.

A [`Problem`](@ref) additionally constrains the admissible versions with user
compat bounds and pins; the answer is the one the same universe with those
versions deleted would give. The requirements-and-keywords form builds one, so
`resolve(data, reqs; compat, pins)` and `resolve(data, Problem(reqs; compat,
pins))` are the same call.

The returned solution is the *layered solution*: requirements are optimized
first, in the priority order induced by `by` (each package maximized to its
best feasible version given the choices already made), then the dependencies
of the chosen versions, layer by layer along the dependency graph. See the
manual's Theory section for the precise semantics and the guarantees it
carries.

`data` may be a `DepsProvider`, a dict of `PkgData`, a dict of `PkgInfo`, or
a `SAT` instance. The `SAT` method also accepts `restore::Bool = true`: when
true the instance's clauses are restored before returning so the instance
can be reused; `restore = false` is faster for single-use instances.

The other methods accept one more keyword:

  * `order`: the *version* preference ordering, as a callable mapping a package
    to a `lt` comparator over its versions ("is preferred to"). The default,
    `nothing`, means the order the universe already lists them in — which for a
    T1 artifact (see [`pkg_info`](@ref Resolver.pkg_info)) is the canonical one
    the registry state fixes. The ordering is not a constraint: it selects among
    the valid solutions rather than changing which ones are valid, so it is a
    parameter here instead of part of the `Problem`, and one artifact serves
    every ordering.

    The universe is built out of *interchangeability classes* — sets of versions
    nothing in the registry can tell apart (see
    [`pkg_info`](@ref Resolver.pkg_info)) — so `order` decides two things: which
    class the answer picks, and which member of that class it names, namely the
    best one the query admits.
"""
function resolve(
    sat  :: SAT{P,V},
    reqs :: SetOrVec{P} = keys(sat.info);
    by   :: Function = identity, # priority ordering
    restore :: Bool = true, # restore the SAT instance's state before returning
) where {P,V}
    sol = resolve_core(sat, reqs; by, restore)
    sol === nothing && return nothing
    # a solution names classes; the version each stands for is `reps`
    return Dict{P,V}(
        p => sat.info[p].versions[sat.reps[p][i]] for (p, i) in sol)
end

# resolve a universe that `prepare_pkg_info` has already ranked & filtered
function resolve_prepared(
    univ :: Universe{P,V},
    prob :: Problem{P};
    by   :: Function = identity, # package ordering
) where {P,V}
    sat = SAT(univ)
    # the instance is single-use, so don't bother restoring its state
    try resolve(sat, prob.reqs; by, restore=false)
    finally
        finalize(sat)
    end
end

# primary entry points: a Problem (requirements + user constraints). each
# builds the T1 artifact for the problem's requirements, then prepares it —
# the T1 info is freshly built and nobody else holds it, so preparation is
# allowed to shrink it in place

function resolve(
    deps :: DepsProvider{P},
    prob :: Problem{P};
    by   :: Function = identity, # package ordering
    order = nothing, # version ordering
) where {P}
    info = pkg_info(deps, prob)
    resolve_prepared(prepare_pkg_info(info, prob, info; order), prob; by)
end

function resolve(
    data :: AbstractDict{P,<:PkgData{P}},
    prob :: Problem{P};
    by   :: Function = identity, # package ordering
    order = nothing, # version ordering
) where {P}
    info = pkg_info(data, prob)
    resolve_prepared(prepare_pkg_info(info, prob, info; order), prob; by)
end

# a caller-supplied info may be a reusable (or cached) T1 artifact, so this
# method prepares into a dict of its own and leaves the argument alone
function resolve(
    info :: AbstractDict{P,PkgInfo{P,V}},
    prob :: Problem{P};
    by   :: Function = identity, # package ordering
    order = nothing, # version ordering
) where {P,V}
    resolve_prepared(prepare_pkg_info(info, prob; order), prob; by)
end

# convenience entry points: bare requirements, with the user constraints (if
# any) given as keywords instead of a `Problem`

resolve(
    deps :: DepsProvider{P},
    reqs :: SetOrVec{P} = deps.packages;
    compat :: AbstractDict{P} = EmptyDict{P,Any}(),
    pins   :: AbstractDict{P} = EmptyDict{P,Any}(),
    by     :: Function = identity, # package ordering
    order  = nothing, # version ordering
) where {P} = resolve(deps, Problem(reqs; compat, pins); by, order)

resolve(
    data :: AbstractDict{P,<:PkgData{P}},
    reqs :: SetOrVec{P} = keys(data);
    compat :: AbstractDict{P} = EmptyDict{P,Any}(),
    pins   :: AbstractDict{P} = EmptyDict{P,Any}(),
    by     :: Function = identity, # package ordering
    order  = nothing, # version ordering
) where {P} = resolve(data, Problem(reqs; compat, pins); by, order)

resolve(
    info :: AbstractDict{P,PkgInfo{P,V}},
    reqs :: SetOrVec{P} = keys(info);
    compat :: AbstractDict{P} = EmptyDict{P,Any}(),
    pins   :: AbstractDict{P} = EmptyDict{P,Any}(),
    by     :: Function = identity, # package ordering
    order  = nothing, # version ordering
) where {P,V} = resolve(info, Problem(reqs; compat, pins); by, order)

"""
    issatisfiable(data, reqs; compat = ..., pins = ...) -> Bool

Compute `resolve(data, reqs; ...) !== nothing` more efficiently, with a single
satisfiability check instead of a full optimizing resolve. Use this when you
want to know whether the requirements can be satisfied but don't need an
actual solution.

`data` may be a `DepsProvider`, a dict of `PkgData`, a dict of `PkgInfo`, or a
`SAT` instance, and a `Problem` may be passed in place of the requirements and
keywords, as with `resolve`. The `by` and `order` keywords are not accepted
since they affect which solution `resolve` picks, not whether one exists. The `SAT` method leaves the instance unchanged, so it can be reused.
"""
function issatisfiable(
    sat  :: SAT{P},
    reqs :: SetOrVec{P} = keys(sat.info),
) where {P}
    # a requirement the instance doesn't know has no installable version —
    # the filter drops such packages — so it cannot be satisfied
    all(p -> haskey(sat.info, p), reqs) || return false
    # the requirements are satisfiable iff this one solve says so; it only
    # assumes, so the instance is left exactly as it was found
    return is_satisfiable(sat, reqs)
end

# check a universe that `prepare_pkg_info` has already ranked & filtered
function issatisfiable_prepared(
    univ :: Universe{P,V},
    prob :: Problem{P},
) where {P,V}
    sat = SAT(univ)
    try issatisfiable(sat, prob.reqs)
    finally
        finalize(sat)
    end
end

# primary entry points, mirroring `resolve`'s: a Problem against each shape of
# package data, each building the T1 artifact and preparing it the same way

function issatisfiable(
    deps :: DepsProvider{P},
    prob :: Problem{P},
) where {P}
    info = pkg_info(deps, prob)
    issatisfiable_prepared(prepare_pkg_info(info, prob, info), prob)
end

function issatisfiable(
    data :: AbstractDict{P,<:PkgData{P}},
    prob :: Problem{P},
) where {P}
    info = pkg_info(data, prob)
    issatisfiable_prepared(prepare_pkg_info(info, prob, info), prob)
end

# a caller-supplied info may be a reusable (or cached) T1 artifact, so this
# method prepares into a dict of its own and leaves the argument alone
issatisfiable(
    info :: AbstractDict{P,PkgInfo{P,V}},
    prob :: Problem{P},
) where {P,V} = issatisfiable_prepared(prepare_pkg_info(info, prob), prob)

# convenience entry points: bare requirements, with the user constraints (if
# any) given as keywords instead of a `Problem`

issatisfiable(
    deps :: DepsProvider{P},
    reqs :: SetOrVec{P} = deps.packages;
    compat :: AbstractDict{P} = EmptyDict{P,Any}(),
    pins   :: AbstractDict{P} = EmptyDict{P,Any}(),
) where {P} = issatisfiable(deps, Problem(reqs; compat, pins))

issatisfiable(
    data :: AbstractDict{P,<:PkgData{P}},
    reqs :: SetOrVec{P} = keys(data);
    compat :: AbstractDict{P} = EmptyDict{P,Any}(),
    pins   :: AbstractDict{P} = EmptyDict{P,Any}(),
) where {P} = issatisfiable(data, Problem(reqs; compat, pins))

issatisfiable(
    info :: AbstractDict{P,PkgInfo{P,V}},
    reqs :: SetOrVec{P} = keys(info);
    compat :: AbstractDict{P} = EmptyDict{P,Any}(),
    pins   :: AbstractDict{P} = EmptyDict{P,Any}(),
) where {P,V} = issatisfiable(info, Problem(reqs; compat, pins))
