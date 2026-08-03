# find the optimal solution determined by the priority ordering `by`:
# a package => version-index dict, or nothing when the requirements are
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
        # optimistic joint probe: assume the best version of every
        # package in the layer that the current model doesn't already
        # have at its best. when jointly feasible — the common case —
        # one solve pins the whole layer, and joint feasibility of the
        # best versions witnesses each of the sequential optimizations
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
            # fix optimized version
            sat_add(sat, p, sol[p])
            sat_add(sat)
        end
        # next optimization set: new dependencies of the chosen versions
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
    resolve(data, prob::Problem; with = ..., without = ..., diagnose = true)
        -> Union{Dict{P,V}, Diagnosis{P,V}, Nothing}
    resolve(data, reqs; by = identity, diagnose = true, compat = ..., pins = ...)

Resolve the requirements `reqs` against the package universe described by
`data`, returning the optimal solution as a dict mapping each needed package
to its chosen version. The solution covers every requirement and is exactly the
dependency closure of the requirements under its own chosen versions.

When the requirements are not jointly satisfiable the result is a
[`Diagnosis`](@ref) — the independent conflicts, the facts that tell each one's
story, and a menu of verified fixes — which `show` renders as a report. Pass
`diagnose = false` to skip computing it and get `nothing` instead, which is all
a caller that only wants the verdict needs and costs nothing to produce.

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

The other methods accept two more keywords:

  * `order`: the *version* preference ordering, as a callable mapping a package
    to a `lt` comparator over its versions ("is preferred to"). The default,
    `nothing`, means the order the universe already lists them in — which for a
    T1 artifact (see [`pkg_info`](@ref Resolver.pkg_info)) is the canonical one
    the registry state fixes. The ordering is not a constraint: it selects among
    the valid solutions rather than changing which ones are valid, so it is a
    parameter here instead of part of the `Problem`, and one artifact serves
    every ordering.
  * `group::Bool = true`: whether to collapse each package's interchangeable
    versions to one representative before solving (see `prepare_pkg_info`).
    Grouping is an optimization and cannot change the answer; `group = false` is
    the escape hatch that says so.

## Goals

`with` and `without` say what the resolve is *for*, beyond satisfying `prob`:

```julia
resolve(info, prob; with = "CUDA")                        # present
resolve(info, prob; with = "DataFrames" => v"1.8.2")      # present, at a version
resolve(info, prob; without = "FillArrays")               # absent
resolve(info, prob; with = ["CUDA" => spec("5"), "MKL"], without = "GLMakie")
```

`with` takes a package, a `pkg => versions` pair (any `in`-supporting set, or a
bare version), or a collection of either; `without` takes a package or a
collection. Together they express a conjunction.

The return contract is unchanged: the optimal solution satisfying `prob` *and*
the goal, a [`Diagnosis`](@ref) of why it cannot be had, or `nothing` under
`diagnose = false` — which with a goal is the cheap feasibility probe. (Cheaper
still, when the versions are not wanted either: [`issatisfiable`](@ref).)

What separates `with = "CUDA"` from adding CUDA to `prob.reqs` is one sentence:
`prob`'s constraints are *negotiable* — a diagnosis may propose relaxing any of
them — and the goal is the *fixed ask*, so "drop requirement CUDA" can never
appear as a fix.

A goal is answered on a sub-instance built from the T1 artifact rather than on
the per-resolve one (see `prepare_goal_info`), because the per-resolve filter
prunes exactly the escape routes a goal needs. When `data` is a provider or a
`PkgData` dict the closure is widened to cover the goal's packages; a
caller-supplied `PkgInfo` artifact is used as it stands, so a goal naming a
package outside it reads as unsatisfiable.
"""
function resolve(
    sat  :: SAT{P,V},
    reqs :: SetOrVec{P} = keys(sat.info);
    by   :: Function = identity, # priority ordering
    restore :: Bool = true, # restore the SAT instance's state before returning
    diagnose :: Bool = true, # explain an unsatisfiable result
) where {P,V}
    # a `with` goal's packages are in the solution by fiat, so the descent has
    # to reach and optimize them the way it does a requirement; otherwise they
    # are true in the model and pruned out of the answer
    sol = resolve_core(sat, descent_roots(sat, reqs); by, restore)
    if sol === nothing
        # the failed check added no clauses, so the instance is still exactly
        # the one that failed: diagnose it in place, before anything else
        # touches it
        diagnose || return nothing
        return Resolver.diagnose(sat, P[p for p in reqs]; by)
    end
    return Dict{P,V}(p => sat.info[p].versions[i] for (p, i) in sol)
end

# resolve a universe that `prepare_pkg_info` has already laid out, collapsed
# & filtered
function resolve_prepared(
    info :: Dict{P,PkgInfo{P,V}},
    prob :: Problem{P};
    by   :: Function = identity, # package ordering
    diagnose :: Bool = true,
    goal :: Union{Nothing,Goal{P}} = nothing,
) where {P,V}
    sat = SAT(info, prob, goal)
    # the instance is single-use, so don't bother restoring its state -- but
    # the diagnosis, if any, runs on it before it is freed
    try resolve(sat, prob.reqs; by, restore=false, diagnose)
    finally
        finalize(sat)
    end
end

# the packages the layered descent optimizes and keeps: the requirements, plus
# whatever a `with` goal insists on having
function descent_roots(sat::SAT{P}, reqs::SetOrVec{P}) where {P}
    goal = sat.goal
    goal === nothing && return reqs
    isempty(goal.with) && return reqs
    sort!(unique!(P[collect(reqs);
                    P[p for (p, _) in goal.with if haskey(sat.info, p)]]))
end

## goal routing
#
# A goal is answered on a different universe than the goal-⊤ resolve is. The
# per-resolve preparation keeps only what an *optimal* solution could use, which
# is exactly what a goal asks to step outside of: `A→B`, `B@3.1→X`, `B@2.3→∅`
# leaves only `B@3.1` reachable, so the prepared instance calls X unavoidable
# while the raw problem avoids it via `B@2.3`. So a goal query runs on a
# sub-instance built from T1 instead -- see `prepare_goal_info`, and the theory
# page's goal-safety section.
#
# The closure has to cover the goal's packages too. A package the universe does
# not have at all stays in the goal rather than being dropped from it: `with` on
# a package that does not exist is unsatisfiable, and saying so beats answering
# a question nobody asked.

goal_universe(deps::DepsProvider{P}, prob, goal; kws...) where {P} =
    goal_universe(pkg_info(deps, goal_reqs(prob, goal, Set{P}(deps.packages))),
                  prob, goal; kws...)

goal_universe(data::AbstractDict{P,<:PkgData{P}}, prob, goal; kws...) where {P} =
    goal_universe(pkg_info(data, goal_reqs(prob, goal, keys(data))),
                  prob, goal; kws...)

goal_universe(info::AbstractDict{P,PkgInfo{P,V}}, prob, goal;
              group = true, order = nothing) where {P,V} =
    prepare_goal_info(info, prob, goal; group, order)

goal_reqs(prob::Problem{P}, goal::Goal{P}, known) where {P} =
    sort!(unique!(P[prob.reqs; P[p for p in goal_packages(goal) if p ∈ known]]))

"""
    issatisfiable(data, prob::Problem; with = ..., without = ...) :: Bool
    issatisfiable(data, reqs; with = ..., without = ..., compat = ..., pins = ...)

Can `prob` be satisfied — and, with a goal, satisfied *while* getting what the
goal asks for? One solve, no optimizing descent, no explanation: the cheapest
question in the API, for the caller that is asking it in a loop.

This is the bottom rung of a three-tier ladder, all three answering the same
question with more or less work:

| call | cost | answer |
|:--|:--|:--|
| `issatisfiable(data, prob; with)` | one solve | `Bool` |
| `resolve(data, prob; with, diagnose = false)` | one solve + descent | the solution, or `nothing` |
| `resolve(data, prob; with)` | + explanation | the solution, or a [`Diagnosis`](@ref) |

Nothing distinguishes them but effort — the verdicts agree — so a UI greying out
an upgrade button uses the first, a caller that wants the versions uses the
second, and only the one that will actually show a report pays for the third.
"""
function issatisfiable(
    data,
    prob :: Problem{P};
    group :: Bool = true,
    order = nothing,
    with = nothing,
    without = nothing,
) where {P}
    goal = make_goal(P, with, without)
    info = goal === nothing ?
        prepare_pkg_info_for(data, prob; group, order) :
        goal_universe(data, prob, goal; group, order)
    all(p -> haskey(info, p), prob.reqs) || return false
    sat = SAT(info, prob, goal)
    try is_satisfiable(sat, prob.reqs)
    finally
        finalize(sat)
    end
end

issatisfiable(
    data,
    reqs :: SetOrVec{P};
    compat :: AbstractDict{P} = EmptyDict{P,Any}(),
    pins   :: AbstractDict{P} = EmptyDict{P,Any}(),
    group  :: Bool = true,
    order  = nothing,
    with = nothing,
    without = nothing,
) where {P} = issatisfiable(data, Problem(reqs; compat, pins);
                            group, order, with, without)

# the goal-⊤ universe, per data flavour -- the same three shapes `resolve` has
prepare_pkg_info_for(deps::DepsProvider{P}, prob; group, order) where {P} =
    (info = pkg_info(deps, prob);
     prepare_pkg_info(info, prob, info; group, order))
prepare_pkg_info_for(data::AbstractDict{P,<:PkgData{P}}, prob;
                     group, order) where {P} =
    (info = pkg_info(data, prob);
     prepare_pkg_info(info, prob, info; group, order))
prepare_pkg_info_for(info::AbstractDict{P,PkgInfo{P,V}}, prob;
                     group, order) where {P,V} =
    prepare_pkg_info(info, prob; group, order)

# primary entry points: a Problem (requirements + user constraints). each
# builds the T1 artifact for the problem's requirements, then prepares it —
# the T1 info is freshly built and nobody else holds it, so preparation is
# allowed to shrink it in place

function resolve(
    deps :: DepsProvider{P},
    prob :: Problem{P};
    by   :: Function = identity, # package ordering
    group :: Bool = true, # collapse interchangeable versions
    order = nothing, # version ordering
    diagnose :: Bool = true,
    with = nothing,
    without = nothing,
) where {P}
    goal = make_goal(P, with, without)
    goal === nothing || return resolve_prepared(
        goal_universe(deps, prob, goal; group, order), prob; by, diagnose, goal)
    info = pkg_info(deps, prob)
    resolve_prepared(prepare_pkg_info(info, prob, info; group, order), prob;
                     by, diagnose)
end

function resolve(
    data :: AbstractDict{P,<:PkgData{P}},
    prob :: Problem{P};
    by   :: Function = identity, # package ordering
    group :: Bool = true, # collapse interchangeable versions
    order = nothing, # version ordering
    diagnose :: Bool = true,
    with = nothing,
    without = nothing,
) where {P}
    goal = make_goal(P, with, without)
    goal === nothing || return resolve_prepared(
        goal_universe(data, prob, goal; group, order), prob; by, diagnose, goal)
    info = pkg_info(data, prob)
    resolve_prepared(prepare_pkg_info(info, prob, info; group, order), prob;
                     by, diagnose)
end

# a caller-supplied info may be a reusable (or cached) T1 artifact, so this
# method prepares into a dict of its own and leaves the argument alone
function resolve(
    info :: AbstractDict{P,PkgInfo{P,V}},
    prob :: Problem{P};
    by   :: Function = identity, # package ordering
    group :: Bool = true, # collapse interchangeable versions
    order = nothing, # version ordering
    diagnose :: Bool = true,
    with = nothing,
    without = nothing,
) where {P,V}
    goal = make_goal(P, with, without)
    goal === nothing || return resolve_prepared(
        goal_universe(info, prob, goal; group, order), prob; by, diagnose, goal)
    resolve_prepared(prepare_pkg_info(info, prob; group, order), prob;
                     by, diagnose)
end

# convenience entry points: bare requirements, with the user constraints (if
# any) given as keywords instead of a `Problem`

resolve(
    deps :: DepsProvider{P},
    reqs :: SetOrVec{P} = deps.packages;
    compat :: AbstractDict{P} = EmptyDict{P,Any}(),
    pins   :: AbstractDict{P} = EmptyDict{P,Any}(),
    by     :: Function = identity, # package ordering
    group  :: Bool = true, # collapse interchangeable versions
    order  = nothing, # version ordering
    diagnose :: Bool = true,
    with = nothing,
    without = nothing,
) where {P} = resolve(deps, Problem(reqs; compat, pins);
                       by, group, order, diagnose, with, without)

resolve(
    data :: AbstractDict{P,<:PkgData{P}},
    reqs :: SetOrVec{P} = keys(data);
    compat :: AbstractDict{P} = EmptyDict{P,Any}(),
    pins   :: AbstractDict{P} = EmptyDict{P,Any}(),
    by     :: Function = identity, # package ordering
    group  :: Bool = true, # collapse interchangeable versions
    order  = nothing, # version ordering
    diagnose :: Bool = true,
    with = nothing,
    without = nothing,
) where {P} = resolve(data, Problem(reqs; compat, pins);
                       by, group, order, diagnose, with, without)

resolve(
    info :: AbstractDict{P,PkgInfo{P,V}},
    reqs :: SetOrVec{P} = keys(info);
    compat :: AbstractDict{P} = EmptyDict{P,Any}(),
    pins   :: AbstractDict{P} = EmptyDict{P,Any}(),
    by     :: Function = identity, # package ordering
    group  :: Bool = true, # collapse interchangeable versions
    order  = nothing, # version ordering
    diagnose :: Bool = true,
    with = nothing,
    without = nothing,
) where {P,V} = resolve(info, Problem(reqs; compat, pins);
                       by, group, order, diagnose, with, without)
