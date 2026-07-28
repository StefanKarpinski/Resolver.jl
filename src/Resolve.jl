# find the optimal solution determined by the priority ordering `by`:
# a package => version-index dict, or nothing when the requirements are
# not jointly satisfiable
function resolve_core(
    sat  :: SAT{P},
    reqs :: SetOrVec{P} = keys(sat.info);
    by   :: Function = identity, # priority ordering
    restore :: Bool = true, # restore the SAT instance's state before returning
) where {P}
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
        for p in sort!(collect(opts); by)
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
    resolve(data, reqs; by = identity) -> Union{Dict{P,V}, Nothing}

Resolve the requirements `reqs` against the package universe described by
`data`, returning the optimal solution as a dict mapping each needed package
to its chosen version, or `nothing` when the requirements are not jointly
satisfiable. The solution covers every requirement and is exactly the
dependency closure of the requirements under its own chosen versions.

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
"""
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

function resolve(
    deps :: DepsProvider{P},
    reqs :: SetOrVec{P} = deps.packages;
    by   :: Function = identity, # package ordering
) where {P}
    info = pkg_info(deps, reqs)
    resolve(info, reqs; by)
end

function resolve(
    data :: AbstractDict{P,<:PkgData{P}},
    reqs :: SetOrVec{P} = keys(data);
    by   :: Function = identity, # package ordering
) where {P}
    info = pkg_info(data, reqs)
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
