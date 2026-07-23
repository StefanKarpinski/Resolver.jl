DEFAULT_MAX_SOLUTIONS::Int = 8

# find the first `max` optimal solutions in lexicographical ordering
function resolve_core(
    sat  :: SAT{P},
    reqs :: SetOrVec{P} = keys(sat.info);
    max  :: Integer = DEFAULT_MAX_SOLUTIONS,
    by   :: Function = identity, # ordering
    restore :: Bool = true, # restore the SAT instance's state before returning
) where {P}
    # solution search data structures
    sol = Dict{P,Int}() # current solution
    sols = typeof(sol)[] # all solutions

    # optimize each package in `opts` to its best feasible version wrt
    # quality, pinning each optimized version as it goes (in priority order)
    function optimize_opts!(opts::Set{P})
        for p in sort!(collect(opts); by)
            # sol[p] == 1 is already optimal: the improvement clause
            # below would be empty, forcing a guaranteed-UNSAT solve
            if sol[p] > 1
                optimize_solution!(sat, sol) do
                    # some strictly better version
                    for i = 1:sol[p]-1
                        sat_add(sat, p, i)
                    end
                    sat_add(sat)
                end
            end
            # fix optimized version
            sat_add(sat, p, sol[p])
            sat_add(sat)
        end
    end

    # next optimization set: new dependencies of the versions chosen in `opts`
    function next_opts(opts::Set{P}, rest::Set{P})
        opts′ = empty(opts)
        for p in opts
            p in keys(sol) || continue
            info_p = sat.info[p]
            i = sol[p]
            for (j, q) in enumerate(info_p.depends)
                info_p.conflicts[i, j] || continue
                q ∈ rest && push!(opts′, q)
            end
        end
        return opts′
    end

    # record `sol` as a solution, minus unreachable packages
    function save_solution!(rest::Set{P})
        for p in rest
            delete!(sol, p)
        end
        push!(sols, copy(sol))
    end

    function find_sat_solutions(
        opts :: Set{P}, # packages to optimize next
        rest :: Set{P}, # other unoptimized packages
    )
        while is_satisfiable(sat)
            extract_solution!(sat, sol)
            @assert opts ⊆ keys(sol)

            with_temp_clauses(sat) do
                optimize_opts!(opts)
                opts′ = next_opts(opts, rest)
                if !isempty(opts′) # recursive search
                    rest′ = setdiff(rest, opts′)
                    find_sat_solutions(opts′, rest′)
                else # nothing left to optimize
                    save_solution!(rest)
                end
            end
            # max if we've hit max number of solutions
            0 < max ≤ length(sols) && break

            # next solution must be undominated by this one
            # require: some package with better version than now
            for p in opts
                # any strictly better versions
                for i = 1:sol[p]-1
                    sat_add(sat, p, i)
                end
            end
            sat_add(sat)
        end
    end

    # specialized descent for max == 1: the enumeration machinery above needs
    # nested temp scopes and a re-solve per dependency layer so it can rewind
    # and search for further solutions; when only the first solution is wanted
    # neither is necessary. With restore = false it runs outside any temp
    # scope, which lets picosat keep its learned clauses for the whole search
    # and is substantially faster on large instances.
    function find_first_solution(
        opts :: Set{P}, # packages to optimize next
        rest :: Set{P}, # other unoptimized packages
    )
        is_satisfiable(sat) || return
        extract_solution!(sat, sol)
        while true
            @assert opts ⊆ keys(sol)
            optimize_opts!(opts)
            opts′ = next_opts(opts, rest)
            isempty(opts′) && break
            setdiff!(rest, opts′)
            opts = opts′
        end
        save_solution!(rest)
    end

    function find_unsat_solutions(
        opts :: Set{P}, # packages to optimize next
        rest :: Set{P}, # other unoptimized packages
    )
        while is_satisfiable(sat)
            extract_solution!(sat, sol)
            @assert opts ⊈ keys(sol)

            with_temp_clauses(sat) do
                # optimize wrt coverage
                for p in sort!(collect(opts); by)
                    sat_assume(sat, p)
                    is_satisfiable(sat) || continue
                    sat_add(sat, p)
                    sat_add(sat)
                end

                # next optimization set: subset of satisfied opts
                extract_solution!(sat, sol)
                opts′ = filter(in(keys(sol)), opts)
                @assert !isempty(opts′)

                # recursive search call
                find_sat_solutions(opts′, rest)
            end
            # max if we've hit max number of solutions
            0 < max ≤ length(sols) && break

            # next solution must be undominated by this one
            # require: some package covered that isn't now
            for p in opts
                p in keys(sol) && continue
                sat_add(sat, p)
            end
            sat_add(sat)
        end
    end

    # start the recursive search...
    function search()
        opts = Set{P}(reqs)
        rest = setdiff(keys(sat.info), opts)
        if is_satisfiable(sat, reqs)
            # force all requirements
            for p in reqs
                sat_add(sat, p)
                sat_add(sat)
            end
            if max == 1
                find_first_solution(opts, rest)
            else
                find_sat_solutions(opts, rest)
            end
        else
            # force some requirement
            for p in reqs
                sat_add(sat, p)
            end
            sat_add(sat)
            find_unsat_solutions(opts, rest)
        end
    end
    if max == 1 && !restore
        # run at the top level: keeping a picosat context open for the whole
        # search measurably degrades solving. The instance is left constrained
        # by the requirements and the solution's version pins, so this mode is
        # only for single-use instances (as in the info-level resolve methods)
        search()
    else
        with_temp_clauses(search, sat)
    end

    # return solutions
    return sols
end

function resolve(
    sat  :: SAT{P,V},
    reqs :: SetOrVec{P} = keys(sat.info);
    max  :: Integer = DEFAULT_MAX_SOLUTIONS,
    by   :: Function = identity, # ordering
    restore :: Bool = true, # restore the SAT instance's state before returning
) where {P,V}
    # generate solutions
    sols = resolve_core(sat, reqs; max, by, restore)

    # sort packages
    pkgs = collect(reqs)
    @assert pkgs !== reqs
    if !isempty(sols)
        pkgs = mapreduce(keys, union!, sols, init=pkgs)
        sort!(pkgs; by)
        sort!(pkgs, by = !in(reqs)) # required ones first
    end

    # sort solutions
    for p in reverse(pkgs)
        sort!(sols, by = sol -> get(sol, p, 0))
    end

    # versions matrix
    vers = Union{V,Nothing}[
        haskey(sol, p) ? sat.info[p].versions[sol[p]] : nothing
        for p in pkgs, sol in sols
    ]

    # return solutions
    return pkgs::Vector{P}, vers::Matrix{Union{V,Nothing}}
end

# convenience entry points

function resolve(
    deps :: DepsProvider{P},
    reqs :: SetOrVec{P} = deps.packages;
    max  :: Integer = DEFAULT_MAX_SOLUTIONS,
    by   :: Function = identity, # package ordering
) where {P}
    info = pkg_info(deps, reqs)
    resolve(info, reqs; max, by)
end

function resolve(
    data :: AbstractDict{P,<:PkgData{P}},
    reqs :: SetOrVec{P} = keys(data);
    max  :: Integer = DEFAULT_MAX_SOLUTIONS,
    by   :: Function = identity, # package ordering
) where {P}
    info = pkg_info(data, reqs)
    resolve(info, reqs; max, by)
end

function resolve(
    info :: Dict{P,PkgInfo{P,V}},
    reqs :: SetOrVec{P} = keys(info);
    max  :: Integer = DEFAULT_MAX_SOLUTIONS,
    by   :: Function = identity, # package ordering
) where {P,V}
    sat = SAT(info)
    # the instance is single-use, so no need to restore its state
    try resolve(sat, reqs; max, by, restore=false)
    finally
        finalize(sat)
    end
end
