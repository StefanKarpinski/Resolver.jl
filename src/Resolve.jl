function resolve(
    sat  :: SAT{P,V},
    reqs :: SetOrVec{P} = keys(sat.info),
) where {P,V}
    # sort reqs for determinism
    reqs = sort!(collect(reqs))

    # solution search data strutures
    sol = Dict{P,Int}() # current solution
    sols = typeof(sol)[] # all solutions

    function find_optimal_solutions(
        opts :: Set{P}, # packages to optimize next
        rest :: Set{P}, # other unoptimized packages
    )
        while is_satisfiable(sat)
            extract_solution!(sat, sol)

            # optimize wrt coverage
            opts ⊆ keys(sol) ||
            optimize_solution!(sat, sol) do
                # clauses: disallow non-improvements
                for p in opts
                    p in keys(sol) || continue
                    sat_add(sat, p)
                    sat_add(sat)
                end
                # clause: require some improvement
                for p in opts
                    p in keys(sol) && continue
                    sat_add(sat, p)
                end
                sat_add(sat)
            end

            # optimize wrt quality
            optimize_solution!(sat, sol) do
                # clauses: disallow non-improvements
                for p in opts
                    p in keys(sol) || continue
                    # allow as-good-or-better versions
                    for i = 1:sol[p]
                        sat_add(sat, p, i)
                    end
                    sat_add(sat)
                end
                # clause: require some improvement
                for p in opts
                    p in keys(sol) || continue
                    # any strictly better versions
                    for i = 1:sol[p]-1
                        sat_add(sat, p, i)
                    end
                end
                sat_add(sat)
            end

            # find next optimization set
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

            if !isempty(opts′)
                # recursion required
                with_temp_clauses(sat) do
                    # clauses: fix already optimized versions
                    for p in opts
                        p in keys(sol) || continue
                        sat_add(sat, p, sol[p])
                        sat_add(sat)
                    end
                    # recursive search call
                    rest′ = setdiff(rest, opts′)
                    find_optimal_solutions(opts′, rest′)
                end
            else # nothing left to optimize
                for p in rest
                    # delete unreachable packages
                    delete!(sol, p)
                end
                # save optimal solution
                push!(sols, copy(sol))
            end

            # clause: require some improvement
            #   unlike the optimization process, this allows other
            #   aspects to get worse; it only forces non-domination

            if opts ⊆ keys(sol) # complete solution
                # preserve coverage
                for p in opts
                    sat_add(sat, p)
                    sat_add(sat)
                end
                # improve quality
                for p in opts
                    # any strictly better versions
                    for i = 1:sol[p]-1
                        sat_add(sat, p, i)
                    end
                end
                sat_add(sat)
            else # incomplete solution
                # improve coverage
                for p in opts
                    p in keys(sol) && continue
                    sat_add(sat, p)
                end
                sat_add(sat)
            end
        end
    end

    with_temp_clauses(sat) do
        if is_satisfiable(sat, reqs)
            # force all requirements
            for p in reqs
                sat_add(sat, p)
            end
            sat_add(sat)
        else
            # force some requirement
            for p in reqs
                sat_add(sat, p)
                sat_add(sat)
            end
        end

        # start recursion
        opts = Set{P}(reqs)
        rest = setdiff(keys(sat.info), opts)
        find_optimal_solutions(opts, rest)
    end

    # sort packages
    if !isempty(sols)
        pkgs = sort!(collect(mapreduce(keys, union, sols, init=reqs)))
        sort!(pkgs, by = !in(reqs)) # required ones first
    else
        pkgs = reqs
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
    filter :: Bool = true,
) where {P}
    info = pkg_info(deps, reqs; filter)
    @timeit "resolve" resolve(info, reqs)
end

function resolve(
    data :: AbstractDict{P,<:PkgData{P}},
    reqs :: SetOrVec{P} = keys(data);
    filter :: Bool = true,
) where {P}
    info = pkg_info(data, reqs; filter)
    @timeit "resolve" resolve(info, reqs)
end

function resolve(
    info :: Dict{P,PkgInfo{P,V}},
    reqs :: SetOrVec{P} = keys(info),
) where {P,V}
    sat = SAT(info)
    pkgs, vers = resolve(sat, reqs)
    finalize(sat)
    return pkgs, vers
end

#=
using PrettyTables, Revise, Resolver
deps = registry_provider();
(str -> begin
reqs = String.(split(str, ','))
pkgs, vers = resolve(deps, reqs)
mv = vec(mapslices(r->maximum(v->something(v, v"0+1-"), r), vers, dims=2))
pretty_table([pkgs vers],
    highlighters = Highlighter(
        (data, i, j) -> j > 1 &&
            something(data[i, j], i ≤ length(reqs) ? v"0+0-" : mv[i]) < mv[i],
        foreground = :red,
    )
)
end)("CEnum,CUDAdrv,HELICS,DocStringExtensions,FunctionalTables")
=#
