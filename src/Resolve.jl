function resolve(
    deps :: DepsProvider{P},
    reqs :: SetOrVec{P} = deps.packages;
    filter :: Bool = true,
) where {P}
    info = load_pkg_info(deps, reqs; filter)
    @timeit "resolve" resolve(info, reqs)
end

function resolve(
    data :: AbstractDict{P,<:PkgData{P}},
    reqs :: SetOrVec{P} = keys(data);
    filter :: Bool = true,
) where {P}
    info = make_pkg_info(data, reqs; filter)
    @timeit "resolve" resolve(info, reqs)
end

function resolve(
    info :: Dict{P,PkgInfo{P,V}},
    reqs :: SetOrVec{P} = keys(info),
) where {P,V}
    # sort stuff for determinism
    reqs = sort!(collect(reqs))
    names = sort!(collect(keys(info)))
    sort!(names, by = !in(reqs))

    # variable indices:
    #   p    var[p]     package p chosen
    #   p@i  var[p]+i   version i of p chosen
    N = 1
    var = Dict{P, Int}()
    for p in names
        var[p] = N
        n_p = length(info[p].versions)
        # 1 varible for package
        # n_p variables for versions
        N += 1 + n_p
    end
    N -= 1 # last used variable

    # instantiate picosat solver
    ps = PicoSAT.init()
    sols = try
        PicoSAT.adjust(ps, N)

        # generate SAT problem
        for p in names
            info_p = info[p]
            n_p = length(info_p.versions)
            v_p = var[p]

            # package implies some version
            #   p => OR_i p@i
            PicoSAT.add(ps, -v_p)
            for i = 1:n_p
                PicoSAT.add(ps, v_p + i)
            end
            PicoSAT.add(ps, 0)

            # version implies its package
            #   p@i => p
            for i = 1:n_p
                PicoSAT.add(ps, -(v_p + i))
                PicoSAT.add(ps, v_p)
                PicoSAT.add(ps, 0)
            end

            # versions are mutually exclusive
            #   !p@i OR !p@j
            for i = 1:n_p-1, j = i+1:n_p
                PicoSAT.add(ps, -(v_p + i))
                PicoSAT.add(ps, -(v_p + j))
                PicoSAT.add(ps, 0)
            end

            # dependencies
            #   ∀ q ∈ deps(p, i): p@i => q
            for i = 1:n_p
                for (j, q) in enumerate(info_p.depends)
                    info_p.conflicts[i, j] || continue
                    PicoSAT.add(ps, -(v_p + i))
                    PicoSAT.add(ps, var[q])
                    PicoSAT.add(ps, 0)
                end
            end

            # conflicts
            #   ∀ q@j ∈ conflicts(p, i): !p@i OR !q@j
            for i = 1:n_p
                for (q, b) in info_p.interacts
                    p < q || continue # conflicts are symmetrical
                    n_q = length(info[q].versions)
                    v_q = var[q]
                    for j = 1:n_q
                        info_p.conflicts[i, b+j] || continue
                        # conflicting versions are mutually exclusive
                        PicoSAT.add(ps, -(v_p + i))
                        PicoSAT.add(ps, -(v_q + j))
                        PicoSAT.add(ps, 0)
                    end
                end
            end
        end

        # check if problem is satisfiable or not
        for p in reqs
            PicoSAT.assume(ps, var[p])
        end
        if PicoSAT.sat(ps) == PicoSAT.SATISFIABLE
            # all reqs must be satisfied
            for p in reqs
                PicoSAT.add(ps, var[p])
                PicoSAT.add(ps, 0)
            end
        else
            # some reqs must be satisfied
            for p in reqs
                PicoSAT.add(ps, var[p])
            end
            PicoSAT.add(ps, 0)
        end

        # solution search data strutures
        sols = Vector{Dict{P,Int}}() # all solutions
        sol = Dict{P,Int}() # current solution
        sats = Set{P}() # satisfied package set

        # helper for extracting a solution
        function extract_solution(opts::Set{P})
            empty!(sol)
            for (p, v_p) in var
                PicoSAT.deref(ps, v_p) < 0 && continue
                for i = 1:length(info[p].versions)
                    if PicoSAT.deref(ps, v_p + i) > 0
                        sol[p] = i
                        break
                    end
                end
            end
            intersect!(copy!(sats, opts), keys(sol))
            return sol
        end

        # helper for optimizing a solution
        function optimize_solution(improve::Function, opts::Set{P})
            done = false
            while !done
                PicoSAT.push(ps)
                improve() # callback adds SAT clauses
                done = PicoSAT.sat(ps) ≠ PicoSAT.SATISFIABLE
                done || extract_solution(opts) # -> sol, sats
                PicoSAT.pop(ps)
            end
        end

        function find_optimal_solutions(
            opts :: Set{P}, # packages to optimize next
            rest :: Set{P}, # other unoptimized packages
        )
            while PicoSAT.sat(ps) == PicoSAT.SATISFIABLE
                extract_solution(opts) # -> sol, sats
                @info "unoptimized solution of size $(length(sol))"

                # optimize wrt coverage
                @info "optimizing coverage..."
                length(sats) < length(opts) &&
                optimize_solution(opts) do
                    # clauses: disallow non-improvements
                    for p in sats
                        PicoSAT.add(ps, var[p])
                        PicoSAT.add(ps, 0)
                    end

                    # clause: require some improvement
                    for p in opts
                        p ∉ sats && PicoSAT.add(ps, var[p])
                    end
                    PicoSAT.add(ps, 0)
                end

                # optimize wrt quality
                @info "optimizing quality..."
                optimize_solution(opts) do
                    # clauses: disallow non-improvements
                    for p in sats
                        v_p = var[p]
                        # allow as-good-or-better versions
                        for i = 1:sol[p]
                            PicoSAT.add(ps, v_p + i)
                        end
                        PicoSAT.add(ps, 0)
                    end

                    # clause: require some improvement
                    for p in sats
                        v_p = var[p]
                        # any strictly better versions
                        for i = 1:sol[p]-1
                            PicoSAT.add(ps, v_p + i)
                        end
                    end
                    PicoSAT.add(ps, 0)
                end

                # find next optimization set
                opts′ = empty(opts)
                for p in sats
                    i = sol[p]
                    info_p = info[p]
                    for (j, q) in enumerate(info_p.depends)
                        info_p.conflicts[i, j] || continue
                        q ∈ rest && push!(opts′, q)
                    end
                end

                if isempty(opts′) # nothing left to optimize
                    # delete unreachable packages
                    for p in rest
                        delete!(sol, p)
                    end

                    # save optimal solution
                    @info "optimized solution of size $(length(sol))"
                    push!(sols, copy(sol))

                else # recursion required
                    PicoSAT.push(ps)

                    # clauses: fix already optimized versions
                    for p in sats
                        PicoSAT.add(ps, var[p] + sol[p])
                        PicoSAT.add(ps, 0)
                    end

                    # recursive search call
                    rest′ = setdiff(rest, opts′)
                    find_optimal_solutions(opts′, rest′)

                    # restore satisfied dependencies
                    intersect!(copy!(sats, opts), keys(sol))

                    # pop version fixing clauses
                    PicoSAT.pop(ps)
                end

                # clause: require some improvement
                #   unlike the optimization process above, this allows other
                #   aspects to get worse; it only forces non-domination
                if length(sats) < length(opts)
                    # incomplete solution--improve coverage
                    for p in opts
                        p ∉ sats && PicoSAT.add(ps, var[p])
                    end
                    PicoSAT.add(ps, 0)
                else # complete solution--
                    # preserve coverage
                    for p in opts
                        PicoSAT.add(ps, var[p])
                        PicoSAT.add(ps, 0)
                    end
                    # improve quality
                    for p in opts
                        v_p = var[p]
                        # any strictly better versions
                        for i = 1:sol[p]-1
                            PicoSAT.add(ps, v_p + i)
                        end
                    end
                    PicoSAT.add(ps, 0)
                end
            end
        end

        # start recursion
        opts = Set{P}(reqs)
        rest = setdiff!(Set{P}(names), opts)
        find_optimal_solutions(opts, rest)

        sols # value of the try block
    finally
        # destroy picosat solver
        PicoSAT.reset(ps)
    end

    # sort packages
    if !isempty(sols)
        pkgs = sort!(collect(mapreduce(keys, union, sols, init=reqs)))
        sort!(pkgs, by = !in(reqs)) # required ones first
    else
        pkgs = reqs
    end

    # sort solutions
    # for p in reverse(pkgs)
    #     sort!(sols, by = sol -> get(sol, p, 0))
    # end

    # versions matrix
    vers = Union{V,Nothing}[
        haskey(sol, p) ? info[p].versions[sol[p]] : nothing
        for p in pkgs, sol in sols
    ]

    # return solutions
    return pkgs::Vector{P}, vers::Matrix{Union{V,Nothing}}
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
