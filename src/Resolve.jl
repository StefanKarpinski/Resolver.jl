function resolve(
    info :: Dict{P, PkgInfo{P,V}},
    reqs :: SetOrVec{P},
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

        # add requirements clauses
        #   ∀ p in reqs: p
        for p in reqs
            PicoSAT.add(ps, var[p])
            PicoSAT.add(ps, 0)
        end

        # helper for finding optimal solutions
        function extract_solution!(sol::Dict{P,Int})
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
            return sol
        end

        # find all optimal solutions
        sols = Vector{Dict{P,Int}}()
        sol = Dict{P,Int}()

        function find_solutions(
            opts :: Set{P}, # packages to optimize next
            rest :: Set{P}, # other unoptimized packages
        )
            while true
                # find some solution
                sat = PicoSAT.sat(ps)
                sat == PicoSAT.SATISFIABLE || break
                extract_solution!(sol)

                # optimize solution wrt opts
                while true
                    PicoSAT.push(ps)

                    # clauses: disallow non-improvements
                    for p in opts
                        v_p, i = var[p], sol[p]
                        # allow as-good-or-better versions
                        for j = 1:i
                            PicoSAT.add(ps, v_p + j)
                        end
                        PicoSAT.add(ps, 0)
                    end

                    # clause: require some improvement
                    for p in opts
                        v_p, i = var[p], sol[p]
                        # any strictly better versions
                        for j = 1:i-1
                            PicoSAT.add(ps, v_p + j)
                        end
                    end
                    PicoSAT.add(ps, 0)

                    # check satisfiability
                    sat′ = PicoSAT.sat(ps)
                    sat′ == PicoSAT.SATISFIABLE && extract_solution!(sol)

                    # pop temporary clauses
                    PicoSAT.pop(ps)

                    sat′ == PicoSAT.SATISFIABLE || break
                end

                # find next optimization set
                opts′ = empty(opts)
                for p in opts
                    info_p, i = info[p], sol[p]
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
                    push!(sols, copy(sol))

                else # recursion required
                    PicoSAT.push(ps)

                    # clauses: fix optimized versions
                    for p in opts
                        PicoSAT.add(ps, var[p] + sol[p])
                        PicoSAT.add(ps, 0)
                    end

                    # recursive search call
                    find_solutions(opts′, setdiff(rest, opts′))

                    # pop version fixing clauses
                    PicoSAT.pop(ps)
                end

                # clause: require some improvement
                for p in opts
                    v_p, i = var[p], sol[p]
                    # any strictly better versions
                    for j = 1:i-1
                        PicoSAT.add(ps, v_p + j)
                    end
                end
                PicoSAT.add(ps, 0)
            end
        end

        # start recursion
        opts = Set{P}(reqs)
        rest = setdiff!(Set{P}(names), opts)
        find_solutions(opts, rest)

        sols # value of the try block
    finally
        # destroy picosat solver
        PicoSAT.reset(ps)
    end

    # sort packages
    if !isempty(sols)
        pkgs = sort!(collect(mapreduce(keys, union, sols)))
        sort!(pkgs, by = !in(reqs)) # required ones first
    else
        pkgs = Vector{P}()
    end

    # sort solutions
    for p in reverse(pkgs)
        sort!(sols, by = sol -> get(sol, p, 0))
    end

    # versions as a matrix
    vers = Union{V, Nothing}[
        haskey(sol, p) ? info[p].versions[sol[p]] : nothing
        for p in pkgs, sol in sols
    ]

    # return solution
    return pkgs, vers
end

#=
(str -> begin
    reqs = String.(split(str, ','))
    info = load_pkg_info(deps, reqs)
    find_reachable!(info, reqs)
    find_redundant!(info)
    shrink_pkg_info!(info)
    pkgs, vers = resolve(info, reqs)
    mv = vec(mapslices(r->maximum(v->something(v, v"0-"), r), vers, dims=2))
    pretty_table([pkgs vers],
        highlighters = Highlighter(
            (data, i, j) -> j > 1 && something(data[i, j], mv[i]) < mv[i],
            foreground = :red,
        )
    )
end)("CEnum,CUDAdrv,HELICS,DocStringExtensions,FunctionalTables")
=#
