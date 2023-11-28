using PicoSAT_jll

const PICOSAT_UNKNOWN = 0
const PICOSAT_SATISFIABLE = 10
const PICOSAT_UNSATISFIABLE = 20

picosat_init() =
    ccall((:picosat_init, libpicosat), Ptr{Cvoid}, ())
picosat_reset(p::Ptr{Cvoid}) =
    ccall((:picosat_reset, libpicosat), Cvoid, (Ptr{Cvoid},), p)
picosat_adjust(p::Ptr{Cvoid}, N::Integer) =
    ccall((:picosat_adjust, libpicosat), Cvoid, (Ptr{Cvoid}, Cint), p, N)
picosat_add(p::Ptr{Cvoid}, lit::Integer) =
    ccall((:picosat_add, libpicosat), Cint, (Ptr{Cvoid}, Cint), p, lit)
picosat_sat(p::Ptr{Cvoid}, limit::Integer = -1) =
    ccall((:picosat_sat, libpicosat), Cint, (Ptr{Cvoid}, Cint), p, limit)
picosat_deref(p::Ptr{Cvoid}, lit::Integer) =
    ccall((:picosat_deref, libpicosat), Cint, (Ptr{Cvoid}, Cint), p, lit)
picosat_push(p::Ptr{Cvoid}) =
    ccall((:picosat_push, libpicosat), Cint, (Ptr{Cvoid},), p)
picosat_pop(p::Ptr{Cvoid}) =
    ccall((:picosat_pop, libpicosat), Cint, (Ptr{Cvoid},), p)

function picosat_print(p::Ptr{Cvoid}, path::AbstractString)
    f = ccall(:fopen, Ptr{Cvoid}, (Cstring, Cstring), path, "w")
    ccall((:picosat_print, libpicosat), Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}), p, f)
    @assert ccall(:fclose, Cint, (Ptr{Cvoid},), f) == 0
end

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
    ps = picosat_init()
    sols = try
        picosat_adjust(ps, N)

        # generate SAT problem
        for p in names
            info_p = info[p]
            n_p = length(info_p.versions)
            v_p = var[p]

            # package implies some version
            #   p => OR_i p@i
            picosat_add(ps, -v_p)
            for i = 1:n_p
                picosat_add(ps, v_p + i)
            end
            picosat_add(ps, 0)

            # version implies its package
            #   p@i => p
            for i = 1:n_p
                picosat_add(ps, -(v_p + i))
                picosat_add(ps, v_p)
                picosat_add(ps, 0)
            end

            # versions are mutually exclusive
            #   !p@i OR !p@j
            for i = 1:n_p-1, j = i+1:n_p
                picosat_add(ps, -(v_p + i))
                picosat_add(ps, -(v_p + j))
                picosat_add(ps, 0)
            end

            # dependencies
            #   ∀ q ∈ deps(p, i): p@i => q
            for i = 1:n_p
                for (j, q) in enumerate(info_p.depends)
                    info_p.conflicts[i, j] || continue
                    picosat_add(ps, -(v_p + i))
                    picosat_add(ps, var[q])
                    picosat_add(ps, 0)
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
                        picosat_add(ps, -(v_p + i))
                        picosat_add(ps, -(v_q + j))
                        picosat_add(ps, 0)
                    end
                end
            end
        end

        # add requirements clauses
        #   ∀ p in reqs: p
        for p in reqs
            picosat_add(ps, var[p])
            picosat_add(ps, 0)
        end

        # helper for finding optimal solutions
        function extract_solution!(sol::Dict{P,Int})
            empty!(sol)
            for (p, v_p) in var
                picosat_deref(ps, v_p) < 0 && continue
                for i = 1:length(info[p].versions)
                    if picosat_deref(ps, v_p + i) > 0
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
                sat = picosat_sat(ps)
                sat == PICOSAT_SATISFIABLE || break
                extract_solution!(sol)

                # optimize solution wrt opts
                while true
                    picosat_push(ps)

                    # clauses: disallow non-improvements
                    for p in opts
                        v_p, i = var[p], sol[p]
                        # allow as-good-or-better versions
                        for j = 1:i
                            picosat_add(ps, v_p + j)
                        end
                        picosat_add(ps, 0)
                    end

                    # clause: require some improvement
                    for p in opts
                        v_p, i = var[p], sol[p]
                        # any strictly better versions
                        for j = 1:i-1
                            picosat_add(ps, v_p + j)
                        end
                    end
                    picosat_add(ps, 0)

                    # check satisfiability
                    sat′ = picosat_sat(ps)
                    sat′ == PICOSAT_SATISFIABLE && extract_solution!(sol)

                    # pop temporary clauses
                    picosat_pop(ps)

                    sat′ == PICOSAT_SATISFIABLE || break
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
                    picosat_push(ps)

                    # clauses: fix optimized versions
                    for p in opts
                        picosat_add(ps, var[p] + sol[p])
                        picosat_add(ps, 0)
                    end

                    # recursive search call
                    find_solutions(opts′, setdiff(rest, opts′))

                    # pop version fixing clauses
                    picosat_pop(ps)
                end

                # clause: require some improvement
                for p in opts
                    v_p, i = var[p], sol[p]
                    # any strictly better versions
                    for j = 1:i-1
                        picosat_add(ps, v_p + j)
                    end
                end
                picosat_add(ps, 0)
            end
        end

        # start recursion
        opts = Set{P}(reqs)
        rest = setdiff!(Set{P}(names), opts)
        find_solutions(opts, rest)

        sols # value of the try block
    finally
        # destroy picosat solver
        picosat_reset(ps)
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
using PrettyTables
pkgs, vers = resolve(info, reqs)
let max_vers = vec(mapslices(r->maximum(v->something(v, v"0-"), r), vers, dims=2))
    pretty_table([pkgs vers],
        highlighters = Highlighter(
            (data, i, j) -> j > 1 && something(data[i, j], max_vers[i]) < max_vers[i],
            foreground = :red,
        )
    )
end
=#

#=
(str -> begin
    reqs = String.(split(str, ','))
    info = load_pkg_info(deps, reqs)
    find_reachable!(info, reqs)
    find_redundant!(info)
    shrink_pkg_info!(info)
    pkgs, vers = resolve(info, reqs)
    max_vers = vec(mapslices(r->maximum(v->something(v, v"0-"), r), vers, dims=2))
    pretty_table([pkgs vers],
        highlighters = Highlighter(
            (data, i, j) -> j > 1 && something(data[i, j], max_vers[i]) < max_vers[i],
            foreground = :red,
        )
    )
end)("CEnum,CUDAdrv,HELICS,DocStringExtensions,FunctionalTables")
=#
