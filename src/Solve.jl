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
    #   p    var[p]        package p chosen
    #   p@i  var[p]+i      version i of p chosen
    #   p!i  var[p]+n_p+i  version i of p conflicted
    N = 1
    var = Dict{P, Int}()
    for p in names
        var[p] = N
        n_p = length(info[p].versions)
        # 1 varible for package
        # n_p variables for versions
        # n_p variables for conflicts
        N += 2n_p + 1
    end
    N -= 1 # last used variable

    # instantiate picosat solver
    ps = picosat_init()
    sols = try
        picosat_adjust(ps, N)

        # track reverse dependencies
        rdeps = Dict{P, Vector{Pair{P,Int}}}()

        # generate SAT problem
        for p in names
            info_p = info[p]
            v_p = var[p]
            n_p = length(info_p.versions)
            x_p = v_p + n_p

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
            for i = 1:n_p, j = 1:i-1
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
                    # save reverse dependency
                    push!(get!(() -> Pair{P,Int}[], rdeps, q), p => i)
                end
            end

            # conflict indicators
            interacts = sort!(collect(keys(info_p.interacts)))
            for i = 1:n_p
                # collect this version's conflicts
                conflicts = Vector{Pair{P,Int}}()
                for q in interacts
                    b = info_p.interacts[q]
                    n_q = length(info[q].versions)
                    for j = 1:n_q
                        info_p.conflicts[i, b+j] || continue
                        push!(conflicts, q => j)
                    end
                end

                # conflict indicator implies some conflicting version
                #   p!i => conflicts(p, i)
                picosat_add(ps, -(x_p + i))
                for (q, j) in conflicts
                    picosat_add(ps, var[q] + j)
                end
                picosat_add(ps, 0)

                # each conflicting version implies conflict indicator
                #   ∀ q@j ∈ conflicts(p, i): q@j => p!i
                for (q, j) in conflicts
                    picosat_add(ps, -(var[q] + j))
                    picosat_add(ps, x_p + i)
                    picosat_add(ps, 0)
                end
            end

            # conflict indicator blocks actual version
            #   p!i => !p@i
            for i = 1:n_p
                picosat_add(ps, -(x_p + i))
                picosat_add(ps, -(v_p + i))
                picosat_add(ps, 0)
            end

            # optimality (package-wise unilateral improvement)
            #   ∀ j < i: p@i => p!j
            #
            # note: it's still possible to get a suboptimal solution in a case where
            # p@i′ > p@i and p@j′ > p@j but p@i conflicts p@j′ and p@j conflicts
            # p@i′ so you have to change both at once in order to actually improve;
            # this rule can't prevent that but we force actual improvement when
            # searching and then we filter out fully dominated solutions
            #
            for i = 1:n_p, j = 1:i-1
                picosat_add(ps, -(v_p + i))
                picosat_add(ps, x_p + j)
                picosat_add(ps, 0)
            end
        end

        # minimality (only packages that something depends on)
        #   p => revdeps(p)
        #
        # note: this can be foiled by dependency loops, e.g. where p@i depends on
        # q@j and vice versa, but nothing external actually depends on either one;
        # this rule can't't prevent that, but we filter these out later
        #
        for (q, rdeps_q) in sort!(collect(rdeps), by=first)
            picosat_add(ps, -var[q])
            for (p, i) in rdeps_q
                picosat_add(ps, -(var[p] + i))
            end
            picosat_add(ps, 0)
        end

        # add requirements clauses
        #   these are the packages that are required
        for p in reqs
            picosat_add(ps, var[p])
            picosat_add(ps, 0)
        end

        function extract_solution!(sol::Dict{P,Int})
            empty!(sol)
            for p in names
                v_p = var[p]
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

        while true
            sat = picosat_sat(ps)
            sat == PICOSAT_SATISFIABLE || break
            extract_solution!(sol)

            # minimize solution
            while true
                # save main problem state
                picosat_push(ps)

                # disallow adding packages
                for p in setdiff(names, keys(sol))
                    picosat_add(ps, -var[p])
                    picosat_add(ps, 0)
                end

                # disallow changing versions
                for (p, i) in sol
                    v_p = var[p]
                    picosat_add(ps, v_p + i)
                    picosat_add(ps, -v_p)
                    picosat_add(ps, 0)
                end

                # require dropping some packages
                for (p, i) in sol
                    picosat_add(ps, -var[p])
                end
                picosat_add(ps, 0)

                # recheck satisfiability
                sat′ = picosat_sat(ps)
                sat′ == PICOSAT_SATISFIABLE && extract_solution!(sol)

                # pop temporary clauses
                picosat_pop(ps)

                # fully minimized
                sat′ == PICOSAT_SATISFIABLE || break
            end

            # drop earlier solutions that are strictly worse
            filter!(sols) do sol′
                for (p, i) in sol
                    j = get(sol′, p, i)
                    j < i && return true # keep
                end
                return false # reject
            end

            # save minimized solution
            push!(sols, copy(sol))

            # next solution has to improve somehow
            for (p, i) in sol
                v_p = var[p]
                for j = 1:i-1
                    picosat_add(ps, var[p] + j)
                end
            end
            picosat_add(ps, 0)
        end

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
