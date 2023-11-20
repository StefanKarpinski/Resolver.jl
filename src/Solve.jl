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

function resolve(
    info :: Dict{P, PkgInfo{P,V}},
    reqs :: SetOrVec{P},
) where {P,V}
    # sort stuff for determinism
    names = sort!(collect(keys(info)))
    reqs = sort!(collect(reqs))

    # variable indices:
    #   var[p]      package p chosen
    #   var[p]+i    version p@i chosen
    N = 0
    var = Dict{P, Int}()
    for p in names
        var[p] = N += 1
        N += length(info[p].versions)
    end

    # instantiate picosat solver
    ps = picosat_init()
    picosat_adjust(ps, N)

    # track reverse dependencies
    rdeps = Dict{P, Vector{Pair{P,Int}}}()

    # generate SAT problem
    for p in names
        info_p = info[p]
        v_p = var[p]
        n_p = length(info_p.versions)

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

        # conflicts
        #   !p@i OR !q@j
        for q in sort!(collect(keys(info_p.interacts)))
            # symmetric, only emit once
            q ≥ p && break
            b = info_p.interacts[q]
            v_q = var[q]
            n_q = length(info[q].versions)
            for i = 1:n_p, j = 1:n_q
                info_p.conflicts[i, b+j] || continue
                picosat_add(ps, -(v_p + i))
                picosat_add(ps, -(v_q + j))
                picosat_add(ps, 0)
            end
        end

        # optimality -- one of these has to hold:
        #   1. package isn't chosen
        #   2. version or better chosen
        #   3. conflicting version chosen
        # that is:
        #   !p OR p@1 ... p@i OR conflicts(p, i)
        for i = 1:n_p
            picosat_add(ps, -v_p)
            for j = 1:i
                picosat_add(ps, v_p + j)
            end
            for q in sort!(collect(keys(info_p.interacts)))
                b = info_p.interacts[q]
                info_q = info[q]
                v_q = var[q]
                for j = 1:length(info_q.versions)
                    info_p.conflicts[i, b+j] || continue
                    picosat_add(ps, v_q + j)
                end
            end
            picosat_add(ps, 0)
        end
    end

    # minimality -- only choose necessary packages
    #   p => revdeps(p)
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

    # find all optimal solutions
    sols = Vector{Dict{P,V}}()
    while true
        sat = picosat_sat(ps)
        sat == PICOSAT_SATISFIABLE || break

        # extract one solution
        sol = Dict{P,V}()
        blk = Int[]
        for p in names
            v_p = var[p]
            picosat_deref(ps, v_p) < 0 && continue
            for (i, ver) in enumerate(info[p].versions)
                if picosat_deref(ps, v_p + i) > 0
                    push!(blk, v_p + i)
                    sol[p] = ver
                    break
                end
            end
        end
        push!(sols, sol)

        # block it for next time
        for v in blk
            picosat_add(ps, -v)
        end
        picosat_add(ps, 0)
    end

    # destroy picosat solver
    picosat_reset(ps)

    # filter out solutions that strictly contain other solutions
    filter!(sols) do sol
        !any(sol ⊇ sol′ for sol′ in sols if sol′ ≠ sol)
    end

    # sort solutions lexicographically
    pkgs = sort!(collect(mapreduce(keys, union, sols)))
    sort!(pkgs, by = !in(reqs)) # required ones first
    for p in reverse(pkgs)
        sort!(sols, by = sol -> get(sol, p, v"0-"), rev=true)
    end

    # versions as a matrix
    vers = Union{V, Nothing}[
        get(sol, p, nothing) for p in pkgs, sol in sols
    ]

    # return solution
    return pkgs, vers
end

#=
using DataFrames
x = vec(mapslices(!allequal, vers, dims=2))
df = DataFrame(pkgs = pkgs[x])
for (i, col) in enumerate(eachslice(vers[x, :], dims=2))
    df[:, "$i"] = col
end
=#
