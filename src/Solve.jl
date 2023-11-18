using PicoSAT_jll

const PICOSAT_UNKNOWN = 0
const PICOSAT_SATISFIABLE = 10
const PICOSAT_UNSATISFIABLE = 20

# TODO: search and replace this when done
const PicoPtr = Ptr{Cvoid}

picosat_init() =
    ccall((:picosat_init, libpicosat), PicoPtr, ())
picosat_reset(p::PicoPtr) =
    ccall((:picosat_reset, libpicosat), Cvoid, (PicoPtr,), p)
picosat_adjust(p::PicoPtr, n::Integer) =
    ccall((:picosat_adjust, libpicosat), Cvoid, (PicoPtr, Cint), p, n)
picosat_add(p::PicoPtr, lit::Integer) =
    ccall((:picosat_add, libpicosat), Cint, (PicoPtr, Cint), p, lit)
picosat_sat(p::PicoPtr, limit::Integer = -1) =
    ccall((:picosat_sat, libpicosat), Cint, (PicoPtr, Cint), p, limit)
picosat_deref(p::PicoPtr, lit::Integer) =
    ccall((:picosat_deref, libpicosat), Cint, (PicoPtr, Cint), p, lit)

function resolve(
    info :: Dict{P, PkgInfo{P,V}},
    reqs :: SetOrVec{P},
) where {P,V}
    # replace args with sorted versions
    infos = sort!(collect(info), by=first)
    reqs = sort!(collect(reqs))

    # compute variable indices
    n = 0 # number of variables
    var = Dict{P, Int}() # variable indices
    for (p, info_p) in infos
        var[p] = n + 1
        n += 1 + length(info_p.versions)
    end

    # instantiate picosat solver
    ps = picosat_init()
    picosat_adjust(ps, n)

    # add requirements clauses
    for p in reqs
        picosat_add(ps, var[p])
        picosat_add(ps, 0)
    end

    # add package version clauses
    for (p, info_p) in infos
        var_p = var[p]
        picosat_add(ps, -var_p)
        for i = 1:length(info_p.versions)
            picosat_add(ps, var_p + i)
        end
        picosat_add(ps, 0)
    end

    # output dependency clauses
    for (p, info_p) in infos
        var_p = var[p]
        for i = 1:length(info_p.versions)
            for (j, q) in enumerate(info_p.depends)
                info_p.conflicts[i, j] || continue
                picosat_add(ps, -(var_p + i))
                picosat_add(ps, var[q])
                picosat_add(ps, 0)
            end
        end
    end

    # output incompatibility clauses
    for (p, info_p) in infos
        var_p = var[p]
        for q in sort!(collect(keys(info_p.interacts)))
            q < p || break # info recorded twice, only emit once
            b = info_p.interacts[q]
            info_q = info[q]
            var_q = var[q]
            for i = 1:length(info_p.versions),
                j = 1:length(info_q.versions)
                info_p.conflicts[i, b+j] || continue
                picosat_add(ps, -(var_p + i))
                picosat_add(ps, -(var_q + j))
                picosat_add(ps, 0)
            end
        end
    end

    # output optimality clauses
    for (p, info_p) in infos
        var_p = var[p]
        for i = 1:length(info_p.versions)
            picosat_add(ps, -var_p)
            picosat_add(ps, var_p + i)
            # can be "excused" if a better version is chosen
            for j = 1:i-1
                picosat_add(ps, var_p + j)
            end
            # or if it conflicts with something else chosen
            for q in sort!(collect(keys(info_p.interacts)))
                b = info_p.interacts[q]
                info_q = info[q]
                var_q = var[q]
                for j = 1:length(info_q.versions)
                    info_p.conflicts[i, b+j] || continue
                    picosat_add(ps, var_q + j)
                end
            end
            picosat_add(ps, 0)
        end
    end

    # find all optimal solutions
    sols = Vector{Dict{P,V}}()
    while true
        sat = picosat_sat(ps)
        sat == PICOSAT_SATISFIABLE || break

        # extract one solution
        sol = Dict{P,V}()
        blk = Int[]
        for (p, info_p) in infos
            var_p = var[p]
            picosat_deref(ps, var_p) < 0 && continue
            for (i, ver) in enumerate(info_p.versions)
                if picosat_deref(ps, var_p + i) > 0
                    push!(blk, var_p + i)
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

    # return solution
    return sols
end
