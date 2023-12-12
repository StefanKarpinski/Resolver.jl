mutable struct SAT{P,V}
    info :: Dict{P,PkgInfo{P,V}}
    pico :: Ptr{Cvoid}
    vars :: Dict{P,Int}
end

function SAT(
    info :: Dict{P,PkgInfo{P,V}},
) where {P,V}
    # sort names for predictability
    names = sort!(collect(keys(info)))

    # variable indices:
    #   p    vars[p]    package p chosen
    #   p@i  vars[p]+i  version i of p chosen
    N = 1
    vars = Dict{P, Int}()
    for p in names
        vars[p] = N
        n_p = length(info[p].versions)
        # 1 varible for package
        # n_p variables for versions
        N += 1 + n_p
    end
    N -= 1 # last used variable

    # instantiate picosat solver
    pico = PicoSAT.init() # TODO: use jl_malloc?
    try # free memory on error
        PicoSAT.adjust(pico, N)

        # generate SAT problem
        for p in names
            info_p = info[p]
            n_p = length(info_p.versions)
            v_p = vars[p]

            # package implies some version
            #   p => OR_i p@i
            PicoSAT.add(pico, -v_p)
            for i = 1:n_p
                PicoSAT.add(pico, v_p + i)
            end
            PicoSAT.add(pico, 0)

            # version implies its package
            #   p@i => p
            for i = 1:n_p
                PicoSAT.add(pico, -(v_p + i))
                PicoSAT.add(pico, v_p)
                PicoSAT.add(pico, 0)
            end

            # versions are mutually exclusive
            #   !p@i OR !p@j
            for i = 1:n_p-1, j = i+1:n_p
                PicoSAT.add(pico, -(v_p + i))
                PicoSAT.add(pico, -(v_p + j))
                PicoSAT.add(pico, 0)
            end

            # dependencies
            #   ∀ q ∈ deps(p, i): p@i => q
            for i = 1:n_p
                for (j, q) in enumerate(info_p.depends)
                    info_p.conflicts[i, j] || continue
                    PicoSAT.add(pico, -(v_p + i))
                    PicoSAT.add(pico, vars[q])
                    PicoSAT.add(pico, 0)
                end
            end

            # conflicts
            #   ∀ q@j ∈ conflicts(p, i): !p@i OR !q@j
            for i = 1:n_p
                for (q, b) in info_p.interacts
                    p < q || continue # conflicts are symmetrical
                    n_q = length(info[q].versions)
                    v_q = vars[q]
                    for j = 1:n_q
                        info_p.conflicts[i, b+j] || continue
                        # conflicting versions are mutually exclusive
                        PicoSAT.add(pico, -(v_p + i))
                        PicoSAT.add(pico, -(v_q + j))
                        PicoSAT.add(pico, 0)
                    end
                end
            end
        end
    catch # on error free picosat solver
        PicoSAT.reset(pico)
        rethrow()
    end
    finalizer(finalize, SAT(info, pico, vars))
end

function finalize(sat::SAT)
    pico = sat.pico
    pico == C_NULL && return
    sat.pico = C_NULL
    PicoSAT.reset(pico)
end

sat_add(sat::SAT, lit::Integer) = PicoSAT.add(sat.pico, lit)
sat_assume(sat::SAT, lit::Integer) = PicoSAT.assume(sat.pico, lit)

is_satisfiable(sat::SAT) = PicoSAT.sat(sat.pico) == PicoSAT.SATISFIABLE

const is_unsatisfiable = !is_satisfiable

function each_solution_ind(f::Function, sat::SAT)
    is_satisfiable(sat) || return
    for (p, v_p) in sat.vars
        PicoSAT.deref(sat.pico, v_p) < 0 && continue
        i = 1
        while true
            if PicoSAT.deref(sat.pico, v_p + i) > 0
                # guaranteed to happen by SAT construction:
                # v_p => v_p + i for some i = 1:n_p
                f(p, i)
                break
            end
            i += 1
        end
    end
end

function get_solution_inds!(sat::SAT{P}, sol::Dict{P,Int}) where {P}
    empty!(sol)
    each_solution_ind(sat) do p, i
        sol[p] = i
    end
    return sol
end

function get_solution(sat::SAT{P,V}) where {P,V}
    sol = Dict{P,V}()
    each_solution_ind(sat) do p, i
        sol[p] = sat.info[p].versions[i]
    end
    return sol
end

function with_temp_clauses(body::Function, sat::SAT)
    PicoSAT.push(sat.pico)
    try body()
    finally
        PicoSAT.pop(sat.pico)
    end
end

function optimize_solution(
    improve :: Function,
    sat     :: SAT,
    extract :: Function,
)
    done::Bool = false
    while !done
        with_temp_clauses(sat) do
            improve() # callback adds SAT clauses
            done = is_unsatisfiable(sat)
            done || extract()
        end
    end
end
