mutable struct SAT{P,V}
    info :: Dict{P,PkgInfo{P,V}}
    pico :: Ptr{Cvoid}
    vars :: Dict{P,Int}
end

function Base.show(io::IO, sat::SAT)
    show(io, typeof(sat))
    p = length(sat.vars)
    v = PicoSAT.var_count(sat.pico)
    c = PicoSAT.clause_count(sat.pico)
    print(io,
        "(packages: ", p,
        ", versions: ", v-p,
        ", clauses: ", c, ")")
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
    vars = Dict{P,Int}()
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

sat_new_variable(sat::SAT) = PicoSAT.inc_max_var(sat.pico)

sat_add_var(sat::SAT, v::Integer) = PicoSAT.add(sat.pico, v)
sat_assume_var(sat::SAT, v::Integer) = PicoSAT.assume(sat.pico, v)

sat_add(sat::SAT) = sat_add_var(sat, 0)
sat_add(sat::SAT{P}, p::P, i::Integer=0) where {P} =
    sat_add_var(sat, sat.vars[p] + i)

sat_assume(sat::SAT{P}, p::P, i::Integer=0) where {P} =
    sat_assume_var(sat, sat.vars[p] + i)
sat_assume(sat::SAT{P}, v::SetOrVec{P}) where {P} =
    foreach(p -> sat_assume(sat, p), v)
sat_assume(sat::SAT{P}, d::Dict{P,<:Integer}) where {P} =
    for (p, i) in d
        sat_assume(sat, p, i)
    end

is_satisfiable(sat::SAT) =
    PicoSAT.sat(sat.pico) == PicoSAT.SATISFIABLE

function is_satisfiable(sat::SAT{P}, reqs::Union{P,SetOrVec{P}}) where {P}
    sat_assume(sat, reqs)
    is_satisfiable(sat)
end

const is_unsatisfiable = !is_satisfiable

function each_solution_index(f::Function, sat::SAT)
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

function extract_solution!(sat::SAT{P}, sol::Dict{P,Int}) where {P}
    empty!(sol)
    each_solution_index(sat) do p, i
        sol[p] = i
    end
    return sol
end

function solution(sat::SAT{P,V}) where {P,V}
    sol = Dict{P,V}()
    each_solution_index(sat) do p, i
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

function optimize_solution!(
    improve :: Function,
    sat     :: SAT{P},
    sol     :: Dict{P,Int},
) where {P}
    done::Bool = false
    while !done
        with_temp_clauses(sat) do
            improve() # callback adds SAT clauses
            done = is_unsatisfiable(sat)
            done || extract_solution!(sat, sol)
        end
    end
end

# fix a given set of versions

function fix_versions(
    sat  :: SAT{P,V},
    pkgs :: AbstractVector{P},
    vers :: AbstractVector{Union{Nothing,V}},
) where {P,V}
    for (p, v) in zip(pkgs, vers)
        v === nothing && continue
        i = findfirst(==(v), sat.info[p].versions)
        i === nothing && throw(ArgumentError("package $p: unknown version $v"))
        sat_add(sat, p, i)
        sat_add(sat)
    end
end

# extra SAT exploration functions

function sat_mus(
    sat  :: SAT{P},
    reqs :: SetOrVec{P} = keys(sat.info),
) where {P}
    sat_assume(sat, reqs)
    is_satisfiable(sat) && return empty(reqs)
    # find initial unsatisfiable set
    mus = Set{P}()
    for p in reqs
        PicoSAT.failed(sat.pico, sat.vars[p]) ≠ 0 && push!(mus, p)
    end
    @assert is_unsatisfiable(sat, mus)
    # try shrinking it
    @label again
    for p in mus
        delete!(mus, p)
        is_unsatisfiable(sat, mus) && @goto again
        push!(mus, p)
    end
    # can't be shrunk
    return mus
end

function sat_mice(
    sat  :: SAT{P},
    reqs :: SetOrVec{P} = keys(sat.info),
) where {P}
    reqs = Set{P}(reqs)
    mice = Set{P}[]
    while true
        mus = sat_mus(sat, reqs)
        isempty(mus) && break
        println(mus)
        push!(mice, mus)
        setdiff!(reqs, mus)
    end
    return mice
end

function sat_humus(
    sat  :: SAT{P},
    reqs :: SetOrVec{P} = keys(sat.info),
) where {P}
    reqs = Set{P}(reqs)
    humus = Set{P}()
    while true
        mus = sat_mus(sat, reqs)
        isempty(mus) && break
        println(mus)
        union!(humus, mus)
        setdiff!(reqs, mus)
    end
    return sort!(collect(humus))
end
