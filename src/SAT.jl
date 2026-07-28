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
        # 1 variable for package
        # n_p variables for versions
        N += 1 + n_p
    end
    N -= 1 # last used variable

    # instantiate picosat solver
    pico = PicoSAT.init() # TODO: use jl_malloc?
    try # free memory on error
        PicoSAT.adjust(pico, N)
        # default unconstrained variables to false: models then carry
        # fewer spuriously-true packages ("junk"), which makes solves
        # faster and improvement steps land on better versions
        PicoSAT.set_global_default_phase(pico, 0)

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

            # conflicts are added below, rectangle-encoded
        end

        # conflicts, encoded as rectangles: within each interacting pair's
        # conflict bitmap, group p's versions by their conflict pattern
        # against q (identical patterns are the norm, since conflicts come
        # from shared compat entries) and split each pattern into maximal
        # runs of q's versions. each (version group) × (run) rectangle
        # becomes one clause ¬x ∨ ¬y over trigger literals meaning "some
        # version in the group/run is chosen": the version literal itself
        # for singletons, the package variable when every version is
        # covered, and otherwise a shared auxiliary variable defined by
        # one implication per covered version. triggers occur only
        # negatively, so the one-directional definitions suffice: a real
        # conflict forces both triggers true and violates the rectangle
        # clause, while every conflict-free assignment extends by setting
        # each trigger to the value of its defining disjunction.
        psets = Dict{Tuple{Int,Vector{Int}},Int}() # (v_p, versions) => var
        qints = Dict{NTuple{3,Int},Int}()          # (v_q, lo, hi)   => var
        pat = UInt64[]
        S = Int[]
        runs = Tuple{Int,Int}[]

        # trigger literal for "some of versions lo:hi of the package with
        # variable base v and n versions is chosen"
        function run_trigger(v::Int, lo::Int, hi::Int, n::Int)
            lo == hi && return v + lo
            lo == 1 && hi == n && return v
            get!(qints, (v, lo, hi)) do
                t = PicoSAT.inc_max_var(pico)
                for j = lo:hi
                    PicoSAT.add(pico, -(v + j))
                    PicoSAT.add(pico, t)
                    PicoSAT.add(pico, 0)
                end
                t
            end
        end

        # trigger literal for "some version in S of the package with
        # variable base v and n versions is chosen"
        function set_trigger(v::Int, S::Vector{Int}, n::Int)
            length(S) == 1 && return v + S[1]
            length(S) == n && return v
            # contiguous sets share the interval triggers
            S[end] - S[1] + 1 == length(S) &&
                return run_trigger(v, S[1], S[end], n)
            get!(() -> begin
                t = PicoSAT.inc_max_var(pico)
                for i in S
                    PicoSAT.add(pico, -(v + i))
                    PicoSAT.add(pico, t)
                    PicoSAT.add(pico, 0)
                end
                t
            end, psets, (v, copy(S)))
        end

        for p in names
            info_p = info[p]
            n_p = length(info_p.versions)
            v_p = vars[p]
            for (q, b) in info_p.interacts
                p < q || continue # conflicts are symmetrical
                info_q = info[q]
                n_q = length(info_q.versions)
                v_q = vars[q]
                Y = info_q.conflicts
                c = info_q.interacts[p]
                W = col_words(Y)
                resize!(pat, W)
                # group p's versions by conflict pattern: the pattern of
                # p@i is the contiguous column c+i of q's matrix
                groups = Dict{Vector{UInt64},Vector{Int}}()
                for i = 1:n_p
                    col_copy!(pat, Y, c + i)
                    clear_rows_above!(pat, n_q)
                    all(iszero, pat) && continue
                    push!(get!(() -> Int[], groups, copy(pat)), i)
                end
                for (pattern, Sg) in groups
                    # maximal runs of the pattern
                    empty!(runs)
                    prev = -2
                    @inbounds for w = 1:W
                        z = pattern[w]
                        while !iszero(z)
                            j = ((w - 1) << 6) + trailing_zeros(z) + 1
                            z &= z - 1
                            if j == prev + 1
                                runs[end] = (runs[end][1], j)
                            else
                                push!(runs, (j, j))
                            end
                            prev = j
                        end
                    end
                    x = set_trigger(v_p, Sg, n_p)
                    for (lo, hi) in runs
                        y = run_trigger(v_q, lo, hi, n_q)
                        PicoSAT.add(pico, -x)
                        PicoSAT.add(pico, -y)
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
sat_add(sat::SAT{P}, p::P, i::Integer=0; not::Bool=false) where {P} =
    sat_add_var(sat, (-1)^not*(sat.vars[p] + i))

sat_assume(sat::SAT{P}, p::P, i::Integer=0; not::Bool=false) where {P} =
    sat_assume_var(sat, (-1)^not*(sat.vars[p] + i))
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

# iterate the current solution's package => version-index assignments by
# dereferencing the solver's assignment without re-solving: only valid right
# after a solve that returned satisfiable, with no clauses added since --
# ensuring that is the caller's responsibility
function each_solution_index(f::Function, sat::SAT)
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

# extract the solver's current solution; same contract as
# each_solution_index: the last solve must have returned satisfiable,
# with no clauses added since
function extract_solution!(sat::SAT{P}, sol::Dict{P,Int}) where {P}
    empty!(sol)
    each_solution_index(sat) do p, i
        sol[p] = i
    end
    return sol
end

function solution(sat::SAT{P,V}) where {P,V}
    sol = Dict{P,V}()
    is_satisfiable(sat) || return sol
    each_solution_index(sat) do p, i
        sol[p] = sat.info[p].versions[i]
    end
    return sol
end

sat_push(sat::SAT) = PicoSAT.push(sat.pico)
sat_pop(sat::SAT) = PicoSAT.pop(sat.pico)

function with_temp_clauses(body::Function, sat::SAT)
    sat_push(sat)
    try body()
    finally
        sat_pop(sat)
    end
end

# improve package p to its best feasible version: repeatedly demand some
# strictly better version until that becomes unsatisfiable, keeping the last
# good model in sol
function optimize_version!(
    sat :: SAT{P},
    sol :: Dict{P,Int},
    p   :: P,
) where {P}
    # cheap first try: after filtering, the best version is feasible more
    # often than not, and probing it by assumption needs no temp clauses
    # at all — one solve replaces the whole improvement loop
    if sol[p] > 1
        sat_assume(sat, p, 1)
        if is_satisfiable(sat)
            extract_solution!(sat, sol)
            @assert sol[p] == 1
        end
    end
    # sol[p] == 1 is already optimal: the improvement clause below would be
    # empty, forcing a guaranteed-UNSAT solve
    while sol[p] > 1
        improved = with_temp_clauses(sat) do
            # some strictly better version of p
            for i = 1:sol[p]-1
                sat_add(sat, p, i)
            end
            sat_add(sat)
            is_satisfiable(sat) || return false
            extract_solution!(sat, sol)
            return true
        end
        improved || break
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
