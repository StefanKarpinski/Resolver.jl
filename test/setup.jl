using Resolver
using Random
using Test

@isdefined(includet) ? includet("tiny_data.jl") : include("tiny_data.jl")
@isdefined(includet) ? includet("registry.jl")  : include("registry.jl")

module TestResolver

using Resolver: resolve, DepsProvider, PkgData, PkgInfo, pkg_data, pkg_info
using Test

export test_resolver

function test_resolver(
    deps :: DepsProvider{P},
    reqs :: AbstractVector{P},
) where {P}
    data = pkg_data(deps, reqs)
    test_resolver(data, reqs)
end

const Π⁺ = 1e6 # enumeration budget for optimality testing

function test_resolver(
    data :: AbstractDict{P,<:PkgData{P,V}},
    reqs :: AbstractVector{P},
) where {P,V}
    # resolve: a single optimal solution, or nothing when unsatisfiable
    sol = resolve(data, reqs)

    if sol === nothing
        # verify that no valid solution covers all the requirements
        pkgs = sort!(collect(keys(data)))
        Π = prod(init=1.0, float(length(data[p].versions)+1) for p in pkgs)
        if Π ≤ Π⁺
            ridx = indexin(collect(reqs), pkgs)
            each_potential_solution(data, pkgs) do s
                is_valid_solution(data, pkgs, s) || return
                @test any(isnothing(s[i]) for i in ridx)
            end
        end
        return nothing
    end

    # the solution covers all requirements with known packages & versions
    @test reqs ⊆ keys(sol)
    @test all(haskey(data, p) && sol[p] ∈ data[p].versions for p in keys(sol))

    # the solution is exactly the dependency closure of the requirements
    reach = Set{P}(reqs)
    queue = collect(reach)
    while !isempty(queue)
        p = pop!(queue)
        haskey(data[p].depends, sol[p]) || continue
        for q in data[p].depends[sol[p]]
            q in reach && continue
            push!(reach, q)
            push!(queue, q)
        end
    end
    @test Set(keys(sol)) == reach

    # vector form for the helpers below
    pkgs = sort!(collect(keys(sol)))
    vers = Union{V,Nothing}[sol[p] for p in pkgs]

    # helper that captures data & pkgs
    is_valid_sol(v) = is_valid_solution(data, pkgs, v)

    # the returned solution is valid
    @test is_valid_sol(vers)

    # it can't be trivially improved: raising any package to a better version
    # breaks validity
    each_trivial_improvement(data, pkgs, vers) do t
        @test !is_valid_sol(t)
    end

    # generate partial order predicate for solutions (may expand pkgs; vers
    # stays as-is and is read via get(...), missing entries default per rank)
    ≤ₛ = make_solution_partial_order!(data, reqs, pkgs)

    # `sol` must be Pareto-maximal: no valid solution strictly dominates it.
    # With the deterministic priority-order descent this pins the unique optimum.
    strictly_dominates(s, w) = (w ≤ₛ s) && !(s ≤ₛ w)

    # estimate how many potential solutions there would be
    Π = prod(init=1.0, float(length(data[p].versions)+1) for p in pkgs)
    if Π ≤ Π⁺
        # @info "optimality testing full data"
        info = data # type unstable but 🤷
    else
        info = pkg_info(data, reqs)
        Π = prod(float(length(ip.versions)+1) for ip in values(info))
        if Π > Π⁺
            # @info "no optimality testing"
            return sol
        end
        # @info "optimality testing filtered info"
    end

    # no valid potential solution strictly dominates the returned one
    each_potential_solution(info, pkgs) do s
        is_valid_sol(s) || return
        @test !strictly_dominates(s, vers)
    end

    return sol
end

# checking validity of a solution
function is_valid_solution(
    data :: AbstractDict{P,<:PkgData{P,V}},
    pkgs :: AbstractVector{P},
    vers :: AbstractVector{<:Union{V,Nothing}},
) where {P,V}
    # check satisfaction of dependencies
    for (i, v) in enumerate(vers)
        v === nothing && continue
        data_p = data[pkgs[i]]
        haskey(data_p.depends, v) || continue
        for q in data_p.depends[v]
            j = something(findfirst(==(q), pkgs), 0)
            w = get(vers, j, nothing)
            isnothing(w) && return false
        end
    end
    # check compatibility of versions
    M = size(vers, 1)
    for i = 1:M, j = 1:M
        v = vers[i]; isnothing(v) && continue
        w = vers[j]; isnothing(w) && continue
        p = pkgs[i]
        q = pkgs[j]
        haskey(data[p].compat, v) &&
        haskey(data[p].compat[v], q) &&
        w ∉ data[p].compat[v][q] && return false
        haskey(data[q].compat, w) &&
        haskey(data[q].compat[w], p) &&
        v ∉ data[q].compat[w][p] && return false
    end
    return true
end

# generate ≤ₛ partial order predicate on solutions
function make_solution_partial_order!(
    data :: AbstractDict{P,<:PkgData{P,V}},
    reqs :: AbstractVector{P},
    pkgs :: AbstractVector{P}, # may be expanded
) where {P,V}
    # version dependencies (may expand pkgs)
    deps = Dict{Tuple{Int,V},Vector{Int}}()
    i = 0
    while (i += 1) ≤ length(pkgs)
        p = pkgs[i]
        for v in data[p].versions
            deps[i,v] = Int[]
            haskey(data[p].depends, v) || continue
            for q in data[p].depends[v]
                j = findfirst(==(q), pkgs)
                if isnothing(j)
                    push!(pkgs, q)
                    j = length(pkgs)
                end
                push!(deps[i,v], j)
            end
        end
    end
    # After this:
    # - pkgs contains all packages that could be needed by any
    #   version of any package in the original pkgs
    # - deps has dependencies for all packages and versions
    # - allows any possible solution to be expressed, not just
    #   the solutions originally presented

    # ranking versions (higher = better)
    ranks = Dict{Tuple{V,Int},Int}()
    for (i, p) in enumerate(pkgs),
        (r, v) in enumerate(data[p].versions)
        ranks[v,i] = -r
    end
    # i: package index
    # d: default rank
    rank(v::V, i::Int, d::Int) = ranks[v,i]
    rank(v::Nothing, i::Int, d::Int) = d # default
    rank(s::AbstractVector{Union{V,Nothing}}, i::Int, d::Int) =
        rank(get(s, i, nothing), i, d)

    # partial order on solutions
    function ≤ₛ(
        s::AbstractVector{Union{V,Nothing}},
        t::AbstractVector{Union{V,Nothing}},
    )
        # necessary package indices
        done = 0 # already compared
        need = indexin(reqs, pkgs)
        # check necessary packages
        while done < length(need)
            # first compare satisfaction of needs
            strict = false
            for k = done+1:length(need)
                i = need[k]
                sᵢ = !isnothing(get(s, i, nothing))
                tᵢ = !isnothing(get(t, i, nothing))
                sᵢ > tᵢ && return false
                sᵢ < tᵢ && (strict = true)
            end
            strict && return true
            # then compare version quality
            for k = done+1:length(need)
                i = need[k]
                # no version = worst = typemin
                sᵢ = rank(s, i, typemin(Int))
                tᵢ = rank(t, i, typemin(Int))
                @assert (sᵢ == typemin(Int)) == (tᵢ == typemin(Int))
                sᵢ > tᵢ && return false
                sᵢ < tᵢ && (strict = true)
            end
            strict && return true
            # find newly necessary packages
            for k = done+1:length(need)
                i = need[k]
                v = get(s, i, nothing)
                @assert v == get(t, i, nothing)
                isnothing(v) || for j in deps[i,v]
                    j ∉ need && push!(need, j)
                end
                done = k
            end
        end
        # check unnecessary packages
        for i = 1:length(pkgs)
            i in need && continue
            # no version = best = typemax
            sᵢ = rank(s, i, typemax(Int))
            tᵢ = rank(t, i, typemax(Int))
            sᵢ > tᵢ && return false
        end
        return true
    end
end

function each_trivial_improvement(
    body :: Function, # callback
    data :: AbstractDict{P,<:PkgData{P,V}},
    pkgs :: AbstractVector{P},
    vers :: AbstractVector{Union{V,Nothing}},
) where {P,V}
    s = Vector{Union{V,Nothing}}(vers) # copy
    for (i, v) in enumerate(vers)
        vers_p = data[pkgs[i]].versions
        r = something(findfirst(==(v), vers_p), 0)
        for (r′, v′) in enumerate(vers_p)
            r′ < r || break
            s[i] = v′
            body(s)
        end
        s[i] = v # return to original value
    end
end

PkgVers{P,V} = Union{PkgData{P,V},PkgInfo{P,V}}

function each_potential_solution(
    body :: Function, # callback
    data :: AbstractDict{P,<:PkgVers{P,V}},
    pkgs :: AbstractVector{P},
) where {P,V}
    L = length(pkgs)
    s = Vector{Union{V,Nothing}}(nothing, L)
    function gen_solutions!(i::Int = 1)
        # call body if solution is complete:
        i ≤ L || (body(s); return)
        # otherwise iterate versions of next package:
        p = pkgs[i]
        vers_p = haskey(data, p) ? data[p].versions : V[]
        for r = 0:length(vers_p)
            v = get(vers_p, r, nothing)
            s[i] = v
            gen_solutions!(i+1)
        end
    end
    gen_solutions!()
end

end # module

using .TestResolver
