using Resolver
using Random
using Test

@isdefined(includet) ? includet("tiny_data.jl") : include("tiny_data.jl")
@isdefined(includet) ? includet("registry.jl")  : include("registry.jl")

module TestResolver

using Resolver: resolve, DepsProvider, PkgData, PkgInfo, pkg_data, pkg_info
using Test

export test_resolver, ref_resolve

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

    # it can't be trivially improved: raising any package to a better
    # version breaks validity or tightness (a valid but unjustified-junk
    # variant is not a legitimate dominator under the front semantics)
    each_trivial_improvement(data, pkgs, vers) do t
        if is_valid_sol(t)
            tsol = Dict{P,V}(
                pkgs[i] => t[i] for i = 1:length(pkgs) if t[i] !== nothing)
            seen = Set{P}(q for q in reqs)
            queue = P[q for q in reqs]
            while !isempty(queue)
                q = pop!(queue)
                haskey(data[q].depends, tsol[q]) || continue
                for r in data[q].depends[tsol[q]]
                    r in seen && continue
                    push!(seen, r)
                    push!(queue, r)
                end
            end
            @test Set(keys(tsol)) != seen
        end
    end

    # optimality: the solution is exactly the one the contract specifies,
    # checked against the brute-force reference when enumeration fits budget
    Π = prod(init=1.0, float(length(data[p].versions)+1) for p in keys(data))
    if Π ≤ Π⁺
        @test sol == ref_resolve(data, reqs)
    end

    return sol
end

# brute-force reference for the front-necessity contract, straight from the
# spec (docs/design/front-necessity.md): enumerate all valid closure-tight
# solutions covering the requirements, form the Pareto front under coverage-
# monotone version dominance, and repeatedly pin the highest-priority package
# present in every remaining front member at the best version among them
function ref_resolve(
    data :: AbstractDict{P,<:PkgData{P,V}},
    reqs :: AbstractVector{P};
    by   :: Function = identity,
) where {P,V}
    pkgs = sort!(collect(keys(data)))
    rank = Dict{P,Dict{V,Int}}(
        p => Dict{V,Int}(v => r for (r, v) in enumerate(data[p].versions))
        for p in pkgs)
    # all valid closure-tight solutions covering the requirements, as
    # package => version dicts
    cands = Dict{P,V}[]
    each_potential_solution(data, pkgs) do s
        is_valid_solution(data, pkgs, s) || return
        sol = Dict{P,V}(
            pkgs[i] => s[i] for i = 1:length(pkgs) if s[i] !== nothing)
        all(haskey(sol, q) for q in reqs) || return
        # tight: exactly the dependency closure of the requirements
        seen = Set{P}(q for q in reqs)
        queue = P[q for q in reqs]
        while !isempty(queue)
            q = pop!(queue)
            haskey(data[q].depends, sol[q]) || continue
            for r in data[q].depends[sol[q]]
                r in seen && continue
                push!(seen, r)
                push!(queue, r)
            end
        end
        Set(keys(sol)) == seen || return
        push!(cands, sol)
    end
    isempty(cands) && return nothing
    # coverage-monotone version dominance: t covers everything s has, at
    # versions at least as good, strictly better somewhere
    function dominates(t, s)
        strict = false
        for (q, v) in s
            haskey(t, q) || return false
            rt = rank[q][t[q]]
            rs = rank[q][v]
            rt > rs && return false
            rt < rs && (strict = true)
        end
        return strict
    end
    front = [s for s in cands if !any(dominates(t, s) for t in cands)]
    # forced descent over the front
    sol = Dict{P,V}()
    order = sort(pkgs; by)
    while true
        i = findfirst(p -> !haskey(sol, p) &&
                           all(haskey(f, p) for f in front), order)
        i === nothing && break
        p = order[i]
        best = data[p].versions[minimum(rank[p][f[p]] for f in front)]
        sol[p] = best
        filter!(f -> f[p] == best, front)
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
