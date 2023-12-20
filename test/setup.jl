using Resolver
using Random
using Test

includet("tiny_data.jl")

function test_resolver(
    deps :: DepsProvider{P},
    reqs :: AbstractVector{P},
) where {P}
    data = Resolver.pkg_data(deps, reqs)
    test_resolver(data, reqs)
end

function test_resolver(
    data :: AbstractDict{P,<:PkgData{P,V}},
    reqs :: AbstractVector{P},
) where {P,V}
    # call resolve
    pkgs, vers = resolve(data, reqs)

    # number of packages & solutions
    M, N = size(vers)
    # @info "$M packages, $N solutions"

    # check basic structure
    @test reqs ⊆ pkgs
    @test allunique(pkgs)
    @test M == length(pkgs)
    @test all(haskey(data, p) for p in pkgs)
    @test all(
        isnothing(vers[i, k]) ||
        vers[i, k] ∈ data[p].versions
        for (i, p) in enumerate(pkgs)
        for k = 1:N
    )

    # helper that captures data, pkgs & vers
    is_valid_sol(vers) = is_valid_solution(data, pkgs, vers)

    # check validity of returned solutions
    for s in eachcol(vers)
        @test is_valid_sol(s)
    end

    # check that solutions can't be trivially improved
    for s in eachcol(vers)
        each_trivial_improvement(data, pkgs, s) do t
            @test !is_valid_sol(t)
        end
    end

    # generate partial order predicate for solutions
    ≤ₛ = make_solution_partial_order!(data, reqs, pkgs)
    # This modifies pkgs so that:
    # - pkgs contains all packages that could be needed by any
    #   version of any package in the original pkgs
    # - allows any possible solution to be expressed, not just
    #   the solutions originally presented

    # check that no solution is dominated by any other
    for s in eachcol(vers), t in eachcol(vers)
        s !== t && @test !(s ≤ₛ t)
    end

    # estimate how many potential solutions there would be
    Π = prod(init=1.0, float(length(data[p].versions)+1) for p in pkgs)
    Π⁺ = 1e6
    if Π ≤ Π⁺
        # @info "optimality testing full data"
        info = data # type unstable but 🤷
    else
        info = Resolver.pkg_info(data, reqs)
        Π = prod(float(length(ip.versions)+1) for ip in values(info))
        if Π > Π⁺
            # @info "no optimality testing"
            return
        end
        # @info "optimality testing filtered info"
    end

    # generate all Π potential solutions
    if isempty(vers)
        # technically, when no solutions, this is dominant
        vers = typeof(vers)(fill(nothing, M, 1))
    end
    each_potential_solution(info, pkgs) do s
        is_valid_sol(s) || return
        # each valid solution is dominated by some returned solution:
        @test any(s ≤ₛ t for t in eachcol(vers))
    end

    return pkgs, vers
end

# checking validity of a solution
function is_valid_solution(
    data :: AbstractDict{P,<:PkgData{P,V}},
    pkgs :: AbstractVector{P},
    vers :: AbstractVector{<:Union{V,Nothing}},
) where {P,V}
    # check satisfiaction of dependencies
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
    # verion dependencies (may expand pkgs)
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

# helpers

function fill_data!(m, n, data, deps, comp)
    for i = 1:m
        data[i] = PkgData(TinyRange(n), deps[i], comp[i])
    end
end
