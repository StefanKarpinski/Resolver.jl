using Test
using Resolver

function test_resolve(
    data :: AbstractDict{P,<:PkgData{P,V}},
    reqs :: AbstractVector{P},
) where {P,V}
    pkgs, vers = resolve(data, reqs)
    test_resolve(data, reqs, pkgs, vers)
end

function test_resolve(
# problem:
    data :: AbstractDict{P,<:PkgData{P,V}},
    reqs :: AbstractVector{P},
# solutions:
    pkgs :: AbstractVector{P},
    vers :: AbstractMatrix{Union{V,Nothing}},
) where {P,V}
    # number of packages & solutions
    M, N = size(vers)

    # check basic structure
    @test reqs ⊆ pkgs
    @test allunique(pkgs)
    @test M == length(pkgs)
    @test all(haskey(data, p) for p in pkgs)
    @test all(
        isnothing(vers[i, s]) ||
        vers[i, s] ∈ data[p].versions
        for (i, p) in enumerate(pkgs)
        for s = 1:N
    )

    # checking validity of a solution
    function is_valid_solution(s::AbstractVector{Union{V,Nothing}})
        # check satisfiaction of dependencies
        for (i, v) in enumerate(s)
            v === nothing && continue
            data_p = data[pkgs[i]]
            haskey(data_p.depends, v) || continue
            for q in data_p.depends[v]
                q ∈ pkgs || return false
                j = findfirst(==(q), pkgs)
                isnothing(s[j]) && return false
            end
        end
        # check compatibility of versions
        for i = 1:M, j = 1:M
            v = s[i]; isnothing(v) && continue
            w = s[j]; isnothing(w) && continue
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

    # check each returned solution
    for k = 1:N
        s = @view vers[:, k]
        @test is_valid_solution(s)
    end

    # verion dependencies (may replace & expand pkgs)
    deps = Dict{Tuple{V,Int},Vector{Int}}()
    i = 0
    while (i += 1) ≤ length(pkgs)
        p = pkgs[i]
        for v in data[p].versions
            deps[v,i] = Int[]
            haskey(data[p].depends, v) || continue
            for q in data[p].depends[v]
                j = findfirst(==(q), pkgs)
                if isnothing(j)
                    if length(pkgs) == M
                        pkgs = copy(pkgs)
                    end
                    push!(pkgs, q)
                    j = length(pkgs)
                end
                push!(deps[v,i], j)
            end
        end
    end
    # pkgs now contains all packages that could be needed by
    # any version of any package in the original pkgs and deps
    # has dependencies for all packages and versions of them

    # ranking versions (higher = better)
    ranks = Dict{Tuple{V,Int},Int}()
    for (i, p) in enumerate(pkgs),
        (r, v) in enumerate(data[p].versions)
        ranks[v,i] = -r
    end
    # i: package index
    # d: default version
    rank(v::V, i::Int, d::V) = ranks[v,i]
    rank(v::Nothing, i::Int, d::V) = d # default
    rank(s::AbstractVector{Union{V,Nothing}}, i::Int, d::V) =
        rank(get(s, i, nothing), i, d)

    # partial order on solutions
    function ≤ₛ(
        s::AbstractVector{Union{V,Nothing}},
        t::AbstractVector{Union{V,Nothing}},
    )
        # set of necessary package indices
        need = Set{Int}(indexin(reqs, pkgs))
        # check necessary packages
        while true
            strict = false
            for i in need
                # no version = worst = typemin(V)
                sᵢ = rank(s, i, typemin(V))
                tᵢ = rank(t, i, typemin(V))
                sᵢ > tᵢ && return false
                sᵢ < tᵢ && (strict = true)
            end
            strict && return true
            n = length(need)
            for i in need
                v = get(s, i, nothing)
                @assert v == get(t, i, nothing)
                union!(need, deps[v,i])
            end
            n == length(need) && break
        end
        # check unnecessary packages
        for i = 1:M
            i in need && continue
            # no version = best = typemax(V)
            sᵢ = rank(s, i, typemax(V))
            tᵢ = rank(t, i, typemax(V))
            sᵢ > tᵢ && return false
        end
        return true
    end
    ≥ₛ(s::Int, t::Int) = ≤ₛ(t, s)

    # check that no solution is dominated by any other
    for k = 1:N, l = 1:N
        k ≠ l || continue
        s = @view vers[:, k]
        t = @view vers[:, l]
        @test !(s ≤ₛ t)
    end
end
