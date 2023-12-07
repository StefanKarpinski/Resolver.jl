using Test
using Random
using Resolver
using Resolver: resolve_core, SetOrVector, DepsProvider, PkgData

function check_resolved(
# problem:
    data :: AbstractDict{P,<:PkgData{P,V}},
    reqs :: AbstractVector{P},
# solutions:
    pkgs :: AbstractVector{P},
    vers :: AbstractMatrix{Union{V,Nothing}},
) where {P,V}
    # number of packages & solutions
    P, S = size(vers)

    # check basic structure
    @test reqs ⊆ pkgs
    @test allunique(pkgs)
    @test P == length(pkgs)
    @test all(haskey(data, p) for p in pkgs)
    @test all(
        isnothing(vers[i, s]) ||
        vers[i, s] ∈ data[p].versions
        for (i, p) in enumerate(pkgs)
        for s = 1:S
    )

    # check validity of each solution
    for s = 1:S
        # check satisfiaction of dependencies
        for i = 1:P
            v = vers[i, s]
            v === nothing && continue
            data_p = data[pkgs[i]]
            haskey(data_p.depends, v) || continue
            for q in data_p.depends[v]
                @test q ∈ pkgs
                j = findfirst(==(q), pkgs)
                @test !isnothing(vers[j, s])
            end
        end
        # check compatibility of versions
        for i = 1:P, j = 1:P
            v = vers[i, s]
            isnothing(v) && continue
            w = vers[j, s]
            isnothing(w) && continue
            p = pkgs[i]
            q = pkgs[j]
            haskey(data[p].compat, v) &&
            haskey(data[p].compat[v], q) &&
            @test w in data[p].compat[v][q]
            haskey(data[q].compat, w) &&
            haskey(data[q].compat[w], p) &&
            @test v in data[q].compat[w][p]
        end
    end

    # matrices for computing solution ordering
    V = zeros(Int, P, S) # matrix of version "goodness" ranks
    D = fill(Int[], P, S) # matrix of dependency vectors
    for (i, p) in enumerate(pkgs)
        versions = data[p].versions
        depends = data[p].depends
        for s = 1:S
            v = vers[i, s]
            isnothing(v) && continue
            k = findfirst(==(v), versions)
            V[i, s] = length(versions) - k + 1
            haskey(depends, v) || continue
            D[i, s] = indexin(depends[v], pkgs)
        end
    end

    # partial order on solutions
    opts = Set{Int}() # optimization set
    # ==ₛ(s::Int, t::Int) = all(V[i,s] == V[i,t] for i=1:P)
    function ≥ₛ(s::Int, t::Int)
        # reset opts = reqs
        union!(empty!(opts), indexin(reqs, pkgs))
        while true
            strict = false
            for i in opts
                V[i,s] < V[i,t] && return false
                V[i,s] > V[i,t] && (strict = true)
            end
            strict && return true
            @assert all(V[i,s] == V[i,t] for i in opts)
            n = length(opts)
            for i in opts
                union!(opts, D[i,s])
            end
            n == length(opts) && return true
        end
    end
    ≤ₛ(s::Int, t::Int) = ≥ₛ(t, s)

end

function find_undominated(opts, done, deps)
    if isempty(deps)
        deps = find new deps
        if isempty(deps)
            return isempty(opts)
        end
    end
    for p in deps
        m = maximum(s[p] for s in opts)
        for i = 1:m-1
            conf[p][i] > 0 && continue
            # update conf for p@i
            opts′ = filter(s->s[p] ≤ i, opts)
            done′ = push!(copy(done), p)
            deps′ = delete!(copy(deps), p)
            if find_undominated(opts′, done′, deps′)
                return true
            end
        end
    end
    return false
end
