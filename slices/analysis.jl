using LinearAlgebra
using ProgressMeter

# all package versions in info
V = [p => v for (p, info_p) in info for v=1:length(info_p.versions)]

# total depenency graph
n = length(V)
D = spzeros(Bool, n, n)
for (j, (p, v)) in enumerate(V)
    info_p = info[p]
    for (i, (q, w)) in enumerate(V)
        k = findfirst(==(q), info_p.depends)
        k === nothing && continue
        info_p.conflicts[v, k] || continue
        D[i, j] = true
    end
end
# D[i,j] iff V[i] satisfies a dependency of V[j]

# possible packages
pkgs = sort!(unique(first.(V)), by=popularity)

# try populating the "cache" and then hit ratio
let
    N = 1000 # iterations
    r = 10 # number of roots
    hit = tot = 0
    cache = Dict{Vector{Int},Int}()
    prog = Progress(N)
    m = length(pkgs)
    for time = 1:N
        reqs = pkgs[rand(m) .< r/m]
        sol = only(resolve_core(sat, reqs; max=1, by=popularity))
        s = indexin(sol, V)
        # compute precompile sets
        Ds = (D[s,s] + I) .> 0
        while true
            Ds′ = Ds^2 .> 0
            Ds′ == Ds && break
            Ds .= Ds′
        end
        # check/populate cache
        for i = 1:size(Ds,2)
            p = findall(Ds[:,i])
            hit += p ∈ keys(cache)
            tot += 1
            cache[p] = time
        end
        # evict from cache (LRU)
        max_cache = 1_000
        if length(cache) > max_cache
            c = sort!(collect(cache), by=reverse)
            for i = 1:length(cache)-max_cache
                delete!(cache, first(c[i]))
            end
        end
        @assert length(cache) ≤ max_cache
        # progress update
        next!(prog, showvalues = [
            ("cache", length(cache)),
            ("hits", hit),
            ("total", tot),
            ("ratio", hit/tot),
        ])
    end
end
