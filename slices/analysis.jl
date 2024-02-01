using LinearAlgebra
using ProgressMeter

# all package versions in info
V = [p => v for (p, info_p) in info for v=1:length(info_p.versions)]
sort!(V, by=popularity)

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

# packages we sample from
pkgs = sort!(unique(first.(vertices)), by=popularity)
cdf = cumsum(Float64[Slices.NREQS[p] for p in pkgs])
cdf ./= cdf[end]
function sample_package()
    x = rand()
    for (i, c) in enumerate(cdf)
        c > x && return pkgs[i]
    end
    # shouldn't happen
    return pkgs[end]
end
sample_packages(r::Integer) =
    sort!(unique(sample_package() for _ = 1:r), by=popularity)

# try populating the "cache" and then hit ratio
function cache_test(
    solve      :: Function;
    sample     :: Integer = 5,
    iterations :: Integer = 2000,
    cache_size :: Integer = 10000,
)
    hit = tot = 0
    cache = Dict{Vector{Int},Int}()
    prog = Progress(iterations)
    prob = sample/length(pkgs)
    for time = 1:iterations
        reqs = sample_packages(sample)
        sol = solve(reqs)
        s = indexin(sol, V)
        # compute precompile sets
        Ds = (D[s,s] + I) .> 0
        while true
            Ds′ = Ds*Ds .> 0
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
        if length(cache) > cache_size
            c = sort!(collect(cache), by=reverse)
            for i = 1:length(cache)-cache_size
                delete!(cache, first(c[i]))
            end
        end
        @assert length(cache) ≤ cache_size
        # progress update
        next!(prog, showvalues = [
            ("cache", length(cache)),
            ("hits", hit),
            ("total", tot),
            ("ratio", hit/tot),
        ])
    end
end

local_resolve(reqs) =
    only(resolve_core(sat, reqs, max=1, by=popularity))

# per-slice package sets
S = [Set(first.(vertices[s])) for s in slices]
DS = typeof(D)[]
let inds = indexin(vertices, V)
    for s in slices
        Ds = (D[inds[s],inds[s]] + I) .> 0
        while true
            Ds′ = Ds^2 .> 0
            Ds′ == Ds && break
            Ds .= Ds′
        end
        push!(DS, Ds)
    end
end

function slice_resolve(reqs::AbstractVector{String})
    for (i, Si) in enumerate(S)
        reqs ⊆ Si || continue
        s = slices[i]
        vers = vertices[s]
        P = first.(vers)
        x = Vector{Int}(indexin(reqs, P))
        v = DS[i][:,x]*fill(true, length(x))
        return vers[v .> 0]
    end
    # fall back to local resolve
    local_resolve(reqs)
end
