# we need to know package groupings so that we know when to stop searching a
# given package because we've already found a solution for it

function resolve!(
    emit       :: Function,
    candidates :: BitVector,
    compatbile :: Vector{BitVector},
    packages   :: Vector{Pair{Int,Int}},
    package    :: Int,
)
    if package ≥ length(packages)
        # emit solution
        return true
    end
    found_any = false
    i = 0
    while true
    end
end

function resolve!(candidates, n = 0)
    if n == n_packages
        # emit a solution
        return true
    end
    i = 0
    found_any = false
    while true
        # TODO: at the last level, we return a solution
        i = find_candidate(candidates, i) # only > i searched
        i === nothing && break
        candidates′ = candidates .& compatible[i] .& .!package[i]
        found = resolve!(candidates′, n+1)
        if found
            found_any = true
            # only: incompatible candidates from other packages
            candidates .&= .!compatible[i] .& .!package[i]
            # could also bump `i` past last package index
        elseif !package[i][i+1]
            # search failed: this package had no versions 
            break
        end
    end
    return found_any
end
