module Resolver

export resolve

# for debugging
dd(A::AbstractArray) = map(reverse∘bitstring, permutedims(A))

function resolve(
    packages  :: AbstractVector{<:AbstractVector{<:Integer}},
    conflicts :: AbstractVector{<:Tuple{Integer,Integer}};
    Block     :: Type{<:Base.BitUnsigned} = UInt,
    sort      :: Bool = true,
)
    # vector of solution vectors
    solutions = Vector{Int}[]

    # counts & sizes
    N = length(packages)                    # number of packages
    M = mapreduce(maximum, max, packages)   # number of versions
    d = 8*sizeof(Block)                     # size of a version block
    m = div(M, d, RoundUp)                  # number of version blocks

    # check packages
    let counts = zeros(Int, M)
        for (p, versions) in enumerate(packages), v in versions
            1 ≤ v ≤ M ||
                throw(ArgumentError("invalid version index: $v"))
            counts[v] += 1
        end
        for v = 1:M
            counts[v] == 1 ||
                throw(ArgumentError("version $v in $(counts[v]) packages"))
        end
    end

    # check conflicts
    for (v1, v2) in conflicts
        1 ≤ v1 ≤ M  || throw(ArgumentError("invalid version index: $v1"))
        1 ≤ v2 ≤ M  || throw(ArgumentError("invalid version index: $v2"))
    end

    # no versions, no solutions
    M == 0 && return solution

    # construct blocked package matrix
    P = zeros(Block, m, M)
        # each column is for a version
        # each row is a block of package bitmask

    # bitmask compressed indicators of same package versions
    for (p, versions) in enumerate(packages)
        for v1 in versions, v2 in versions
            b, s = divrem(v2-1, d)
            P[b+1, v1] |= 1 << s
        end
    end

    # construct blocked conflicts matrix
    X = copy(P) # same package versions are incompatible
        # each column is for a version
        # each row is a block of conflict bitmask

    # invert the packages matrix
    @. P = ~P
    # turn off the extra bits in the last block of each column
    for j = 1:M
        P[m, j] &= typemax(Block) >> mod(-M, d)
    end

    # explicit conflicts are also incompatible
    for (v1, v2) in conflicts
        # conflicts are symmetrized
        b1, s1 = divrem(v1-1, d)
        b2, s2 = divrem(v2-1, d)
        X[b1+1, v2] |= 1 << s1
        X[b2+1, v1] |= 1 << s2
    end

    # allocate iteration candidates matrix
    C = zeros(Block, m, N)
        # each column is for a recursion level
        # each row is a block of candidate bitmask
    # turn all candidates on in first column
    for i = 1:m-1
        C[i, 1] = typemax(Block)
    end
    # except the extra bits in the last block
    let s = mod(-M, d)
        C[m, 1] = typemax(Block) << s >> s
    end

    # allocate recursion candidates matrix
    R = zeros(Block, m, N-1)
        # each column is for a recursion level
        # each row is a block of candidate bitmask

    # allocate selections vector
    S = zeros(Int, N)

    function search!(r::Int=1)
        # copy recursion candidates from iteration candidates
        if r < N
            for i = 1:m
                R[i, r] = C[i, r]
            end
        end
        # initial block & shift values
        b = s = 0
        while true
            # look for the next candidate
            let c = C[b+1, r]
                if c >> s == 0 # no more candidates in current block
                    b += 1 # next block
                    while b < m
                        c = C[b+1, r] # candidates bitmask
                        c != 0 && break # there's one in this block
                        b += 1
                    end
                end
                # shift for next candidate
                s += trailing_zeros(c >> s)
            end
            v = b*d + s + 1
            v ≤ M || break
            # record candidate version
            S[r] = v
            # recurse or save solution
            if r < N
                for i = 1:m
                    x = X[i, v]
                    # next iteration: only conflicts
                    C[i, r] &= x
                    # next recursion: skip conflicts
                    C[i, r+1] = R[i, r] & ~x
                end
                # do recursive search
                search!(r+1)
            else # complete solution
                # record solution
                push!(solutions, copy(S))
                # stop search for each found package
                for r′ = 1:N, i = 1:m
                    C[i, r′] &= P[i, S[r′]]
                end
                break
            end
            # next candidate
            s += 1
        end
    end
    search!()

    # return sorted vector of sorted solutions
    sort && sort!(map(sort!, solutions))
    return solutions
end

end # module
