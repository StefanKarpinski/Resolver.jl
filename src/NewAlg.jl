function resolve(
    packages  :: AbstractVector{<:AbstractVector{<:Integer}},
    conflicts :: AbstractVector{<:Tuple{Integer,Integer}};
    Block     :: Type{<:Base.BitUnsigned} = UInt,
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

    # construct blocked packages matrix
    P = zeros(Block, m, M)
        # each column is for a version
        # each row is a block of package bitmask
    for (p, versions) in enumerate(packages)
        for v1 in versions, v2 in versions
            b, s = divrem(v2-1, d)
            P[b+1, v1] |= 1 << s
        end
    end

    # construct blocked conflicts matrix
    X = zeros(Block, m, M)
        # each column is for a version
        # each row is a block of conflict bitmask
    for (v1, v2) in conflicts
        # conflicts are symmetrized
        b1, s1 = divrem(v1-1, d)
        b2, s2 = divrem(v2-1, d)
        X[b1+1, v2] |= 1 << s1
        X[b2+1, v1] |= 1 << s2
    end

    # remove conflicts within a package
    # (this shouldn't happen; we could error instead)
    @. X &= ~P

    # allocate recursion candidates matrix
    R = zeros(Block, m, N)
        # each column is for a recursion level
        # each row is a block of candidate bitmask
    # turn all candidates on in first column
    for i = 1:m-1
        R[i, 1] = typemax(Block)
    end
    # except the extra bits in the last block
    let s = mod(-M, d)
        R[m, 1] = typemax(Block) << s >> s
    end

    # allocate iteration candidates matrix 
    C = zeros(Block, m, N)
        # each column is for a recursion level
        # each row is a block of candidate bitmask

    # allocate selections vector
    S = zeros(Int, N)

    function search!(r::Int=1)
        b = s = 0
        found = false
        # copy iteration candidates from recursion candidates
        @. C[:, r] = R[:, r]
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
            @assert C[b+1, r] & (1 << s) != 0
            # record candidate version
            S[r] = v
            # recurse or save solution
            if r < N
                # next recursion: skip this package and conflicts
                @. R[:, r+1] = R[:, r] & ~P[:, v] & ~X[:, v]
                # do recursive search
                if search!(r+1)
                    # next iteration: only conflicting
                    @. C[:, r] &= X[:, v]
                    found = true
                else
                    # next iteration: only conflicting or same package
                    @. C[:, r] &= X[:, v] | P[:, v]
                end
            else # record complete solution
                push!(solutions, copy(S))
                return true
            end
            # next candidate
            s += 1
        end
        return found
    end
    search!()

    return solutions
end

# debugging...

packages = [1:2, 3:6, 7:10]
conflicts = [(1,3), (1,7), (3,7), (2, 10)]
Block = UInt8

dd(A::AbstractArray) = map(reverse∘bitstring, permutedims(A))
