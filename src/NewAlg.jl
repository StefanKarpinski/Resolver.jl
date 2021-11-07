function resolve(
    packages  :: AbstractVector{<:AbstractVector{<:Integer}},
    conflicts :: AbstractVector{Tuple{Integer,Integer}};
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

    # allocate candidates matrix
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

    # allocate selections vector
    S = zeros(Int, N)

    function search!(r::Int=1)
        found = false
        b = s = 0
        while true
            # look for the next candidate
            let c = C[b+1, r]
                if c >> s == 0
                    b += 1
                    while b < m
                        c = C[b+1, r]
                        c != 0 && break
                        b += 1
                    end
                end
                s = trailing_zeros(c)
            end
            v = b*d + s + 1
            v ≤ M || break
            # viable candidate found
            S[r] = v
            # recurse or record
            if r < N
                # compute iterative & recursive candidates
                for i = 1:m
                    c = C[i, r]
                    p = P[i, v]
                    x = X[i, v]
                    # itertion: skip non-conflicting versions
                    C[i, r] = c & x
                    # recursion: skip same package & conflicting versions
                    C[i, r+1] = c & ~p & ~x
                end
                # do recursive search
                if search!(r+1)
                    # skip worse versions of the same package
                    for i = 1:m
                        c = C[i, r]
                        p = P[i, v]
                        # itertion: skip same package versions
                        C[i, r] = c & ~p
                    end
                    found = true
                end
            else # record complete solution
                push!(solutions, copy(S))
                found = true
            end
            # next candidate
            s += 1
        end
        return found
    end
    search!()

    return solutions
end
