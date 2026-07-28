function filter_pkg_info!(
    info :: Dict{P, PkgInfo{P,V}},
    reqs :: SetOrVec{P} = keys(info),
) where {P,V}
    mark_reachable!(info, reqs)
    mark_necessary!(info)
    drop_unmarked!(info)
end

"""
    find_reachable(info, reqs) :: Dict{P, Int}

This function finds a minimal "reachable" subset of package and versions that
could appear in pareto-optimal solutions to version resolution for the given set
of required "root" packages, using the following recursive logic:

- P in reqs => P[1] reachable
- P[i] reachable & P[i] depends on D => D[1] reachable
- P[i] reachable & P[i] conflicts w. reachable => P[i+1] reachable
- D[end] conflicts w. reachable & P[i] depends on D => P[i+1] reachable

The function returns a dictionary mapping packages to the maximum version index
of that package that could be reached in an optimal solution. If a package
cannot appear in an optimal solution, it will not appear in this dictionary.
"""
function find_reachable(
    info :: Dict{P, PkgInfo{P,V}},
    reqs :: SetOrVec{P} = keys(info),
) where {P,V}
    # intern packages as integer indices so that the fixpoint loop below
    # hashes no package names; all derived tables are indexed by pkg id
    pkgs = sort!(collect(keys(info)))
    N = length(pkgs)
    ix = Dict{P,Int}(p => i for (i, p) in enumerate(pkgs))
    infos = Vector{PkgInfo{P,V}}(undef, N)
    nvers = Vector{Int}(undef, N)
    deps = Vector{Vector{Int}}(undef, N)     # interned depends, same order
    partners = Vector{Vector{Tuple{Int,Int}}}(undef, N)
    for p = 1:N
        info_p = info[pkgs[p]]
        infos[p] = info_p
        nvers[p] = length(info_p.versions)
        deps[p] = Int[ix[q] for q in info_p.depends]
        # (partner id, offset of p's version block in the partner's matrix)
        prt = Tuple{Int,Int}[
            (ix[q], info[q].interacts[pkgs[p]])
            for q in keys(info_p.interacts)
        ]
        partners[p] = sort!(prt)
    end

    # meaning (both map packages to version indices):
    #   - reach tracks fully processed reachable versions
    #   - queue tracks newly reachable versions not yet processed
    reach = zeros(Int, N)
    queue = zeros(Int, N)
    stack = Int[] # packages with queued work
    instack = falses(N)

    function enqueue(p::Int, j::Int)
        queue[p] = j
        if !instack[p]
            instack[p] = true
            push!(stack, p)
        end
    end

    # add next active version of p *after* i to the queue
    # do nothing if there's already version > i in reach/queue
    function next(p::Int, i::Int)
        reach[p] > i && return false
        queue[p] > i && return false
        m = nvers[p]
        X = infos[p].conflicts
        # first active version after i (the flag column is the active set)
        j = col_min_from(X, size(X, 2), i + 1, m)
        # j == 0: we're out of versions (i.e. saturated; see below)
        enqueue(p, j == 0 ? m + 1 : j)
        return true
    end

    rdeps = [Dict{Int,Int}() for _ = 1:N] # reverse dependency map
    # rdeps[p][q] == k means
    #   "k is latest reachable version of q that depends on p"

    for p in reqs
        enqueue(ix[p], 1)
    end

    # notation:
    #   - p, q: packages
    #   - info_p = infos[p]
    #   - indices of versions of package p: i, j
    #   - indices of versions of package q: k
    while !isempty(stack)
        # get unprocessed package + version
        p = pop!(stack)
        instack[p] = false
        i = queue[p]
        queue[p] = 0
        info_p = infos[p]
        m = nvers[p]
        # check for saturation
        #   p saturated means: conflicts can force p to be uninstallable
        #   saturation represented by i > nvers[p]
        if i > m
            # p has become saturated
            for (q, k) in rdeps[p]
                # q@k depends on p, therefore
                # q@k conflicts with p being uninstallable
                # p being saturated means that can happen
                next(q, k)
            end
        end
        # process each newly reachable version of p
        for j = reach[p]+1:min(i, m)
            # dependencies
            for (k, q) in enumerate(deps[p])
                info_p.conflicts[j, k] || continue
                # p@j depends on q
                rdeps_q = rdeps[q]
                rdeps_q[p] = max(get(rdeps_q, p, 0), j)
                next(q, 0) # q can be required
                # check if q is saturated:
                if reach[q] > nvers[q]
                    # p@j depends on q, therefore
                    # p@j conflicts with q being uninstallable
                    # q being saturated means that can happen
                    next(p, j)
                end
            end
            # find p@j's conflicts with reached versions of each partner:
            # conflicts are stored symmetrically, so scan the partner-side
            # column Y[k, c+j] == X[j, b+k], which is contiguous
            for (q, c) in partners[p]
                r = min(reach[q], nvers[q])
                r > 0 || continue
                k = col_max_upto(infos[q].conflicts, c + j, r)
                k == 0 && continue
                # p@j conflicts with q@k (the highest such k; pushing past
                # it subsumes pushing past any lower conflicting version)
                next(p, j)
                next(q, k)
            end
        end
        # update the reach map
        reach[p] = i
    end

    return Dict{P,Int}(pkgs[p] => reach[p] for p = 1:N if reach[p] > 0)
end

function mark_reachable!(
    info :: Dict{P, PkgInfo{P,V}},
    reqs :: SetOrVec{P} = keys(info),
) where {P,V}
    reach = find_reachable(info, reqs)
    for (p, info_p) in info
        r = min(get(reach, p, 0), length(info_p.versions))
        info_p.conflicts[1:r, end] .= true
        info_p.conflicts[r+1:end, end] .= false
    end
    return reach
end

function mark_necessary!(
    info :: Dict{P, PkgInfo{P,V}},
) where {P,V}
    # intern packages as integer indices so that the work loop below
    # hashes no package names (sorted so the processing order — and with
    # it any order-dependent tie — is deterministic)
    pkgs = sort!(collect(keys(info)))
    N = length(pkgs)
    ix = Dict{P,Int}(p => i for (i, p) in enumerate(pkgs))
    infos = Vector{PkgInfo{P,V}}(undef, N)
    nvers = Vector{Int}(undef, N)
    partners = Vector{Vector{NTuple{3,Int}}}(undef, N)
    for p = 1:N
        info_p = info[pkgs[p]]
        infos[p] = info_p
        nvers[p] = length(info_p.versions)
        # (partner id, partner's version block offset in p's matrix,
        #  p's version block offset in the partner's matrix)
        prt = NTuple{3,Int}[
            (ix[q], b, info[q].interacts[pkgs[p]])
            for (q, b) in info_p.interacts
        ]
        partners[p] = sort!(prt)
    end
    # some work buffers
    A = UInt64[]        # active version mask
    D = UInt64[]        # per-version domination candidate masks
    T = UInt64[]        # mask of versions still tracked by the sweep
    R = Int[]           # redundant indices vector
    # initialize active column flags
    for p = 1:N
        info_p = infos[p]
        X = info_p.conflicts
        m = nvers[p]
        n = size(X, 2) - 1
        resize!(A, col_words(X))
        col_copy!(A, X, n + 1)
        clear_rows_above!(A, m)
        for j = 1:n
            # we don't have to look at columns that have no
            # conflicts for any active versions of package
            X[m+1, j] = col_intersects(X, j, A)
        end
        # conflict columns for inactive partner versions are also
        # irrelevant: no model over the active versions can include
        # the partner, so the conflict cannot constrain anything
        for (q, b, c) in partners[p]
            Y = infos[q].conflicts
            for k = 1:nvers[q]
                X[m+1, b+k] &= Y[k, end]
            end
        end
    end
    # process every package once, cheapest first (stable sort breaks ties
    # by name); packages woken by a partner's deletions are appended to
    # the queue unless already pending
    queue = sort!(collect(1:N), by = p -> length(partners[p]))
    inwork = trues(N)
    head = 1
    while head ≤ length(queue)
        p = queue[head]
        head += 1
        inwork[p] || continue
        inwork[p] = false
        # get conflicts & dimensions
        info_p = infos[p]
        X = info_p.conflicts
        m = nvers[p]
        m > 1 || continue # unique version cannot be redundant
        n = size(X, 2) - 1
        W = col_words(X)
        # active version mask
        resize!(A, W)
        col_copy!(A, X, n + 1)
        clear_rows_above!(A, m)
        nact = 0
        @inbounds for w = 1:W
            nact += count_ones(A[w])
        end
        nact > 1 || continue # unique version cannot be redundant
        # find redundant versions: i < j and X[i, k] => X[j, k] for all
        # active k means an earlier version is strictly more compatible,
        # therefore i will always be chosen instead of j. the candidates
        # dominated by i are the active versions worse than i that have a
        # constraint in every column i does. sweep the active columns
        # once, restricting the candidates of every version with a bit in
        # the column to the column's row set; versions whose candidate
        # set empties drop out of the sweep. a version dominated only by
        # dominated versions is also dominated by their (transitively
        # live) dominators, so batching over the initial active set finds
        # exactly the sequentially dominated versions.
        resize!(D, m * W)
        resize!(T, W)
        copyto!(T, A)
        @inbounds for w = 1:W
            c = A[w]
            while !iszero(c)
                i = ((w - 1) << 6) + trailing_zeros(c) + 1
                c &= c - 1
                # candidates of i: active versions worse than i
                o = (i - 1) * W
                for w′ = 1:W
                    D[o + w′] = A[w′]
                end
                wt = i >> 6
                for w′ = 1:min(wt, W)
                    D[o + w′] = 0
                end
                r = i & 63
                if r != 0 && wt < W
                    D[o + wt + 1] &= ~((UInt64(1) << r) - 1)
                end
                live = UInt64(0)
                for w′ = 1:W
                    live |= D[o + w′]
                end
                iszero(live) && (T[w] &= ~(UInt64(1) << ((i - 1) & 63)))
            end
        end
        for k = 1:n
            X[m+1, k] || continue # inactive column
            base = (k - 1) * W
            @inbounds for w = 1:W
                c = X.chunks[base + w] & T[w]
                while !iszero(c)
                    i = ((w - 1) << 6) + trailing_zeros(c) + 1
                    c &= c - 1
                    o = (i - 1) * W
                    live = UInt64(0)
                    for w′ = 1:W
                        live |= (D[o + w′] &= X.chunks[base + w′])
                    end
                    iszero(live) && (T[w] &= ~(UInt64(1) << ((i - 1) & 63)))
                end
            end
        end
        # union of all candidate sets = the dominated versions
        fill!(T, UInt64(0)) # reuse as the union accumulator
        @inbounds for w = 1:W
            c = A[w]
            while !iszero(c)
                i = ((w - 1) << 6) + trailing_zeros(c) + 1
                c &= c - 1
                o = (i - 1) * W
                for w′ = 1:W
                    T[w′] |= D[o + w′]
                end
            end
        end
        empty!(R)
        @inbounds for w = 1:W
            c = T[w]
            while !iszero(c)
                push!(R, ((w - 1) << 6) + trailing_zeros(c) + 1)
                c &= c - 1
            end
        end
        isempty(R) && continue
        # deactivate redundant versions
        X[R, end] .= false
        for (q, b, c) in partners[p]
            infos[q].conflicts[nvers[q] + 1, c .+ R] .= false
            if !inwork[q] # can create new redundancies
                inwork[q] = true
                push!(queue, q)
            end
        end
    end
end

function drop_unmarked!(
    info′ :: Dict{P, <: PkgInfo{P}},
    info  :: Dict{P, <: PkgInfo{P}} = info′,
) where {P}
    # first pass: per package, compute the kept versions and kept columns
    # from the (definitive) version flags — all of it before any package
    # is rebuilt, since kept columns are read off the partners' matrices
    masks = Dict{P, Tuple{BitVector, BitVector}}()
    A = UInt64[] # active version mask buffer
    for (p, info_p) in info
        X = info_p.conflicts
        m = length(info_p.versions)
        n = size(X, 2) - 1
        W = col_words(X)
        # active versions
        I = X[1:m, end]
        resize!(A, W)
        col_copy!(A, X, n + 1)
        clear_rows_above!(A, m)
        # kept columns: dependency columns some active version depends
        # on, and conflict columns of active partner versions — a
        # partner's column block must stay aligned with its version list,
        # so kept partner versions keep their columns even if empty
        K = falses(n)
        for k = 1:length(info_p.depends)
            K[k] = col_intersects(X, k, A)
        end
        for (q, b) in info_p.interacts
            Y = info[q].conflicts
            copy_col_bits!(K.chunks, b, Y, size(Y, 2), length(info[q].versions))
        end
        masks[p] = (I, K)
    end
    # save original version counts (needed if info′ === info)
    N = Dict{P,Int}(
        p => length(info_p.versions)
        for (p, info_p) in info
    )
    # second pass: shrink each PkgInfo
    for (p, info_p) in info
        # abbreviate components
        V = info_p.versions
        D = info_p.depends
        T = info_p.interacts
        X = info_p.conflicts
        m = length(V)
        n = size(X, 2) - 1
        W = col_words(X)
        I, K = masks[p]
        m′ = count(I)
        # delete if no active versions
        if m′ == 0
            delete!(info′, p)
            continue
        end
        # keep as is if everything is active (with the column flags
        # normalized, as rebuilding would leave them)
        if m′ == m && all(K)
            X[m+1, 1:n] .= true
            if info !== info′
                V′ = copy(V)
                D′ = copy(D)
                T′ = copy(T)
                X′ = copy(X)
                info′[p] = PkgInfo(V′, D′, T′, X′)
            end
            continue
        end
        # compute shrunken components
        V′ = V[I]
        D′ = D[K[1:length(D)]]
        T′ = Dict{P,Int}()
        b′ = length(D′)
        for (q, b) in sort!(collect(T), by=last)
            n′ = count(K[b .+ (1:N[q])])
            if n′ > 0
                T′[q] = b′
                b′ += n′
            end
        end
        R′ = padded_rows(m′)
        W′ = R′ >> 6
        X′ = falses(R′, b′ + 1)
        src = X.chunks
        dst = X′.chunks
        # kept-version mask words for the gather below
        resize!(A, W)
        col_copy!(A, X, n + 1)
        clear_rows_above!(A, m)
        prefix = findlast(I) == m′ # kept versions are 1:m′
        j′ = 0
        for j = 1:n
            K[j] || continue
            j′ += 1
            sb = (j - 1) * W
            db = (j′ - 1) * W′
            if prefix
                # kept versions are a prefix: straight word copy
                nw = (m′ + 63) >> 6
                @inbounds for w = 1:nw
                    dst[db + w] = src[sb + w]
                end
                r = m′ & 63
                r != 0 && @inbounds (dst[db + nw] &= (UInt64(1) << r) - 1)
            else
                # gather the kept versions' bits
                acc = UInt64(0)
                na = 0
                dw = db + 1
                @inbounds for w = 1:W
                    avail = A[w]
                    iszero(avail) && continue
                    v = src[sb + w]
                    while !iszero(avail)
                        acc |= ((v >> trailing_zeros(avail)) & 1) << na
                        avail &= avail - 1
                        na += 1
                        if na == 64
                            dst[dw] = acc
                            dw += 1
                            acc = UInt64(0)
                            na = 0
                        end
                    end
                end
                na > 0 && @inbounds (dst[dw] = acc)
            end
        end
        @assert j′ == b′
        X′[1:m′, end] .= true  # kept versions are active
        X′[m′+1, 1:b′] .= true # kept columns are active
        # assign new struct into info′
        info′[p] = PkgInfo(V′, D′, T′, X′)
    end
    return info′
end

function drop_excluded!(
    info :: Dict{P, PkgInfo{P,V}},
    drop :: SetOrVec{P},
) where {P,V}
    drop isa AbstractSet || (drop = Set{P}(drop))
    dirty = false
    for (p, info_p) in info
        # deactivate all versions of dropped packages
        if p ∈ drop
            info_p.conflicts[:, end] .= false
            dirty = true
            continue
        end
    end
    while dirty
        # deactivate versions that depend on dropped packages
        dirty = false
        for (p, info_p) in info
            X = info_p.conflicts
            for (k, q) in enumerate(info_p.depends)
                Y = info[q].conflicts
                any(Y[i, end] for i = 1:length(info[q].versions)) && continue
                # no active versions of q, delete p@i that depend on it
                for i = 1:length(info_p.versions)
                    info_p.conflicts[i, end] || continue
                    info_p.conflicts[i, k] || continue
                    # p@i depends on q, so drop p@i too
                    info_p.conflicts[i, end] = false
                    dirty = true
                end
            end
        end
    end
    drop_unmarked!(info)
end

function check_info_structure(
    info :: Dict{P, PkgInfo{P,V}},
) where {P,V}
    for (p, info_p) in info
        for q in info_p.depends
            q in keys(info) ||
                return :depends, p, q
        end
        interacts = sort!(collect(info_p.interacts), by=last)
        for (i, (q, b)) in enumerate(interacts)
            if i == 1
                b == length(info_p.depends) ||
                    return :conflicts, p, q
            end
            if i < length(interacts)
                b + length(info[q].versions) == interacts[i+1][2] ||
                    return :conflicts, p, q
            end
        end
    end
end
