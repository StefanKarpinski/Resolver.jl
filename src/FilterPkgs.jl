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
    # meaning (both map packages to version indices):
    #   - reach tracks fully processed reachable versions
    #   - queue tracks newly reachable versions not yet processed
    reach = Dict{P,Int}(p => 0 for p in reqs)
    queue = Dict{P,Int}(p => 1 for p in reqs)

    # add next active version of p *after* i to the queue
    # do nothing if there's already version > i in reach/queue
    function next(p::P, i::Int)
        get(reach, p, 0) > i && return false
        get(queue, p, 0) > i && return false
        info_p = info[p]
        m = length(info_p.versions)
        X = info_p.conflicts
        # first active version after i (the flag column is the active set)
        j = col_min_from(X, size(X, 2), i + 1, m)
        # j == 0: we're out of versions (i.e. saturated; see below)
        queue[p] = j == 0 ? m + 1 : j
        return true
    end

    rdeps = Dict{P,Dict{P,Int}}() # reverse dependency map
    # p => q => k means
    #   "k is latest reachable version of q that depends on p"
    # add new reverse dependency:
    function rdep(p::P, q::P, k::Int)
        rdeps_p = get!(() -> valtype(rdeps)(), rdeps, p)
        rdeps_p[q] = max(get(rdeps_p, q, 0), k)
    end

    # notation:
    #   - p, q: packages
    #   - info_p = info[p]
    #   - indices of versions of package p: i, j
    #   - indices of versions of package q: k
    while !isempty(queue)
        # get unprocessed package + version
        p, i = pop!(queue)
        info_p = info[p]
        # check for saturation
        #   p saturated means: conflicts can force p to be uninstallable
        #   saturation represented by i > length(info_p.versions)
        if i > length(info_p.versions)
            # p has become saturated
            haskey(rdeps, p) && for (q, k) in rdeps[p]
                # q@k depends on p, therefore
                # q@k conflicts with p being uninstallable
                # p being saturated means that can happen
                next(q, k)
            end
        end
        # process each newly reachable version of p
        for j = get(reach, p, 0)+1:min(i, length(info_p.versions))
            # dependencies
            for (k, q) in enumerate(info_p.depends)
                info_p.conflicts[j, k] || continue
                rdep(q, p, j) # p@j depends on q
                next(q, 0) # q can be required
                # check if q is saturated:
                if get(reach, q, 0) > length(info[q].versions)
                    # p@j depends on q, therefore
                    # p@j conflicts with q being uninstallable
                    # q being saturated means that can happen
                    next(p, j)
                end
            end
            # find p@j's conflicts with reached versions of each partner:
            # conflicts are stored symmetrically, so scan the partner-side
            # column Y[k, c+j] == X[j, b+k], which is contiguous
            for q in keys(info_p.interacts)
                r = min(get(reach, q, 0), length(info[q].versions))
                r > 0 || continue
                info_q = info[q]
                k = col_max_upto(info_q.conflicts, info_q.interacts[p] + j, r)
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

    return reach
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
    work = copy(keys(info))
    names = sort!(collect(work))
    # some work buffers
    A = UInt64[]        # active version mask
    C = UInt64[]        # domination candidates mask
    J = Int[]           # active versions vector
    R = Int[]           # redundant indices vector
    cols = Vector{Int}[] # per-version active constraint columns
    # initialize active column flags
    for (p, info_p) in info
        X = info_p.conflicts
        m = length(info_p.versions)
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
        for (q, b) in info_p.interacts
            Y = info[q].conflicts
            for k = 1:length(info[q].versions)
                X[m+1, b+k] &= Y[k, end]
            end
        end
    end
    sort!(names, by = p -> length(info[p].interacts))
    for p in Iterators.cycle(names)
        isempty(work) && break
        p in work || continue
        delete!(work, p)
        # get conflicts & dimensions
        info_p = info[p]
        X = info_p.conflicts
        m = length(info_p.versions)
        m > 1 || continue # unique version cannot be redundant
        n = size(X, 2) - 1
        W = col_words(X)
        # active version mask & indices
        resize!(A, W)
        col_copy!(A, X, n + 1)
        clear_rows_above!(A, m)
        empty!(J)
        @inbounds for w = 1:W
            c = A[w]
            while !iszero(c)
                push!(J, ((w - 1) << 6) + trailing_zeros(c) + 1)
                c &= c - 1
            end
        end
        length(J) > 1 || continue # unique version cannot be redundant
        # collect every active version's active constraint columns in one
        # sweep over the columns, each of which is a contiguous bit span
        while length(cols) < m
            push!(cols, Int[])
        end
        for j in J
            empty!(cols[j])
        end
        for k = 1:n
            X[m+1, k] || continue # inactive column
            base = (k - 1) * W
            @inbounds for w = 1:W
                c = X.chunks[base + w] & A[w]
                while !iszero(c)
                    i = ((w - 1) << 6) + trailing_zeros(c) + 1
                    push!(cols[i], k)
                    c &= c - 1
                end
            end
        end
        # find redundant versions: i < j and X[i, k] => X[j, k] for all
        # active k means an earlier version is strictly more compatible,
        # therefore i will always be chosen instead of j. the candidates
        # dominated by i are the active versions worse than i that have a
        # constraint in every column i does — the word-parallel
        # intersection of i's constraint columns
        empty!(R)
        resize!(C, W)
        for i in J
            getbit(A, i) || continue # dominated versions can't dominate
            copyto!(C, A)
            clear_rows_upto!(C, i)
            live = false
            @inbounds for w = 1:W
                live |= !iszero(C[w])
            end
            live || continue
            ok = true
            for k in cols[i]
                col_and!(C, X, k) && continue
                ok = false
                break
            end
            ok || continue
            # remove dominated versions from the active set & record them
            @inbounds for w = 1:W
                c = C[w]
                while !iszero(c)
                    push!(R, ((w - 1) << 6) + trailing_zeros(c) + 1)
                    c &= c - 1
                end
                A[w] &= ~C[w]
            end
        end
        isempty(R) && continue
        sort!(R)
        # deactivate redundant versions
        X[R, end] .= false
        for q in keys(info_p.interacts)
            b = info[q].interacts[p]
            info[q].conflicts[length(info[q].versions) + 1, b .+ R] .= false
            push!(work, q) # can create new redundancies
        end
    end
end

function drop_unmarked!(
    info′ :: Dict{P, <: PkgInfo{P}},
    info  :: Dict{P, <: PkgInfo{P}} = info′,
) where {P}
    # info[p].conflicts version-flag column bits definitive
    # info[p].conflicts column-flag row bits made to match
    for (p, info_p) in info
        X = info_p.conflicts
        m = length(info_p.versions)
        X[m+1, :] .= false
        # dependency column is active if:
        # - some active version has that dependency
        for i = 1:m
            X[i, end] || continue # skip inactive versions
            for (k, q) in enumerate(info_p.depends)
                X[m+1, k] |= X[i, k]
            end
        end
        # conflict column is active if:
        # - the conflicting version is active
        for (q, b) in info_p.interacts
            Y = info[q].conflicts
            for j = 1:length(info[q].versions)
                X[m+1, b + j] = Y[j, end]
            end
        end
    end
    # save original version counts (needed if info′ === info)
    N = Dict{P,Int}(
        p => length(info_p.versions)
        for (p, info_p) in info
    )
    # go through again and shrink each PkgInfo
    for (p, info_p) in info
        # abbreviate components
        V = info_p.versions
        D = info_p.depends
        T = info_p.interacts
        X = info_p.conflicts
        m = length(V)
        n = size(X, 2) - 1
        # active version & column masks
        I = X[1:m, end]
        K = X[m+1, 1:n]
        # delete if no active versions
        if !any(I)
            delete!(info′, p)
            continue
        end
        # copy if everything is active
        if all(I) && all(K)
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
        m′ = length(V′)
        X′ = falses(padded_rows(m′), b′ + 1)
        rows = findall(I)
        j′ = 0
        for j = 1:n
            K[j] || continue
            j′ += 1
            for (i′, i) in enumerate(rows)
                X′[i′, j′] = X[i, j]
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
