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
        n = length(info_p.versions)
        for j = i+1:n
            info_p.conflicts[j, end] || continue # inactive
            queue[p] = j
            return true
        end
        # we're out of versions (i.e. saturated; see below)
        queue[p] = n+1
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
            # find all p@j's conflicts
            for (q, b) in info_p.interacts
                for k = 1:min(get(reach, q, 0), length(info[q].versions))
                    info_p.conflicts[j, b+k] || continue
                    # p@j conflicts with q@k
                    next(p, j)
                    next(q, k)
                end
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
        r = get(reach, p, 0)
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
    # initialize active column flags
    for (p, info_p) in info
        X = info_p.conflicts
        m = size(X, 1)-1
        n = size(X, 2)-1
        for j = 1:n
            # we don't have to look at columns that have no
            # conflicts for any active versions of package
            X[end, j] = any(X[i, j] for i = 1:m if X[i, end])
        end
    end
    # some work vectors
    J = Int[] # active versions vector
    K = Int[] # active conflicts vector
    R = Int[] # redundant indices vector
    sort!(names, by = p -> length(info[p].interacts))
    for p in Iterators.cycle(names)
        isempty(work) && break
        p in work || continue
        delete!(work, p)
        # get conflicts & dimensions
        X = info[p].conflicts
        m = size(X, 1) - 1
        m > 1 || continue # unique version cannot be redundant
        n = size(X, 2) - 1
        # active indices
        append!(empty!(J), (j for j = 1:m if X[j, end]))
        length(J) > 1 || continue # unique version cannot be redundant
        append!(empty!(K), (k for k = 1:n if X[end, k]))
        # find redundant versions
        empty!(R)
        for j in J
            # don't combine loops--it changes what break & continue do
            for i in J
                i < j || break
                i ∈ R && continue
                all(!X[i, k] | X[j, k] for k in K) || continue
                # an earlier version is strictly more compatible
                # i.e. i < j and X[i, k] => X[j, k] for all k
                # therefore i will always be chosen instead of j
                push!(R, j)
                break
            end
        end
        isempty(R) && continue
        # deactivate redundant versions
        X[R, end] .= false
        for q in keys(info[p].interacts)
            b = info[q].interacts[p]
            info[q].conflicts[end, b .+ R] .= false
            push!(work, q) # can create new redundancies
        end
    end
end

function drop_unmarked!(
    info′ :: Dict{P, <: PkgInfo{P}},
    info  :: Dict{P, <: PkgInfo{P}} = info′,
) where {P}
    # info[p].conflicts[:, end] bits definitive (row flags)
    # info[p].conflicts[end, :] bits made to match (column flags)
    for (p, info_p) in info
        X = info_p.conflicts
        X[end, :] .= false
        # dependency column is active if:
        # - some active version has that dependency
        for i = 1:size(X,1)-1
            X[i, end] || continue # skip inactive versions
            for (k, q) in enumerate(info_p.depends)
                X[end, k] |= X[i, k]
            end
        end
        # conflict column is active if:
        # - the conflicting version is active
        for (q, b) in info_p.interacts
            Y = info[q].conflicts
            for j = 1:size(Y,1)-1
                X[end, b + j] = Y[j, end]
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
        # active version masks
        I = X[1:end-1, end]
        K = X[end, 1:end-1]
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
        X′ = X[[I; true], [K; true]]
        @assert size(X′, 1) == length(V′) + 1
        @assert size(X′, 2) == b′ + 1
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
                any(Y[i, end] for i=1:size(Y,1)-1) && continue
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
