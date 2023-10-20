"""
    find_reachable(info, init) :: Dict{P, Int}

This function finds a minimal "reachable" subset of package and versions that
could appear in pareto-optimal solutions to version resolution for the given set
of required "root" packages, using the following recursive logic:

- P in reqs => P[1] reachable
- P[i] reachable & P[i] depends on D => D[1] reachable
- P[i] reachable & P[i] conflicts w. reachable => P[i+1] reachable

The function returns a dictionary mapping packages to the maximum version index
of that package that could be reached in an optimal solution. If a pacakge
cannot appear in an optimal solution, it will not appear in this dictionary.
"""
function find_reachable(
    info :: Dict{P, PkgInfo{P,V}},
    reqs :: SetOrVec{P},
) where {P,V}
    reach = Dict{P,Int}(p => 0 for p in reqs)
    queue = Dict{P,Int}(p => 1 for p in reqs)
    # meaning:
    #   - both map packages to version indices
    #   - reach tracks fully processed reachable versions
    #   - queue tracks newly reachable versions not yet processed

    # add to the work queue unless already processed
    function enqueue(p::P, i::Int)
        i > get(reach, p, 0) || return false
        i > get(queue, p, 0) || return false
        queue[p] = min(i, length(info[p].versions)+1)
        return true
    end

    rdeps = Dict{P, Dict{P,Int}}() # reverse dependency map
    # p => q => k means k is latest version of q that depends on p

    # add new reverse dependency
    function add_rdep(p::P, q::P, k::Int)
        rdeps_q = get!(() -> valtype(rdeps)(), rdeps, q)
        rdeps_q[p] = max(get(rdeps_q, p, 0), k)
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
                enqueue(q, k+1)
            end
        end
        # process each newly reachable version of p
        for j = get(reach, p, 0)+1:min(i, length(info_p.versions))
            # dependencies
            for (k, q) in enumerate(info_p.depends)
                info_p.conflicts[j, k] || continue
                # p@j depends on q
                enqueue(q, 1) # q@1 is reachable
                add_rdep(q, p, j)
                # check if q is saturated:
                if get(reach, q, 0) > length(info[q].versions)
                    # p@j conflicts with q being uninstallable
                    # q being saturated means that can happen
                    enqueue(p, j+1)
                end
            end
            # find all p@j's conflicts
            for (q, b) in info_p.interacts
                for k = 1:min(get(reach, q, 0), length(info[q].versions))
                    info_p.conflicts[j, b+k] || continue
                    # p@j conflicts with q@k
                    enqueue(p, j+1)
                    enqueue(q, k+1)
                end
            end
        end
        # update the reach map
        reach[p] = i
    end

    return reach
end
