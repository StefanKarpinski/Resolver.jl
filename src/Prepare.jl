# Type Parameters:
#  P = package type (String)
#  V = version type (VersionNumber)
#  S = version set type (VersionSpec)

struct PkgInfo{P,V,S}
    versions :: Vector{V}
    depends  :: Dict{V, Vector{P}}
    compat   :: Dict{V, Dict{P, S}}
end

struct DepsProvider{P,V,S,F<:Function}
    packages :: Vector{P}
    provider :: F
end

const SetOrVec{T} = Union{AbstractSet{T}, AbstractVector{T}}

function DepsProvider{P,V,S}(
    provider :: Function,
    packages :: SetOrVec{P},
) where {P,V,S}
    packages = sort!(P[p for p in packages])
    DepsProvider{P,V,S,typeof(provider)}(packages, provider)
end

(deps::DepsProvider{P,V,S,F})(pkg::P) where {P,V,S,F<:Function} =
    deps.provider(pkg) :: PkgInfo{P,V,S}

# load all package info that might be needed for a set of requirements
# given no requirements, load all packages the provider knows about

function load_packages(
    deps :: DepsProvider{P,V,S},
    reqs :: SetOrVec{P} = deps.packages,
) where {P,V,S}
    pkgs = Dict{P, PkgInfo{P,V,S}}()
    work = Set(reqs)
    while !isempty(work)
        pkg = pop!(work)
        @assert pkg ∉ keys(pkgs)
        info = pkgs[pkg] = deps(pkg)
        for pkgs′ in values(info.depends), pkg′ in pkgs′
            pkg′ in keys(pkgs) && continue
            push!(work, pkg′)
        end
    end
    return pkgs
end

# two packages interact if there is some conflict between them

function find_interacts(
    pkgs :: Dict{P, PkgInfo{P,V,S}},
) where {P,V,S}
    interacts = Dict{P,Vector{P}}(p => P[] for p in keys(pkgs))
    for (pkg₁, info₁) in pkgs
        interact₁ = interacts[pkg₁]
        for ver₁ in info₁.versions
            ver₁ in keys(info₁.compat) || continue
            compat₁ = info₁.compat[ver₁]
            for (pkg₂, spec₁) in compat₁
                pkg₂ in interact₁ && continue
                interact₂ = interacts[pkg₂]
                for ver₂ in pkgs[pkg₂].versions
                    if ver₂ ∉ spec₁
                        push!(interact₁, pkg₂)
                        push!(interact₂, pkg₁)
                        break
                    else
                        compat₂ = pkgs[pkg₂].compat
                        pkg₁ in keys(compat₂) || continue
                        spec₂ = compat₂[pkg₁]
                        if ver₁ ∉ spec₂
                            push!(interact₁, pkg₂)
                            push!(interact₂, pkg₁)
                            break
                        end
                    end
                end
            end
        end
    end
    filter!(interacts) do (pkg, ix)
        !isempty(ix)
    end
    foreach(sort!, values(interacts))
    return interacts
end

# filter down to versions that resolve might actually pick

"""
    find_reachable(deps, reqs) :: Dict{P, Int}

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
    pkgs :: Dict{P, PkgInfo{P,V,S}},
    reqs :: Vector{P},
) where {P,V,S}
    interacts = find_interacts(pkgs)

    queue = Dict{P,Int}(p => 1 for p in reqs)
    reach = Dict{P,Int}(p => 0 for p in reqs)
    saturated = Set{P}()

    function add_queue!(p::P, k::Int)
        get(reach, p, 0) ≥ k && return false
        get(queue, p, 0) ≥ k && return false
        queue[p] = k
        return true
    end

    # empty deps & compat
    deps_∅ = Vector{P}()
    comp_∅ = Dict{P,S}()
    intx_∅ = deps_∅ # same structure

    while !isempty(queue)
        # get unprocessed package + version
        p, k = pop!(queue)
        j = get(reach, p, 0)
        # look up some stuff about p
        intx = get(interacts, p, intx_∅)
        info_p = pkgs[p]
        vers_p = info_p.versions
        deps_p = info_p.depends
        comp_p = info_p.compat
        # check for saturation
        if k > length(vers_p)
            push!(saturated, p)
            for (q, j) in reach
                info_q = pkgs[q]
                vers_q = info_q.versions
                1 ≤ j ≤ length(vers_q) || continue
                w = vers_q[j]
                if p in get(info_q.depends, w, deps_∅)
                    # q@j conflicts with p being uninstallable
                    # p in saturated means that can happen
                    add_queue!(q, j+1)
                end
            end
        end
        # main work loop
        for i = j+1:min(k, length(vers_p))
            v = vers_p[i]
            deps_pv = get(deps_p, v, deps_∅)
            comp_pv = get(comp_p, v, comp_∅)
            # dependencies
            for q in deps_pv
                add_queue!(q, 1)
                if q in saturated
                    # p@i conflicts with q being uninstallable
                    # q in saturated means that can happen
                    add_queue!(p, i+1)
                end
            end
            # conflicts
            for q in intx
                info_q = pkgs[q]
                vers_q = info_q.versions
                comp_q = info_q.compat
                l = get(reach, q, 0)
                for (m, w) in enumerate(vers_q)
                    m ≤ l || break # only consider reachable
                    comp_qw = get(comp_q, w, comp_∅)
                    v ∈ keys(comp_p) &&
                    q ∈ keys(comp_pv) &&
                    w ∉ comp_pv[q] ||
                    w ∈ keys(comp_q) &&
                    p ∈ keys(comp_qw) &&
                    v ∉ comp_qw[p] || continue
                    # v & w have a conflict
                    add_queue!(p, i+1)
                    add_queue!(q, m+1)
                end
            end
        end
        # update the reach map
        reach[p] = k
    end

    # TODO: if reach[p] > length(vers_p) then if the highest reachable
    # version of a package depends on p, we need to add it's successor

    return reach
end

function filter_reachable!(
    pkgs  :: Dict{P, PkgInfo{P,V,S}},
    reach :: Dict{P, Int},
) where {P,V,S}
    for pkg in sort!(collect(keys(reach)))
        k = reach[pkg]
        info = pkgs[pkg]
        vers = info.versions
        for j = k+1:length(vers)
            ver = vers[j]
            delete!(info.depends, ver)
            delete!(info.compat, ver)
        end
        k < length(vers) && resize!(vers, k)
    end
    filter!(pkgs) do (pkg, info)
        pkg in keys(reach) && !isempty(info.versions)
    end
end

function filter_reachable!(
    pkgs :: Dict{P, PkgInfo{P,V,S}},
    reqs :: Vector{P},
) where {P,V,S}
    reach = find_reachable(pkgs, reqs)
    filter_reachable!(pkgs, reach)
end

# compute the set of conflicts between package versions

struct PkgEntry{P,V}
    versions  :: Vector{V}
    depends   :: Vector{P}
    interacts :: Dict{P, Int}
    conflicts :: BitMatrix
end

function PkgEntry(
    pkg    :: P,
    pkgs   :: Dict{P, PkgInfo{P,V,S}},
    intx   :: Vector{P},
    active :: Bool,
) where {P,V,S}
    # shorter name
    p = pkg
    # look up some stuff about p
    info_p = pkgs[p]
    vers_p = info_p.versions
    deps_p = info_p.depends
    comp_p = info_p.compat
    # collect all dependency packages
    dx = P[]
    for (v, deps) in deps_p
        union!(dx, deps)
    end
    sort!(dx)
    ix = Dict{P, Int}(p => 0 for p in intx)
    # compute interactions matrix (m × n)
    m = length(vers_p) # per package version
    n = length(dx) + # per dependency package + interacting package version
        sum(init=0, length(pkgs[q].versions) for q in intx)
    X = falses(m + active, n + active) # conflicts & actives
    # empty deps & compat
    deps_∅ = Vector{P}()
    comp_∅ = Dict{P,S}()
    # main work loop
    for (i, v) in enumerate(vers_p)
        deps_pv = get(deps_p, v, deps_∅)
        comp_pv = get(comp_p, v, comp_∅)
        # dependencies
        for (j, q) in enumerate(dx)
            X[i, j] = q ∈ deps_pv
        end
        # conflicts
        b = length(dx)
        for q in intx
            ix[q] = b
            info_q = pkgs[q]
            vers_q = info_q.versions
            comp_q = info_q.compat
            for (j, w) in enumerate(vers_q)
                comp_qw = get(comp_q, w, comp_∅)
                X[i, b + j] =
                    v ∈ keys(comp_p) &&
                    q ∈ keys(comp_pv) &&
                    w ∉ comp_pv[q] ||
                    w ∈ keys(comp_q) &&
                    p ∈ keys(comp_qw) &&
                    v ∉ comp_qw[p]
            end
            b += length(vers_q)
        end
    end
    # initially all versions are active
    if active
        X[1:m, end] .= true
        # columns can start inactive if they have no conflicts
        for j = 1:n
            X[end, j] = any(X[i, j] for i = 1:m)
        end
    end
    return PkgEntry(vers_p, dx, ix, X)
end

function pkg_entries(
    pkgs :: Dict{P,PkgInfo{P,V,S}},
    active :: Bool = false,
) where {P,V,S}
    ∅ = P[] # empty vector (reused)
    interacts = find_interacts(pkgs)
    Dict{P, PkgEntry{P,V}}(
        p => PkgEntry(p, pkgs, get(interacts, p, ∅), active)
        for p in keys(pkgs)
    )
end

# eliminate versions that can never be chosen

function filter_redundant!(
    pkgs :: Dict{P, PkgInfo{P,V,S}},
    data :: Dict{P, PkgEntry{P,V}} = pkg_entries(pkgs, true),
) where {P,V,S}
    work = copy(keys(data))
    names = sort!(collect(work))
    # some work vectors
    J = Int[] # active versions vector
    K = Int[] # active conflicts vector
    R = Int[] # redundant indices vector
    sort!(names, by = p -> length(data[p].interacts))
    for p in Iterators.cycle(names)
        isempty(work) && break
        p in work || continue
        delete!(work, p)
        # get conflicts & dimensions
        X = data[p].conflicts
        m = size(X, 1) - 1
        m > 1 || continue # unique version cannot be reundant
        n = size(X, 2) - 1
        # active indices
        append!(empty!(J), (j for j = 1:m if X[j, end]))
        length(J) > 1 || continue # unique version cannot be reundant
        append!(empty!(K), (k for k = 1:n if X[end, k]))
        # find redundant versions
        empty!(R)
        for j in J
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
        for q in keys(data[p].interacts)
            b = data[q].interacts[p]
            data[q].conflicts[end, b .+ R] .= false
            push!(work, q) # can create new redundancies
        end
    end
    for (p, info) in pkgs
        X = data[p].conflicts
        m = length(info.versions)
        # find all the redundant versions of p
        append!(empty!(R), (j for j = 1:m if !X[j, end]))
        for (i, v) in enumerate(info.versions)
            i in R || continue
            delete!(info.depends, v)
            delete!(info.compat, v)
        end
        deleteat!(info.versions, R)
    end
end

# primary entry-point for loading package data

function get_pkg_data(
    deps :: DepsProvider{P,V,S},
    reqs :: SetOrVec{P} = deps.packages,
) where {P,V,S}
    pkgs = load_packages(deps, reqs)
    filter_reachable!(pkgs, reqs)
    filter_redundant!(pkgs)
    pkg_entries(pkgs)
end

using ArgTools

const magic = "\xfapkg data v1\0"

function write_strs(io::IO, v::Vector)
    write(io, UInt16(length(v)))
    for x in v
        s = x isa String ? x : string(x)
        write(io, UInt8(ncodeunits(s)))
        write(io, s)
    end
end

function read_strs(io::IO)
    n = read(io, UInt16)
    v = Vector{String}(undef, n)
    for i = 1:n
        l = read(io, UInt8)
        v[i] = String(read(io, l))
    end
    return v
end

function write_strs(io::IO, inds::Dict{<:Any,<:Integer}, v::Vector)
    write(io, UInt16(length(v)))
    for x in v
        s = x isa String ? x : string(x)
        write(io, UInt16(inds[s]))
    end
end

function read_strs(io::IO, vals::Vector{String})
    n = read(io, UInt16)
    v = Vector{String}(undef, n)
    for i = 1:n
        j = read(io, UInt16)
        v[i] = vals[j]
    end
    return v
end

function write_bits(io::IO, X::BitMatrix)
    l = length(X)
    write(io, UInt32(l))
    bytes = @view(reinterpret(UInt8, X.chunks)[1:cld(l,8)])
    write(io, bytes)
end

function read_bits(io::IO, m::Int)
    l = read(io, UInt32)
    n, r = divrem(l, m)
    iszero(r) || error("unexpected bit matrix length: $l not divisible by $m")
    X = BitMatrix(undef, m, n)
    read!(io, @view(reinterpret(UInt8, X.chunks)[1:cld(l,8)]))
    return X
end

function save_pkg_data(
    data :: Dict{P, PkgEntry{P,V}},
    out  :: Union{ArgWrite, Nothing} = nothing,
) where {P,V}
    arg_write(out) do out
        write(out, magic)
        pv = sort!(collect(keys(data)))
        pm = Dict{P,Int}(p => i for (i, p) in enumerate(pv))
        write_strs(out, pv)
        for (p, d) in data
            write_strs(out, d.versions)
            write_strs(out, pm, d.depends)
            write_strs(out, pm, sort!(collect(keys(d.interacts))))
            write_bits(out, d.conflicts)
        end
    end
end

# variables for each package:
#  - one for the package
#  - one for each version
#
# clauses for each package:
#  - one/zero if the package is required or not
#  - one specifying package versions
#  - one for each conflict
#  - optimality clauses, one per version

function gen_sat(
    out  :: Union{ArgWrite, Nothing},
    data :: Dict{P, PkgEntry{P}},
    reqs :: Vector{P},
) where {P}
    arg_write(out) do out
        # compute & output header
        names = sort!(collect(keys(data)))
        var = Dict{P, Int}() # variable indices
        v = 0 # number of variables
        x = 0 # number of conflicts
        for p in names
            var[p] = v + 1
            v += length(data[p].versions) + 1
            x += sum(data[p].conflicts)
            x += sum(@view(data[p].conflicts[:, 1:length(data[p].depends)]))
        end
        # conflicts are double-counted
        @assert iseven(x)
        x >>= 1
        # number of clauses
        c = length(reqs) + v + x
        println(out, "p cnf $v $c")
        println(out)
        # output requirements clauses
        for p in sort(reqs)
            println(out, "$(var[p]) 0")
        end
        println(out)
        # output package version clauses
        for p in names
            print(out, "-$(var[p]) ")
            for i = 1:length(data[p].versions)
                print(out, "$(var[p]+i) ")
            end
            println(out, "0")
        end
        println(out)
        # output dependency clauses
        for p in names
            for i = 1:length(data[p].versions)
                for (j, q) in enumerate(data[p].depends)
                    data[p].conflicts[i, j] || continue
                    println(out, "-$(var[p]+i) $(var[q]) 0")
                end
            end
        end
        println(out)
        # output incompatibility clauses
        for p in names
            for q in sort!(collect(keys(data[p].interacts)))
                q < p || break
                b = data[p].interacts[q]
                for i = 1:length(data[p].versions),
                    j = 1:length(data[q].versions)
                    data[p].conflicts[i, b+j] || continue
                    println(out, "-$(var[p]+i) -$(var[q]+j) 0")
                end
            end
        end
        println(out)
        # output optimality clauses
        for p in names
            for i = 1:length(data[p].versions)
                print(out, "-$(var[p]) $(var[p]+i) ")
                # can be "excused" if a better version is chosen
                for j = 1:i-1
                    print(out, "$(var[p]+j) ")
                end
                # or if it conflicts with something else chosen
                for q in sort!(collect(keys(data[p].interacts)))
                    b = data[p].interacts[q]
                    for j = 1:length(data[q].versions)
                        data[p].conflicts[i, b+j] || continue
                        print(out, "$(var[q]+j) ")
                    end
                end
                println(out, "0")
            end
        end
    end
end

function gen_sat(
    data :: Dict{P, PkgEntry{P}},
    reqs :: Vector{P},
) where {P}
    gen_sat(nothing, data, reqs)
end
