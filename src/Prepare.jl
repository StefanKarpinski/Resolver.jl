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

function Base.:(==)(a::PkgEntry, b::PkgEntry)
    a.versions  == b.versions  &&
    a.depends   == b.depends   &&
    a.interacts == b.interacts &&
    a.conflicts == b.conflicts
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

# write LEB128 integer value
function write_int(io::IO, n::Integer)
    @assert n ≥ 0
    more = true
    while more
        b = (n % UInt8) & 0x7f
        n >>= 7
        more = !iszero(n)
        b |= UInt8(more) << 7
        write(io, b)
    end
end

# read LEB128 integer value
function read_int(io::IO, ::Type{T} = Int) where {T<:Integer}
    n::T = s = 0
    while true
        b = read(io, UInt8)
        n |= T(b & 0x7f) << s
        (b % Int8) ≥ 0 && break
        s += 7
    end
    return n
end

function write_vals(io::IO, v::Vector)
    write_int(io, length(v))
    for x in v
        s = x isa String ? x : string(x)
        write_int(io, ncodeunits(s))
        write(io, s)
    end
end

function read_vals(io::IO, ::Type{T}) where {T}
    n = read_int(io, Int)
    v = Vector{T}(undef, n)
    for i = 1:n
        l = read_int(io)
        s = String(read(io, l))
        v[i] = T === String ? s : parse(T, s)
    end
    return v
end

function write_vals(io::IO, inds::Dict{<:Any,<:Integer}, v::Vector)
    write_int(io, length(v))
    for x in v
        s = x isa String ? x : string(x)
        write_int(io, inds[s])
    end
end

function read_vals(io::IO, vals::Vector{T}) where {T}
    n = read_int(io)
    v = Vector{T}(undef, n)
    for i = 1:n
        j = read_int(io)
        v[i] = vals[j]
    end
    return v
end

function write_bits(io::IO, X::BitMatrix)
    l = length(X)
    write_int(io, l)
    bytes = @view(reinterpret(UInt8, X.chunks)[1:cld(l,8)])
    write(io, bytes)
end

function read_bits(io::IO, m::Int)
    l = read_int(io)
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
        write_vals(out, pv)
        for p in pv
            d = data[p]
            write_vals(out, d.versions)
            write_vals(out, pm, d.depends)
            write_vals(out, pm, sort!(collect(keys(d.interacts))))
            write_bits(out, d.conflicts)
        end
    end
end

function load_pkg_data(
    in :: ArgRead,
    :: Type{P} = String,
    :: Type{V} = VersionNumber,
) where {P,V}
    arg_read(in) do in
        m = String(read(in, ncodeunits(magic)))
        m == magic || error("""
            Unexpected package data format, magic doesn't match:
                $(repr(m)) != $(repr(magic))
            """)
        pv = read_vals(in, P)
        data = Dict{P, PkgEntry{P,V}}()
        for p in pv
            vers = read_vals(in, V)
            deps = read_vals(in, pv)
            intx = read_vals(in, pv)
            conf = read_bits(in, length(vers))
            ix = Dict{P, Int}(q => 0 for q in intx)
            data[p] = PkgEntry(vers, deps, ix, conf)
        end
        # compute interacts dict values
        for (p, d) in data
            b = length(d.depends)
            for q in sort!(collect(keys(d.interacts)))
                d.interacts[q] = b
                b += length(data[q].versions)
            end
        end
        return data
    end
end
