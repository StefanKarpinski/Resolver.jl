module Resolver

export graph, parse_deps, resolve

const Version = Tuple{String,Int}

function parse_pkg(pkg::AbstractString)::String
    occursin(r"^[a-z]+$"i, pkg) ? pkg : error("invalid package name: $pkg")
end

function parse_ver(ver::AbstractString)::Version
    m = match(r"^([a-z]+)([1-9][0-9]*)$"i, ver)
    m === nothing && error("invalid package version: $ver")
    String(m.captures[1]), parse(Int, m.captures[2])
end

function parse_list(f::Function, list::AbstractString)
    map(f, split(strip(list), r"\s*,\s*", keepempty=false))
end
parse_pkgs(list::AbstractString) = parse_list(parse_pkg, list)
parse_vers(list::AbstractString) = parse_list(parse_ver, list)

function parse_deps(
    requires::Vector{String},
    dependencies::Dict{String,String},
    conflicts::Vector{Tuple{String,String}} = Tuple{String,String}[],
)
    packages = String[]
    versions = Dict{String,Vector{Int}}()
    excludes = Set{Tuple{Version,Version}}()

    function add_version!((pkg, ver)::Version)
        if pkg ∉ packages
            push!(packages, pkg)
            versions[pkg] = [0]
        end
        ver in versions[pkg] || push!(versions[pkg], ver)
    end

    function add_conflict!(v1::Version, v2::Version)
        add_version!(v1)
        add_version!(v2)
        (v1, v2) in excludes || push!(excludes, (v1, v2))
        (v2, v1) in excludes || push!(excludes, (v2, v1))
    end

    function add_dependency!(ver::Version, pkg::String)
        add_conflict!(ver, (pkg, 0))
    end

    for (vers, deps) in dependencies
        for ver in parse_vers(vers)
            add_version!(ver)
            for pkg in parse_pkgs(deps)
                add_dependency!(ver, pkg)
            end
        end
    end

    for (vers1, vers2) in conflicts
        for ver1 in parse_vers(vers1),
            ver2 in parse_vers(vers2)
            add_conflict!(ver1, ver2)
        end
    end

    sort!(packages, by = pkg -> (pkg ∉ requires, pkg))
    for (pkg, vers) in versions
        if pkg in requires
            sort!(vers, by = v -> v ≠ 0 ? v : typemax(v))
        else
            sort!(vers)
        end
    end

    return packages, versions, excludes
end

function parse_deps(
    requires::String,
    dependencies::Dict{String,String},
    conflicts::Vector{Tuple{String,String}} = Tuple{String,String}[],
)
    parse_deps(parse_pkgs(requires), dependencies, conflicts)
end

function resolve(
    requires::String,
    dependencies::Dict{String,String},
    conflicts::Vector{Tuple{String,String}} = Tuple{String,String}[],
)
    pkgs, vers, excl = parse_deps(requires, dependencies, conflicts)

    # TODO: resolve
end

end # module
