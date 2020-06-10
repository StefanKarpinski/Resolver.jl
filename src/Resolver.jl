module Resolver

function parse_pkg(pkg::AbstractString)
    occursin(r"^[a-z]+$"i, pkg) ? pkg : error("invalid package name: $pkg")
end

function parse_deps(deps::AbstractString)
    map(parse_pkg, split(strip(deps), r"\s*,\s*", keepempty=false))
end

function parse_ver(ver::AbstractString)
    m = match(r"^([a-z]+)([1-9][0-9]*)$"i, ver)
    m === nothing && error("invalid package version: $ver")
    String(m.captures[1]), parse(Int, m.captures[2])
end

function graph(
    deps::Dict{String,String},
    conflicts::Vector{Tuple{String,String}},
)
    nodes = sort!(unique!(map(parse_ver, collect(keys(deps)))))
    edges = Tuple{Tuple{String,Int},Tuple{String,Int}}[
        (parse_ver(a), parse_ver(b)) for (a, b) in conflicts
    ]
    for (v, d) in deps
        ver = parse_ver(v)
        for dep in parse_deps(d)
            push!(edges, (ver, (dep, 0)))
        end
    end
    sort!(union!(edges, map(reverse, edges)))
    sort!(union!(nodes, map(first, edges)))
    # construct the neighbor lists
    n = length(nodes)
    N = [[i] for i = 1:n]
    for (src, dst) in edges
        i = findfirst(==(src), nodes)
        j = findfirst(==(dst), nodes)
        j in N[i] || push!(N[i], j)
        i in N[j] || push!(N[j], i)
    end
    foreach(sort!, N)
    return nodes, N
end

end # module
