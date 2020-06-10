module Resolver

export graph

const Node = Tuple{String,Int}

function parse_pkg(pkg::AbstractString)::String
    occursin(r"^[a-z]+$"i, pkg) ? pkg : error("invalid package name: $pkg")
end

function parse_ver(ver::AbstractString)::Node
    m = match(r"^([a-z]+)([1-9][0-9]*)$"i, ver)
    m === nothing && error("invalid package version: $ver")
    String(m.captures[1]), parse(Int, m.captures[2])
end

function parse_list(f::Function, list::AbstractString)
    map(f, split(strip(list), r"\s*,\s*", keepempty=false))
end
parse_pkgs(list::AbstractString) = parse_list(parse_pkg, list)
parse_vers(list::AbstractString) = parse_list(parse_ver, list)

function add_node!(nodes::Vector{Node}, node::Node)
    node in nodes && return
    not = (node[1], 0)
    not in nodes || push!(nodes, not)
    push!(nodes, node)
    return
end

function add_edge!(
    nodes::Vector{Node},
    edges::Vector{Tuple{Node,Node}},
    edge::Tuple{Node,Node},
)
    edge in edges && return
    add_node!(nodes, edge[1])
    add_node!(nodes, edge[2])
    push!(edges, edge)
    push!(edges, reverse(edge))
    return
end

function add_edge!(
    nodes::Vector{Node},
    edges::Vector{Tuple{Node,Node}},
    edge::Pair{Node,String},
)
    add_edge!(nodes, edges, (edge[1], (edge[2], 0)))
end

function graph(
    dependencies::Dict{String,String},
    conflicts::Vector{Tuple{String,String}} = Tuple{String,String}[],
)
    # construct the graph
    nodes = Node[]
    edges = Tuple{Node,Node}[]
    for (vers, deps) in dependencies
        for ver in parse_vers(vers),
            pkg in parse_pkgs(deps)
            add_edge!(nodes, edges, ver => pkg)
        end
    end
    for (vers1, vers2) in conflicts
        for v1 in parse_vers(vers1),
            v2 in parse_vers(vers2)
            add_edge!(nodes, edges, (v1, v2))
        end
    end
    sort!(nodes)
    sort!(edges)

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
