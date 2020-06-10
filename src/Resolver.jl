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
    node in nodes || push!(nodes, node)
    not = (node[1], 0)
    not in nodes || push!(nodes, not)
end

function add_edge!(
    nodes::Vector{Node},
    edges::Vector{Tuple{Node,Node}},
    (a, b)::Tuple{Node,Node},
)
    (a, b) in edges && return
    add_node!(nodes, a)
    add_node!(nodes, b)
    push!(edges, (a, b))
    push!(edges, (b, a))
end

function add_edge!(
    nodes::Vector{Node},
    edges::Vector{Tuple{Node,Node}},
    (v, p)::Pair{Node,String},
)
    add_edge!(nodes, edges, (v, (p, 0)))
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
    sort!(nodes, by = v -> (v[1], -v[2]))
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

function resolve(
    V::Vector{Node},
    N::Vector{Vector{Int}},
)
    pkgs = parse_pkgs(pkg_list)
    # by(v::Node) = (v[1], v[2] == 0 && v[1] ∉ pkgs ? typemin(Int) : -v[2])
    # p = sortperm(V; by)
    # q = invperm(p)
    # S = solutions(V[p], [q[I] for I in N[p]])
    S = solutions(V, N)
end

const \ = setdiff

function solutions(V::Vector{Node}, N::Vector{Vector{Int}})
    S = Vector{Int}[]
    function BronKerbosch(R::Vector{Int}, P::Vector{Int}, X::Vector{Int})
        if isempty(P) && isempty(X)
            push!(S, R)
            return true
        end
        found = false
        for v in P
            if BronKerbosch(R ∪ [v], P \ N[v], X \ N[v])
                filter!(w -> w[1] == v[1], P)
                found = true
            else
                filter!(!=(v), P)
            end
            push!(X, v)
        end
        return found
    end
    BronKerbosch(Int[], collect(1:length(V)), Int[])
    return S
end

end # module
