module Resolver

export graph, solutions, resolve

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
    node′ = (node[1], 0)
    node′ in nodes || push!(nodes, node′)
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
    requires::Vector{String},
    dependencies::Dict{String,String},
    conflicts::Vector{Tuple{String,String}} = Tuple{String,String}[],
)
    # construct the graph
    nodes = Node[]
    edges = Tuple{Node,Node}[]
    for (vers, deps) in dependencies
        for ver in parse_vers(vers)
            add_node!(nodes, ver)
            for pkg in parse_pkgs(deps)
                add_edge!(nodes, edges, ver => pkg)
            end
        end
    end
    for (vers1, vers2) in conflicts
        for v1 in parse_vers(vers1),
            v2 in parse_vers(vers2)
            add_edge!(nodes, edges, (v1, v2))
        end
    end
    by(v) = (v[1], v[2] == 0 && v[1] ∉ requires ? typemin(Int) : -v[2])
    sort!(nodes; by)
    sort!(edges)

    # construct the neighbor lists & graph matrix
    n = length(nodes)
    N = [[i] for i = 1:n]
    G = zeros(Bool, n, n)
    for (src, dst) in edges
        i = findfirst(==(src), nodes)
        j = findfirst(==(dst), nodes)
        j in N[i] || push!(N[i], j)
        i in N[j] || push!(N[j], i)
        G[i,j] = true
    end
    foreach(sort!, N)

    return nodes, N, G
end

function graph(
    dependencies::Dict{String,String},
    conflicts::Vector{Tuple{String,String}} = Tuple{String,String}[],
)
    graph(String[], dependencies, conflicts)
end

const \ = setdiff

function solutions(
    V::Vector{Node},
    N::Vector{Vector{Int}},
)
    S = Vector{Int}[]

    # find dominating independent kernels
    function BronKerbosch(R::Vector{Int}, P::Vector{Int}, X::Vector{Int})
        # @show R, P, X
        if isempty(P) && isempty(X)
            push!(S, R)
            return true
        end
        found = false
        while !isempty(P)
            v = popfirst!(P)
            R′, P′, X′ = R ∪ [v], P \ N[v], X \ N[v]
            at_least_one = all(
                any(V[w][1] == p for w in R′) ||
                any(V[w][1] == p for w in P′)
            for p in unique(V[w][1] for w in N[v]))
            if at_least_one && BronKerbosch(R′, P′, X′)
                # don't consider sub-optimal solutions
                # filter!(P) do w
                #     V[w][1] != V[v][1] || V[w][2] == 0
                # end
                found = true
            end
            # require at least one version of each package
            p = V[v][1]
            any(V[w][1] == p for w in R) ||
            any(V[w][1] == p for w in P) || break
            push!(X, v)
        end
        return found
    end
    BronKerbosch(Int[], collect(1:length(V)), Int[])

    # expand each kernel maximal independent set
    # function expand!(s::Vector{Int})
    #     nodes = collect(1:length(V))
    #     for v in s
    #         setdiff!(nodes, N[v])
    #     end
    #     while !isempty(nodes)
    #         v = pop!(nodes)
    #         push!(s, v)
    #         setdiff!(nodes, N[v])
    #     end
    #     return nodes
    # end
    # foreach(expand!, S)

    # sort and take best versions
    for soln in S
        sort!(soln)
        unique!(v -> V[v][1], soln)
    end

    [[V[j] for j in J] for J in S]
end

function resolve(
    requires::AbstractString,
    dependencies::Dict{String,String},
    conflicts::Vector{Tuple{String,String}} = Tuple{String,String}[],
)
    V, N = graph(parse_pkgs(requires), dependencies, conflicts)
    solutions(V, N)
end

end # module
