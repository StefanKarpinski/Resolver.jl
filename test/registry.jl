using Pkg: depots1
using Pkg.Registry: RegistryInstance, init_package_info!
using Pkg.Types: stdlibs
using Pkg.Versions: VersionSpec

const excludes = push!(Set(first.(values(stdlibs()))), "julia")
const reg_path = joinpath(depots1(), "registries", "General.toml")
const reg_inst = RegistryInstance(reg_path)
const reg_dict = Dict(p.name => p for p in values(reg_inst.pkg_entries) if p.name ∉ excludes)

dp = DepsProvider{String, VersionNumber, VersionSpec}(keys(reg_dict)) do pkg::String
    info = init_package_info!(reg_dict[pkg])
    vers = sort!(collect(keys(info.version_info)), rev=true)
    deps = Dict(v => String[] for v in vers)
    comp = Dict(v => Dict{String,VersionSpec}() for v in vers)
    # scan versions and populate deps & compat data
    for v in vers
        for (r, d) in info.deps
            v in r && union!(deps[v], keys(d))
        end
        for (r, c) in info.compat
            v in r && mergewith!(intersect, comp[v], c)
        end
    end
    foreach(sort!, values(deps))
    # scrub out excluded deps (stdlibs, julia itself)
    for d in values(deps)
        setdiff!(d, excludes)
    end
    for c in values(comp), x in excludes
        delete!(c, x)
    end
    # deduplicate data structures to save memory
    for i = 1:length(vers)-1, j = i+1:length(vers)
        v, w = vers[i], vers[j]
        deps[v] == deps[w] && (deps[v] = deps[w])
        comp[v] == comp[w] && (comp[v] = comp[w])
    end
    # return resolver PkgData structure
    PkgData{String, VersionNumber, VersionSpec}(vers, deps, comp)
end

using GraphModularDecomposition
using GraphModularDecomposition.StrongModuleTrees

reqs = String.(split("Druid,HTTP", ','))
data = get_pkg_data(dp, reqs)
nodes, G = to_graph(data)
parts = let d = Dict{String,Vector{Int}}()
    for (i, (p, v)) in enumerate(nodes)
        push!(get!(()->valtype(d)(), d, p), i)
    end
    sort!(collect(values(d)))
end
L = [v == v"0-" ? "!$p" : "$p/$v" for (p, v) in nodes]

H = complete_graph(G, parts)
D = H .!= G
# E = sort!([(i, j) for (i, j) in zip(findnz(D)...) if i < j])

S = sort!(StrongModuleTree(G))
T = sort!(StrongModuleTree(H))

# get a prime node
P = S[1]
v = map(first_leaf, P)
c = Tuple{Int,Int}[]
for i in v, j in v
    i < j || continue
    r = setdiff(v, (i, j))
    x = G[r, i] .!= G[r, j]
    y = D[r, i] .| D[r, j]
    z = count(@. x & !y)
    # println("$i, $j => $z")
    z == 0 && println("$i, $j")
    push!(c, (i, j))
end

x = map(first_leaf, S[1])
y = setdiff(1:length(nodes), leaves(S[1]))
G1 = G[x,x] # prime module quotient subgraph



# filter_reachable!(all_pkg_entries, all_names)
# filter_redundant!(all_pkg_entries)

#=
# uninstallable pair with the least total versions
reqs = ["ClassicalOrthogonalPolynomials", "PoincareInvariants"]
pkg_entries = deepcopy(all_pkg_entries)
filter_reachable!(pkg_entries, reqs)
filter_redundant!(pkg_entries)
solve(pkg_entries, reqs)
=#

#=
for line in eachline("test/pkg_pairs.csv")
    p, q = split(line, ',')
    reqs = String[p, q]
    println(repr(reqs))
    pkg_entries = deepcopy(all_pkg_entries)
    filter_reachable!(pkg_entries, reqs)
    filter_redundant!(pkg_entries)
    !solve(pkg_entries, reqs) && break
end
=#

#=
uninstallable = String[]
all_pkg_entries = load_pkg_entries(dp, all_names)
filter_reachable!(all_pkg_entries, all_names)
filter_redundant!(all_pkg_entries)

for pkg in all_names
    println("[[[ $pkg ]]]")
    reqs = [pkg]
    pkg_entries = deepcopy(all_pkg_entries)
    filter_reachable!(pkg_entries, reqs)
    filter_redundant!(pkg_entries)
    solve(pkg_entries, reqs) && continue
    push!(uninstallable, pkg)
end
=#

#=
all_pkg_entries = load_pkg_entries(dp, all_names)
filter_reachable!(all_pkg_entries, all_names)
filter_redundant!(all_pkg_entries)
all_ix = find_interacts(all_pkg_entries)

const pairs = Tuple{String,String}[]
for p in all_names, q in get(all_ix, p, String[])
    p < q || continue
    reqs = [p, q]
    pkg_entries = load_pkg_entries(dp, reqs)
    filter_reachable!(pkg_entries, reqs)
    filter_redundant!(pkg_entries)
    all(length(pkg_entries[p].versions) > 1 for p in reqs) || continue
    push!(pairs, (p, q))
    @show reqs
end
=#
