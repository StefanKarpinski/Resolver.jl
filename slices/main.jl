include("../test/registry.jl")
include("Slices.jl"); using .Slices

using SparseArrays
using Resolver:
    Resolver,
    resolve_core,
    sat_assume,
    sat_add,
    is_satisfiable,
    extract_solution!,
    optimize_solution!,
    with_temp_clauses

## top-level code

@info "Loading pkg info..."
info = Resolver.pkg_info(registry.provider())
sat = Resolver.SAT(info)
best = compute_best_versions(sat)

# select best versions of most popular packages
top = 1000
vertices = select_popular_packages(best, top)
@assert length(vertices) ≥ top

# pare down the SAT problem to only popular packages and deps
Resolver.filter_pkg_info!(info, map(first, vertices))
sat = Resolver.SAT(info)

G = conflicts(vertices, sat)
while true
    global vertices, colors, slices, G
    @show length(vertices)
    # compute colors and expand to slices
    colors = color_sat_dsatur(vertices, sat, G)
    slices = expand_slices(vertices, sat, colors)
    # if slices are self-contained, we're done:
    vertices′ = slices_support(vertices, sat, slices)
    vertices′ == vertices && break
    # expand old conflict graph to new vertices
    n = length(vertices′)
    P = Vector{Int}(indexin(vertices, vertices′))
    I, J = findnz(G)
    G = sparse(P[I], P[J], true, n, n)
    # combine with explicit conflicts
    G .|= explicit_conflicts(vertices′, info)
    # expand the conflict graph...
    C = color_sat_dsatur(vertices′, sat, G)
    S = expand_slices(vertices′, sat, C)
    G .|= implicit_conflicts(vertices′, sat, S, P)
    # replace global vertices
    vertices = vertices′
end

# check self-containedness of slices
for slice in slices
    sol = solve_slice(vertices, sat, slice)
    @assert sol ⊆ vertices[slice]
end
