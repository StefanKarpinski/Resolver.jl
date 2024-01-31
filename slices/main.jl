include("../test/registry.jl")
include("Slices.jl"); using .Slices

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

# select top most popular packages
top = 1024
vertices₀ = select_popular_packages(best, top)
@assert length(vertices₀) ≥ top

# pare down the SAT problem to only popular packages and deps
Resolver.filter_pkg_info!(info, map(first, vertices₀))
sat = Resolver.SAT(info)

colors₁ = color_sat_dsatur(vertices₀, sat)

# find all packages needed for optimal solutions
vertices₁ = mapreduce(union!, colors₁) do color
    with_temp_clauses(sat) do
        fixed = vertices₀[color]
        for (p, v) in fixed
            sat_add(sat, p, v)
            sat_add(sat)
        end
        sol = only(resolve_core(sat, first.(fixed); max=1, by=popularity))
        sort!(collect(sol))
    end
end
sort!(vertices₁)
sort!(vertices₁, by = popularity)
@assert first.(vertices₀) ==
    unique(first.(vertices₁))[1:length(vertices₀)]

colors₂ = color_sat_dsatur(vertices₁, sat)

# we have complete colors, grow to slices
vertices = vertices₁
colors = colors₂

# expand slices
slices = expand_slices(vertices, sat, colors)

for slice in slices
    with_temp_clauses(sat) do
        fixed = vertices₁[slice]
        for (p, v) in fixed
            sat_add(sat, p, v)
            sat_add(sat)
        end
        sol = only(resolve_core(sat, first.(fixed), max=1, by=popularity))
        @assert sol ⊆ vertices # fails
        @assert sol ∩ vertices ⊆ fixed # also fails
    end
end
