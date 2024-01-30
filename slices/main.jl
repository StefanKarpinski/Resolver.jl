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
packages₀ = select_popular_packages(best, top)
@assert length(packages₀) ≥ top

# pare down the SAT problem to only popular packages and deps
Resolver.filter_pkg_info!(info, map(first, packages₀))
sat = Resolver.SAT(info)

colors₁ = color_sat_dsatur(packages₀, sat)

# find all packages needed for optimal solutions
packages₁ = mapreduce(union!, colors₁) do color
    with_temp_clauses(sat) do
        fixed = packages₀[color]
        for (p, v) in fixed
            sat_add(sat, p, v)
            sat_add(sat)
        end
        sols = resolve_core(sat, first.(fixed), max=0)
        mapreduce(collect, union!, sols)
    end
end
sort!(packages₁)
sort!(packages₁, by = popularity)
@assert first.(packages₀) ==
    unique(first.(packages₁))[1:length(packages₀)]

colors₂ = color_sat_dsatur(packages₁, sat)

# we have complete colors, grow to slices
packages = packages₁
colors = colors₂

# expand slices
slices = expand_slices(packages, sat, colors)

for slice in slices
    with_temp_clauses(sat) do
        fixed = packages₁[slice]
        for (p, v) in fixed
            sat_add(sat, p, v)
            sat_add(sat)
        end
        sols = resolve_core(sat, first.(fixed))
        @assert length(sols) == 1
        sol = only(sols)
        @assert sol ⊆ packages # fails
        @assert sol ∩ packages ⊆ fixed # also fails
    end
end
