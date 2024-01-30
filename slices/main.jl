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

G₀ = explicit_conflicts(packages₀, info)
colors₀ = color_sat_dsatur(packages₀, sat, G₀)
slices₀ = expand_slices(packages₀, sat, colors₀)

G₁ = implicit_conflicts(packages₀, sat, slices₀)
colors₁ = color_sat_dsatur(packages₀, sat, G₁)

# find all packages needed for optimal solutions
packages₁ = mapreduce(union!, colors₁) do color
    with_temp_clauses(sat) do
        reqs = packages₀[color]
        for (p, v) in reqs
            sat_add(sat, p, v)
            sat_add(sat)
        end
        mapreduce(collect, union!, resolve_core(sat, first.(reqs)))
    end
end
sort!(packages₁)
sort!(packages₁, by = popularity)
@assert first.(packages₀) ==
    unique(first.(packages₁))[1:length(packages₀)]

G₂ = explicit_conflicts(packages₁, info)
colors₂ = color_sat_dsatur(packages₁, sat, G₂)
slices₂ = expand_slices(packages₁, sat, colors₂)

G₃ = implicit_conflicts(packages₁, sat, slices₂)
colors₃ = color_sat_dsatur(packages₁, sat, G₃)

packages = packages₁
colors = colors₃

vers = Dict{String,Vector{Int}}()
for color in colors
    with_temp_clauses(sat) do
        reqs = packages[color]
        for (p, v) in reqs
            sat_add(sat, p, v)
            sat_add(sat)
        end
        sols = resolve_core(sat, first.(reqs))
        for sol in sols, (p, v) in sol
            V = get!(()->Int[], vers, p)
            v in V || push!(V, v)
        end
    end
end
foreach(sort!, values(vers))

#=
# turn slices into solutions
sols = [slice_dict(packages, best, slice) for slice in slices]
for sol in sols
    with_temp_clauses(sat) do
        for (p, v) in sol
            sat_add(sat, p, v)
            sat_add(sat)
        end
        for p in packages
            p in keys(sol) && continue
            sat_assume(sat, p)
            is_satisfiable(sat) || continue

            # require some version of p
            sat_add(sat, p)
            sat_add(sat)
            extract_solution!(sat, sol)

            # optimize p's version
            best[p] < sol[p] &&
            optimize_solution!(sat, sol) do
                for i = 1:sol[p]-1
                    sat_add(sat, p, i)
                end
                sat_add(sat)
            end

            # fix p's version
            sat_add(sat, p, sol[p])
            sat_add(sat)
        end
    end
end

vers = Dict{String,Vector{Int}}()
for sol in sols
    for (p, v) in sol
        V = get!(()->Int[], vers, p)
        v in V || push!(V, v)
    end
end
map(sort!, values(vers))
=#
