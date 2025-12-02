# Minimal example demonstrating the yanked version issue in Resolver.jl
# See https://github.com/StefanKarpinski/Resolver.jl/issues/15
#
# This script demonstrates that Resolver.jl may select yanked versions
# of packages, which should be avoided since yanked versions are typically
# broken or have been superseded.

using Pkg
using Pkg.Registry: RegistryInstance, init_package_info!, isyanked, reachable_registries
using Pkg.Versions: VersionSpec

using Resolver: DepsProvider, PkgData, resolve

# Get registry instance
const REG = reachable_registries()[1]

# Helper function to find packages with yanked versions
function find_yanked_packages()
    yanked_pkgs = Dict{String, Vector{VersionNumber}}()
    for (uuid, entry) in REG.pkgs
        info = init_package_info!(entry)
        for v in keys(info.version_info)
            if isyanked(info, v)
                if !haskey(yanked_pkgs, entry.name)
                    yanked_pkgs[entry.name] = VersionNumber[]
                end
                push!(yanked_pkgs[entry.name], v)
            end
        end
    end
    return yanked_pkgs
end

# Create a DepsProvider that does NOT filter out yanked versions
# (this mimics the current behavior of Resolver.jl)
function provider_without_yank_filter(pkg_names::Vector{String})
    reg_dict = Dict(p.name => p for p in values(REG.pkgs) if p.name in pkg_names)

    DepsProvider(keys(reg_dict)) do pkg::AbstractString
        info = init_package_info!(reg_dict[pkg])
        # Note: We do NOT filter out yanked versions here - this is the bug
        vers = sort!(collect(keys(info.version_info)), rev=true)
        deps = Dict(v => String[] for v in vers)
        comp = Dict(v => Dict{String, VersionSpec}() for v in vers)

        for v in vers
            for (r, d) in info.deps
                v in r && union!(deps[v], keys(d))
            end
            for (r, c) in info.compat
                v in r && mergewith!(∩, comp[v], c)
            end
        end
        foreach(sort!, values(deps))

        # Filter deps to only include packages we're considering
        for d in values(deps)
            filter!(in(pkg_names), d)
        end
        for c in values(comp)
            for k in collect(keys(c))
                k in pkg_names || delete!(c, k)
            end
        end

        PkgData(vers, deps, comp)
    end
end

# Create a minimal test case using Compat which has yanked version 4.0.0
function demonstrate_yanked_issue()
    println("=" ^ 60)
    println("Demonstrating yanked version issue in Resolver.jl")
    println("=" ^ 60)
    println()

    # First, show which versions of Compat are yanked
    println("Step 1: Identify yanked versions of Compat")
    println("-" ^ 40)

    for (uuid, entry) in REG.pkgs
        entry.name == "Compat" || continue
        info = init_package_info!(entry)
        vers = sort(collect(keys(info.version_info)))

        yanked_vers = VersionNumber[]
        for v in vers
            if isyanked(info, v)
                push!(yanked_vers, v)
            end
        end

        println("Yanked versions of Compat: $yanked_vers")
        println()

        # Find version just before and after 4.0.0
        v400_idx = findfirst(==(v"4.0.0"), vers)
        if v400_idx !== nothing
            println("Versions around 4.0.0:")
            for i in max(1, v400_idx-2):min(length(vers), v400_idx+2)
                yanked_str = isyanked(info, vers[i]) ? " [YANKED]" : ""
                println("  $(vers[i])$yanked_str")
            end
        end
    end
    println()

    # Step 2: Create a scenario where resolver might pick a yanked version
    println("Step 2: Test resolution with compat constraint requiring yanked version")
    println("-" ^ 40)

    # Create a minimal package setup where someone depends on Compat ~4.0
    # This should ideally resolve to 4.1.0+ but NOT 4.0.0 since it's yanked

    # For this test, we'll use a synthetic package "TestPkg" that depends on Compat
    # with a compat bound of "~4" (meaning 4.0.0 - 4.x.x)

    pkg_names = ["Compat"]
    deps = provider_without_yank_filter(pkg_names)

    # Get data and inspect versions
    data = Dict{String, PkgData{String,VersionNumber,VersionSpec}}()
    for pkg in pkg_names
        data[pkg] = deps.provider(pkg)
        println("Available versions for $pkg (first 10): ", data[pkg].versions[1:min(10, end)])
    end
    println()

    # Resolve - Resolver will pick the "best" (newest) versions
    println("Step 3: Running resolve...")
    println("-" ^ 40)

    pkgs, vers = resolve(data, ["Compat"])

    println("Resolution result:")
    for (i, pkg) in enumerate(pkgs)
        solutions = [vers[i, j] for j in 1:size(vers, 2) if vers[i, j] !== nothing]
        println("  $pkg => $solutions")

        # Check if any selected version is yanked
        for (uuid, entry) in REG.pkgs
            entry.name == pkg || continue
            info = init_package_info!(entry)
            for sol in solutions
                if isyanked(info, sol)
                    println("  ⚠️  WARNING: Selected version $sol is YANKED!")
                end
            end
        end
    end
    println()

    # Step 4: Create a more specific test case
    println("Step 4: Force resolution to pick yanked version via compat")
    println("-" ^ 40)

    # Create a synthetic scenario: a package "TestPkg" with compat Compat = "4.0"
    # This forces the resolver to pick exactly v4.0.0, which is yanked

    test_data = Dict{String, PkgData{String,VersionNumber,VersionSpec}}()

    # TestPkg v1.0.0 depends on Compat with compat "4.0" (only 4.0.x)
    test_data["TestPkg"] = PkgData(
        [v"1.0.0"],                              # versions
        Dict(v"1.0.0" => ["Compat"]),           # depends
        Dict(v"1.0.0" => Dict("Compat" => VersionSpec("4.0")))  # compat
    )

    # Add Compat data (without yank filtering)
    test_data["Compat"] = data["Compat"]

    pkgs2, vers2 = resolve(test_data, ["TestPkg"])

    println("Resolution result for TestPkg depending on Compat ~4.0:")
    for (i, pkg) in enumerate(pkgs2)
        solutions = [vers2[i, j] for j in 1:size(vers2, 2) if vers2[i, j] !== nothing]
        println("  $pkg => $solutions")

        if pkg == "Compat"
            for (uuid, entry) in REG.pkgs
                entry.name == "Compat" || continue
                info = init_package_info!(entry)
                for sol in solutions
                    if isyanked(info, sol)
                        println("  ** PROBLEM: Resolver selected yanked version Compat@", sol, "!")
                        println("     This demonstrates the bug - yanked versions should be excluded.")
                    end
                end
            end
        end
    end
    println()

    println("=" ^ 60)
    println("CONCLUSION:")
    println("Resolver.jl currently does not filter out yanked versions.")
    println("The fix should add yank filtering in the DepsProvider creation,")
    println("similar to how it filters prereleases.")
    println("=" ^ 60)
end

# Run the demonstration
demonstrate_yanked_issue()
