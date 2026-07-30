#!/usr/bin/env julia
#
# Regression tests for the `bin/` resolver tooling (bin/Registries.jl).
#
# Run from the repo root with the bin/ environment:
#
#     julia --project=bin bin/test/runtests.jl
#
# These are integration tests: they need the General registry installed and,
# on a cold cache, network access to fetch the list of released Julia versions
# (bin/julia_versions.json). They are intentionally *not* part of the package
# test suite (test/runtests.jl), which is offline and exercises only src/.

using Test
import Base: UUID
import Pkg.Registry: JULIA_UUID, PkgEntry, reachable_registries
import Pkg.Versions: VersionSpec
import Resolver

include(joinpath(@__DIR__, "..", "Registries.jl"))

# UUIDs of the packages involved in issue #24.
const COMPILER_SUPPORT_LIBRARIES_JLL = UUID("e66e0078-7015-5450-92f7-15fbd957f2ae")
const LINEAR_ALGEBRA = UUID("37e2e46d-f89d-539d-b4ee-838fcccc9c8e")
const JSON = UUID("682c06a0-de6a-54ab-a142-c8b1cf79cde6")
const STATISTICS = UUID("10745b16-79ce-11e8-11f9-7d13ad32a3b2")

# Load packages from the installed registries (mirrors bin/resolve.jl).
const packages = Dict{UUID,Vector{PkgEntry}}()
for reg in reachable_registries()
    for (uuid, entry) in reg.pkgs
        push!(get!(()->PkgEntry[], packages, uuid), entry)
    end
end

make_provider(julia::VersionSpec) = registry_provider(
    packages;
    julia_versions = julia,
    allow_pre = Dict(UUID(0) => false),
)

# Is the requirement set resolvable for the given Julia spec?
function resolves(reqs::Vector{UUID}; julia::VersionSpec)
    reg = make_provider(julia)
    info = Resolver.pkg_info(reg, reqs)
    Resolver.resolve(info, reqs) !== nothing
end

@testset "bin/Registries.jl" begin
    # Issue #24: Julia 1.10.8 bundles CompilerSupportLibraries_jll 1.1.1, but the
    # General registry's entry for 1.1.x declares `julia = "1.11.0-1"`. The
    # provider must widen 1.1.1's `julia` compat to admit the Julia that bundles
    # it -- while *keeping* the registry bound on versions that aren't bundled by
    # 1.10.8, so those stay correctly incompatible (a bundled stdlib version can
    # also be a plain resolvable dependency on another Julia, where the registry
    # bound must still govern). Julia 1.10.8 is a frozen release, so its bundled
    # CompilerSupportLibraries_jll (1.1.1) never changes.
    @testset "stdlib julia-compat is widened, not cleared" begin
        reg = make_provider(VersionSpec("1.10.8"))
        pd = Resolver.pkg_data(reg, [COMPILER_SUPPORT_LIBRARIES_JLL])[COMPILER_SUPPORT_LIBRARIES_JLL]
        compatible(v) = (spec = get(pd.compat[v], JULIA_UUID, nothing);
                         spec === nothing || v"1.10.8" ∈ spec)
        # the bundled version (1.1.1) is admitted on 1.10.8 (this is the bug fix)
        bundled = [v for v in keys(pd.compat) if (v.major, v.minor, v.patch) == (1, 1, 1)]
        @test !isempty(bundled)
        @test all(compatible, bundled)
        # but the registry bound survives elsewhere: some later version, neither
        # bundled by 1.10.8 nor registry-compatible with it, must stay excluded
        # (guards against clearing the bound outright)
        @test any(!compatible(v) for v in keys(pd.compat))
    end

    # End-to-end: the minimal reproducers from issue #24 must resolve.
    @testset "BLAS stack resolves on Julia 1.10" begin
        # Smallest reproducer: depend directly on the culprit stdlib jll.
        @test resolves([COMPILER_SUPPORT_LIBRARIES_JLL, JULIA_UUID]; julia = VersionSpec("1.10.8"))
        @test resolves([COMPILER_SUPPORT_LIBRARIES_JLL, JULIA_UUID]; julia = VersionSpec("1.10"))
        # Realistic reproducer: LinearAlgebra pulls in the same stack transitively.
        @test resolves([LINEAR_ALGEBRA, JULIA_UUID]; julia = VersionSpec("1.10.8"))
    end

    # A project's compat bounds are user constraints, so the provider no longer
    # applies them: they reach the resolver as part of a `Problem`, which
    # forbids the excluded versions by clause rather than deleting them.
    @testset "project compat travels in the Problem" begin
        reg = make_provider(VersionSpec("1.10"))
        prob = Resolver.Problem([JSON, JULIA_UUID];
            compat = Dict(JSON => VersionSpec("0.20")))
        info = Resolver.pkg_info(reg, prob)
        sol = Resolver.resolve(info, prob)
        @test sol !== nothing
        @test sol[JSON] ∈ VersionSpec("0.20")
        # the provider still offers the newer versions and the filter keeps
        # them: they are constrained away, not deleted
        @test info[JSON].versions[1] ∉ VersionSpec("0.20")
    end

    # The provider offers a bundled stdlib version whatever the registries say,
    # and it used to apply project compat to the registry versions only --
    # before patching the bundled ones back in. `bundled_versions` is what lets
    # bin/resolve.jl widen its Problem's bounds to reproduce that exactly.
    # Julia 1.10.8 is frozen, so what it bundles never changes.
    @testset "bundled stdlib versions" begin
        bundled = bundled_versions(VersionSpec("1.10.8"))
        @test haskey(bundled, STATISTICS)
        @test v"1.10.0" in bundled[STATISTICS]
        @test !haskey(bundled, JSON) # not a stdlib
    end
end
