#!/usr/bin/env julia
#
# Regression tests for the `bin/` resolver tooling (bin/Registries.jl, plus
# option handling exercised through bin/resolve.jl in a subprocess).
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
const DELIMITED_FILES = UUID("8bb1440f-4735-579b-a4ab-409b98df4dab")

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

# Resolve a project with the given dependencies, `[compat]` body and command
# line flags, returning uuid => version -- or `nothing` when the requirements
# are unsatisfiable. This drives the real `bin/resolve.jl` in a subprocess,
# since that is where the project file is read, where the option precedence
# lives, and where the resolved versions are reported. `--print-versions`
# (rather than manifest generation) so it works on any host Julia, prerelease
# included.
const RESOLVE_JL = normpath(joinpath(@__DIR__, "..", "resolve.jl"))
const BIN_PROJECT = normpath(joinpath(@__DIR__, ".."))

function resolve_versions(
    compat :: AbstractString,
    flags  :: Vector{String} = String[];
    deps   :: Vector{Pair{String,UUID}} = ["JSON" => JSON],
)
    dir = mktempdir()
    open(joinpath(dir, "Project.toml"), "w") do io
        println(io, "[deps]")
        for (name, uuid) in deps
            println(io, "$name = \"$uuid\"")
        end
        isempty(compat) && return
        println(io, "\n[compat]")
        println(io, compat)
    end
    out = IOBuffer()
    err = IOBuffer()
    julia = Base.julia_cmd()[1]
    cmd = `$julia --project=$BIN_PROJECT $RESOLVE_JL $dir --print-versions $flags`
    if !success(pipeline(cmd; stdout = out, stderr = err))
        # tell "no solution" apart from "the script broke"
        occursin("Unsatisfiable", String(take!(err))) ||
            error("failed: $cmd")
        return nothing
    end
    vers = Dict{UUID,VersionNumber}()
    for line in eachline(seekstart(out))
        m = match(r"^(\S{36})\s+\S+\s+(\S+)", line)
        isnothing(m) && continue
        vers[UUID(m[1])] = VersionNumber(m[2])
    end
    return vers
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
    # so `bundled_versions` answers "which versions of this stdlib can I get on
    # these Julias at all" -- the version sets the couplings below are stated
    # in terms of. Julia 1.10.8 is frozen, so what it bundles never changes.
    @testset "bundled stdlib versions" begin
        bundled = bundled_versions(VersionSpec("1.10.8"))
        @test haskey(bundled, STATISTICS)
        @test v"1.10.0" in bundled[STATISTICS]
        @test !haskey(bundled, JSON) # not a stdlib
    end

    # The Julia versions to resolve for come from `--julia` if given, otherwise
    # from the project's own `[compat] julia` bound, otherwise from the `1`
    # default. Unlike every other compat entry, the `julia` bound selects a
    # version *universe* rather than constraining one: it decides which Julias
    # exist to be resolved among, and with them which stdlib versions are
    # bundled and pinned.
    @testset "julia compat as the default julia bound" begin
        # the newest release the `1` default admits: what resolving with
        # neither a flag nor a project bound has always picked
        newest = maximum(v for v in Registries.JULIA_VERSIONS
                         if isempty(v.prerelease) && v ∈ VersionSpec("1"))

        # neither a flag nor a bound: unchanged
        plain = resolve_versions("")
        @test plain[JULIA_UUID] == newest

        # the project's bound supplies the default, and here it changes the
        # answer -- which it silently failed to do before
        bound = resolve_versions("julia = \"~1.10\"")
        @test bound[JULIA_UUID] ∈ VersionSpec("1.10")
        @test bound[JULIA_UUID] ≠ plain[JULIA_UUID]

        # `--julia` overrides the bound outright rather than intersecting with
        # it: the answer is exactly the flag's answer on a project with no
        # bound at all
        flag = resolve_versions("julia = \"~1.10\"", ["--julia=1.9"])
        @test flag[JULIA_UUID] ∈ VersionSpec("1.9")
        @test flag == resolve_versions("", ["--julia=1.9"])
    end

    # Julia bundles a version of every stdlib, but it only *pins* some of them:
    # `Pkg.Types.UPGRADABLE_STDLIBS` are bundled unpinned, so their registry
    # versions compete like any other package's and a user bound on one is
    # enforced strictly. We resolve over Julia versions too, so "pinned" here
    # means pinned per candidate Julia -- each candidate carries its own pins,
    # and an upgradable stdlib carries none from any of them.
    @testset "upgradable stdlibs resolve like packages" begin
        @test STATISTICS ∈ UPGRADABLE_STDLIBS_UUIDS
        @test DELIMITED_FILES ∈ UPGRADABLE_STDLIBS_UUIDS
        @test LINEAR_ALGEBRA ∉ UPGRADABLE_STDLIBS_UUIDS

        stat = ["Statistics" => STATISTICS]

        # Julia 1.10 bundles Statistics 1.10.0, but the registry's 1.11.0 and
        # 1.11.1 declare `julia = "1.9.0 - 1"` / `"1.9.4 - 1"`, so they are
        # installable there: a bound demanding 1.11 gets a registry version
        # instead of being pinned to the bundled one. (Also covers the
        # reporting path, which used to substitute the bundled version for any
        # resolved stdlib.)
        up = resolve_versions("Statistics = \"1.11\"", ["--julia=1.10"]; deps = stat)
        @test !isnothing(up)
        @test up[STATISTICS] ∈ VersionSpec("1.11")
        @test up[JULIA_UUID] ∈ VersionSpec("1.10")

        # no Statistics is 2.x, bundled or registered, so the bound is
        # unsatisfiable rather than silently widened away
        @test isnothing(
            resolve_versions("Statistics = \"2\"", ["--julia=1.10"]; deps = stat))

        # the bundled version remains a candidate: on Julia 1.8 every
        # registered Statistics is ruled out by its own `julia` compat, so the
        # answer is whatever 1.8 bundles
        old = resolve_versions("", ["--julia=1.8"]; deps = stat)
        @test !isnothing(old)
        @test old[STATISTICS] ∈ bundled_versions(VersionSpec("1.8"))[STATISTICS]
    end

    # A bound on a stdlib is not inert. A Julia version is compatible with
    # exactly the stdlib versions it bundles -- its compat pins say so -- and a
    # version that exists only as a bundled stdlib is installable only on the
    # Julias that ship it. So a bound on a stdlib constrains the Julia choice
    # through those couplings, and the resolver steers Julia to satisfy it or
    # reports the requirements unsatisfiable.
    @testset "stdlib bounds steer the julia choice" begin
        linalg = ["LinearAlgebra" => LINEAR_ALGEBRA]
        stat = ["Statistics" => STATISTICS]

        # a bound on a pinned stdlib pushes Julia forward: only a Julia that
        # bundles a LinearAlgebra ≥ 1.11 can satisfy it
        steered = resolve_versions("LinearAlgebra = \"1.11\""; deps = linalg)
        @test !isnothing(steered)
        @test steered[LINEAR_ALGEBRA] ≥ v"1.11"
        @test steered[JULIA_UUID] ≥ v"1.11"

        # ... and when the Julia universe is restricted to versions that bundle
        # no such LinearAlgebra, the bound is unsatisfiable rather than inert
        @test isnothing(resolve_versions("LinearAlgebra = \"1.11\"",
                                        ["--julia=1.10"]; deps = linalg))

        # no Julia bundles a LinearAlgebra 99.x, and there is no registry copy
        # of it either
        @test isnothing(resolve_versions("LinearAlgebra = \"99\""; deps = linalg))

        # a bundled-only version is admitted just by the Julias that bundle it,
        # never by any other: Statistics 1.10.0 exists only as Julia 1.10's
        # bundled copy (the registry starts at 1.11.0), so bounding Statistics
        # to the 1.10 series steers Julia to 1.10 as well -- it must not pair
        # 1.10.0 with a newer Julia that ships something else
        tilde = resolve_versions("Statistics = \"~1.10\""; deps = stat)
        @test !isnothing(tilde)
        @test tilde[STATISTICS] ∈ bundled_versions(VersionSpec("1.10"))[STATISTICS]
        @test tilde[JULIA_UUID] ∈ VersionSpec("1.10")
    end
end
