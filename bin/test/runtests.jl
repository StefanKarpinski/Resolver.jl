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
# LLVM_full_jll has prerelease versions in the General registry (11.0.0-rc5 &c),
# and Wine_jll has *only* prerelease versions -- so it must not resolve at all
# unless prereleases are admitted.
const LLVM_FULL_JLL = UUID("a3ccf953-465e-511d-b87f-60a6490c289d")
const WINE_JLL = UUID("9fae3aff-8997-5dd1-9b84-5d0cc5e0bffa")
# Compat 4.0.0 is yanked from General (the julia-downgrade-compat case, see the
# "yanked packages" testset in the package suite), and libsodium_jll's *newest*
# registered version is too.
const COMPAT = UUID("34da2185-b29b-5c13-b0c7-acf172513d20")
const LIBSODIUM_JLL = UUID("a9144af2-ca23-56d9-984f-0d03f7b5ccf8")

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
)

# `--allow-pre`'s dictionary: the zero uuid holds the default
no_pre = Dict(UUID(0) => false)
all_pre = Dict(UUID(0) => true)

# The old baked path, for the knobs that used to be applied by deleting versions
# from the provider's output: `prob`'s excluded versions gone from the data, and
# with them their deps & compat entries (which `pkg_info` reads by version).
function bake(data::AbstractDict{P}, prob::Resolver.Problem{P}) where {P}
    baked = empty(data)
    for (p, d) in data
        vers = [v for v in d.versions if !Resolver.is_excluded(prob, p, v)]
        baked[p] = Resolver.PkgData(vers,
            typeof(d.depends)(v => d.depends[v] for v in vers if haskey(d.depends, v)),
            typeof(d.compat)(v => d.compat[v] for v in vers if haskey(d.compat, v)))
    end
    return baked
end

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

    # The version preference ordering is a resolve parameter, not part of the
    # package universe: the provider always hands over the canonical order
    # (newest first), and a query's ordering is `resolve`'s `order` argument. The
    # oracle is the old path, where the ordering was baked into the provider's
    # version vectors -- so sorting the data that way and resolving canonically
    # must give exactly what passing the comparator gives, at registry scale.
    @testset "the version ordering is a resolve parameter" begin
        reg = make_provider(VersionSpec("1.10"))
        reqs = sort([JSON, COMPILER_SUPPORT_LIBRARIES_JLL, JULIA_UUID])
        prob = Resolver.Problem(reqs)
        data = Resolver.pkg_data(reg, reqs)
        # the provider's own order is canonical, so nothing needs permuting
        info = Resolver.pkg_info(data, reqs)
        @test Resolver.version_permutations(info, _ -> >) === nothing
        # oldest first, then a couple of the real --max-minor / --min-minor
        # comparators the options build
        minor(u, v) = Base.thisminor(u) ≠ Base.thisminor(v) ?
            Base.thisminor(u) > Base.thisminor(v) : u < v
        rminor(u, v) = Base.thisminor(u) ≠ Base.thisminor(v) ?
            Base.thisminor(u) < Base.thisminor(v) : u > v
        for lt in (<, minor, rminor)
            order = _ -> lt
            baked = Dict(u => Resolver.PkgData(
                            sort(collect(d.versions); lt), d.depends, d.compat)
                         for (u, d) in data)
            @test Resolver.resolve(data, prob; order) ==
                  Resolver.resolve(baked, prob)
            # ... and the one T1 artifact answers under any of them
            @test Resolver.resolve(info, prob; order) ==
                  Resolver.resolve(baked, prob)
        end
        # a per-package ordering: only julia reversed, everything else canonical
        order = u -> u == JULIA_UUID ? (<) : (>)
        baked = Dict(u => Resolver.PkgData(
                        sort(collect(d.versions); lt = order(u)),
                        d.depends, d.compat) for (u, d) in data)
        @test Resolver.resolve(info, prob; order) == Resolver.resolve(baked, prob)
        @test Resolver.resolve(info, prob; order)[JULIA_UUID] <
              Resolver.resolve(info, prob)[JULIA_UUID]
    end

    # Prerelease admission is a query constraint, not a property of the package
    # universe: the versions exist either way, and `--allow-pre` says which of
    # them a query accepts. The oracle is the old path, where the provider
    # deleted the prereleases it was not told to keep -- so resolving against the
    # data with them dropped must give exactly what constraining them gives, at
    # registry scale and for every admission setting.
    @testset "prerelease admission is a resolve-time constraint" begin
        reg = make_provider(VersionSpec("1.10"))
        reqs = sort([LLVM_FULL_JLL, WINE_JLL, JSON, JULIA_UUID])
        data = Resolver.pkg_data(reg, reqs)
        # the provider really does offer prereleases now
        @test any(!isempty(v.prerelease) for v in data[LLVM_FULL_JLL].versions)
        # and the one artifact keeps them, for every query to constrain as it likes
        info = Resolver.pkg_info(reg, reqs)
        @test any(!isempty(v.prerelease) for v in info[LLVM_FULL_JLL].versions)

        for allow_pre in (no_pre, all_pre,
                          Dict(UUID(0) => false, LLVM_FULL_JLL => true))
            excludes = [prerelease_exclusion(allow_pre)]
            prob = Resolver.Problem(reqs; excludes)
            old = bake(data, prob) # the versions deleted, as the provider used to
            @test Resolver.resolve(data, prob) == Resolver.resolve(old, reqs)
            @test Resolver.resolve(info, prob) == Resolver.resolve(old, reqs)
        end

        # and it is not inert: every Wine_jll in the registry is a prerelease, so
        # requiring it is unsatisfiable unless prereleases are admitted -- per
        # package or globally, the flag's two shapes
        wine = ["Wine_jll" => WINE_JLL]
        @test isnothing(resolve_versions("", ["--julia=1.10"]; deps = wine))
        pre = resolve_versions("", ["--julia=1.10", "--allow-pre=Wine_jll"];
                               deps = wine)
        @test !isnothing(pre)
        @test !isempty(pre[WINE_JLL].prerelease)
        @test pre == resolve_versions("", ["--julia=1.10", "--allow-pre"];
                                      deps = wine)
    end

    # Yankedness is the same shape of fact as prerelease-ness -- a property of a
    # version that a query may or may not accept -- so it is modeled the same
    # way, as an exclusion kind, even though nothing yet offers the user a way to
    # turn it off. The oracle is the old path, where the provider deleted yanked
    # versions outright.
    @testset "yanked versions are constrained, not deleted" begin
        reg = make_provider(VersionSpec("1.10"))
        reqs = sort([COMPAT, LIBSODIUM_JLL, JULIA_UUID])
        data = Resolver.pkg_data(reg, reqs)
        yanked = yanked_exclusion(packages)
        forbids = last(yanked)

        # `is_yanked` reads the registries, not the version number: Compat 4.0.0
        # is yanked from General and libsodium_jll's newest version is too
        @test forbids(COMPAT, v"4.0.0")
        @test !forbids(COMPAT, v"4.1.0")
        @test forbids(LIBSODIUM_JLL, v"1.0.22+0")
        @test !forbids(LIBSODIUM_JLL, v"1.0.21+0")
        @test !forbids(JSON, v"0.21.4")
        @test !forbids(JULIA_UUID, v"1.10.0") # not a registered package at all
        @test !forbids(COMPAT, v"99.0.0")     # no such version

        # the provider offers them and the one artifact keeps them
        @test v"4.0.0" in data[COMPAT].versions
        info = Resolver.pkg_info(reg, reqs)
        @test v"4.0.0" in info[COMPAT].versions

        prob = Resolver.Problem(reqs; excludes = [yanked])
        old = bake(data, prob) # the versions deleted, as the provider used to
        @test Resolver.resolve(data, prob) == Resolver.resolve(old, reqs)
        @test Resolver.resolve(info, prob) == Resolver.resolve(old, reqs)
        # ... together with the prerelease kind, and under a reversed ordering
        both = Resolver.Problem(reqs;
            excludes = [yanked, prerelease_exclusion(no_pre)])
        @test Resolver.resolve(info, both) ==
              Resolver.resolve(bake(data, both), reqs)
        order = u -> (<)
        @test Resolver.resolve(info, both; order) ==
              Resolver.resolve(bake(data, both), reqs; order)

        # The selector map names the kind, so a diagnostic could one day say
        # "the only version that satisfies this is yanked". libsodium_jll's
        # yanked version is its *newest*, so nothing dominates it and the filter
        # has to keep it -- it is ruled out by clause. Compat's is not, so
        # redundancy elimination deletes it, which is the filter subsuming the
        # deletion the provider used to do.
        work = Resolver.prepare_pkg_info(info, prob)
        sat = Resolver.SAT(work, prob)
        try
            @test (:yanked, LIBSODIUM_JLL) in values(sat.sels)
            @test v"1.0.22+0" in work[LIBSODIUM_JLL].versions
            @test v"4.0.0" ∉ work[COMPAT].versions
        finally
            Resolver.finalize(sat)
        end

        # end to end: libsodium_jll's newest registered version is yanked, so the
        # resolve must land on the one below it
        sol = resolve_versions("", ["--julia=1.10"];
                               deps = ["libsodium_jll" => LIBSODIUM_JLL])
        @test !isnothing(sol)
        @test sol[LIBSODIUM_JLL] == v"1.0.21+0"
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
