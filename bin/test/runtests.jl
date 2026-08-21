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
# (bin/julia_versions.json) and, for the manifest test, to download the packages
# it resolves. They are intentionally *not* part of the package test suite
# (test/runtests.jl), which is offline and exercises only src/.

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
# Issue #54: Python_jll's Compat.toml bounds LibMPDec_jll over `3 - 3.10`, but
# its Deps.toml only lists LibMPDec_jll from 3.10 on -- so the 3.8.x versions
# carry a compat entry for a package that is not one of their dependencies.
const PYTHON_JLL = UUID("93d3a430-8e7c-50da-8e8d-3dfcfb3baf05")
const LIBMPDEC_JLL = UUID("7106de7a-f406-5ef1-84f7-3345f7341bd2")

# MbedTLS_jll is a stdlib of Julia 1.10 but not of 1.11 and later, so it tells
# the target Julia's stdlib set apart from the host's.
const MBEDTLS_JLL = UUID("c8ffd9c3-330d-5841-b78e-0817d7145fa1")

# Load packages from the installed registries (mirrors bin/resolve.jl).
const packages = Dict{UUID,Vector{PkgEntry}}()
for reg in reachable_registries()
    for (uuid, entry) in reg.pkgs
        push!(get!(()->PkgEntry[], packages, uuid), entry)
    end
end

# The one package universe: the only parameter `registry_provider` has left is
# `allow_yanked`, and the default universe -- the one every query but
# `--allow-yanked` resolves against -- is built without the versions the
# registries have struck.
const yanked_dropped = Dict{UUID,Set{VersionNumber}}()
const reg = registry_provider(packages; yanked = yanked_dropped)

# `--allow-pre`'s dictionary: the zero uuid holds the default
no_pre = Dict(UUID(0) => false)
all_pre = Dict(UUID(0) => true)

# The `Problem` bin/resolve.jl builds: the Julia bound as an ordinary compat entry
# on `julia`, the project's other bounds beside it, and the admission kinds.
function make_problem(
    reqs      :: Vector{UUID};
    julia     :: VersionSpec,
    compat    :: AbstractDict = Dict{UUID,VersionSpec}(),
    allow_pre :: Dict{UUID,Bool} = no_pre,
)
    c = Dict{UUID,Any}(compat)
    c[JULIA_UUID] = julia
    Resolver.Problem(reqs;
        compat = c, prerelease = prerelease_exclusion(allow_pre))
end

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
    prob = make_problem(reqs; julia)
    info = Resolver.pkg_info(reg, prob)
    verdict = Resolver.resolve(info, prob; diagnose = false) !== nothing
    # this is exactly what `issatisfiable` answers on its own, so it must agree
    # with the descent's verdict here -- on real registry data, with the Julia
    # bound and the admission kinds in force
    @test Resolver.issatisfiable(info, prob) == verdict
    return verdict
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
    err    :: IOBuffer = IOBuffer(), # the script's stderr, for tests that read it
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
    julia = Base.julia_cmd()[1]
    cmd = `$julia --project=$BIN_PROJECT $RESOLVE_JL $dir --print-versions $flags`
    if !success(pipeline(cmd; stdout = out, stderr = err))
        # tell "no solution" apart from "the script broke", putting the message
        # back so a caller that passed `err` in can read it too
        msg = String(take!(err))
        print(err, msg)
        # when it broke, say what it said -- the command alone explains nothing
        occursin("Unsatisfiable", msg) || error("failed: $cmd\n$msg")
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

# Resolve a project with the given dependencies for the given Julia, returning
# the generated manifest -- or the script's stderr when it failed. Unlike
# `resolve_versions` this drives manifest generation, which is a different code
# path in the script: it builds a Pkg `Context`, downloads the resolved packages
# and reads their project files for extension info.
function resolve_manifest(
    julia_version :: AbstractString;
    deps          :: Vector{Pair{String,UUID}},
)
    dir = mktempdir()
    open(joinpath(dir, "Project.toml"), "w") do io
        println(io, "[deps]")
        for (name, uuid) in deps
            println(io, "$name = \"$uuid\"")
        end
    end
    out, err = IOBuffer(), IOBuffer()
    julia = Base.julia_cmd()[1]
    cmd = `$julia --project=$BIN_PROJECT $RESOLVE_JL $dir --print-manifest --julia=$julia_version`
    success(pipeline(cmd; stdout = out, stderr = err)) || return String(take!(err))
    return String(take!(out))
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

    # Issue #54: a registry compat entry may name a package that is not a
    # dependency of the version the entry covers -- Python_jll bounds
    # LibMPDec_jll over `3 - 3.10` while only listing it as a dependency from
    # 3.10 on, so every 3.8.x build names a non-dependency. Pkg treats such a
    # bound as inert for that version, and so must we.
    @testset "compat naming a non-dependency is inert" begin
        pd = Resolver.pkg_data(reg, [PYTHON_JLL])[PYTHON_JLL]
        # the versions the malformed entry covers are still offered ...
        old = [v for v in pd.versions if Base.thispatch(v) < v"3.10"]
        @test !isempty(old)
        # ... with neither a dependency on LibMPDec_jll nor a bound on it
        @test all(v -> LIBMPDEC_JLL ∉ pd.depends[v], old)
        @test all(v -> !haskey(pd.compat[v], LIBMPDEC_JLL), old)
        # while the versions that *do* depend on it keep the bound
        new = [v for v in pd.versions if LIBMPDEC_JLL in pd.depends[v]]
        @test !isempty(new)
        @test all(v -> haskey(pd.compat[v], LIBMPDEC_JLL), new)
    end

    # The two Pkg representations reach that entry by different routes: the
    # <=1.13 one keys compat by name and resolves through the version's own deps
    # map, where the name is absent (a `KeyError` before the fix); the >=1.14 one
    # keys compat by uuid, resolved through a name table Pkg builds from every
    # range at once, where it is present and the entry silently survives. So the
    # drop has to be explicit on both, or the same registry yields different
    # package data depending on which Julia parsed it. Whichever Julia runs this
    # suite exercises only one of the two, so compare them head to head.
    @testset "both registry representations agree" begin
        dep = UUID("00000000-0000-0000-0000-0000000000d1")
        non = UUID("00000000-0000-0000-0000-0000000000e2")
        # the version's dependencies, in each key space
        names = Dict("Dep" => dep, "julia" => JULIA_UUID)   # <=1.13
        uuids = [dep, JULIA_UUID]                            # >=1.14
        # the same compat entry both ways: a bound on a dependency, a bound on
        # julia, and a bound on a package this version does not depend on
        by_name = Dict("Dep" => VersionSpec("1"),
                       "julia" => VersionSpec("1.6.0-1"),
                       "NotADep" => VersionSpec("2"))
        by_uuid = Dict(dep => VersionSpec("1"),
                       JULIA_UUID => VersionSpec("1.6.0-1"),
                       non => VersionSpec("2"))
        from_names = Dict(Registries.compat_uuid_pairs(by_name, names, uuids))
        from_uuids = Dict(Registries.compat_uuid_pairs(by_uuid, names, uuids))
        @test from_names == from_uuids
        # and what they agree on is: the non-dependency dropped, the rest kept
        @test from_names == Dict(dep => VersionSpec("1"),
                                 JULIA_UUID => VersionSpec("1.6.0-1"))
    end

    # A project's compat bounds are user constraints, so the provider no longer
    # applies them: they reach the resolver as part of a `Problem`, which
    # forbids the excluded versions by clause rather than deleting them.
    @testset "project compat travels in the Problem" begin
        prob = make_problem([JSON, JULIA_UUID]; julia = VersionSpec("1.10"),
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
        reqs = sort([JSON, COMPILER_SUPPORT_LIBRARIES_JLL, JULIA_UUID])
        prob = make_problem(reqs; julia = VersionSpec("1.10"))
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
        reqs = sort([LLVM_FULL_JLL, WINE_JLL, JSON, JULIA_UUID])
        data = Resolver.pkg_data(reg, reqs)
        # the provider really does offer prereleases now
        @test any(!isempty(v.prerelease) for v in data[LLVM_FULL_JLL].versions)
        # and the one artifact keeps them, for every query to constrain as it likes
        info = Resolver.pkg_info(reg, reqs)
        @test any(!isempty(v.prerelease) for v in info[LLVM_FULL_JLL].versions)

        for allow_pre in (no_pre, all_pre,
                          Dict(UUID(0) => false, LLVM_FULL_JLL => true))
            prob = Resolver.Problem(reqs;
                prerelease = prerelease_exclusion(allow_pre))
            old = bake(data, prob) # the versions deleted, as the provider used to
            @test Resolver.resolve(data, prob; diagnose = false) ==
                  Resolver.resolve(old, reqs; diagnose = false)
            @test Resolver.resolve(info, prob; diagnose = false) ==
                  Resolver.resolve(old, reqs; diagnose = false)
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

    # Yankedness is *not* the same shape of fact as prerelease-ness: a yanked
    # version is one the registry has struck, so it is not on offer at all. The
    # provider filters those versions out, which is what makes "the resolver
    # never produces a yanked version and never suggests one" true by
    # construction rather than by a constraint every query must remember to
    # carry. `--allow-yanked` is the one query that asks for a different
    # universe.
    @testset "yanked versions are filtered out of the universe" begin
        reqs = sort([COMPAT, LIBSODIUM_JLL, JULIA_UUID])

        # `yanked_versions` reads the registries, not the version number: Compat
        # 4.0.0 is yanked from General and libsodium_jll's newest version is too
        @test v"4.0.0" in yanked_versions(packages, COMPAT)
        @test v"4.1.0" ∉ yanked_versions(packages, COMPAT)
        @test v"99.0.0" ∉ yanked_versions(packages, COMPAT) # no such version
        @test v"1.0.22+0" in yanked_versions(packages, LIBSODIUM_JLL)
        @test v"1.0.21+0" ∉ yanked_versions(packages, LIBSODIUM_JLL)
        @test v"0.21.4" ∉ yanked_versions(packages, JSON)
        # not a registered package at all, so nothing of it is struck
        @test isempty(yanked_versions(packages, JULIA_UUID))

        # The multi-registry rule: a version yanked by *any* entry that carries
        # it is yanked, because one registry cannot undo another's yank -- which
        # is often a security action, so a registry that still lists the version
        # must not launder it back into the universe. No two installed registries
        # carry the same package with different yank flags (General is the only
        # one carrying any of these), so the rule cannot be exercised end to end
        # on real data; test it where it is stated, over the per-entry sets.
        none = Set{VersionNumber}()
        @test Registries.struck_versions([Set([v"1.0.0"]), none]) == Set([v"1.0.0"])
        @test Registries.struck_versions([none, Set([v"1.0.0"])]) == Set([v"1.0.0"])
        @test Registries.struck_versions([Set([v"1.0.0"]), Set([v"1.0.0"])]) ==
              Set([v"1.0.0"])
        @test isempty(Registries.struck_versions([none, none]))
        # each entry contributes its own, and only versions some entry yanks
        # appear at all -- which is why a version no registry has is never struck
        @test Registries.struck_versions([Set([v"1.0.0"]), Set([v"2.0.0"])]) ==
              Set([v"1.0.0", v"2.0.0"])
        @test isempty(Registries.struck_versions(Set{VersionNumber}[]))

        # the provider does not offer them, so the one artifact does not have
        # them either -- and the versions beside them are untouched
        data = Resolver.pkg_data(reg, reqs)
        info = Resolver.pkg_info(reg, reqs)
        @test v"4.0.0" ∉ data[COMPAT].versions
        @test v"4.1.0" in data[COMPAT].versions
        @test v"4.0.0" ∉ info[COMPAT].versions
        @test v"1.0.22+0" ∉ data[LIBSODIUM_JLL].versions
        @test v"1.0.21+0" in data[LIBSODIUM_JLL].versions
        # ... and the provider says what it dropped, which is what lets
        # bin/resolve.jl name the versions a bound admits that no longer exist
        @test v"4.0.0" in yanked_dropped[COMPAT]
        @test v"1.0.22+0" in yanked_dropped[LIBSODIUM_JLL]
        @test !haskey(yanked_dropped, JULIA_UUID)

        # `--allow-yanked` keeps them: same universe otherwise
        kept = registry_provider(packages; allow_yanked = Dict(UUID(0) => true))
        kept_data = Resolver.pkg_data(kept, reqs)
        @test v"4.0.0" in kept_data[COMPAT].versions
        @test v"1.0.22+0" in kept_data[LIBSODIUM_JLL].versions
        # per package: only the named package gets its struck versions back
        one = Resolver.pkg_data(
            registry_provider(packages;
                allow_yanked = Dict(UUID(0) => false, COMPAT => true)), reqs)
        @test v"4.0.0" in one[COMPAT].versions
        @test v"1.0.22+0" ∉ one[LIBSODIUM_JLL].versions

        # The mechanism this replaced: the unfiltered universe plus an exclusion
        # kind forbidding exactly the struck versions. Filtering must give the
        # same answers -- that is the whole claim of the change -- for the plain
        # problem, alongside the prerelease kind, and under a reversed ordering.
        struck = Dict{UUID,Set{VersionNumber}}() # per package, memoized
        yanked_of(u::UUID) = get!(() -> yanked_versions(packages, u), struck, u)
        # a kind reaches every package in the universe, not just the required
        # ones, so this has to answer for any uuid it is handed
        yanked = (u::UUID, v::VersionNumber) -> v in yanked_of(u)
        for pre in ((;), (; prerelease = prerelease_exclusion(no_pre)))
            old = Resolver.Problem(reqs; yanked, pre...)
            new = Resolver.Problem(reqs; pre...)
            @test Resolver.resolve(reg, new; diagnose = false) ==
                  Resolver.resolve(kept, old; diagnose = false)
            @test Resolver.resolve(info, new; diagnose = false) ==
                  Resolver.resolve(kept, old; diagnose = false)
            order = u -> (<)
            @test Resolver.resolve(info, new; order, diagnose = false) ==
                  Resolver.resolve(kept, old; order, diagnose = false)
        end

        # end to end: libsodium_jll's newest registered version is yanked, so
        # the default resolve lands on the one below it
        libsodium = ["libsodium_jll" => LIBSODIUM_JLL]
        sol = resolve_versions("", ["--julia=1.10"]; deps = libsodium)
        @test !isnothing(sol)
        @test sol[LIBSODIUM_JLL] == v"1.0.21+0"

        # ... and `--allow-yanked` takes it, bare or naming the package
        yes = resolve_versions("", ["--julia=1.10", "--allow-yanked"];
                               deps = libsodium)
        @test !isnothing(yes)
        @test yes[LIBSODIUM_JLL] == v"1.0.22+0"
        @test yes == resolve_versions(
            "", ["--julia=1.10", "--allow-yanked=libsodium_jll"]; deps = libsodium)

        # ... while naming some *other* package leaves it filtered
        other = resolve_versions("", ["--julia=1.10", "--allow-yanked=$COMPAT"];
                                 deps = libsodium)
        @test !isnothing(other)
        @test other[LIBSODIUM_JLL] == v"1.0.21+0"

        # The one failure the filtered universe cannot explain for itself: a
        # bound whose only admissible versions are struck ones. They are plainly
        # there in the registry, so bare "Unsatisfiable" would read as a lie --
        # the note names them and names the flag that brings them back.
        compat_pkg = ["Compat" => COMPAT]
        err = IOBuffer()
        @test isnothing(resolve_versions("Compat = \"=4.0.0\"", ["--julia=1.10"];
                                         deps = compat_pkg, err))
        msg = String(take!(err))
        @test occursin("Unsatisfiable", msg)
        @test occursin("compat on Compat", msg)
        @test occursin("yanked (4.0.0)", msg)
        @test occursin("--allow-yanked", msg)

        # ... and the flag really does accept them
        sol = resolve_versions("Compat = \"=4.0.0\"",
            ["--julia=1.10", "--allow-yanked"]; deps = compat_pkg)
        @test !isnothing(sol)
        @test sol[COMPAT] == v"4.0.0"

        # no note when a bound admits versions that survived: then the yank is
        # not what makes the requirements unsatisfiable
        err = IOBuffer()
        @test isnothing(resolve_versions("Compat = \"99\"", ["--julia=1.10"];
                                         deps = compat_pkg, err))
        @test !occursin("--allow-yanked", String(take!(err)))
    end

    # HistoricalStdlibVersions gives one stdlib version number to entries with
    # different dependency sets across Julia versions -- the provider tells them
    # apart with a synthetic build number -- so a bundled entry has to be matched
    # to the Julias that bundle *it*, not to every Julia that bundles something
    # with the same version number. Julia 1.11 and 1.12 both bundle a Markdown
    # 1.11.0, and only 1.12's depends on JuliaSyntaxHighlighting.
    @testset "a bundled stdlib entry belongs to its own Julias" begin
        MARKDOWN = UUID("d6f4376e-aef5-505a-96c1-9c027394607a")
        JULIA_SYNTAX_HIGHLIGHTING = UUID("ac6e5ff7-fb65-4e79-a425-ec3bc9c03011")
        pd = Resolver.pkg_data(reg, [MARKDOWN])[MARKDOWN]
        # several entries share the version number, told apart by the build
        entries = [v for v in pd.versions if Base.thispatch(v) == v"1.11.0"]
        @test length(entries) ≥ 2
        with = [v for v in entries if JULIA_SYNTAX_HIGHLIGHTING in pd.depends[v]]
        without = setdiff(entries, with)
        @test !isempty(with) && !isempty(without)
        # each is admitted by the Julias that ship it, and by no others: the
        # entry Julia 1.12 bundles is the one with the extra dependency
        @test any(v -> v"1.12.0" ∈ pd.compat[v][JULIA_UUID], with)
        @test all(v -> v"1.12.0" ∉ pd.compat[v][JULIA_UUID], without)
        @test any(v -> v"1.11.0" ∈ pd.compat[v][JULIA_UUID], without)
        @test all(v -> v"1.11.0" ∉ pd.compat[v][JULIA_UUID], with)
        # ... so a resolve on 1.12 gets 1.12's Markdown, dependencies and all
        sol = resolve_versions("", ["--julia=1.12.0"];
                               deps = ["Markdown" => MARKDOWN])
        @test !isnothing(sol)
        @test sol[MARKDOWN] == v"1.11.0"
        @test haskey(sol, JULIA_SYNTAX_HIGHLIGHTING)
        sol = resolve_versions("", ["--julia=1.11.0"];
                               deps = ["Markdown" => MARKDOWN])
        @test !isnothing(sol)
        @test sol[MARKDOWN] == v"1.11.0"
        @test !haskey(sol, JULIA_SYNTAX_HIGHLIGHTING)
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

    # The Julia universe is not a query parameter either: the provider offers
    # *every* Julia version, along with every stdlib version any of them bundles,
    # and the `--julia` / project bound is an ordinary compat entry on `julia`. The
    # stdlib <-> Julia couplings then do the work a restricted universe used to:
    # a Julia the bound rules out cannot be chosen, so neither can a version only
    # that Julia bundles.
    @testset "the julia universe is the whole universe" begin
        pd = Resolver.pkg_data(reg, [JULIA_UUID])[JULIA_UUID]
        # every Julia there is, 0.x and the prerelease included -- the admission
        # kinds and the bound are what narrow it
        @test length(pd.versions) == length(Registries.JULIA_VERSIONS)
        @test any(v -> v.major == 0, pd.versions)
        @test any(v -> !isempty(v.prerelease), pd.versions)

        # ... and each of them still pins the stdlibs it bundles, so the pins are
        # per candidate Julia as before (LinearAlgebra is versioned with Julia)
        @test v"1.10.8" ∈ pd.compat[v"1.10.8"][LINEAR_ALGEBRA]
        @test v"1.11.0" ∉ pd.compat[v"1.10.8"][LINEAR_ALGEBRA]
        @test !haskey(pd.compat[v"1.10.8"], STATISTICS) # upgradable: unpinned

        # A version that exists only as a bundled stdlib is admitted by exactly
        # the Julias that bundle it, drawn from the whole universe -- so a bound
        # elsewhere that only such a version satisfies steers the Julia choice,
        # and a Julia bound that excludes those Julias is unsatisfiable rather
        # than inert. Statistics 1.10.0 ships only with Julia 1.10 (the registry
        # starts at 1.11.0).
        sd = Resolver.pkg_data(reg, [STATISTICS])[STATISTICS]
        @test v"1.10.0" in sd.versions
        julias = sd.compat[v"1.10.0"][JULIA_UUID]
        @test v"1.10.8" ∈ julias
        @test v"1.11.0" ∉ julias
        @test v"1.9.0" ∉ julias
        # a registry version keeps its own bound (Statistics 1.11.1 declares
        # `julia = "1.9.4 - 1"`, so nothing older) and gains the Julias that
        # bundle it (1.11.x, whose pins would otherwise be unsatisfiable there)
        reg_julias = sd.compat[v"1.11.1"][JULIA_UUID]
        @test v"1.9.0" ∉ reg_julias
        @test v"1.11.6" ∈ reg_julias

        # the whole point: one artifact, several Julia bounds, each answer the
        # one a from-scratch resolve of that bound gives
        reqs = sort([LINEAR_ALGEBRA, JSON, JULIA_UUID])
        info = Resolver.pkg_info(reg, reqs)
        for julia in (VersionSpec("1.9"), VersionSpec("1.10"), VersionSpec("1.11"),
                      VersionSpec("1"), VersionSpec("1.10.8"))
            prob = make_problem(reqs; julia)
            sol = Resolver.resolve(info, prob)
            @test !isnothing(sol)
            @test sol[JULIA_UUID] ∈ julia
            @test isempty(sol[JULIA_UUID].prerelease)
            @test sol == Resolver.resolve(reg, prob)
        end
        # and a bound no Julia satisfies is unsatisfiable, not silently widened
        @test isnothing(Resolver.resolve(info,
            make_problem(reqs; julia = VersionSpec("99")); diagnose = false))
        @test !Resolver.issatisfiable(info,
            make_problem(reqs; julia = VersionSpec("99")))
    end

    # The Julia versions to resolve for come from `--julia` if given, otherwise
    # from the project's own `[compat] julia` bound, otherwise from the `1`
    # default -- the one bound whose reach goes beyond its own package, since the
    # stdlib couplings propagate it.
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

    # An unsatisfiable project gets a report rather than a verdict: which
    # requirements cannot be satisfied, which constraint is responsible, and
    # what to change. The exit status is unchanged -- nonzero, which is what
    # `resolve_versions` reads as "no solution".
    @testset "an unsatisfiable project gets a report" begin
        err = IOBuffer()
        @test isnothing(resolve_versions("LinearAlgebra = \"99\"";
            deps = ["LinearAlgebra" => LINEAR_ALGEBRA], err))
        msg = String(take!(err))
        @test occursin("Unsatisfiable", msg)
        @test occursin("1 conflict", msg)
        # the requirement, the package your bound emptied, and the imperative —
        # by name, since a uuid is not what the reader knows the package as
        @test occursin("cannot be satisfied", msg)
        @test occursin("you require LinearAlgebra", msg)
        @test occursin("no version of LinearAlgebra is available", msg)
        # (the report is filled to a line width, so the clause may be broken)
        @test occursin("are excluded by", msg)
        @test occursin("Fix it by any one of:", msg)
        @test occursin("relax your compat on LinearAlgebra", msg)
        @test !occursin(string(LINEAR_ALGEBRA), msg)
        # and nothing of how the answer was found
        for word in ("class", "literal", "assum", "clause", "solver")
            @test !occursin(word, msg)
        end
    end

    # Manifest generation resolves for the target Julia but then hands the
    # result to Pkg, which needs to be told which Julia that was: the packages
    # to download and the stdlib set to reconcile them against are the target's,
    # not the host's. Left to default, Pkg assumed the host's version and looked
    # for a source path for every package the *host* does not bundle -- so
    # MbedTLS_jll, a stdlib on 1.10 but a downloadable package from 1.11 on,
    # broke the build with "could not find source path" on any host past 1.10.
    # (On a 1.10 host there is no version gap and this passes trivially.)
    @testset "a manifest is generated for the target Julia" begin
        manifest = resolve_manifest("1.10"; deps = ["MbedTLS_jll" => MBEDTLS_JLL])
        @test occursin("julia_version = \"1.10", manifest)
        @test occursin("MbedTLS_jll", manifest)
    end
end
