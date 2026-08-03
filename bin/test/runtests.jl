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
# Issue #54: Python_jll's Compat.toml bounds LibMPDec_jll over `3 - 3.10`, but
# its Deps.toml only lists LibMPDec_jll from 3.10 on -- so the 3.8.x versions
# carry a compat entry for a package that is not one of their dependencies.
const PYTHON_JLL = UUID("93d3a430-8e7c-50da-8e8d-3dfcfb3baf05")
const LIBMPDEC_JLL = UUID("7106de7a-f406-5ef1-84f7-3345f7341bd2")

# Load packages from the installed registries (mirrors bin/resolve.jl).
const packages = Dict{UUID,Vector{PkgEntry}}()
for reg in reachable_registries()
    for (uuid, entry) in reg.pkgs
        push!(get!(()->PkgEntry[], packages, uuid), entry)
    end
end

# The one package universe. `allow_yanked` is the only parameter left, and it is
# not a query knob: a yanked version is withdrawn, so the default universe simply
# does not contain one, and `yanked` collects what was left out.
const yanked = Dict{UUID,Vector{VersionNumber}}()
const reg = registry_provider(packages; yanked)

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
    Resolver.Problem(reqs; compat = c,
        excludes = [prerelease_exclusion(allow_pre)])
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

# `Resolver.resolve`, with the diagnosis of an unsatisfiable result skipped and
# the result mapped back to `nothing`. These tests compare resolves against one
# another and against baked-data references, which wants a value that compares
# by content -- a `Diagnosis` is freshly built each time, and two of those are
# not the same object. The diagnostics have their own tests below.
bare_resolve(args...; kws...) = Resolver.resolve(args...; diagnose = false, kws...)

# Is the requirement set resolvable for the given Julia spec?
function resolves(reqs::Vector{UUID}; julia::VersionSpec)
    prob = make_problem(reqs; julia)
    info = Resolver.pkg_info(reg, prob)
    bare_resolve(info, prob) !== nothing
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
        occursin("Unresolvable", String(take!(err))) ||
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
        sol = bare_resolve(info, prob)
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
            @test bare_resolve(data, prob; order) ==
                  bare_resolve(baked, prob)
            # ... and the one T1 artifact answers under any of them
            @test bare_resolve(info, prob; order) ==
                  bare_resolve(baked, prob)
        end
        # a per-package ordering: only julia reversed, everything else canonical
        order = u -> u == JULIA_UUID ? (<) : (>)
        baked = Dict(u => Resolver.PkgData(
                        sort(collect(d.versions); lt = order(u)),
                        d.depends, d.compat) for (u, d) in data)
        @test bare_resolve(info, prob; order) == bare_resolve(baked, prob)
        @test bare_resolve(info, prob; order)[JULIA_UUID] <
              bare_resolve(info, prob)[JULIA_UUID]
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
            excludes = [prerelease_exclusion(allow_pre)]
            prob = Resolver.Problem(reqs; excludes)
            old = bake(data, prob) # the versions deleted, as the provider used to
            @test bare_resolve(data, prob) == bare_resolve(old, reqs)
            @test bare_resolve(info, prob) == bare_resolve(old, reqs)
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

    # Unlike prerelease-ness, yankedness is not a query knob: a yanked version is
    # withdrawn, and the resolver must neither produce one nor propose one. Both
    # hold by construction once the version is not in the universe, so the
    # provider filters yanked versions out -- and records what it dropped, since
    # "the versions your compat admits are yanked" is a story the filtering would
    # otherwise lose.
    @testset "yanked versions are pre-filtered from the universe" begin
        reqs = sort([COMPAT, LIBSODIUM_JLL, JULIA_UUID])

        # `is_yanked` reads the registries, not the version number: Compat 4.0.0
        # is yanked from General and libsodium_jll's newest version is too
        @test is_yanked(packages, COMPAT, v"4.0.0")
        @test !is_yanked(packages, COMPAT, v"4.1.0")
        @test is_yanked(packages, LIBSODIUM_JLL, v"1.0.22+0")
        @test !is_yanked(packages, LIBSODIUM_JLL, v"1.0.21+0")
        @test !is_yanked(packages, JSON, v"0.21.4")
        @test !is_yanked(packages, JULIA_UUID, v"1.10.0") # not registered at all
        @test !is_yanked(packages, COMPAT, v"99.0.0")     # no such version

        # so the provider does not offer them, and the T1 artifact has none
        data = Resolver.pkg_data(reg, reqs)
        @test v"4.0.0" ∉ data[COMPAT].versions
        @test v"1.0.22+0" ∉ data[LIBSODIUM_JLL].versions
        info = Resolver.pkg_info(reg, reqs)
        @test v"4.0.0" ∉ info[COMPAT].versions
        @test v"1.0.22+0" ∉ info[LIBSODIUM_JLL].versions

        # ... and what it dropped is recorded, as explanation data
        @test v"4.0.0" in yanked[COMPAT]
        @test v"1.0.22+0" in yanked[LIBSODIUM_JLL]
        @test !haskey(yanked, JSON)

        # `--allow-yanked` is not a constraint being lifted, it is a different
        # universe: the same provider built with it has the versions back
        loose = registry_provider(packages;
                                  allow_yanked = Dict(UUID(0) => true))
        @test v"4.0.0" in Resolver.pkg_data(loose, [COMPAT])[COMPAT].versions
        # and per package, exactly as `--allow-yanked=Compat` parses
        one = registry_provider(packages; allow_yanked = Dict(COMPAT => true))
        @test v"4.0.0" in Resolver.pkg_data(one, [COMPAT])[COMPAT].versions
        @test v"1.0.22+0" ∉
              Resolver.pkg_data(one, [LIBSODIUM_JLL])[LIBSODIUM_JLL].versions

        # the filtered universe is exactly the old constrained one: resolving
        # with yankedness as an exclusion kind over the unfiltered data gives
        # what resolving the filtered data gives, prereleases and orderings and
        # all
        raw = Resolver.pkg_data(loose, reqs)
        excl = Pair{Symbol,Any}[
            :yanked => (u, v) -> !is_bundled(u, v) && is_yanked(packages, u, v)]
        prob = Resolver.Problem(reqs; excludes = excl)
        @test bare_resolve(raw, prob) == bare_resolve(data, reqs)
        both = Resolver.Problem(reqs;
            excludes = [excl[1], prerelease_exclusion(no_pre)])
        plain = Resolver.Problem(reqs; excludes = [prerelease_exclusion(no_pre)])
        @test bare_resolve(raw, both) == bare_resolve(data, plain)
        order = u -> (<)
        @test bare_resolve(raw, both; order) ==
              bare_resolve(data, plain; order)

        # end to end: libsodium_jll's newest registered version is yanked, so the
        # resolve must land on the one below it
        sol = resolve_versions("", ["--julia=1.10"];
                               deps = ["libsodium_jll" => LIBSODIUM_JLL])
        @test !isnothing(sol)
        @test sol[LIBSODIUM_JLL] == v"1.0.21+0"
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
            sol = bare_resolve(info, prob)
            @test !isnothing(sol)
            @test sol[JULIA_UUID] ∈ julia
            @test isempty(sol[JULIA_UUID].prerelease)
            @test sol == bare_resolve(reg, prob)
        end
        # and a bound no Julia satisfies is unsatisfiable, not silently widened
        @test isnothing(bare_resolve(info,
            make_problem(reqs; julia = VersionSpec("99"))))
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
end

@testset "unsat diagnostics (bin/resolve.jl)" begin
    # bin/resolve.jl diagnoses an unsatisfiable resolve automatically: the
    # report comes back from `resolve` itself, so all the script does is
    # substitute package names for uuids and print it. BulkSMS (HTTP 0.x) and
    # AnthropicClient (HTTP 1.x) require disjoint majors of HTTP.
    julia = Base.julia_cmd()[1]
    dir = mktempdir()
    write(joinpath(dir, "Project.toml"), "[deps]\n")

    err = IOBuffer()
    cmd = `$julia --project=$BIN_PROJECT $RESOLVE_JL $dir
           --extra-deps=BulkSMS,AnthropicClient --julia=1.11`
    ok = success(pipeline(cmd; stderr = err))
    report = String(take!(err))
    @test !ok                            # unsatisfiable -> nonzero exit
    @test occursin("Unresolvable", report)
    @test occursin("cannot be satisfied together", report)
    # uuids are substituted for names, and the contested dependency is named
    @test occursin("BulkSMS", report)
    @test occursin("AnthropicClient", report)
    @test occursin("HTTP", report)
    @test !occursin(r"[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}",
                    report)
    @test occursin("Verified fixes:", report)
    @test occursin("drop requirement", report)
    @test occursin("Upstream fixes:", report)

    # a project compat bound nothing can satisfy is blamed on the bound, and
    # relaxing it is offered as a fix
    dir2 = mktempdir()
    write(joinpath(dir2, "Project.toml"), """
        [deps]
        JSON = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"

        [compat]
        JSON = "0.10"
        """)
    err2 = IOBuffer()
    cmd2 = `$julia --project=$BIN_PROJECT $RESOLVE_JL $dir2 --julia=1.11`
    ok2 = success(pipeline(cmd2; stderr = err2))
    report2 = String(take!(err2))
    @test !ok2
    @test occursin("your compat restricts JSON", report2)
    @test occursin("relax your compat on JSON", report2)

    # the admission knobs are constraints like any other, so they get facts
    # and fixes like any other. Wine_jll has never had a non-prerelease
    # release, so requiring it fails on the prerelease knob alone -- and
    # `--allow-pre` is exactly the fix the report names.
    err3 = IOBuffer()
    cmd3 = `$julia --project=$BIN_PROJECT $RESOLVE_JL $dir
            --extra-deps=Wine_jll --julia=1.11`
    ok3 = success(pipeline(cmd3; stderr = err3))
    report3 = String(take!(err3))
    @test !ok3
    @test occursin("every version of Wine_jll is a prerelease", report3)
    @test occursin("allow prereleases of Wine_jll", report3)
    # ... and passing it really does resolve, to the version the fix promised
    m = match(r"allow prereleases of Wine_jll\n\s*→ allows: ([^\n]*)", report3)
    @test m !== nothing
    @test occursin("Wine_jll", m[1])
    vers = resolve_versions("", ["--allow-pre", "--julia=1.11"];
                            deps = ["Wine_jll" => WINE_JLL])
    @test vers !== nothing && haskey(vers, WINE_JLL)
    @test occursin(string(vers[WINE_JLL]), m[1])

    # Yanked versions are not in the universe at all, so no fix can ever propose
    # accepting one -- by construction, with no policy flag anywhere. What that
    # costs is the story: a bound admitting only yanked versions admits nothing,
    # and "no versions" is not the answer. The provider kept the set it dropped,
    # so the reporting layer says what really happened. libsodium_jll's newest
    # release is yanked, and `=1.0.22` admits only it.
    dir4 = mktempdir()
    write(joinpath(dir4, "Project.toml"), """
        [deps]
        libsodium_jll = "$LIBSODIUM_JLL"

        [compat]
        libsodium_jll = "=1.0.22"
        """)
    err4 = IOBuffer()
    cmd4 = `$julia --project=$BIN_PROJECT $RESOLVE_JL $dir4 --julia=1.11`
    ok4 = success(pipeline(cmd4; stderr = err4))
    report4 = String(take!(err4))
    @test !ok4
    @test occursin(
        r"the versions your compat on libsodium_jll admits are yanked \(1\.0\.22\S*\)",
        report4)
    @test !occursin("allow the yanked", report4)
    # what it does offer is relaxing the bound, whose witness is then a version
    # that is not yanked
    @test occursin("relax your compat on libsodium_jll", report4)
    @test occursin(r"→ allows:.*libsodium_jll 1\.0\.21", report4)

    # ... and `--allow-yanked` is there for the user to reach for themselves:
    # not a constraint lifted, a universe rebuilt with the versions back in
    vers4 = resolve_versions("libsodium_jll = \"=1.0.22\"",
                             ["--allow-yanked", "--julia=1.11"];
                             deps = ["libsodium_jll" => LIBSODIUM_JLL])
    @test vers4 !== nothing
    @test vers4[LIBSODIUM_JLL] == v"1.0.22+0"
    # the per-package shape of the flag works the same way, and is not global:
    # allowing another package's yanked versions leaves this one unresolvable
    @test resolve_versions("libsodium_jll = \"=1.0.22\"",
                           ["--allow-yanked=libsodium_jll", "--julia=1.11"];
                           deps = ["libsodium_jll" => LIBSODIUM_JLL]) == vers4
    @test isnothing(
        resolve_versions("libsodium_jll = \"=1.0.22\"",
                         ["--allow-yanked=Compat", "--julia=1.11"];
                         deps = ["libsodium_jll" => LIBSODIUM_JLL]))

    # julia is modelled as a dependency edge, never a requirement: it is in
    # every solution because packages need it, so no fix can offer to drop it
    # and no witness list names it unless a story blamed it
    dir5 = mktempdir()
    write(joinpath(dir5, "Project.toml"), """
        [deps]
        Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"

        [extras]
        ColorTypes = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"

        [compat]
        ColorTypes = "0.9"
        """)
    err5 = IOBuffer()
    cmd5 = `$julia --project=$BIN_PROJECT $RESOLVE_JL $dir5 --julia=1.11`
    ok5 = success(pipeline(cmd5; stderr = err5))
    report5 = String(take!(err5))
    @test !ok5
    # the right fix, named first, at versions from this decade
    @test occursin("1. relax your compat on ColorTypes", report5)
    @test occursin(r"1\. relax your compat on ColorTypes\n\s*→ allows:[^\n]*Makie 0\.2",
                   report5)
    # julia appears in the story, since its bound really is part of why ...
    @test occursin("restricts julia", report5)
    # ... but is never something a fix offers to drop
    @test !occursin("drop requirement julia", report5)
end
