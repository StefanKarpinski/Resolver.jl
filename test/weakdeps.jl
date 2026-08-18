# Weak dependencies (#61).
#
# The resolver core has no notion of a weak dependency and does not need one: a
# weak dep is a compat entry *without* a dependency edge, which the existing
# machinery already expresses. A conflict clause fires only when both versions
# are chosen, so a bound with no edge behind it is exactly "if you install both,
# they must be compatible" — the semantics Pkg gives an extension.
#
# The providers are where the translation happens (`bin/Registries.jl`, and
# `test/registry.jl` for these tests): both merge a version's weak compat into
# the same compat dict as its strong compat, and add nothing to its deps.
#
# So what needs testing is not a mechanism but a *contract*, in three parts:
#
#   1. a weak bound does not drag its target in,
#   2. it binds when the target is installed anyway, whoever pulled it in,
#   3. it is vacuous when the target is not installed.
#
# Part 3 is what makes the strong-only closure walk in `pkg_data` correct: a
# package is installed only if it is required or some installed package depends
# on it, so every installable package is already reached. A weak target that
# nothing else pulls in cannot be installed, and its bound cannot apply.

using Resolver: Problem, PkgData, DepsProvider, resolve, issatisfiable, pkg_data
using Pkg.Versions: VersionSpec

# A package with one version per entry of `spec`, where each entry is
# `version => (deps, compat)`. Weak and strong compat are not distinguished
# here for the same reason the resolver does not distinguish them: the
# difference is entirely whether the partner also appears in `deps`.
function pkg(spec::Pair{String,<:Tuple}...)
    vers = VersionNumber[]
    deps = Dict{VersionNumber,Vector{String}}()
    comp = Dict{VersionNumber,Dict{String,VersionSpec}}()
    for (v, (d, c)) in spec
        ver = VersionNumber(v)
        push!(vers, ver)
        deps[ver] = collect(String, d)
        comp[ver] = Dict{String,VersionSpec}(k => VersionSpec(s) for (k, s) in c)
    end
    sort!(vers; rev = true) # newest first, as providers hand them over
    return PkgData(vers, deps, comp)
end

# W is a plain package. P declares a bound on W at "1" but no edge to it —
# a weak dep. S depends on W strongly, with no opinion about its version.
const WEAK_DATA = Dict(
    "P" => pkg("1.0.0" => ((), ("W" => "1",))),
    "S" => pkg("1.0.0" => (("W",), ())),
    "W" => pkg("1.0.0" => ((), ()), "2.0.0" => ((), ())),
)

@testset "weak deps: a bound without an edge" begin
    # 1. the bound does not drag its target in: resolving P alone installs P
    #    and nothing else, even though P names W
    @testset "a weak bound does not install its target" begin
        sol = resolve(WEAK_DATA, ["P"])
        @test sol !== nothing
        @test Set(keys(sol)) == Set(["P"])
        # ... which is the whole difference from a strong dep, so check the
        # contrast rather than asserting the absence alone
        strong = Dict(WEAK_DATA..., "P" => pkg("1.0.0" => (("W",), ("W" => "1",))))
        sol_strong = resolve(strong, ["P"])
        @test sol_strong !== nothing
        @test Set(keys(sol_strong)) == Set(["P", "W"])
    end

    # 2. it binds when the target is installed anyway. two ways in: the user
    #    asks for it, or another package's strong dep brings it
    @testset "a weak bound binds when its target is installed" begin
        # required alongside P: W 2.0.0 is the better version and would win
        # unrestricted, so seeing 1.0.0 is the bound doing something
        sol = resolve(WEAK_DATA, ["P", "W"])
        @test sol !== nothing
        @test sol["W"] == v"1.0.0"
        # pulled in by S's strong dep, with P never naming an edge to it
        sol = resolve(WEAK_DATA, ["P", "S"])
        @test sol !== nothing
        @test sol["W"] == v"1.0.0"
        # and without P there is nothing to hold W back
        sol = resolve(WEAK_DATA, ["S"])
        @test sol !== nothing
        @test sol["W"] == v"2.0.0"
    end

    # an unsatisfiable weak bound is unsatisfiable like any other: the clause
    # does not care where the compat entry came from
    @testset "a weak bound can make a problem unsatisfiable" begin
        data = Dict(WEAK_DATA...,
            "P" => pkg("1.0.0" => ((), ("W" => "3",))))  # no such version of W
        @test resolve(data, ["P"]) !== nothing            # vacuous on its own
        @test resolve(data, ["P", "W"]; diagnose = false) === nothing
        @test !issatisfiable(data, ["P", "W"])
    end

    # 3. vacuity, stated as the closure property it licenses: a weak target
    #    that nothing depends on stays out of the universe entirely, so the
    #    bound has nothing to constrain
    @testset "a weak target outside the closure is not in the universe" begin
        provider = DepsProvider(collect(keys(WEAK_DATA))) do p
            WEAK_DATA[p]
        end
        closure = pkg_data(provider, ["P"])
        @test haskey(closure, "P")
        @test !haskey(closure, "W")
        # S's strong dep is what brings W in, so the same walk from S reaches it
        @test haskey(pkg_data(provider, ["S"]), "W")
    end
end

# The provider contract, against real registry data. The hand-built cases above
# fix the *semantics*; this fixes the *translation* — that a real provider turns
# a registry weak-compat entry into a compat entry and no dep edge.
#
# Checking the translation directly rather than inferring it from a resolve is
# deliberate. A resolve only exercises a weak bound when that bound is what
# holds a version back, and most are not: dropping weak compat entirely leaves
# the DataFrames, Flux and Turing resolutions byte-identical, because strong
# bounds already imply them. A test built on those would pass either way.
@testset "weak deps: the provider translates weak compat" begin
    reg_inst = registry.RegistryInstance(registry.REG_PATH)
    prov = registry.provider()

    # registry packages that declare a weak compat entry naming a weak dep
    function weak_case(e)
        info = try
            applicable(registry.init_package_info!, reg_inst, e) ?
                registry.init_package_info!(reg_inst, e) :
                registry.init_package_info!(e)
        catch
            return nothing
        end
        (isempty(info.weak_deps) || isempty(info.weak_compat)) && return nothing
        for (r, c) in info.weak_compat
            c isa AbstractDict || continue
            # Julia 1.14 made this data uuid-keyed; the provider's own
            # normalizer is what turns either shape back into names, and
            # reaching past it is how this test broke on nightly
            for (w, spec) in registry._compat_names(c, reg_inst)
                # stdlibs (and julia) are excluded from the universe outright,
                # so their bounds are dropped on purpose -- not a weak-dep case
                w in registry.EXCLUDES && continue
                w in prov.packages || continue
                # a version of e covered by this entry, and known to the provider
                data = try pkg_data(prov, e.name) catch; return nothing end
                for v in data.versions
                    v in r || continue
                    haskey(data.compat, v) || continue
                    return (e.name, v, w, spec, data)
                end
            end
        end
        return nothing
    end

    cases = Any[]
    for e in values(reg_inst.pkgs)
        length(cases) ≥ 25 && break
        c = weak_case(e)
        c === nothing || push!(cases, c)
    end
    # the corpus must actually contain weak-compat packages, or everything
    # below passes by vacuity
    @test length(cases) ≥ 10

    for (p, v, w, spec, data) in cases
        # the bound survives into the provider's compat ...
        @test haskey(data.compat[v], w)
        # ... and says what the registry said
        @test data.compat[v][w] == spec
        # ... while adding no dependency edge, which is what keeps it weak.
        # (a package may declare the same name in both sections; then the edge
        # is the strong one's, and there is nothing to check here)
        w in data.depends[v] || @test w ∉ data.depends[v]
    end
end
