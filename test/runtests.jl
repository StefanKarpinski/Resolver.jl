include("setup.jl")

@testset "tiny tests, complete" begin
    for m = 1:2, n = 1:2
        make_deps, make_comp, data, d, c = tiny_data_makers(m, n)
        # all possible dependency patterns
        for deps_bits = 0:2^d-1
            deps = make_deps(deps_bits)
            all(deps[i].bits ≥ deps[i+1].bits for i=1:m-1) || continue
            # all possible compatibility patterns
            for comp_bits = 0:2^c-1
                comp = make_comp(comp_bits)
                fill_data!(m, n, deps, comp, data)
                # all possible requirements sets
                for reqs_bits = 0:2^m-1
                    reqs = make_reqs(reqs_bits)
                    test_resolver(data, reqs)
                end
            end
        end
    end
end

@testset "small tests, semi-full" begin
    Random.seed!(rand(RandomDevice(), UInt64))
    for m = 2:3, n = 2:3
        m == n && continue # fully tested or too large
        make_deps, make_comp, data, d, c = tiny_data_makers(m, n)
        # all dependency patterns + random compatibility pattern
        for deps_bits = 0:2^d-1
            deps = make_deps(deps_bits)
            comp = make_comp(randbits(c))
            fill_data!(m, n, deps, comp, data)
            for reqs_bits = 1:2^m-1
                # @show m, n, deps_bits, reqs_bits
                # println("comp = Comp(", comp.bits, ")")
                reqs = make_reqs(reqs_bits)
                test_resolver(data, reqs)
            end
        end
        # all compatibility patterns + random dependency pattern
        for comp_bits = 0:2^c-1
            deps = make_deps(randbits(d))
            comp = make_comp(comp_bits)
            fill_data!(m, n, deps, comp, data)
            for reqs_bits = 1:2^m-1
                # @show m, n, comp_bits, reqs_bits
                # println("deps = Deps(", deps.bits, ")")
                reqs = make_reqs(reqs_bits)
                test_resolver(data, reqs)
            end
        end
    end
end

@testset "small tests, random" begin
    Random.seed!(rand(RandomDevice(), UInt64))
    for m = 2:5, n = 2:5
        16 < (m*n)^2 ≤ 128 || continue
        make_deps, make_comp, data, d, c = tiny_data_makers(m, n)
        for _ = 1:1000
            # random dependencies, compatibility and requirements
            deps = make_deps(randbits(d))
            comp = make_comp(randbits(c))
            reqs = make_reqs(rand(1:2^m-1))
            fill_data!(m, n, deps, comp, data)
            test_resolver(data, reqs)
        end
    end
end

@testset "medium tests, adversarial" begin
    Random.seed!(rand(RandomDevice(), UInt64))
    for m = 2:5, n = 1:5
        (m*n)^2 ≤ 128 || continue
        make_deps, make_comp, data, d, c, bit = tiny_data_makers(m, n)
        for _ = 1:100
            deps = make_deps(0) # start with no dependencies
            comp = make_comp(0) # start with no conflicts
            reqs = make_reqs(rand(1:2^m-1)) # random requirements set
            # iteratively break the solution until unsolvable
            while true
                fill_data!(m, n, deps, comp, data)
                sol = test_resolver(data, reqs)
                sol === nothing && break
                # pick a resolved package
                p = rand(collect(keys(sol)))
                v = sol[p]
                # pick a different package (it may be unresolved)
                q = rand(1:m-1)
                q += q ≥ p
                w = get(sol, q, nothing)
                # make resolved versions of p & q incompatible
                if isnothing(w)
                    # add a dependency p@v => q
                    b = bit(p, v, q)
                    @assert iszero(deps.bits & b)
                    deps = typeof(deps)(deps.bits | b)
                else
                    # add incompatibility p@v ⊼ q@w
                    b = bit(p, v, q, w)
                    @assert iszero(comp.bits & b)
                    comp = typeof(comp)(comp.bits | b)
                end
            end
        end
    end
end

@testset "resolve: single solution" begin
    for ex in tiny_data.examples
        data, reqs = ex.data, ex.reqs
        sol = resolve(data, reqs)
        # a package => version dict covering the requirements, or nothing
        if sol !== nothing
            @test sol isa Dict{Int,Int}
            @test reqs ⊆ keys(sol)
        end
        # deterministic: resolving again, and with the requirements supplied
        # in a different order, yields the same solution
        @test resolve(data, reqs) == sol
        @test resolve(data, reverse(collect(reqs))) == sol
    end
end

@testset "resolve: brute-force reference" begin
    Random.seed!(rand(RandomDevice(), UInt64))
    hi = p -> p  # default priority: lower package id first
    lo = p -> -p # reversed priority
    for m = 1:3, n = 1:3
        make_deps, make_comp, data, d, c = tiny_data_makers(m, n)
        for _ = 1:100
            deps = make_deps(randbits(d))
            comp = make_comp(randbits(c))
            fill_data!(m, n, deps, comp, data)
            for reqs_bits = 0:2^m-1, by in (hi, lo)
                reqs = collect(make_reqs(reqs_bits))
                @test resolve(data, reqs; by) == ref_resolve(data, reqs; by)
            end
        end
    end
end

using Resolver: PkgData
@testset "consistency validation" begin
    # Test missing dependency
    data = Dict(
        :PkgA => PkgData([:v1], Dict(:v1 => [:MissingDep, :PkgB]), Dict{Symbol,Dict{Symbol,Any}}()),
        :PkgB => PkgData([:v1], Dict{Symbol,Vector{Symbol}}(), Dict{Symbol,Dict{Symbol,Any}}())
    )
    @test_throws ArgumentError("Package PkgA depends on MissingDep, but MissingDep is not available in the package data") resolve(data, [:PkgA])

    # Test that compatibility constraints with missing packages are allowed
    data = Dict(
        :PkgA => PkgData([:v1], Dict{Symbol,Vector{Symbol}}(), Dict(:v1 => Dict(:MissingCompat => Any[]))),
        :PkgB => PkgData([:v1], Dict{Symbol,Vector{Symbol}}(), Dict{Symbol,Dict{Symbol,Any}}())
    )
    sol = resolve(data, [:PkgA])
    @test sol !== nothing && haskey(sol, :PkgA)

    # Test missing required package
    data = Dict(
        :PkgA => PkgData([:v1], Dict{Symbol,Vector{Symbol}}(), Dict{Symbol,Dict{Symbol,Any}}()),
        :PkgB => PkgData([:v1], Dict{Symbol,Vector{Symbol}}(), Dict{Symbol,Dict{Symbol,Any}}())
    )
    @test_throws ArgumentError("Required package MissingReq is not available in the package data") resolve(data, [:MissingReq])
end

using Resolver: pkg_info, save_pkg_info_file, load_pkg_info_file
@testset "pkg info files" begin
    # roundtrip with String package names & integer versions
    sdata = Dict(
        "A" => PkgData([2, 1], Dict(2 => ["B"], 1 => ["B"]),
                       Dict(2 => Dict("B" => [1]))),
        "B" => PkgData([2, 1], Dict{Int,Vector{String}}(),
                       Dict{Int,Dict{String,Vector{Int}}}()),
    )
    for filter in (false, true)
        info = pkg_info(sdata, ["A"]; filter)
        path = save_pkg_info_file(info)
        @test load_pkg_info_file(path) == info
        rm(path)
    end
    # roundtrip with Symbol package names & Symbol versions
    ydata = Dict(
        :A => PkgData([:v2, :v1], Dict(:v2 => [:B], :v1 => [:B]),
                      Dict(:v2 => Dict(:B => [:v1]))),
        :B => PkgData([:v2, :v1], Dict{Symbol,Vector{Symbol}}(),
                      Dict{Symbol,Dict{Symbol,Vector{Symbol}}}()),
    )
    for filter in (false, true)
        info = pkg_info(ydata, [:A]; filter)
        path = save_pkg_info_file(info)
        @test load_pkg_info_file(path) == info
        rm(path)
    end
end

@testset "zero-version packages" begin
    # a package with no versions can never be installed: requiring it makes
    # the problem unsatisfiable, versions depending on it are uninstallable,
    # and the filter must not leave references to it behind
    nodeps = Dict{Symbol,Vector{Symbol}}()
    nocomp = Dict{Symbol,Dict{Symbol,Any}}()

    # a zero-version package required directly
    data = Dict(
        :A => PkgData(Symbol[], nodeps, nocomp),
        :B => PkgData([:v1], nodeps, nocomp),
    )
    @test resolve(data, [:A]) === nothing
    @test resolve(data, [:A, :B]) === nothing
    @test resolve(data, [:B]) == Dict(:B => :v1)
    @test isempty(resolve(data, Symbol[]))
    @test !haskey(pkg_info(data, [:A]), :A)
    test_resolver(data, [:A])
    test_resolver(data, [:B])

    # A's only version depends on the zero-version package B
    data = Dict(
        :A => PkgData([:v1], Dict(:v1 => [:B]), nocomp),
        :B => PkgData(Symbol[], nodeps, nocomp),
    )
    @test resolve(data, [:A]) === nothing
    @test isempty(resolve(data, Symbol[]))
    @test isempty(pkg_info(data, [:A]))
    @test isempty(pkg_info(data, Symbol[]))
    # unfiltered info keeps the package (and its dependents) as they are
    info = pkg_info(data, [:A]; filter = false)
    @test isempty(info[:B].versions)
    @test info[:A].depends == [:B]
    test_resolver(data, [:A])

    # uninstallability propagates transitively
    data = Dict(
        :A => PkgData([:v1], Dict(:v1 => [:B]), nocomp),
        :B => PkgData(Symbol[], nodeps, nocomp),
        :C => PkgData([:v1], Dict(:v1 => [:A]), nocomp),
    )
    @test resolve(data, [:C]) === nothing
    @test isempty(pkg_info(data, [:C]))
    test_resolver(data, [:C])

    # A@v2 depends on the zero-version package B, A@v1 doesn't: degrade
    data = Dict(
        :A => PkgData([:v2, :v1], Dict(:v2 => [:B]), nocomp),
        :B => PkgData(Symbol[], nodeps, nocomp),
        :C => PkgData([:v1], Dict(:v1 => [:A]), nocomp),
    )
    @test resolve(data, [:A]) == Dict(:A => :v1)
    @test resolve(data, [:C]) == Dict(:C => :v1, :A => :v1)
    @test isempty(resolve(data, Symbol[]))
    info = pkg_info(data, [:A])
    @test !haskey(info, :B)
    @test info[:A].versions == [:v1]
    @test isempty(info[:A].depends)
    test_resolver(data, [:A])
    test_resolver(data, [:C])

    # resolving an *unpruned* info reaches the same answers: dropping the
    # version-less package leaves `depends` naming a package that is gone, and
    # arc consistency has to read that as uninstallable too
    data = Dict(
        :A => PkgData([:v1], Dict(:v1 => [:B]), nocomp),
        :B => PkgData(Symbol[], nodeps, nocomp),
    )
    info = pkg_info(data, [:A]; filter = false)
    @test resolve(info, [:A]) === nothing
    @test resolve(info, [:A]; group = false) === nothing
    data = Dict(
        :A => PkgData([:v2, :v1], Dict(:v2 => [:B]), nocomp),
        :B => PkgData(Symbol[], nodeps, nocomp),
        :C => PkgData([:v1], Dict(:v1 => [:A]), nocomp),
    )
    info = pkg_info(data, [:C]; filter = false)
    @test resolve(info, [:C]) == Dict(:C => :v1, :A => :v1)
    @test resolve(info, [:C]; group = false) == Dict(:C => :v1, :A => :v1)
end

include("unsat_cores.jl")
include("problem.jl")
include("satisfiable.jl")
include("classes.jl")
include("ordering.jl")

@testset "registry resolve" begin
    rp = registry.provider()
    test_resolver(rp, ["JSON"])
    test_resolver(rp, ["DataFrames"])
    test_resolver(rp, ["DataFrames", "JSON"])
    test_resolver(rp, ["DifferentialEquations"])
    test_resolver(rp, ["DifferentialEquations", "JSON"])
    test_resolver(rp, ["DifferentialEquations", "JSON", "DataFrames"])
    # test some details
    sol = resolve(rp, ["JSON"])
    @test sol isa Dict{String,VersionNumber}
    @test haskey(sol, "JSON")
    # test corner case (empty)
    sol = resolve(rp, String[])
    @test sol isa Dict{String,VersionNumber}
    @test isempty(sol)
    # the cheap verdict on a real registry, straight from the provider
    @test issatisfiable(rp, ["JSON"])
    @test issatisfiable(rp, ["DataFrames", "JSON"])
    @test issatisfiable(rp, String[])
    # ... including through the convenience keywords: no version of JSON can
    # satisfy a pin at a version that doesn't exist
    @test !issatisfiable(rp, ["JSON"]; pins = Dict("JSON" => v"0.0.0"))
    @test resolve(rp, ["JSON"]; pins = Dict("JSON" => v"0.0.0")) === nothing
end

@testset "yanked packages" begin
    # Compat.jl v4.0.0 is yanked from the General registry, see
    # https://github.com/JuliaRegistries/General/blob/1cb12c7cf0c4ce32b24daa8373c18541d787bae2/C/Compat/Versions.toml#L198
    # This leads to issues with Resolver.jl used from
    # https://github.com/julia-actions/julia-downgrade-compat
    # To test this, we use the same commands as the
    # `julia-downgrade-compat` action to resolve a project
    # depending on Compat with "4.0" compatibility requirement.
    #
    # The struck version is not in the package universe at all: bin/Registries.jl
    # filters yanked versions out at the provider, so the resolve below cannot
    # produce v4.0.0 or offer it as an alternative (`--allow-yanked` is what asks
    # for a universe that keeps them). What this testset checks is the outcome --
    # the version the action ends up with -- which is the same whichever
    # mechanism enforces it; the mechanism itself is covered in bin/test/.

    julia = Base.julia_cmd()[1]
    project = pkgdir(Resolver, "bin")
    script = joinpath(project, "resolve.jl")
    # Resolve for the running Julia's release series -- but on a prerelease
    # (e.g. nightly) that series has no registered releases to resolve against,
    # so fall back to the LTS series (a released, ≥1.9 target that yields a
    # modern-format manifest).
    julia_version = isempty(VERSION.prerelease) ?
        string(Int(VERSION.major), ".", VERSION.minor) : "1.10"

    # Since Pkg.test runs the tests in a temporary directory, we write the
    # Project.toml file in a temporary directory as well.
    # dir = joinpath(@__DIR__, "ProjectWithYankedDependency") TODO
    dir = mktempdir()
    open(joinpath(dir, "Project.toml"), "w") do io
        write(io, """
[deps]
Compat = "34da2185-b29b-5c13-b0c7-acf172513d20"
Pkg = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
UUIDs = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[compat]
Compat = "4"
""")
    end

    run(`$julia --project=$project -e 'import Pkg; Pkg.instantiate()'`)
    # In both cases we assert `success` rather than `@test_nowarn`: resolving
    # legitimately prints Pkg's "Installed ..." progress to stderr, which
    # `@test_nowarn` would flag as a warning on a clean depot.
    if isempty(VERSION.prerelease)
        # Released Julia: exercise the full manifest-generation path, then
        # confirm the yanked Compat v4.0.0 was not selected (the second command
        # exits nonzero if it was).
        @test success(`$julia --project=$project $script $dir --min=@deps --julia=$julia_version`)
        @test success(`$julia --project=$dir -e '
            using Pkg, UUIDs
            deps = Pkg.dependencies()
            pkg = deps[UUID("34da2185-b29b-5c13-b0c7-acf172513d20")]
            compat_version = pkg.version
            exit(compat_version == v"4.0.0" ? 1 : 0)'`)
    else
        # Prerelease (e.g. nightly): writing a manifest needs the host and
        # target Julia to match, which they can't when the host has no matching
        # registered release. Read the resolved version directly instead.
        out = IOBuffer()
        @test success(pipeline(`$julia --project=$project $script $dir --print-versions --min=@deps --julia=$julia_version`; stdout=out))
        m = match(r"34da2185-b29b-5c13-b0c7-acf172513d20 +\S+ +(\S+)", String(take!(out)))
        @test m !== nothing
        @test VersionNumber(m.captures[1]) != v"4.0.0"
    end
end

@testset "bin/ tooling regression tests" begin
    # The bin/ resolver tooling has its own integration tests (they need the
    # General registry, so they live in bin/test/ and run in the bin/
    # environment rather than as part of this package's offline suite). Run
    # them here as a subprocess so CI exercises them too.
    julia = Base.julia_cmd()[1]
    project = pkgdir(Resolver, "bin")
    tests = joinpath(project, "test", "runtests.jl")
    run(`$julia --project=$project -e 'import Pkg; Pkg.instantiate()'`)
    @test success(`$julia --project=$project $tests`)
end

include("workspaces.jl")
