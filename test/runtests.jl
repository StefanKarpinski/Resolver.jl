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
            # iteratively pick a solution and break it until unsolvable
            while true
                fill_data!(m, n, deps, comp, data)
                pkgs, vers = test_resolver(data, reqs)
                all(isnothing, vers) && break
                # pick a random solution
                k = rand(1:size(vers,2))
                # there must be some non-nothing version
                @test any(vers[i, k] !== nothing for i=1:size(vers,1))
                # pick a package with non-nothing version
                i = rand(1:size(vers,1))
                while vers[i, k] === nothing
                    i = rand(1:size(vers,1))
                end
                p = pkgs[i]
                v = vers[i, k]
                # pick a different package (version can be nothing)
                q = rand(1:m-1)
                q += q ≥ p
                j = something(findfirst(==(q), pkgs), 0)
                w = get(vers, (j, k), nothing)
                # make resolved versions of p & q incompatible
                if isnothing(w)
                    # add a dependency p@v => q
                    b = bit(p, v, q)
                    @assert iszero(deps.bits & b)
                    deps = typeof(deps)(deps.bits | b)
                else
                    # add incompatibility p@v ⊼ p@w
                    w = vers[j, k]
                    b = bit(p, v, q, w)
                    @assert iszero(comp.bits & b)
                    comp = typeof(comp)(comp.bits | b)
                end
            end
        end
    end
end

@testset "resolve: max keyword" begin
    for ex in tiny_data.examples
        data, reqs = ex.data, ex.reqs
        pkgs, vers = test_resolver(data, reqs)
        nsol = size(vers,2)
        for max = 1:nsol
            pkgs, vers = resolve(data, reqs; max)
            @test max == size(vers,2)
        end
        pkgs, vers = resolve(data, reqs; max=-1)
        @test nsol == size(vers,2)
        pkgs, vers = resolve(data, reqs; max=0)
        @test nsol == size(vers,2)
        pkgs, vers = resolve(data, reqs; max=nsol+1)
        @test nsol == size(vers,2)
        pkgs, vers = resolve(data, reqs; max=typemax(Int))
        @test nsol == size(vers,2)
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
    pkgs, vers = resolve(data, [:PkgA])
    @test :PkgA in pkgs

    # Test missing required package
    data = Dict(
        :PkgA => PkgData([:v1], Dict{Symbol,Vector{Symbol}}(), Dict{Symbol,Dict{Symbol,Any}}()),
        :PkgB => PkgData([:v1], Dict{Symbol,Vector{Symbol}}(), Dict{Symbol,Dict{Symbol,Any}}())
    )
    @test_throws ArgumentError("Required package MissingReq is not available in the package data") resolve(data, [:MissingReq])
end

@testset "registry resolve" begin
    rp = registry.provider()
    test_resolver(rp, ["JSON"])
    test_resolver(rp, ["DataFrames"])
    test_resolver(rp, ["DataFrames", "JSON"])
    test_resolver(rp, ["DifferentialEquations"])
    test_resolver(rp, ["DifferentialEquations", "JSON"])
    test_resolver(rp, ["DifferentialEquations", "JSON", "DataFrames"])
    # test some details
    pkgs, vers = resolve(rp, ["JSON"])
    @test pkgs isa Vector{String}
    @test vers isa Matrix{Union{Nothing,VersionNumber}}
    @test all(!isnothing, vers)
    # test corner case (empty)
    pkgs, vers = resolve(rp, String[])
    @test pkgs isa Vector{String}
    @test vers isa Matrix{Union{Nothing,VersionNumber}}
    @test isempty(pkgs)
    @test isempty(vers)
end

@testset "yanked packages" begin
    # Compat.jl v4.0.0 is yanked from the General registry, see
    # https://github.com/JuliaRegistries/General/blob/1cb12c7cf0c4ce32b24daa8373c18541d787bae2/C/Compat/Versions.toml#L198
    # This leads to issues with Resolver.jl used from
    # https://github.com/julia-actions/julia-downgrade-compat
    # To test this, we use the same commands as the
    # `julia-downgrade-compat` action to resolve a project
    # depending on Compat with "4.0" compatibility requirement.

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
