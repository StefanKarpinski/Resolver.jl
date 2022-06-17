include("setup.jl")

@testset "core solver" begin
    @testset "basic example" begin
        packages = [fill("A", 2); fill("B", 4); fill("C", 4)]
        conflicts = [(1, 3), (1, 7), (3, 7), (2, 10)]
        solutions = [[1, 4, 8], [2, 3, 8], [2, 4, 7]]
        @test solutions == resolve_brute_force(packages, conflicts)
        @test solutions == resolve_brute_force(packages, Set(conflicts))
        @test solutions == resolve_core(packages, conflicts)
        @test solutions == resolve_core(packages, Set(conflicts))
        for relax = 1:5
            solutions′ = resolve_brute_force(packages, conflicts; relax)
            resolved′ = resolve_core(packages, conflicts; relax)
            @test solutions′ == resolved′
        end
        @testset "permutations" for _ = 1:100
            # permute the versions, conflicts and solutions
            packages′ = shuffle(packages)
            p = sortperm(packages′)
            conflicts′ = [(p[c[1]], p[c[2]]) for c in conflicts]
            solutions′ = [p[solution] for solution in solutions]
            foreach(sort!, solutions′)
            sort!(solutions′)
            @test solutions′ == resolve_brute_force(packages′, conflicts′)
            @test solutions′ == resolve_core(packages′, conflicts′)
        end
        @testset "unordered package names" begin
            packages = [fill(1+0im, 2); fill(2+0im, 4); fill(3+0im, 4)]
            @test solutions == resolve_brute_force(packages, conflicts)
            @test solutions == resolve_core(packages, conflicts)
        end
    end

    @testset "comprehensive tests" begin
        Random.seed!(0x53f3ce0656b85450bbb52b15fc58853f)
        count = zeros(Int, 17)
        for M = 2:5, # number of packages
            V = 1:5  # number of versions per package
            T = V^2*(M*(M-1)÷2)
            T ≤ 128 || continue
            packages = shuffle!([p for p=1:M for _=1:V])
            p = sortperm(packages)
            for C in (T ≤ 12 ? (0:2^T-1) : [randu128(3-(k%5)%3) for k=1:2^12])
                G = gen_conflicts(M, V, C)
                conflicts = Tuple{Int,Int}[(p[c[1]], p[c[2]]) for c in G]
                @assert length(conflicts) == count_ones(C % UInt128(2)^T)
                for relax = 0:M^2
                    solutions = resolve_brute_force(packages, conflicts; relax)
                    resolved = resolve_core(packages, conflicts; relax)
                    @test resolved isa Vector{Vector{Int}}
                    @test solutions == resolved
                    count[length(solutions)+1] += 1
                end
            end
        end
        # specific to chosen seed (chosen to have some large solution sets)
        @test count == [
            12693, 605486, 17731, 11964, 5792, 2073,
            777, 341, 163, 65, 23, 11, 9, 1, 0, 0, 1,
        ]
    end

    @testset "edge cases" begin
        @testset "0 pkgs" begin
            packages = []
            conflicts = Tuple{Int,Int}[]
            @test resolve_core(packages, conflicts) == [[]]
        end
        @testset "1 pkgs" begin
            for V = 1:3
                packages = fill("Only", V)
                conflicts = Tuple{Int,Int}[]
                @test resolve_core(packages, conflicts) == [[1]]
            end
        end
    end
end

@testset "resolver API" begin
    @testset "example: 0 pkgs" begin
        # also tests that pacakge and version types flow through
        for (P, V) in [(String, VersionNumber), (Real, Complex{Int})]
            versions = Dict{P,Vector{V}}()
            required = P[]
            dp = make_deps(versions)
            resolved = Resolver.resolve(dp, required)
            @test resolved == [[]]
            @test eltype(resolved) == Vector{Pair{P,V}}
        end
    end

    @testset "example: 1 pkgs" begin
        versions = Dict("A" => 1:2)
        required = ["A"]
        dp = make_deps(versions)
        resolved = Resolver.resolve(dp, required)
        @test resolved == [["A" => 1]]
    end

    @testset "example: 1 pkgs (VersionNumbers)" begin
        versions = Dict("A" => [v"2", v"1"])
        required = ["A"]
        dp = make_deps(versions)
        resolved = Resolver.resolve(dp, required)
        @test resolved == [["A" => v"2"]]
    end

    @testset "example: 2 pkgs, 0 conflicts" begin
        versions = Dict(
            "A" => 1:2,
            "B" => 1:2,
        )
        required = ["A"]
        dp = make_deps(versions)
        resolved = Resolver.resolve(dp, required)
        @test resolved == [["A" => 1]]
    end

    @testset "example: 2 pkgs, 1 conflicts" begin
        versions = Dict(
            "A" => 1:2,
            "B" => 1:2,
        )
        conflicts = Dict(
            ("A", "B") => [(1, 1)]
        )
        required = ["A"]
        dp = make_deps(versions; conflicts)
        resolved = Resolver.resolve(dp, required)
        @test resolved == [["A" => 1]]
    end

    @testset "example: 2 pkgs, 2 deps" begin
        versions = Dict(
            "A" => 1:2,
            "B" => 1:2,
        )
        deps = Dict(
            ("A" => 1) => ["B"],
            ("A" => 2) => ["B"],
        )
        required = ["A"]
        dp = make_deps(versions; deps, conflicts)
        resolved = Resolver.resolve(dp, required)
        @test resolved == [["A" => 1, "B" => 1]]
    end

    @testset "example: 2 pkgs, 2 deps, 1 conflicts" begin
        versions = Dict(
            "A" => 1:2,
            "B" => 1:2,
        )
        deps = Dict(
            ("A" => 1) => ["B"],
            ("A" => 2) => ["B"],
        )
        conflicts = Dict(
            ("A", "B") => [(1, 1)]
        )
        required = ["A"]
        dp = make_deps(versions; deps, conflicts)
        resolved = Resolver.resolve(dp, required)
        @test resolved == [
            ["A" => 1, "B" => 2],
            ["A" => 2, "B" => 1],
        ]
    end

    @testset "example: 2 pkgs, 2 conflicts, 2 reqs" begin
        # package type must be sortable
        # version type can be anything
        versions = Dict(
            "A" => 1:2,
            "B" => 1:2
        )
        conflicts = Dict(
            ("A", "B") => [(1, 1)],
            ("B", "A") => [(2, 2)],
        )
        required = ["B", "A"]
        dp = make_deps(versions; conflicts)
        resolved = Resolver.resolve(dp, required)
        @test resolved == [
            ["A" => 1, "B" => 2],
            ["A" => 2, "B" => 1],
        ]
    end

    @testset "example: 2 pkgs, 2 conflicts, 2 reqs (weird types)" begin
        # package type must be sortable
        # version type can be anything
        versions = Dict(
            :A => [1+0im, 2+0im],
            :B => [1+0im, 2+0im]
        )
        conflicts = Dict(
            (:A, :B) => [(1+0im, 1+0im)],
            (:B, :A) => [(2+0im, 2+0im)],
        )
        required = [:B, :A]
        dp = make_deps(versions; conflicts)
        resolved = Resolver.resolve(dp, required)
        @test resolved == [
            [:A => 1+0im, :B => 2+0im],
            [:A => 2+0im, :B => 1+0im],
        ]
    end

    @testset "example: 2 pkgs, 3 deps, 5 conflicts" begin
        versions = Dict(
            "A" => 1:4,
            "B" => 1:3,
        )
        deps = Dict(
            ("A" => 1) => ["B"],
            ("A" => 2) => ["B"],
            ("B" => 3) => ["A"],
        )
        conflicts = Dict(
            ("A", "B") => [
                (1, 1)
                (2, 1)
                (3, 1)
                (4, 1)
                (4, 3)
            ]
        )
        required = ["A"]
        dp = make_deps(versions; deps, conflicts)
        resolved = Resolver.resolve(dp, required)
        @test resolved == [
            ["A" => 1, "B" => 2],
            ["A" => 3],
        ]
    end

    @testset "example: 4 pkgs, 2 deps, 2 conflicts" begin
        versions = Dict(
            "A" => 1:2,
            "B" => 1:2,
            "C" => 1:1,
            "D" => 1:1,
        )
        deps = Dict(
            ("A" => 1) => ["B", "C", "D"],
            ("A" => 2) => ["B", "C", "D"],
        )
        conflicts = Dict(
            ("A", "C") => [(1, 1)],
            ("B", "D") => [(1, 1)],
        )
        required = ["A"]
        dp = make_deps(versions; deps, conflicts)
        resolved = Resolver.resolve(dp, required)
        @test resolved == [
            ["A" => 2, "B" => 2, "C" => 1, "D" => 1],
        ]
    end

    @testset "example: 3 pkgs, 2 deps, 2 conflicts (incomplete)" begin
        versions = Dict(
            "A" => 1:1,
            "B" => 1:1,
            "C" => 1:2,
        )
        deps = Dict(
            ("A" => 1) => ["C"],
            ("B" => 1) => ["C"],
        )
        conflicts = Dict(
            ("A", "C") => [(1, 1)],
            ("B", "C") => [(1, 2)],
        )
        required = ["A", "B"]
        dp = make_deps(versions; deps, conflicts)
        resolved = Resolver.resolve(dp, required)
        @test resolved == [
            ["A" => 1, "C" => 2],
            ["B" => 1, "C" => 1],
        ]
    end

    @testset "example: 3 pkgs, 2 deps, 1 conflicts (incompatible)" begin
        versions = Dict(
            "A" => 1:1,
            "B" => 1:1,
            "C" => 1:1,
        )
        deps = Dict(
            ("A" => 1) => ["B"],
            ("B" => 1) => ["C"],
        )
        conflicts = Dict(
            ("A", "C") => [(1, 1)],
        )
        required = ["A"]
        dp = make_deps(versions; deps, conflicts)
        resolved = Resolver.resolve(dp, required)
        @test resolved == [["A" => 1, "B" => 1, "C" => 1]]
    end

    @testset "example: 3 pkgs, 2 deps, 1 conflicts (incompatible)" begin
        versions = Dict(
            "A" => 1:1,
            "B" => 1:1,
            "C" => 1:1,
        )
        deps = Dict(
            ("A" => 1) => ["B", "C"],
        )
        conflicts = Dict(
            ("B", "C") => [(1, 1)],
        )
        required = ["A"]
        dp = make_deps(versions; deps, conflicts)
        resolved = Resolver.resolve(dp, required)
        @test resolved == [["A" => 1, "B" => 1, "C" => 1]]
    end

    @testset "example: 3 pkgs, 2 deps, 3 conflicts (incompatible)" begin
        versions = Dict(
            "A" => 1:1,
            "B" => 1:1,
            "C" => 1:1,
        )
        deps = Dict(
            ("A" => 1) => ["B", "C"],
        )
        conflicts = Dict(
            ("A", "B") => [(1, 1)],
            ("A", "C") => [(1, 1)],
            ("B", "C") => [(1, 1)],
        )
        required = ["A"]
        dp = make_deps(versions; deps, conflicts)
        resolved = Resolver.resolve(dp, required)
        @test resolved == [["A" => 1, "B" => 1, "C" => 1]]
    end

    @testset "example: 3 pkgs, 2 deps, 2 conflicts (incompatible)" begin
        versions = Dict(
            "A" => 1:2,
            "B" => 1:1,
            "C" => 1:1,
        )
        deps = Dict(
            ("A" => 1) => ["B"],
            ("A" => 2) => ["C"],
        )
        conflicts = Dict(
            ("A", "B") => [(1, 1)],
            ("A", "C") => [(2, 1)],
        )
        required = ["A"]
        dp = make_deps(versions; deps, conflicts)
        resolved = Resolver.resolve(dp, required)
        @test resolved == [
            ["A" => 1, "B" => 1],
            ["A" => 2, "C" => 1],
        ]
    end
end
