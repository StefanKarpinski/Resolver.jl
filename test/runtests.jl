include("setup.jl")

@testset "core resolver" begin
    @testset "basic example" begin
        packages = [fill(1, 2); fill(2, 4); fill(3, 4)]
        conflicts = [(1,3), (1,7), (3,7), (2,10)]
        solutions = [[1, 4, 8], [2, 3, 8], [2, 4, 7]]
        @test solutions == resolve_brute_force(packages, conflicts)
        @test solutions == resolve_brute_force(packages, Set(conflicts))
        @test solutions == resolve_core(packages, conflicts)
        @test solutions == resolve_core(packages, Set(conflicts))
        @testset "permutations" for _ = 1:100
            # permute the versions, conflicts and solutions
            versions′ = shuffle(packages)
            p = sortperm(versions′)
            conflicts′ = [(p[c[1]], p[c[2]]) for c in conflicts]
            solutions′ = [p[solution] for solution in solutions]
            foreach(sort!, solutions′)
            sort!(solutions′)
            @test solutions′ == resolve_brute_force(versions′, conflicts′)
            @test solutions′ == resolve_core(versions′, conflicts′)
        end
    end

    @testset "comprehensive tests" begin
        Random.seed!(0x53f3ce0656b85450bbb52b15fc58853f)
        count = zeros(Int, 10)
        for M = 2:5, # number of packages
            V = 1:5  # number of versions per package
            T = V^2*(M*(M-1)÷2)
            T ≤ 128 || continue
            packages = shuffle!([p for _=1:V for p=1:M])
            p = sortperm(packages)
            for C in (T ≤ 12 ? (0:2^T-1) : [randu128(3-(k%5)%3) for k=1:2^12])
                G = gen_conflicts(M, V, C)
                conflicts = Tuple{Int,Int}[(p[c[1]], p[c[2]]) for c in G]
                @assert length(conflicts) == count_ones(C % UInt128(2)^T)
                solutions = resolve_brute_force(packages, conflicts)
                resolved = resolve_core(packages, conflicts)
                @test resolved isa Vector{Vector{Int}}
                @test solutions == resolved
                count[length(solutions)+1] += 1
            end
        end
        # specific to chosen seed (chosen to have some large solution sets)
        @test count == [6563, 22938, 11289, 4032, 1343, 382, 104, 25, 4, 2]
    end
end

@testset "resolver API" begin
	compat((p1, v1), (p2, v2)) = !(p1 == "A" && p2 == "B" && (
        v1 in (v"3", v"1")  && v2 == v"2"         ||
        v1 == v"2"          && v2 == v"2"         ||
        v1 == v"0"          && v2 in (v"2", v"0")
    ))
	versions = Dict(
        "A" => [v"3", v"2", v"1", v"0"],
        "B" => [v"2", v"1", v"0"],
    )
	required = ["A"]
	deps = Dict(
        ("A" => v"3") => ["B"],
        ("A" => v"2") => ["B"],
        ("B" => v"1") => ["A"],
    )
	solutions = [
        ["A" => v"3.0.0", "B" => v"1.0.0"],
        ["A" => v"1.0.0"],
    ]
    resolved = Resolver.resolve(compat, versions, required, deps)
    @test typeof(resolved) == typeof(solutions)
    @test resolved == solutions
end
