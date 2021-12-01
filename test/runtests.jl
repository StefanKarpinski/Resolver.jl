include("setup.jl")

@testset "core solver" begin
    @testset "basic example" begin
        packages = [1:2, 3:6, 7:10]
        conflicts = [(1,3), (1,7), (3,7), (2, 10)]
        solutions = [[1, 4, 8], [3, 2, 8], [7, 2, 4]]
        @test solutions == resolve(packages, conflicts, sort=false)
        @test solutions == resolve(packages, conflicts, Block=UInt8, sort=false)
        sort!(map(sort!, solutions))
        @test solutions == resolve(packages, conflicts)
        @test solutions == resolve(packages, conflicts, Block=UInt8)
        @test solutions == resolve_brute_force(packages, conflicts)
    end

    @testset "comprehensive tests" begin
        Random.seed!(0xed3c5f374cf4319a)
        for N = 2:5, # number of packages
            V = 1:5  # number of versions per package
            packages = [i*V+1:(i+1)*V for i=0:N-1]
            T = V^2*(N*(N-1)÷2)
            T ≤ 128 || continue
            for C in (T ≤ 12 ? (0:2^T-1) : [rand(UInt128) for _ = 1:2^12])
                conflicts = gen_conflicts(N, V, C)
                @assert length(conflicts) == count_ones(C % UInt128(2)^T)
                solutions = resolve_brute_force(packages, conflicts)
                @test solutions == resolve(packages, conflicts)
                @test solutions == resolve(packages, conflicts; Block = UInt8)
                resolved = resolve(packages, conflicts)
                # if solutions != resolve(packages, conflicts)
                #     resolved = resolve(packages, conflicts)
                #     @show packages conflicts solutions resolved
                #     return
                # end
            end
        end
    end
end
