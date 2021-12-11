include("setup.jl")

@testset "core solver" begin
    @testset "basic example" begin
        packages = [1:2, 3:6, 7:10]
        conflicts = [(1,3), (1,7), (3,7), (2,10)]
        solutions = [[1, 4, 8], [2, 3, 8], [2, 4, 7]]
        @test solutions == resolve_brute_force(packages, conflicts)
        @test solutions == resolve_core(packages, conflicts)
    end

    @testset "comprehensive tests" begin
        Random.seed!(0x53f3ce0656b85450bbb52b15fc58853f)
        count = zeros(Int, 12)
        for M = 2:5, # number of packages
            V = 1:5  # number of versions per package
            packages = [i*V+1:(i+1)*V for i=0:M-1]
            T = V^2*(M*(M-1)÷2)
            T ≤ 128 || continue
            for C in (T ≤ 12 ? (0:2^T-1) : [randu128(3-(k%5)%3) for k=1:2^12])
                conflicts = gen_conflicts(M, V, C)
                @assert length(conflicts) == count_ones(C % UInt128(2)^T)
                solutions = resolve_brute_force(packages, conflicts)
                resolved = resolve_core(packages, conflicts)
                @test resolved isa Vector{Vector{Int}}
                @test solutions == resolved
                count[length(solutions)+1] += 1
            end
        end
        # specific to chosen seed (chosen to have some large solution sets)
        @test count == [6599, 22763, 11355, 4055, 1375, 407, 98, 22, 7, 0, 1, 0]
    end
end
