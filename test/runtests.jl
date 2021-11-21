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
        for N = 2:5, # number of packages
            V = 1:5  # number of versions per package
            packages = [i*V+1:(i+1)*V for i=0:N-1]
            if N + V <= 5
                for C = 0:2^((N*(N-1)÷2)*V^2)-1
                    conflicts = Tuple{Int,Int}[]
                    s = trailing_zeros(C)
                    while (C >> s) ≠ 0
                        b = s += trailing_zeros(C >> s)
                        b, v₂ = divrem(b, V)
                        b, v₁ = divrem(b, V)
                        p₂, p₁ = divrem(b, N-1)
                        p₂ += 1 + p₁
                        v₁ += 1 + p₁*V
                        v₂ += 1 + p₂*V
                        push!(conflicts, (v₁, v₂))
                        s += 1
                    end
                    @show conflicts
                    solutions = resolve_brute_force(packages, conflicts)
                    @test solutions ⊆ resolve(packages, conflicts)
                    @test solutions ⊆ resolve(packages, conflicts; Block = UInt8)
                end
            else
                # try random conflicts
                conflicts = Tuple{Int,Int}[]
                for _ = 1:2^12

                end
            end
        end
    end
end
