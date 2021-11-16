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
                    for v₁ = 1:V, v₂ = 1:V, p₁ = 1:N-1, p₂ = p₁+1:N
                        s  = v₂-1; s *= V
                        s += v₁-1; s *= V
                        s += p₁-1; s *= N-1
                        s += p₂-p₁-1
                        isodd(C >> s) || continue
                        push!(conflicts, ((p₁-1)*V + v₁, (p₂-1)*V + v₂))
                    end
                    @show conflicts
                    solutions = resolve_brute_force(packages, conflicts)
                    @test solutions == resolve(packages, conflicts)
                    @test solutions == resolve(packages, conflicts; Block = UInt8)
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
