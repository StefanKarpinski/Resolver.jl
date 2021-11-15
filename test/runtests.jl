include("setup.jl")

@testset "core solver" begin
    @testset "basic example" begin
        packages = [1:2, 3:6, 7:10]
        conflicts = [(1,3), (1,7), (3,7), (2, 10)]
        solutions = [[1, 4, 8], [3, 2, 8], [7, 2, 4]]
        @test solutions == resolve(packages, conflicts)
        @test solutions == resolve(packages, conflicts, Block=UInt8)
    end


end
