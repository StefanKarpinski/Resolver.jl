using Resolver
using Test

@testset "simple, 1-conflict example" begin
    required = ["A"]
    versions = Dict(
        "A" => [v"2", v"1"],
        "B" => [v"2", v"1"],
    )
    dependencies = Dict(
        PkgVer("A", v"1") => ["B"],
        PkgVer("A", v"2") => ["B"],
    )
    conflicts = Set([
        (PkgVer("A", v"2"), PkgVer("B", v"2")),
    ])
    @test resolve(required, versions, dependencies, conflicts) == [
        ["A" => v"2", "B" => v"1"],
        ["A" => v"1", "B" => v"2"],
    ]
end

@testset "3-conflict example" begin
    required = ["A"]
    versions = Dict(
        "A" => [v"2", v"1"],
        "B" => [v"2", v"1"],
        "C" => [v"2", v"1"],
        "D" => [v"2", v"1"],
        "E" => [v"2", v"1"],
    )
    dependencies = Dict(
        PkgVer("A", v"1") => ["B"],
        PkgVer("A", v"2") => ["B", "D"],
        PkgVer("B", v"1") => ["C"],
        PkgVer("B", v"2") => ["C"],
        PkgVer("C", v"2") => ["E"],
    )
    conflicts = Set([
        (PkgVer("A", v"2"), PkgVer("B", v"2")),
        (PkgVer("A", v"2"), PkgVer("C", v"2")),
        (PkgVer("B", v"2"), PkgVer("C", v"2")),
    ])
    @test resolve(required, versions, dependencies, conflicts) == [
        ["A" => v"2", "B" => v"1", "C" => v"1", "D" => v"2"],
        ["A" => v"1", "B" => v"2", "C" => v"1"],
        ["A" => v"1", "B" => v"1", "C" => v"2", "E" => v"2"],
    ]
end

#=
A B

A[1] B
A[2] B C
A[3] B C
A[4] B D

B[1] C
B[2] C D E
B[3] A D
B[4] A C D

C[1] A
C[2] D
C[3] B D
C[4] D

D[1] C
D[2]
D[3] B
D[4] E

E[1]
E[2]
=#

#=
m1 = a2 b c
m2 = a2 b c d
m3 = a1 b

a1: B
a2: B C
=#

