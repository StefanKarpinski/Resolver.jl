using Resolver
using Test

@testset "0-conflict examples" begin
    deps = Dict("A1,B1,B2" => "", "A2" => "B")
    resolve("", deps) # == [""]
    resolve("A", deps) # == ["A2,B2"]
    resolve("B", deps) # == ["B2"]
    resolve("A,B", deps) # == ["A2,B2"]
    deps = Dict("A1,B1" => "", "A2" => "B", "B2" => "A")
    resolve("A", deps) # == ["A2,B2"]
    resolve("B", deps) # == ["A2,B2"]
    resolve("A,B", deps) # == ["A2,B2"]
end

@testset "1-conflict examples" begin
    deps = Dict("A1,A2" => "B", "B1,B2" => "")
    conflicts = [("A2", "B2")]
    resolve("A", deps, conflicts) # == ["A2,B1", "A1,B2"]
    resolve("B", deps, conflicts) # == ["B2"]
    resolve("A,B", deps, conflicts) # == ["A2,B1", "A1,B2"]
    deps = Dict("A2" => "B", "B1" => "A", "A1,B2" => "")
    resolve("A", deps, conflicts) # == ["A2,B1", "A1"]
    resolve("B", deps, conflicts) # == ["B2"]
    resolve("A,B", deps, conflicts) # == ["A2,B1", "A1,B2"]
end

@testset "3-conflict example" begin
    deps = Dict(
        "A1" => "B",
        "A2" => "B,D",
        "B1,B2" => "C",
        "C2" => "E",
        "C1,D1,D2,E1,E2" => "",
    )
    conflicts = [("A2", "B2,C2"), ("B2","C2")]
    resolve("A", deps, conflicts) # == ["A2,B1,C1,D2", "A1,B1,C2,E2"]
end

@testset "2-conflict overcommit example" begin
    deps = Dict("A1,A2" => "B,C,D", "B1,B2,C1,D1" => "")
    conflicts = [("A2", "C1"), ("B2", "D1")]
    resolve("A", deps, conflicts) # == ["A1,B1,C1,D1"]
end
