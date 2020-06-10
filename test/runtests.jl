using Resolver
using Test

@testset "0-conflict examples" begin
    deps = Dict("A1,B1,B2" => "", "A2" => "B")
    V, N = graph(deps)
    # resolve(V, N, "") == [""]
    # resolve(V, N, "A") == ["A2,B2"]
    # resolve(V, N, "B") == ["B2"]
    # resolve(V, N, "A,B") == ["A2,B2"]
    deps = Dict("A1,B1" => "", "A2" => "B", "B2" => "A")
    V, N = graph(deps)
    # resolve(V, N, "A") == ["A2,B2"]
    # resolve(V, N, "B") == ["A2,B2"]
    # resolve(V, N, "A,B") == ["A2,B2"]
end

@testset "1-conflict examples" begin
    deps = Dict("A1,A2" => "B", "B1,B2" => "")
    conflicts = [("A2", "B2")]
    V, N = graph(deps, conflicts)
    # resolve(V, N, "A") == ["A2,B1", "A1,B2"]
    # resolve(V, N, "B") == ["B2"]
    # resolve(V, N, "A,B") == ["A2,B1", "A1,B2"]
    deps = Dict("A2" => "B", "B1" => "A", "A1,B2" => "")
    V, N = graph(deps, conflicts)
    # resolve(V, N, "A") == ["A2,B1", "A1"]
    # resolve(V, N, "B") == ["B2"]
    # resolve(V, N, "A,B") == ["A2,B1", "A1,B2"]
end

@testset "3-conflict example" begin
    deps = Dict(
        "A1" => "B",
        "A2" => "B,D",
        "B1,B2" => "C",
        "C1" => "",
        "C2" => "E",
    )
    conflicts = [("A2", "B2,C2"), ("B2","C2")]
    V, N = graph(deps, conflicts)
    # resolve(V, N, "A") == ["A2,B1,C1,D2", "A1,B1,C2,E2"]
end

@testset "2-conflict overcommit example" begin
    deps = Dict("A1,A2" => "B,C,D", "B1,B2,C1,D1" => "")
    conflicts = [("A2", "C1"), ("B2", "D1")]
    V, N = graph(deps, conflicts)
    # resolve(V, N, "A") == ["A1,B1,C1,D1"]
end
