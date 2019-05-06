using Resolver
using Test

@testset "0-conflict examples" begin
    versions = Dict(
        "A" => [v"2", v"1"],
        "B" => [v"2", v"1"],
    )
    dependencies = Dict(
        ("A" => v"1") => String[],
        ("A" => v"2") => ["B"],
    )
    @test resolve(String[], versions, dependencies) |> isempty
    @test resolve(["A"], versions, dependencies) == [
        ["A" => v"2", "B" => v"2"],
    ]
    @test resolve(["B"], versions, dependencies) == [
        ["B" => v"2"],
    ]
    @test resolve(["A", "B"], versions, dependencies) == [
        ["A" => v"2", "B" => v"2"],
    ]
    dependencies["B" => v"2"] = ["A"]
    @test resolve(["A"], versions, dependencies) == [
        ["A" => v"2", "B" => v"2"],
    ]
    @test resolve(["B"], versions, dependencies) == [
        ["B" => v"2", "A" => v"2"],
    ]
    @test resolve(["A", "B"], versions, dependencies) == [
        ["A" => v"2", "B" => v"2"],
    ]
end

@testset "1-conflict examples" begin
    versions = Dict(
        "A" => [v"2", v"1"],
        "B" => [v"2", v"1"],
    )
    dependencies = Dict(
        ("A" => v"1") => ["B"],
        ("A" => v"2") => ["B"],
    )
    conflicts = Set([
        ("A" => v"2", "B" => v"2"),
    ])
    @test resolve(["A"], versions, dependencies, conflicts) == [
        ["A" => v"2", "B" => v"1"],
        ["A" => v"1", "B" => v"2"],
    ]
    @test resolve(["B"], versions, dependencies, conflicts) == [
        ["B" => v"2"],
    ]
    @test resolve(["A", "B"], versions, dependencies, conflicts) == [
        ["A" => v"2", "B" => v"1"],
        ["A" => v"1", "B" => v"2"],
    ]
    dependencies = Dict(
        ("A" => v"2") => ["B"],
        ("B" => v"1") => ["A"],
    )
    @test_broken resolve(["A"], versions, dependencies, conflicts) == [
        ["A" => v"2", "B" => v"1"],
        ["A" => v"1"],
    ]
    @test resolve(["B"], versions, dependencies, conflicts) == [
        ["B" => v"2"],
    ]
    @test resolve(["A", "B"], versions, dependencies, conflicts) == [
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
        ("A" => v"1") => ["B"],
        ("A" => v"2") => ["B", "D"],
        ("B" => v"1") => ["C"],
        ("B" => v"2") => ["C"],
        ("C" => v"2") => ["E"],
    )
    conflicts = Set([
        ("A" => v"2", "B" => v"2"),
        ("A" => v"2", "C" => v"2"),
        ("B" => v"2", "C" => v"2"),
    ])
    @test resolve(required, versions, dependencies, conflicts) == [
        ["A" => v"2", "B" => v"1", "C" => v"1", "D" => v"2"],
        ["A" => v"1", "B" => v"2", "C" => v"1"],
        ["A" => v"1", "B" => v"1", "C" => v"2", "E" => v"2"],
    ]
end
