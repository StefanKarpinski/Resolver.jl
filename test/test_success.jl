# Tests for the @test_success macro (test/TestSuccess.jl).

using Test

@isdefined(TestSuccess) || include(joinpath(@__DIR__, "TestSuccess.jl"))
using .TestSuccess

module TestSuccessTests

using Test
using ..TestSuccess

export collect_results, only_result, message, child

# A testset that keeps its results instead of reporting them, so that a failing
# test can be inspected as a value rather than failing the run.
mutable struct Collector <: Test.AbstractTestSet
    results :: Vector{Test.Result}
end
Collector(::AbstractString; kws...) = Collector(Test.Result[])
Test.record(ts::Collector, res::Test.Result) = (push!(ts.results, res); res)
Test.finish(ts::Collector) = ts # deliberately not recorded in the parent

collect_results(body::Function) =
    (@testset Collector "collected" begin body() end).results
only_result(body::Function) = only(collect_results(body))

# what a failure looks like in the log
message(res::Test.Result) = sprint(show, res)

# a child process that runs `code`, portably (no shell)
child(code::AbstractString) = `$(Base.julia_cmd()[1]) --startup-file=no -e $code`

end # module

using .TestSuccessTests

@testset "@test_success" begin
    @testset "a command that exits 0 passes, quietly" begin
        res, echoed = mktemp() do log, io
            res = redirect_stdout(io) do
                only_result() do
                    @test_success child("print(\"loud\"); flush(stdout)")
                end
            end
            close(io)
            res, read(log, String)
        end
        @test res isa Test.Pass
        @test isempty(echoed) # the child's output is not echoed
    end

    @testset "a nonzero exit fails, reporting stdout" begin
        res = only_result() do
            @test_success child("print(\"spoken on stdout\"); exit(3)")
        end
        @test res isa Test.Fail
        msg = message(res)
        @test occursin("Exit code: 3", msg)
        @test occursin("spoken on stdout", msg)
    end

    @testset "a nonzero exit fails, reporting stderr" begin
        res = only_result() do
            @test_success child("print(stderr, \"spoken on stderr\"); exit(4)")
        end
        @test res isa Test.Fail
        msg = message(res)
        @test occursin("Exit code: 4", msg)
        @test occursin("spoken on stderr", msg)
    end

    @testset "the report shows the command as run, not as written" begin
        code = "exit(5)"
        res = only_result() do
            @test_success child(code)
        end
        msg = message(res)
        @test occursin("Expression: child(code)", msg) # the source expression
        @test occursin("exit(5)", msg)                 # what it expanded to
    end

    @testset "both streams, in the order the command wrote them" begin
        res = only_result() do
            @test_success child("""
                print(stdout, "first\\n");  flush(stdout)
                print(stderr, "second\\n"); flush(stderr)
                print(stdout, "third\\n");  flush(stdout)
                exit(1)
                """)
        end
        msg = message(res)
        at(s) = first(something(findfirst(s, msg)))
        @test at("first") < at("second") < at("third")
    end

    @testset "output past the limit is truncated, from the front" begin
        res = only_result() do
            @test_success child("""
                for i = 1:2000; println("line \$i"); end
                exit(1)
                """) limit = 256
        end
        @test res isa Test.Fail
        msg = message(res)
        @test occursin(r"\d+ truncated", msg)
        @test occursin("| line 2000", msg)  # the tail is what is kept
        @test !occursin("| line 1\n", msg)  # the head is what is dropped
    end

    @testset "a command that cannot start fails, and does not throw" begin
        res = only_result() do
            @test_success `no-such-program-9f3a1c`
        end
        @test res isa Test.Fail
        msg = message(res)
        @test occursin("could not run", msg)
        @test occursin("no-such-program-9f3a1c", msg)
    end
end
