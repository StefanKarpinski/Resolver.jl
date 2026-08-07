# A `@test` for subprocesses. Nothing here knows about Resolver: it uses Base
# and Test only, so it can be lifted out of this repo unchanged. Inside Test
# itself, `test_success` below would build its Pass/Fail the way the other test
# macros do; from outside, it sticks to constructor forms that have been stable
# since Julia 1.9.

"""
Defines [`@test_success`](@ref), a `@test` for commands whose failures report
the command, how it exited and what it printed.
"""
module TestSuccess

using Test

export @test_success

# how many bytes of a failed command's output to report, by default
const OUTPUT_LIMIT = 8 * 1024

# Test aligns the labels of a failure report in a 14-column gutter
# ("  Expression: "); the lines we add line up with them.
label(name::AbstractString) = lpad(name * ":", 13) * " "
const INDENT = " "^14

"""
    @test_success cmd
    @test_success cmd limit = 8192

Test that `cmd` runs and exits with status 0 — the same contract as
`@test success(cmd)`, but a failure says why: it reports the command as it was
actually run, how it exited, and what it printed.

Output is captured rather than shown, so a command that succeeds is silent.
Standard output and standard error are captured together, interleaved in the
order the command wrote them. A failure reports the last `limit` bytes of that
— the end being where an explanation usually is — and says how many bytes it
dropped to fit. A stream `cmd` redirects for itself goes where `cmd` sends it,
and only the other one is captured.

A command that cannot be started at all — a misspelled program, say — fails the
test and reports why, rather than throwing.

The test is recorded in the enclosing testset like any other, counting as one
`Pass` or `Fail` in its summary.
"""
macro test_success(cmd, kws...)
    limit = OUTPUT_LIMIT
    for kw in kws
        Meta.isexpr(kw, :(=), 2) && kw.args[1] === :limit ||
            throw(ArgumentError("@test_success: expected `limit = n`, got `$kw`"))
        limit = kw.args[2]
    end
    quote
        test_success($(esc(cmd)), $(string(cmd)), $(esc(limit)),
                     $(QuoteNode(__source__)))
    end
end

# Run `cmd` with its output captured; return `nothing` if it exited 0, and
# otherwise the report of what happened instead.
#
# Both streams are pointed at one file opened in append mode, so the kernel
# interleaves them in write order exactly as a shared terminal would, and
# however much the command writes is bounded by disk rather than by memory.
# Only the part of it that gets reported is ever read back in.
function capture(cmd, limit::Integer)
    path, io = mktemp()
    close(io) # the command opens the file itself; we only wanted the name
    try
        redirected = pipeline(cmd; stdout = path, stderr = path, append = true)
        proc = try
            run(redirected, wait = false)
        catch err
            err isa Base.IOError || rethrow()
            return string("could not run ", cmd, "\n",
                          label("Error"), sprint(showerror, err))
        end
        wait(proc)
        success(proc) && return nothing
        return string(cmd, "\n",
                      label("Exit code"), status(proc), "\n",
                      label("Output"), report_output(path, limit))
    finally
        rm(path, force = true)
    end
end

status(proc) = proc.termsignal != 0 ?
    string("killed by signal ", proc.termsignal) : string(proc.exitcode)

# The tail of the captured output, each line indented and marked as the
# command's, headed by a note of how much was dropped to fit in `limit`.
function report_output(path::AbstractString, limit::Integer)
    total = filesize(path)
    total == 0 && return "(none)"
    dropped = max(total - limit, 0)
    bytes = open(path) do io
        seek(io, dropped)
        read(io)
    end
    if dropped > 0
        # start at a line boundary: the cut lands mid-line, and as easily
        # mid-character, so the leading partial line is worth less than the
        # confusion it causes
        nl = findfirst(==(UInt8('\n')), bytes)
        if nl !== nothing && nl < length(bytes)
            dropped += nl
            bytes = bytes[nl+1:end]
        end
    end
    header = dropped == 0 ? "($total bytes)" :
        "(last $(total - dropped) of $total bytes, $dropped truncated)"
    lines = split(chomp(String(bytes)), '\n')
    return header * join("\n" * INDENT * "| " * rstrip(line, '\r') for line in lines)
end

function test_success(cmd, expr::AbstractString, limit::Integer,
                      source::LineNumberNode)
    report = capture(cmd, limit)
    result = report === nothing ?
        Test.Pass(:test, expr, nothing, true, source, false) :
        Test.Fail(:test, expr, report, "false", nothing, source, false)
    Test.record(Test.get_testset(), result)
end

end # module
