# Does a proof account for every version it splits on? (issue #97)
#
# A line that splits a package into cases has to cover it: if the cases leave a
# gap, the versions in the gap are ones the reader is never told the fate of,
# and the argument does not close. Overlap is harmless -- two cases may both
# rule out a version -- but a gap is not.
#
# Under resolution this is not a property to be arranged, it is what the rule
# does. Resolving on a package intersects that package's literal, so the
# versions ruled out accumulate as a union and the derivation reaches ⊥ only
# when nothing is left. A gap means it never gets there.
#
# So what is checked here is that no printed line is false: a gap makes the
# elimination claim more than its inputs give, and a claim more than they give
# is a statement the universe does not support. That covers the old question
# and every other way a line could be wrong.

using Resolver: SAT, Problem, Diagnosis, resolve, pkg_info, prepare_pkg_info
using Resolver.Clauses: Clauses, Clause, packages, isbottom, literal
using Resolver.Diagnostics: Diagnostics, Conflict, Line, clause_versions
@isdefined(ProofCheck) || include(joinpath(@__DIR__, "proof_check.jl"))
using .ProofCheck

function coverage_gaps(data, prob::Problem)
    d = resolve(data, prob)
    d isa Diagnosis || return d, String[]
    sat = SAT(prepare_pkg_info(pkg_info(data, prob), prob))
    return d, reduce(vcat, [lines_are_true(sat, prob, printed_lines(c))
                            for c in d.conflicts]; init = String[])
end

@testset "chain coverage: every case split covers" begin
    for (data, prob) in (
        (varied_bound, Problem([:A]; compat = Dict(:C => [:c1]))),
        (two_hops,     Problem([:A]; compat = Dict(:C => [:c1]))),
        (weak_bound,   Problem([:P, :S]; compat = Dict(:W => [:w2]))),
        (shared_bound, Problem([:A, :B];
            compat = Dict(:A => [:a2], :B => [:b2], :C => [:c1]))),
        (late_speaker, Problem([:A]; compat = Dict(:C => [:c1]))),
        (whole_run,    Problem([:A]; compat = Dict(:B => [:b1]))),
        (narrowed_run, Problem([:A]; compat = Dict(:A => [:a2, :a1],
                                                   :B => [:b1]))),
        (two_conflicts, Problem([:A, :B, :E, :F])),
        (conflict_path, Problem([:A, :B, :C, :D])),
    )
        _, gaps = coverage_gaps(data, prob)
        @test isempty(gaps) || (@show gaps; false)
    end

    # ... and so do the registry-scale ones the diagnosis tests already use
    rp = registry.provider()
    for prob in (
        Problem(["DataFrames"]; compat = Dict(
            "DataFrames" => VersionSpec("1"),
            "PrettyTables" => VersionSpec("0.9"))),
        Problem(["Plots"]; compat = Dict(
            "Plots" => VersionSpec("1.40"),
            "RecipesPipeline" => VersionSpec("0.5"))),
        Problem(["Flux"]; compat = Dict(
            "Flux" => VersionSpec("0.16.11"),
            "Tables" => VersionSpec("0.2"))),
    )
        _, gaps = coverage_gaps(rp, prob)
        @test isempty(gaps) || (@show gaps; false)
    end

    # The case #100 was filed over. A walk used to narrow every package it
    # forced by every *other* forced package's bounds, so the run a line spoke
    # for was smaller than the bound the line before it stated, and
    # StructArrays 0.6.11–0.6.18 was accounted for nowhere. Resolution cannot
    # produce that: what a resolvent says about a package is exactly the
    # intersection of what its parents said.
    d, gaps = coverage_gaps(rp, Problem(["Zygote"]; compat = Dict(
        "Zygote" => VersionSpec("0.7.12"), "TableTraits" => VersionSpec("0.3"))))
    @test isempty(gaps) || (@show gaps; false)
    c = only(d.conflicts)
    @test c.reqs == ["Zygote"]
    named = Set(p for l in c.lines for p in packages(l.clause))
    @test "Zygote" in named && "TableTraits" in named
    # the contradiction is reached, and it is not one of the printed lines
    @test !any(l -> isbottom(l.clause), c.lines)
end

# Eliminating a package has to account for every version of it. Pairs are not
# enough to do that: a literal is a *set*, so three clauses can leave a package
# nothing to be while every two of them still agree on something -- and the
# pairwise resolvents then say strictly less than what they came from, which is
# a gap by another name.
@testset "chain coverage: eliminating covers the package" begin
    P = String
    # `A` in one of `vs` (of three versions), or `q` in `{1}`
    say(vs, q) = Line{P}(Clauses.clause([
        "A" => literal(3, vs), q => literal(2, [1])]), P[], false)
    # no version of A is in all three, but every two of them overlap
    cs = Line{P}[say([1, 2], "X"), say([2, 3], "Y"), say([1, 3], "Z")]
    out = Diagnostics.project(cs, Set{P}(["X", "Y", "Z"]))
    @test out !== nothing
    # what survives has to say that one of the three is installed; pairwise
    # elimination cannot reach it, and dropping A without it would be a gap
    @test any(l -> Set(packages(l.clause)) == Set(["X", "Y", "Z"]), out)
    @test !any(l -> "A" in packages(l.clause), out)
end
