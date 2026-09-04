# The statements a proof is made of (src/Clauses.jl).
#
# Packages are letters and versions are integers throughout: with integers for
# both, `1 1 requires 2 1` is unreadable and a wrong assertion looks right.

using Resolver.Clauses
using Resolver.Clauses: clause, Lit, absent, present, version_order

const L = Clauses.literal
# three packages: A with versions 10,20,30; B with 5,6; C with 1,2
const VERS = Dict(1 => [10, 20, 30], 2 => [5, 6], 3 => [1, 2])
vers(p) = VERS[p]
say(c) = Clauses.clause_phrase(c, vers, Clauses.letters)

@testset "clauses: normal form" begin
    # a package covering everything makes the clause say nothing at all
    @test clause([1 => L(3, [1,2,3]; absent = true)]) === nothing
    @test clause([1 => L(3, [1]), 2 => L(2, [1,2]; absent = true)]) === nothing
    # a disjunct nothing can satisfy is not carried around
    c = clause([1 => L(3, [1]), 2 => L(2, Int[])])
    @test packages(c) == [1]
    # one entry per package, in package order, so equal clauses are ==
    @test clause([2 => L(2, [1]), 1 => L(3, [1])]) ==
          clause([1 => L(3, [1]), 2 => L(2, [1])])
    # the empty clause is the contradiction
    @test isbottom(clause(Pair{Int,Lit}[]))
    @test !isbottom(clause([1 => L(3, [1])]))
end

@testset "clauses: one rule" begin
    # intersect on the shared package, union on the rest
    c = clause([1 => L(3, [1,2]), 2 => L(2, [1])])
    d = clause([1 => L(3, [2,3]), 3 => L(2, [1])])
    r = resolve_on(c, d, 1)
    @test r[1] == L(3, [2])
    @test r[2] == L(2, [1])
    @test r[3] == L(2, [1])

    # an intersection no smaller than a parent is that parent, weakened: no
    # progress, so the rule declines it
    @test resolve_on(clause([1 => L(3, [1,2])]), clause([1 => L(3, [1,2,3])]), 1) ===
          nothing
    # nothing to resolve on
    @test resolve_on(clause([1 => L(3, [1])]), clause([2 => L(2, [1])]), 1) === nothing
    # disjoint literals on one package leave nothing: the contradiction
    @test isbottom(resolve_on(clause([1 => L(3, [1])]), clause([1 => L(3, [2,3])]), 1))
    # a resolvent that would be a tautology is not a statement
    # the two halves of B cover it between them, so what comes out says nothing
    @test resolve_on(clause([1 => L(3, [1,2]), 2 => L(2, [1]; absent = true)]),
                     clause([1 => L(3, [2,3]), 2 => L(2, [2])]), 1) === nothing
end

@testset "clauses: direction is not part of the statement" begin
    # `A@10 requires B@5` and its contrapositive are one object, so a rule
    # written for either reaches both
    fwd = clause([1 => L(3, [1], true; absent = true), 2 => L(2, [1])])
    bwd = clause([2 => L(2, [1]), 1 => L(3, [1], true; absent = true)])
    @test fwd == bwd
    # ... and which way it is said is chosen when it is said
    @test say(fwd) == "A 10 requires B 5"
    @test Clauses.clause_phrase(fwd, vers, Clauses.letters; subject = 1) ==
          "B 6 or absent constrains A ≥20"
end

@testset "clauses: subsumption" begin
    @test subsumes(clause([1 => L(3, [1])]), clause([1 => L(3, [1,2])]))
    @test !subsumes(clause([1 => L(3, [1,2])]), clause([1 => L(3, [1])]))
    # a package the other does not mention is a literal it cannot cover
    @test !subsumes(clause([1 => L(3, [1])]), clause([2 => L(2, [1])]))
end

@testset "clauses: said in English" begin
    # whether the verb forces is whether ⊥ is in the consequent, which is data
    # rather than a flag a rule has to carry
    @test say(clause([1 => L(3, [1], true; absent = true), 2 => L(2, [1])])) ==
          "A 10 requires B 5"
    @test say(clause([1 => L(3, [1], true; absent = true),
                      2 => L(2, [1]; absent = true)])) == "A 10 constrains B 5"
    # a clause about one package is a verdict, said as what it rules out
    @test say(clause([1 => L(3, [1,2], true; absent = true)])) ==
          "A ≤20 cannot be installed"
    @test say(clause([1 => L(3, Int[]; absent = true)])) ==
          "no version of A can be installed"
    # ... and without ⊥ it demands the package instead
    @test say(clause([1 => L(3, [1,2,3])])) == "A must be installed"
    @test say(clause([1 => L(3, [1,2])])) == "A must be installed at ≤20"
    # more than two literals is a conjunctive antecedent
    @test say(clause([1 => L(3, [1], true; absent = true),
                      2 => L(2, [1], true; absent = true),
                      3 => L(2, [2])])) == "A 10 and B 5 together require C 2"
    # ranges read low to high, and separate runs are listed
    @test say(clause([1 => L(3, [1,3], true; absent = true)])) ==
          "A 10, 30 cannot be installed"
    @test say(clause([1 => L(3, [2], true; absent = true)])) ==
          "A 20 cannot be installed"
end

@testset "clauses: packages as letters" begin
    @test Clauses.letters(1) == "A"
    @test Clauses.letters(26) == "Z"
    @test Clauses.letters(27) == "AA"
    @test Clauses.letters("Flux") == "Flux"
end
