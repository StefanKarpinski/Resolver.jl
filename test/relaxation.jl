# The substrate a diagnosis is assembled from: the raw literals for the two
# things it puts to the solver, the lift that makes every class of the universe
# choosable again, and the attribution that says why a class is empty. Nothing
# here is a diagnosis, and nothing outside this file calls any of it.
#
# One property ties the three together. A constraint forbids *versions*, and a
# class is empty when every member of it is forbidden; emptiness is what the
# solver sees, as `reps[p][c] == 0` and one forbidding unit clause. So:
#
#     reps[p][c] == 0  ⟺  every member of class c has an excluding source
#
# Attribution is a statement about constraints and versions, deactivation is a
# statement about the instance, and the sweeps below are that they are the same
# statement — checked over the tiny data grids crossed with randomized compat
# bounds, pins, an admission kind, and every combination of the three.

using Resolver: PkgData, Problem, SAT, PicoSAT, pkg_info, nclasses, NoKinds,
    rank_pkg_info, finalize, is_satisfiable, is_excluded,
    exclusion_sources, class_exclusions, installed_lit, forbidden_lit,
    with_classes_relaxed, with_temp_clauses, sat_assume_var, sat_push

const NO_DEPS = Dict{Symbol,Vector{Symbol}}()
const NO_COMP = Dict{Symbol,Dict{Symbol,Vector{Symbol}}}()

## references, none of which asks the solver anything

# the sources excluding a version, recomputed straight from the problem
function ref_sources(prob::Problem{P}, p::P, v) where {P}
    srcs = Symbol[]
    haskey(prob.compat, p) && v ∉ prob.compat[p] && push!(srcs, :compat)
    haskey(prob.pins, p) && v ≠ prob.pins[p] && push!(srcs, :pin)
    for (kind, forbids) in prob.excludes
        forbids(p, v) && push!(srcs, kind)
    end
    return srcs
end

# every source of `prob` that could name a version of `p` at all
sources_of(prob::Problem{P}, p::P) where {P} = Symbol[
    (haskey(prob.compat, p) ? (:compat,) : ())...,
    (haskey(prob.pins, p) ? (:pin,) : ())...,
    first.(prob.excludes)...]

# `prob` with one of its sources kept and the rest dropped, so that what a
# single source excludes can be asked of `is_excluded` directly
one_source(prob::Problem, src::Symbol) =
    src == :compat ? Problem(prob.reqs; compat = prob.compat) :
    src == :pin    ? Problem(prob.reqs; pins = prob.pins) :
                     Problem(prob.reqs;
                         excludes = filter(kf -> first(kf) == src, prob.excludes))

# a random constraint set of the requested shape over `m` packages with `n`
# versions each. `random_constraints` (problem.jl) draws the compat bounds and
# pins, including entries for a package that is not in the data; the admission
# kind forbids a random set of (package, version) pairs, reaching every package
# at once the way a kind does
function random_problem(reqs, m::Int, n::Int, shape::Symbol)
    compat, pins = random_constraints(m, n)
    shape in (:compat, :all) || (compat = Dict{Int,Vector{Int}}())
    shape in (:pins,   :all) || (pins   = Dict{Int,Int}())
    banned = Set((p, v) for p = 1:m+1, v = 1:n if rand(Bool))
    excludes = shape in (:kind, :all) ?
        Pair{Symbol,Any}[:ban => (p, v) -> (p, v) in banned] : NoKinds
    return Problem(reqs; compat, pins, excludes)
end

const SHAPES = (:none, :compat, :pins, :kind, :all)

## instance-level helpers

# the instance's answers to a fixed battery of queries, every solve made under
# `assume` as well: satisfiable on its own, with each package required, with
# each class demanded, and with each class forbidden while its package is
# required. Two states of one instance that answer this identically are
# indistinguishable to anything that asks about packages and classes
function verdicts(sat::SAT{P}, assume::AbstractVector{Int} = Int[]) where {P}
    out = Bool[]
    function ask(lits::Int...)
        for l in assume
            sat_assume_var(sat, l)
        end
        for l in lits
            sat_assume_var(sat, l)
        end
        push!(out, is_satisfiable(sat))
    end
    ask()
    for p in sort!(collect(keys(sat.vars)))
        ask(installed_lit(sat, p))
        for c = 1:nclasses(sat.info[p])
            ask(-forbidden_lit(sat, p, c))
            ask(installed_lit(sat, p), forbidden_lit(sat, p, c))
        end
    end
    return out
end

# the classes this query has emptied, sorted — which is the order their
# literals are in, since variables are laid out in package-name order
empty_classes(univ) = sort!(
    [(p, c) for (p, rs) in univ.reps for c in eachindex(rs) if iszero(rs[c])])

## attribution

@testset "attribution: every source that excludes a version, and only those" begin
    prob = Problem([:A];
        compat = Dict(:A => [:v1]),
        pins   = Dict(:A => :v2),
        excludes = [:pre => (p, v) -> v == :v3, :odd => (p, v) -> p == :B])
    # each version of :A is excluded by the sources that name it, in source
    # order, and by no others
    @test exclusion_sources(prob, :A, :v1) == [:pin]
    @test exclusion_sources(prob, :A, :v2) == [:compat]
    @test exclusion_sources(prob, :A, :v3) == [:compat, :pin, :pre]
    # compat and pins name their packages; a kind reaches every package
    @test isempty(exclusion_sources(prob, :C, :v1))
    @test exclusion_sources(prob, :C, :v3) == [:pre]
    @test exclusion_sources(prob, :B, :v1) == [:odd]
    @test exclusion_sources(prob, :B, :v3) == [:pre, :odd]
    # an unconstrained problem excludes nothing at all
    @test isempty(exclusion_sources(Problem([:A]), :A, :v1))
    # ... and a source that excludes nothing names nothing
    @test isempty(exclusion_sources(Problem([:A];
        compat = Dict(:A => [:v1, :v2])), :A, :v1))
end

@testset "attribution: why a class is empty" begin
    # nothing in the registry tells :A's :v2 and :v1 apart, so they are one
    # class of two members, and two sources empty it between them: the compat
    # bound takes one member, the admission kind the other
    data = Dict(
        :A => PkgData([:v2, :v1], Dict(:v2 => [:B], :v1 => [:B]), NO_COMP),
        :B => PkgData([:w1], NO_DEPS, NO_COMP),
    )
    prob = Problem([:A]; compat = Dict(:A => [:v1]),
                         excludes = [:pre => (p, v) -> v == :v1])
    info = pkg_info(data, prob)
    @test info[:A].members == [[1, 2]]
    @test class_exclusions(prob, :A, info[:A], 1) ==
        [:v2 => [:compat], :v1 => [:pre]]
    # ... which is exactly the class the query has emptied
    @test rank_pkg_info(info, prob).reps[:A] == [0]
    # :B is not constrained, so its class holds a version nothing excludes
    @test class_exclusions(prob, :B, info[:B], 1) == [:w1 => Symbol[]]
    @test rank_pkg_info(info, prob).reps[:B] == [1]
end

@testset "a class is empty iff every member of it is excluded" begin
    Random.seed!(rand(RandomDevice(), UInt64))
    empties = classes = shared = partial = 0
    for (m, n) in ((2, 2), (2, 3), (3, 2), (3, 3), (2, 4), (4, 2), (2, 5))
        make_deps, make_comp, data, d, c = tiny_data_makers(m, n)
        for _ = 1:12
            fill_data!(m, n, make_deps(randbits(d)), make_comp(randbits(c)), data)
            reqs = collect(make_reqs(rand(1:2^m-1)))
            info = pkg_info(data, keys(data); filter = false)
            for shape in SHAPES
                prob = random_problem(reqs, m, n, shape)
                # the unfiltered universe, so that every class is still there
                univ = rank_pkg_info(info, prob)
                for (p, info_p) in univ.info, k = 1:nclasses(info_p)
                    ex = class_exclusions(prob, p, info_p, k)
                    # the pairs are the class's members, in version order
                    @test first.(ex) == info_p.versions[info_p.members[k]]
                    # the property itself
                    empty = all(!isempty(srcs) for (_, srcs) in ex)
                    @test iszero(univ.reps[p][k]) == empty
                    empties += empty
                    classes += 1
                    if length(ex) > 1
                        # a class of several members is where the property has
                        # something to say: emptying one takes every member
                        shared += empty
                        partial += !empty && any(!isempty(s) for (_, s) in ex)
                    end
                    # and a class that is not empty stands at a member nothing
                    # excludes
                    r = univ.reps[p][k]
                    if !iszero(r)
                        @test isempty(exclusion_sources(prob, p, info_p.versions[r]))
                    end
                    for (v, srcs) in ex
                        # attribution agrees with a recomputation from the
                        # problem, and with `is_excluded` on whether there is
                        # anything to name
                        @test srcs == ref_sources(prob, p, v)
                        @test isempty(srcs) == !is_excluded(prob, p, v)
                        # every source named really does exclude the version on
                        # its own, and every source that does is named
                        for src in sources_of(prob, p)
                            @test (src in srcs) ==
                                is_excluded(one_source(prob, src), p, v)
                        end
                    end
                end
            end
        end
    end
    # the sweep really did produce empty classes, and non-empty ones — and
    # both of the cases a class of several members can be in: emptied between
    # its members, and standing despite some of them being excluded
    @test 0 < empties < classes
    @test shared > 0
    @test partial > 0
end

## the literals and the lift

@testset "the frame's literals are the ones the accessors name" begin
    Random.seed!(rand(RandomDevice(), UInt64))
    deactivated = 0
    for (m, n) in ((2, 2), (2, 3), (3, 2), (3, 3))
        make_deps, make_comp, data, d, c = tiny_data_makers(m, n)
        for _ = 1:6
            fill_data!(m, n, make_deps(randbits(d)), make_comp(randbits(c)), data)
            reqs = collect(make_reqs(rand(1:2^m-1)))
            info = pkg_info(data, keys(data); filter = false)
            for shape in SHAPES
                prob = random_problem(reqs, m, n, shape)
                univ = rank_pkg_info(info, prob)
                sat = SAT(univ)
                try
                    @test all(installed_lit(sat, p) == v for (p, v) in sat.vars)
                    # the frame forbids exactly the classes the query emptied,
                    # each by the literal the accessor names
                    lits = sort!(Int[forbidden_lit(sat, p, k)
                                     for (p, k) in empty_classes(univ)]; by = abs)
                    @test sat.deact == lits
                    deactivated += length(lits)
                finally
                    finalize(sat)
                end
            end
        end
    end
    @test deactivated > 0
end

@testset "the frame lift round-trips" begin
    Random.seed!(rand(RandomDevice(), UInt64))
    lifted = freed = 0
    for (m, n) in ((2, 2), (2, 3), (3, 2), (3, 3))
        make_deps, make_comp, data, d, c = tiny_data_makers(m, n)
        for _ = 1:4
            fill_data!(m, n, make_deps(randbits(d)), make_comp(randbits(c)), data)
            reqs = collect(make_reqs(rand(1:2^m-1)))
            info = pkg_info(data, keys(data); filter = false)
            for shape in SHAPES
                prob = random_problem(reqs, m, n, shape)
                sat = SAT(rank_pkg_info(info, prob))
                try
                    before = verdicts(sat)
                    lifted += !isempty(sat.deact)
                    with_classes_relaxed(sat) do
                        # relaxed, the instance can only answer yes more often
                        relaxed = verdicts(sat)
                        @test all(before .≤ relaxed)
                        freed += relaxed ≠ before
                        # every deactivation re-imposed jointly, by assumption,
                        # states what the frame's clauses stated
                        @test verdicts(sat, sat.deact) == before
                        # ... and each one singly is weaker than all of them
                        # and stronger than none
                        for l in sat.deact
                            single = verdicts(sat, [l])
                            @test all(before .≤ single .≤ relaxed)
                        end
                    end
                    # the frame is back, so the instance answers as it did
                    @test verdicts(sat) == before
                    # and it lifts again, as many times as asked
                    with_classes_relaxed(sat) do
                        @test verdicts(sat, sat.deact) == before
                    end
                    @test verdicts(sat) == before
                finally
                    finalize(sat)
                end
            end
        end
    end
    # the sweep really did lift frames, and lifting really did make queries
    # answerable that the query's constraints had ruled out
    @test 0 < freed ≤ lifted
end

@testset "the lift is scoped, and selective inside" begin
    # :P's two classes are emptied by two different sources: the compat bound
    # takes the class holding :v3 and :v1, the admission kind the one holding
    # :v2 (which is a class of its own because it conflicts with :Q@:w2)
    data = Dict(
        :P => PkgData([:v3, :v2, :v1], NO_DEPS, Dict(:v2 => Dict(:Q => [:w1]))),
        :Q => PkgData([:w2, :w1], NO_DEPS, NO_COMP),
    )
    prob = Problem([:P, :Q];
        compat = Dict(:P => [:v2]),
        excludes = [:pre => (p, v) -> v == :v2])
    info = pkg_info(data, prob; filter = false)
    univ = rank_pkg_info(info, prob)
    @test univ.reps[:P] == [0, 0]
    @test class_exclusions(prob, :P, univ.info[:P], 1) ==
        [:v3 => [:compat], :v1 => [:compat]]
    @test class_exclusions(prob, :P, univ.info[:P], 2) == [:v2 => [:pre]]

    sat = SAT(univ)
    try
        before = verdicts(sat)
        @test !is_satisfiable(sat, [:P])
        @test with_classes_relaxed(() -> is_satisfiable(sat, [:P]), sat)
        with_classes_relaxed(sat) do
            # re-imposing what the compat emptied leaves the class the kind
            # emptied choosable, and the other way round: one literal at a time
            # is a different question from the frame
            sat_assume_var(sat, forbidden_lit(sat, :P, 1))
            @test is_satisfiable(sat, [:P])
            @test PicoSAT.deref(sat.pico, -forbidden_lit(sat, :P, 2)) > 0
            sat_assume_var(sat, forbidden_lit(sat, :P, 2))
            @test is_satisfiable(sat, [:P])
            @test PicoSAT.deref(sat.pico, -forbidden_lit(sat, :P, 1)) > 0
            # both at once is the frame, and :P is unsatisfiable again
            sat_assume_var(sat, forbidden_lit(sat, :P, 1))
            sat_assume_var(sat, forbidden_lit(sat, :P, 2))
            @test !is_satisfiable(sat, [:P])
        end
        @test verdicts(sat) == before

        # the lift hands back what the body returns, and puts the frame back
        # even when the body throws
        @test with_classes_relaxed(() -> :done, sat) === :done
        @test_throws ErrorException with_classes_relaxed(sat) do
            error("body threw")
        end
        @test verdicts(sat) == before
        @test !is_satisfiable(sat, [:P])
    finally
        finalize(sat)
    end

    # an unconstrained query deactivates nothing, so there is no frame to lift
    # and the body runs against the instance as it stands
    plain = Problem([:P, :Q])
    sat = SAT(rank_pkg_info(info, plain))
    try
        @test isempty(sat.deact)
        before = verdicts(sat)
        @test with_classes_relaxed(() -> verdicts(sat), sat) == before
        @test verdicts(sat) == before
    finally
        finalize(sat)
    end
end

@testset "push frames have to balance" begin
    # A helper that opens a frame is relying on the frames below it being the
    # ones it thinks they are, since `sat_pop` retracts whatever was pushed
    # last. `sat.depth` is what lets it say so, and what these check is that it
    # does: a body that pushes without popping leaves the instance a frame
    # deeper than the helper found it, which is a wrong `sat_pop` somewhere
    # later and nothing at all at the call that caused it.
    data = Dict(
        :P => PkgData([:v2, :v1], NO_DEPS, Dict(:v2 => Dict(:Q => [:w1]))),
        :Q => PkgData([:w2, :w1], NO_DEPS, NO_COMP),
    )
    prob = Problem([:P, :Q]; compat = Dict(:P => [:v1]))
    info = pkg_info(data, prob; filter = false)

    for (name, prob) in ((:constrained, prob), (:plain, Problem([:P, :Q])))
        sat = SAT(rank_pkg_info(info, prob))
        try
            # the frame the query's deactivations sit in, and nothing else
            @test sat.depth == (isempty(sat.deact) ? 0 : 1)
            before = verdicts(sat)
            # a balanced body is what both helpers are for, and leaves no trace
            @test with_temp_clauses(() -> :ok, sat) === :ok
            @test with_classes_relaxed(() -> :ok, sat) === :ok
            @test sat.depth == (isempty(sat.deact) ? 0 : 1)
            @test verdicts(sat) == before
        finally
            finalize(sat)
        end
        # an unbalanced one is caught by each of them, at the call that did it
        sat = SAT(rank_pkg_info(info, prob))
        try
            @test_throws AssertionError with_temp_clauses(() -> sat_push(sat), sat)
        finally
            finalize(sat)
        end
        sat = SAT(rank_pkg_info(info, prob))
        try
            @test_throws AssertionError with_classes_relaxed(() -> sat_push(sat), sat)
        finally
            finalize(sat)
        end
    end

    # ... and lifting the deactivations from inside a temp-clause frame would
    # pop that frame instead of theirs, which is caught on the way in
    sat = SAT(rank_pkg_info(info, prob))
    try
        @test !isempty(sat.deact) # there is a frame to get wrong
        @test_throws AssertionError with_temp_clauses(sat) do
            with_classes_relaxed(() -> nothing, sat)
        end
    finally
        finalize(sat)
    end
end

@testset "the substrate stays internal" begin
    # substrate, reached by qualified name: `Resolver` exports none of it
    for name in (:installed_lit, :forbidden_lit, :with_classes_relaxed,
                 :push_deactivations!, :exclusion_sources, :class_exclusions)
        @test isdefined(Resolver, name)
        @test name ∉ names(Resolver)
    end
end
