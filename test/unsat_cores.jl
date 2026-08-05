# The `Resolver.UnsatCores` submodule: minimal unsatisfiable subsets and
# minimal correction sets over assumption literals (src/UnsatCores.jl).
#
# The oracle throughout is brute-force enumeration: on instances small enough,
# every potential solution is enumerated and checked for validity, and a set of
# assumption literals counts as satisfiable iff some valid solution satisfies
# all of them. That is an independent reference for "satisfiable assuming X" —
# it never asks the SAT solver anything — so it pins the contracts the
# primitives claim: minimality of cores, maximality of satisfiable subsets, the
# complement relationship, and completeness of the enumeration.

# the API under test is exactly what the submodule exports — three functions,
# taken unqualified the way the diagnostics machinery will take them
using Resolver.UnsatCores

# instance plumbing, from Resolver proper: none of this is part of the API above
using Resolver: PkgData, pkg_info, SAT, finalize, is_satisfiable,
    is_unsatisfiable, sat_new_variable, sat_add, sat_add_var, with_temp_clauses

# white-box, deliberately: the single-repair pair is unexported — the enumeration
# is built out of it — and reaching past the export list is how its contracts get
# pinned directly rather than through the enumeration that wraps it
using Resolver.UnsatCores: sat_mss, sat_mcs

# the module speaks only literals, so where these tests want to talk about
# packages they do the translation themselves: "package p is installed" is p's
# package variable in the instance's variable map
pkg_lits(sat::SAT{P}, pkgs::AbstractVector{P}) where {P} =
    Int[sat.vars[p] for p in pkgs]

function pkg_names(sat::SAT{P}, pkgs::AbstractVector{P}, sub) where {P}
    L = pkg_lits(sat, pkgs)
    return P[pkgs[findfirst(==(l), L)] for l in sub]
end

## brute-force oracle

# every valid solution of `data` over `pkgs`, as version-or-nothing vectors
function valid_solutions(
    data :: AbstractDict{P,<:PkgData{P,V}},
    pkgs :: AbstractVector{P},
) where {P,V}
    sols = Vector{Union{V,Nothing}}[]
    TestResolver.each_potential_solution(data, pkgs) do s
        TestResolver.is_valid_solution(data, pkgs, s) && push!(sols, copy(s))
    end
    return sols
end

# a candidate assumption literal, paired with what it means as a predicate on a
# brute-force solution vector
struct Cand
    lit  :: Int
    pred :: Function
end

# "package p is installed" for every package, plus "not p's first class" for
# every package with classes, so that negative literals — the shape the
# diagnostics rebuild will use — get exercised too. A class variable stands for
# every member of its class at once, which is what the predicate has to say
function candidates(
    sat  :: SAT{P},
    pkgs :: AbstractVector{P};
    versions :: Bool = true,
) where {P}
    cands = Cand[]
    for (i, p) in enumerate(pkgs)
        push!(cands, Cand(sat.vars[p], s -> s[i] !== nothing))
    end
    versions || return cands
    for (i, p) in enumerate(pkgs)
        info_p = sat.info[p]
        isempty(info_p.members) && continue
        mem = Set(info_p.versions[j] for j in info_p.members[1])
        push!(cands, Cand(-(sat.vars[p] + 1), s -> s[i] ∉ mem))
    end
    return cands
end

lits(cands::Vector{Cand}) = Int[c.lit for c in cands]

# is the instance satisfiable assuming all of `cands`?
oracle_sat(sols, cands::Vector{Cand}) =
    any(sol -> all(c -> c.pred(sol), cands), sols)

# every minimal correction set, by exhaustive enumeration of the subsets of the
# candidate set: the complements of the maximal satisfiable ones
function brute_force_mcses(sols, cands::Vector{Cand})
    n = length(cands)
    n ≤ 16 || error("$n candidates is too many to enumerate")
    bit(m, i) = m & (1 << (i-1)) ≠ 0
    ok = Bool[oracle_sat(sols, Cand[cands[i] for i = 1:n if bit(m, i)])
              for m = 0:(1<<n)-1]
    mcses = Set{Set{Int}}()
    for m = 0:(1<<n)-1
        ok[m+1] || continue
        # maximal: no candidate can be added without becoming unsatisfiable
        all(!ok[(m | 1 << (i-1)) + 1] for i = 1:n if !bit(m, i)) || continue
        push!(mcses, Set{Int}(cands[i].lit for i = 1:n if !bit(m, i)))
    end
    return mcses
end

# is `sub` a subsequence of `all`?
function is_subsequence(sub, all)
    i = 0
    for x in sub
        j = findnext(==(x), all, i + 1)
        j === nothing && return false
        i = j
    end
    return true
end

## the full property check, run on every instance below

# the primitives' answers on `L`, in a fixed call order. two fresh instances of
# the same problem, asked the same questions in the same order, must transcribe
# identically -- that is the determinism the primitives promise. what they do
# *not* promise is idempotence within one instance: the MUS seed is the solver's
# failed-assumption core, which depends on the solver's accumulated state, so
# asking twice can legitimately land on a different (equally minimal) core. see
# the checks below for which answers are state-independent and which aren't.
transcript(sat::SAT, L::Vector{Int}) = (
    mus   = sat_mus(sat, L),
    cover = sat_disjoint_muses(sat, L),
    mss   = sat_mss(sat, L),
    mcs   = sat_mcs(sat, L),
    mcses = sat_mcses(sat, L),
)

function check_mus_mcs(sat::SAT, info, sols, cands::Vector{Cand})
    L = lits(cands)
    by = Dict{Int,Cand}(c.lit => c for c in cands)
    pick(sub) = Cand[by[l] for l in sub]
    all_sat = oracle_sat(sols, cands)

    # a MUS is unsatisfiable, and satisfiable with any one element dropped
    function check_is_mus(mus, order = L)
        @test is_subsequence(mus, order)
        @test all_sat == isempty(mus)
        all_sat && return
        @test !oracle_sat(sols, pick(mus))
        for l in mus
            @test oracle_sat(sols, pick(filter(≠(l), mus)))
        end
    end

    t = transcript(sat, L)
    check_is_mus(t.mus)

    # the disjoint cover: each a MUS, pairwise disjoint, remainder satisfiable
    @test all_sat == isempty(t.cover)
    covered = Int[]
    for m in t.cover
        check_is_mus(m)
        append!(covered, m)
    end
    @test allunique(covered)
    @test oracle_sat(sols, pick(setdiff(L, covered)))

    # an MSS is satisfiable and admits no further candidate
    @test t.mss !== nothing
    @test is_subsequence(t.mss, L)
    @test oracle_sat(sols, pick(t.mss))
    for l in setdiff(L, t.mss)
        @test !oracle_sat(sols, pick([t.mss; l]))
    end
    @test all_sat == (t.mss == L)

    # the MCS is its complement, and no proper subset of it is a repair
    @test t.mcs == Int[l for l in L if l ∉ Set(t.mss)]
    @test oracle_sat(sols, pick(setdiff(L, t.mcs)))
    for l in t.mcs
        @test !oracle_sat(sols, pick(setdiff(L, setdiff(t.mcs, [l]))))
    end

    # enumeration: complete, undominated, each exactly once
    @test all(m -> is_subsequence(m, L), t.mcses)
    @test allunique(t.mcses)
    @test Set(Set.(t.mcses)) == brute_force_mcses(sols, cands)
    @test length(t.mcses) == length(Set(Set.(t.mcses)))
    @test Set(t.mcs) ∈ Set.(t.mcses) # the greedy one is one of them

    # a fresh instance of the same problem answers identically
    sat′ = SAT(info)
    try
        @test transcript(sat′, L) == t
    finally
        finalize(sat′)
    end

    # asked again on *this* instance, the greedy answers are unchanged -- they
    # are fixed by the candidate order alone, since a candidate is kept exactly
    # when keeping it is satisfiable -- and the enumeration is still complete,
    # while the MUS only has to be a MUS again
    @test sat_mss(sat, L) == t.mss
    @test sat_mcs(sat, L) == t.mcs
    @test Set(Set.(sat_mcses(sat, L))) == Set(Set.(t.mcses))
    check_is_mus(sat_mus(sat, L))

    # `limit` truncates the enumeration: which correction sets come back depends
    # on discovery order, but they are always that many distinct valid ones
    for k in unique([0, 1, length(t.mcses) ÷ 2, length(t.mcses)])
        some = sat_mcses(sat, L; limit = k)
        @test length(some) == min(k, length(t.mcses))
        @test allunique(Set.(some))
        @test Set.(some) ⊆ Set.(t.mcses)
    end

    # reversing the candidate order changes which witness is picked, but not
    # what the answers have to satisfy: MUS/MSS/MCS remain valid and the
    # enumeration remains complete
    R = reverse(L)
    check_is_mus(sat_mus(sat, R), R)
    rmss = sat_mss(sat, R)
    @test oracle_sat(sols, pick(rmss))
    for l in setdiff(R, rmss)
        @test !oracle_sat(sols, pick([rmss; l]))
    end
    @test Set(Set.(sat_mcses(sat, R))) == Set(Set.(t.mcses))
end

# build the SAT instance for `data` with every package present & unpruned, run
# every check against brute force, and clean up
function test_mus_mcs(
    data :: AbstractDict{P};
    versions :: Bool = true,
) where {P}
    pkgs = sort!(collect(keys(data)))
    info = pkg_info(data, pkgs; filter = false)
    @test Set(keys(info)) == Set(pkgs)
    sat = SAT(info)
    try
        sols = valid_solutions(data, pkgs)
        @test !isempty(sols) # the empty solution is always valid
        check_mus_mcs(sat, info, sols, candidates(sat, pkgs; versions))
        return sat
    finally
        finalize(sat)
    end
end

## hand-built instances

const nodeps = Dict{Symbol,Vector{Symbol}}()
const nocomp = Dict{Symbol,Dict{Symbol,Any}}()

# A needs C@v1, B needs C@v2: A & B conflict, D is unrelated
const conflict_data = Dict(
    :A => PkgData([:v1], Dict(:v1 => [:C]), Dict(:v1 => Dict(:C => [:v1]))),
    :B => PkgData([:v1], Dict(:v1 => [:C]), Dict(:v1 => Dict(:C => [:v2]))),
    :C => PkgData([:v2, :v1], nodeps, nocomp),
    :D => PkgData([:v1], nodeps, nocomp),
)

# two independent conflicts: A vs B (over C) and E vs F (over G)
const two_conflicts_data = Dict(
    :A => PkgData([:v1], Dict(:v1 => [:C]), Dict(:v1 => Dict(:C => [:v1]))),
    :B => PkgData([:v1], Dict(:v1 => [:C]), Dict(:v1 => Dict(:C => [:v2]))),
    :C => PkgData([:v2, :v1], nodeps, nocomp),
    :E => PkgData([:v1], Dict(:v1 => [:G]), Dict(:v1 => Dict(:G => [:v1]))),
    :F => PkgData([:v1], Dict(:v1 => [:G]), Dict(:v1 => Dict(:G => [:v2]))),
    :G => PkgData([:v2, :v1], nodeps, nocomp),
)

# nothing conflicts with anything
const no_conflict_data = Dict(
    :A => PkgData([:v2, :v1], Dict(:v1 => [:B]), nocomp),
    :B => PkgData([:v1], nodeps, nocomp),
)

@testset "unsat cores: the API is the submodule's export list" begin
    api = [:sat_mus, :sat_disjoint_muses, :sat_mcses]
    @test Set(names(Resolver.UnsatCores)) == Set([:UnsatCores; api])
    # the machinery underneath — including the single-repair pair the enumeration
    # is built out of — is defined but stays inside
    for internal in (:sat_solve_assuming, :failed_assumptions, :sat_grow!,
                     :sat_mss, :sat_mcs, :sat_mus_core, :mus_shrink_linear)
        @test isdefined(Resolver.UnsatCores, internal)
        @test internal ∉ names(Resolver.UnsatCores)
    end
    # and `Resolver` re-exports none of it
    @test isempty(intersect(names(Resolver), api))
    # nothing survives of what was dropped: the package-keyed conveniences and
    # their package/literal translation, the union of a cover, and the `_lits`
    # suffix that only existed to keep the conveniences from colliding
    for gone in (:sat_mus_union, :sat_pkg_lits, :sat_lit_pkgs, :sat_mus_lits,
                 :sat_disjoint_muses_lits, :sat_mus_union_lits, :sat_mss_lits,
                 :sat_mcs_lits, :sat_mcses_lits)
        @test !isdefined(Resolver.UnsatCores, gone)
        @test !isdefined(Resolver, gone)
    end
end

@testset "unsat cores: hand-built instances" begin
    for data in (conflict_data, two_conflicts_data, no_conflict_data)
        test_mus_mcs(data)
        test_mus_mcs(data; versions = false)
    end
end

@testset "unsat cores: candidate order is the preference knob" begin
    info = pkg_info(conflict_data; filter = false)
    sat = SAT(info)
    try
        # ask about packages by assuming each one's "installed" literal, and read
        # the answers back as package names
        ask(f, pkgs) = pkg_names(sat, pkgs, f(sat, pkg_lits(sat, pkgs)))
        # A & B are the conflict; whichever comes last is what gets dropped
        # (white-box: the single-repair pair)
        @test ask(sat_mcs, [:A, :B, :C, :D]) == [:B]
        @test ask(sat_mcs, [:B, :A, :C, :D]) == [:A]
        @test ask(sat_mss, [:A, :B, :C, :D]) == [:A, :C, :D]
        @test ask(sat_mss, [:B, :A, :C, :D]) == [:B, :C, :D]
        # the conflict is the same set either way, listed in candidate order
        @test ask(sat_mus, [:A, :B, :C, :D]) == [:A, :B]
        @test ask(sat_mus, [:D, :C, :B, :A]) == [:B, :A]
        # both repairs, and only those, come out of the enumeration
        pkgs = [:A, :B, :C, :D]
        @test Set(Set(pkg_names(sat, pkgs, m))
                  for m in sat_mcses(sat, pkg_lits(sat, pkgs))) ==
            Set([Set([:A]), Set([:B])])
        # requiring neither A nor B is fine, so there is nothing to report
        @test isempty(ask(sat_mus, [:C, :D]))
        @test isempty(sat_disjoint_muses(sat, pkg_lits(sat, [:C, :D])))
        @test ask(sat_mss, [:C, :D]) == [:C, :D]
        @test isempty(ask(sat_mcs, [:C, :D]))
        @test sat_mcses(sat, pkg_lits(sat, [:C, :D])) == [Int[]]
    finally
        finalize(sat)
    end
end

@testset "unsat cores: two independent conflicts" begin
    info = pkg_info(two_conflicts_data; filter = false)
    sat = SAT(info)
    try
        pkgs = [:A, :B, :C, :E, :F, :G]
        L = pkg_lits(sat, pkgs)
        names_of(sub) = pkg_names(sat, pkgs, sub)
        # a single conflict covers one of them -- which one is up to the solver's
        # failed-assumption core -- while the disjoint cover finds both
        @test Set(names_of(sat_mus(sat, L))) ∈
            (Set([:A, :B]), Set([:E, :F]))
        cover = sat_disjoint_muses(sat, L)
        @test Set(Set(names_of(m)) for m in cover) ==
            Set([Set([:A, :B]), Set([:E, :F])])
        # repairing needs one package from each conflict: 2 × 2 = 4 ways
        mcses = sat_mcses(sat, L)
        @test length(mcses) == 4
        @test Set(Set(names_of(m)) for m in mcses) == Set(
            Set([x, y]) for x in (:A, :B), y in (:E, :F))
        # ... and a single grown repair keeps the earlier of each conflicting
        # pair, dropping the later one: pure candidate order, no solver state
        # (white-box: the single-repair pair)
        @test names_of(sat_mcs(sat, L)) == [:B, :F]
        @test pkg_names(sat, reverse(pkgs), sat_mcs(sat, reverse(L))) ==
            [:E, :A]
    finally
        finalize(sat)
    end
end

@testset "unsat cores: unsatisfiable instance" begin
    # with no assumption to blame, the empty set is the MUS and there is no
    # correction set at all: the primitives have to say so rather than lie
    info = pkg_info(conflict_data; filter = false)
    sat = SAT(info)
    try
        L = pkg_lits(sat, [:A, :B, :C, :D])
        with_temp_clauses(sat) do
            x = sat_new_variable(sat)
            sat_add_var(sat, x); sat_add(sat)
            sat_add_var(sat, -x); sat_add(sat)
            @test is_unsatisfiable(sat)
            @test isempty(sat_mus(sat, L))
            @test isempty(sat_disjoint_muses(sat, L))
            @test isempty(sat_mcses(sat, L))
            # white-box: no repair exists, and the pair says so rather than
            # inventing an empty one
            @test sat_mss(sat, L) === nothing
            @test sat_mcs(sat, L) === nothing
        end
        # popping the contradiction restores the instance
        @test is_satisfiable(sat)
        @test sat_mus(sat, L) == first(sat_disjoint_muses(sat, L))
    finally
        finalize(sat)
    end
end

@testset "unsat cores: empty candidate set" begin
    info = pkg_info(conflict_data; filter = false)
    sat = SAT(info)
    try
        @test isempty(sat_mus(sat, Int[]))
        @test isempty(sat_disjoint_muses(sat, Int[]))
        @test sat_mcses(sat, Int[]) == [Int[]]
        @test sat_mss(sat, Int[]) == Int[] # white-box
        @test sat_mcs(sat, Int[]) == Int[] # white-box
    finally
        finalize(sat)
    end
end

## randomized sweep over tiny instances

@testset "unsat cores: brute-force reference on tiny data" begin
    Random.seed!(rand(RandomDevice(), UInt64))
    for m = 2:3, n = 2:3
        make_deps, make_comp, data, d, c = tiny_data_makers(m, n)
        for _ = 1:15
            deps = make_deps(randbits(d))
            comp = make_comp(randbits(c))
            fill_data!(m, n, deps, comp, data)
            # `data` is reused across iterations, so copy it
            test_mus_mcs(Dict(p => data[p] for p = 1:m))
        end
    end
    # a wider instance, package literals only (candidate sets have to stay
    # small enough to enumerate every subset)
    for _ = 1:5
        m, n = 5, 2
        make_deps, make_comp, data, d, c = tiny_data_makers(m, n)
        fill_data!(m, n, make_deps(randbits(d)), make_comp(randbits(c)), data)
        test_mus_mcs(Dict(p => data[p] for p = 1:m); versions = false)
    end
end
