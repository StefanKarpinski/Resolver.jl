# Is a conflict's proof a proof?
#
# Two things make it one, and each is checked here.
#
#   * every printed line is true of the universe the query left;
#   * the lines cannot all hold at once, so they entail the verdict.
#
# Together those say the report is a set of true statements that contradict,
# which is the whole of what proving unsatisfiability is.
#
# The lines are flat: none is derived from another on the page, so there is no
# step between them to check. What each one *is* derived from is the registry,
# and that is what the first question asks -- put to the instance the failed
# resolve ran on, so the query's own constraints are in force. A line that a
# withheld argument reaches through other packages is judged the same way as
# one read straight off the registry: is it true there.
#
# The only thing assumed beyond the lines themselves is what "installed at a
# version" means: a package takes one of its versions or none at all. A reader
# assumes that much; assuming more would be checking the universe rather than
# the page.

module ProofCheck

using Resolver: SAT, PicoSAT, nclasses, Problem, exclusion_kinds,
    installed_lit, forbidden_lit, sat_assume_var, sat_solve
using Resolver.Clauses: Clause, packages, isbottom, subsumes, clause_phrase,
    resolve_on, resolve_raw, resolve_all
using Resolver.Clauses: absent, present, Lit
using Resolver.Clauses: Clauses
using Resolver.Diagnostics: Diagnostics, Conflict, Line, clause_versions,
    clauses_satisfiable

export chain_is_a_proof, chain_hole, proof_problems, lines_are_true,
    printed_lines, names_what_it_uses, proofs_stand_alone

# The lines, taken as all the reader has, cannot hold together.
chain_is_a_proof(sat::SAT{P,V}, chain::Vector{Clause{P}}) where {P,V} =
    !clauses_satisfiable(sat, chain)

# ... and when they can, what an assignment satisfying all of them looks like
function chain_hole(sat::SAT{P,V}, chain::Vector{Clause{P}}) where {P,V}
    ps = P[]
    for c in chain, p in packages(c)
        p in ps || push!(ps, p)
    end
    pico = PicoSAT.init()
    try
        var = Dict{P,Vector{Int}}(); n = 0
        for p in ps
            var[p] = Int[(n += 1) for _ = 1:(length(clause_versions(sat, p)) + 1)]
        end
        PicoSAT.adjust(pico, max(n, 1))
        add(ls) = (for l in ls; PicoSAT.add(pico, l) end; PicoSAT.add(pico, 0))
        for p in ps
            add(var[p])
            for i in eachindex(var[p]), j = i+1:length(var[p])
                add([-var[p][i], -var[p][j]])
            end
        end
        for c in chain
            add(Int[var[p][i] for (p, m) in c.lits for i in eachindex(m) if m[i]])
        end
        PicoSAT.sat(pico) == PicoSAT.UNSATISFIABLE && return nothing
        out = Pair{P,Any}[]
        for p in ps
            vs = clause_versions(sat, p)
            k = findfirst(i -> PicoSAT.deref(pico, var[p][i]) > 0,
                          eachindex(var[p]))
            push!(out, p => (k === nothing || k > length(vs) ? nothing : vs[k]))
        end
        return out
    finally
        PicoSAT.reset(pico)
    end
end

# Can these clauses hold together, over these packages' versions? Asked of the
# clauses alone -- no registry, no query -- because entailment between
# statements is a question about the statements.
function clauses_hold(versions::Dict{P,Vector{V}}, cs::Vector{Clause{P}}) where {P,V}
    ps = P[]
    for c in cs, p in packages(c)
        p in ps || push!(ps, p)
    end
    isempty(cs) && return true
    pico = PicoSAT.init()
    try
        var = Dict{P,Vector{Int}}(); n = 0
        for p in ps
            var[p] = Int[(n += 1) for _ = 1:(length(versions[p]) + 1)]
        end
        PicoSAT.adjust(pico, max(n, 1))
        add(ls) = (for l in ls; PicoSAT.add(pico, l) end; PicoSAT.add(pico, 0))
        for p in ps                     # exactly one version, or none
            add(var[p])
            for i in eachindex(var[p]), j = i+1:length(var[p])
                add([-var[p][i], -var[p][j]])
            end
        end
        for c in cs
            add(Int[var[p][i] for (p, m) in c.lits for i in eachindex(m) if m[i]])
        end
        return PicoSAT.sat(pico) != PicoSAT.UNSATISFIABLE
    finally
        PicoSAT.reset(pico)
    end
end

# Is a line true of the universe this query left?
#
# Asked of the instance the failed resolve ran on, so the query's own
# constraints are in force: a line can be true because the compat removed the
# counterexample, and the registry on its own would call it false. What is
# asserted is the denial of the line -- every package outside what its literal
# allows -- and a model of that is a counterexample.
function untrue_witness(sat::SAT{P,V}, prob::Problem{P},
                        c::Clause{P}) where {P,V}
    isbottom(c) && return nothing     # ⊥ is never true; it is what is derived
    # A clause no literal of which admits ⊥ demands that something be
    # installed, and nothing in a registry demands that -- only the query
    # does. It is a premise, or derived from one, and either way it is not a
    # claim about the universe for this to judge.
    any(absent(m) for (_, m) in c.lits) || return nothing
    assume = Int[]
    # The universe as *this query* left it. A line can be true because the
    # compat removed the counterexample, so asking without the query's own
    # constraints in force calls it false and is asking about a different
    # registry than the one the report is about.
    for (q, _) in c.lits
        haskey(sat.reps, q) || continue
        for (cl, r) in enumerate(sat.reps[q])
            iszero(r) && push!(assume, forbidden_lit(sat, q, cl))
        end
    end
    for (p, m) in c.lits
        haskey(sat.info, p) && haskey(sat.vars, p) || return nothing
        info = sat.info[p]
        mem = clause_versions(sat, p)
        inm = Set{V}(mem[i] for i in eachindex(m) if i ≤ length(mem) && m[i])
        # a class is ruled out only if every version of it is
        for cl = 1:nclasses(info)
            vs = V[info.versions[j] for j in eachindex(info.versions)
                   if info.classes[j] == cl]
            isempty(vs) && continue
            all(v -> v in inm, vs) || continue
            push!(assume, forbidden_lit(sat, p, cl))
        end
        # the literal admits absence, so a counterexample has it installed
        absent(m) && push!(assume, installed_lit(sat, p))
    end
    # A line with nothing assumable about it -- every class of every package
    # it names left standing -- is one this cannot refute either way
    isempty(assume) && return nothing
    for l in assume
        sat_assume_var(sat, l)
    end
    sat_solve(sat) || return nothing
    out = Pair{P,Any}[]
    for (q, _) in c.lits
        info = sat.info[q]
        k = findfirst(cl -> PicoSAT.deref(sat.pico, sat.vars[q] + cl) > 0,
                      1:nclasses(info))
        push!(out, q => (k === nothing ? nothing :
            info.versions[findfirst(==(k), info.classes)]))
    end
    # The instance chooses a *class*, and a query's exclusions are finer than
    # one -- a class can hold an excluded version beside an admitted one. So a
    # model is not automatically a counterexample: read it back and keep it
    # only if it really does falsify every literal. Reporting one that does
    # not is how an oracle cries wolf.
    for (q, m) in c.lits
        v = last(out[findfirst(x -> first(x) === q, out)])
        mem = clause_versions(sat, q)
        if v === nothing
            absent(m) && return nothing        # ⊥ satisfies this literal
        else
            # ... and a version this query rules out is no counterexample: the
            # instance picks a class, and a class can hold an excluded version
            # beside an admitted one, so the exclusion has to be re-checked
            # here rather than assumed away
            isempty(exclusion_kinds(prob, q, v)) || return nothing
            i = findfirst(==(v), mem)
            i === nothing && return nothing
            m[i] && return nothing             # the chosen version satisfies it
        end
    end
    return out
end

function lines_are_true(sat::SAT{P,V}, prob::Problem{P}, clauses) where {P,V}
    bad = String[]
    vers(p) = clause_versions(sat, p)
    for c in clauses
        w = untrue_witness(sat, prob, c)
        w === nothing && continue
        push!(bad, "not true of the universe: " * clause_phrase(c, vers) *
              " -- " * join(["$p=$(v === nothing ? "none" : v)" for (p, v) in w], ", "))
    end
    return bad
end

# every line a report prints
printed_lines(c::Conflict{P,V}) where {P,V} = Clause{P}[l.clause for l in c.lines]

# Does the report use what it names?
#
# A conflict names the requirements it answers for, and it answers for them by
# way of the reasons they belong to -- so a requirement it names and no line
# mentions is one whose part in the story was never shown. The reader is told
# `A and B cannot both be satisfied` and then given an argument in which B does
# not appear, which reads as a mistake and is one: it means a reason went
# unproved, or was proved by leaning on a premise some other reason brought.
#
# Cheap to check and worth checking every time. It has caught the same fault
# twice, in different guises, where reading a sample of reports did not.
function names_what_it_uses(c::Conflict{P,V}) where {P,V}
    isempty(c.lines) && return String[]   # nothing to prove, nothing to name
    used = Set{P}(p for l in c.lines for p in packages(l.clause))
    return String["names $r and no line mentions it" for r in c.reqs
                  if !(r in used)]
end

# Does each proof stand on its own?
#
# A conflict answers for several reasons and sets each out as a proof of its
# own. The union of them contradicting is not enough: where one reason is
# already impossible the union stays impossible whatever the others say, so a
# check over the union would pass a proof that borrows its conclusion from the
# proof beside it. Each is asked separately -- its own lines, and the query's
# facts, which the reader has on the page either way.
function proofs_stand_alone(sat::SAT{P,V}, c::Conflict{P,V}) where {P,V}
    bad = String[]
    given = Clause{P}[l.clause for l in c.lines if l.given]
    for n in unique!(Int[l.proof for l in c.lines if !l.given])
        mine = Clause{P}[l.clause for l in c.lines if !l.given && l.proof == n]
        isempty(mine) && continue
        clauses_satisfiable(sat, Clause{P}[given; mine]) || continue
        push!(bad, "proof $n does not contradict on its own: " *
              join([clause_phrase(x, q -> clause_versions(sat, q)) for x in mine], "; "))
    end
    return bad
end

function proof_problems(sat::SAT{P,V}, prob::Problem{P},
                        c::Conflict{P,V}) where {P,V}
    cs = printed_lines(c)
    bad = names_what_it_uses(c)
    append!(bad, proofs_stand_alone(sat, c))
    append!(bad, lines_are_true(sat, prob, cs))
    isempty(cs) || chain_is_a_proof(sat, cs) ||
        push!(bad, "the printed lines can all hold at once: " *
                   string(chain_hole(sat, cs)))
    return bad
end

end # module
