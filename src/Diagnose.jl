# UNSAT diagnostics: explain why a resolution problem is unsatisfiable and
# list verified fixes. Built on a separate, diagnosis-only SAT instance in
# which every blameable constraint group is guarded by a selector variable,
# toggled via picosat assumptions.

# facts — the vocabulary of explanations and fixes

abstract type Fact end

struct Requirement{P} <: Fact
    pkg :: P
end

struct UserCompat{P,V} <: Fact
    pkg     :: P
    allowed :: Vector{V}   # versions of pkg (from data) the user's compat admits
end

struct Hold{P,V} <: Fact
    pkg     :: P
    version :: V
end

struct Bound{P,V} <: Fact
    pkg      :: P          # whose compat declaration this is
    versions :: Vector{V}  # versions of pkg sharing this declaration
    dep      :: P          # the target package
    allowed  :: Vector{V}  # versions of dep NOT excluded by the declaration
end

fact_pkg(f::Fact) = f.pkg

# value equality for facts (immutable, compared field-wise)
Base.:(==)(a::F, b::F) where {F<:Fact} =
    all(getfield(a, i) == getfield(b, i) for i = 1:nfields(a))
function Base.hash(a::Fact, h::UInt)
    for i = 1:nfields(a)
        h = hash(getfield(a, i), h)
    end
    hash(typeof(a), h)
end

struct Fix{P,V}
    actions  :: Vector{Fact}  # Requirement=drop, UserCompat=relax, Hold=unhold
    solution :: Dict{P,V}     # verified resulting versions
end

struct UpstreamFix{P,V}
    bound    :: Bound{P,V}
    solution :: Dict{P,V}
end

struct Conflict{P,V}
    reqs     :: Vector{P}               # requirement-level MUS (cluster)
    chain    :: Vector{Fact}            # group-MUS, rendered in story order
    upstream :: Vector{UpstreamFix{P,V}}
end

struct Diagnosis{P,V}
    reqs      :: Vector{P}
    conflicts :: Vector{Conflict{P,V}}
    fixes     :: Vector{Fix{P,V}}       # global: each verified to fix everything
    versions  :: Dict{P,Vector{V}}      # package => full version list, for rendering
end

# diagnostic SAT instance: one selector variable per blameable group

mutable struct DiagSAT{P,V,D<:AbstractDict{P,<:PkgData{P,V}}}
    data  :: D
    pico  :: Ptr{Cvoid}
    vars  :: Dict{P,Int}                 # p => base var (p@i is vars[p]+i)
    reqs  :: Vector{P}                   # requirements, sorted (validated ⊆ data)
    Rs    :: Vector{Tuple{Int,Requirement{P}}}
    Cs    :: Vector{Tuple{Int,UserCompat{P,V}}}
    Hs    :: Vector{Tuple{Int,Hold{P,V}}}
    Bs    :: Vector{Tuple{Int,Bound{P,V}}}
    facts :: Dict{Int,Fact}              # selector var => fact
end

function DiagSAT(
    data   :: AbstractDict{P,<:PkgData{P,V}},
    reqs   :: SetOrVec{P};
    compat :: AbstractDict{P} = Dict{P,Vector{V}}(),
    holds  :: AbstractDict{P,V} = Dict{P,V}(),
) where {P,V}
    # sort names for predictability
    names = sort!(collect(keys(data)))

    # variable indices:
    #   p    vars[p]    package p chosen
    #   p@i  vars[p]+i  version i of p chosen
    N = 1
    vars = Dict{P,Int}()
    for p in names
        vars[p] = N
        N += 1 + length(data[p].versions)
    end
    N -= 1 # last used variable

    pico = PicoSAT.init()
    try # free memory on error
        PicoSAT.adjust(pico, N)

        # hard (unguarded) clauses — mirror SAT() minus symmetric compat conflicts
        for p in names
            data_p = data[p]
            v_p = vars[p]
            n_p = length(data_p.versions)

            # package implies some version:  p => OR_i p@i
            PicoSAT.add(pico, -v_p)
            for i = 1:n_p
                PicoSAT.add(pico, v_p + i)
            end
            PicoSAT.add(pico, 0)

            # version implies its package:  p@i => p
            for i = 1:n_p
                PicoSAT.add(pico, -(v_p + i))
                PicoSAT.add(pico, v_p)
                PicoSAT.add(pico, 0)
            end

            # versions are mutually exclusive:  !p@i OR !p@j
            for i = 1:n_p-1, j = i+1:n_p
                PicoSAT.add(pico, -(v_p + i))
                PicoSAT.add(pico, -(v_p + j))
                PicoSAT.add(pico, 0)
            end

            # dependency edges (existence, not blameable):  p@i => q
            vidx = Dict{V,Int}(v => i for (i, v) in enumerate(data_p.versions))
            for (v, deps_pv) in data_p.depends
                haskey(vidx, v) || continue
                i = vidx[v]
                for q in deps_pv
                    (q == p || q ∉ keys(data)) && continue
                    PicoSAT.add(pico, -(v_p + i))
                    PicoSAT.add(pico, vars[q])
                    PicoSAT.add(pico, 0)
                end
            end
        end

        # selector-guarded group clauses: each group's clauses get ¬s prepended
        # so the group is active iff its selector s is assumed true.

        # R(p) — requirement:  (¬s, p)
        for p in reqs
            p in keys(data) || throw(ArgumentError(
                "Required package $p is not available in the package data"))
        end
        reqs_sorted = sort!(P[p for p in reqs])
        Rs = Tuple{Int,Requirement{P}}[]
        for p in reqs_sorted
            s = PicoSAT.inc_max_var(pico)
            PicoSAT.add(pico, -s)
            PicoSAT.add(pico, vars[p])
            PicoSAT.add(pico, 0)
            push!(Rs, (s, Requirement{P}(p)))
        end

        # C(p) — user compat/pin:  (¬s, ¬p@i) for each excluded version i
        Cs = Tuple{Int,UserCompat{P,V}}[]
        for p in sort!(collect(keys(compat)))
            p in keys(data) || continue
            set = compat[p]
            versions_p = data[p].versions
            excluded = Int[i for (i, v) in enumerate(versions_p) if v ∉ set]
            isempty(excluded) && continue   # no version excluded → no group
            s = PicoSAT.inc_max_var(pico)
            v_p = vars[p]
            for i in excluded
                PicoSAT.add(pico, -s)
                PicoSAT.add(pico, -(v_p + i))
                PicoSAT.add(pico, 0)
            end
            allowed = V[v for v in versions_p if v ∈ set]
            push!(Cs, (s, UserCompat{P,V}(p, allowed)))
        end

        # H(p) — manifest hold:  (¬s, ¬p@i) for each version i ≠ held version
        Hs = Tuple{Int,Hold{P,V}}[]
        for p in sort!(collect(keys(holds)))
            p in keys(data) || continue
            ver = holds[p]
            versions_p = data[p].versions
            excluded = Int[i for (i, v) in enumerate(versions_p) if v != ver]
            isempty(excluded) && continue   # no version excluded → no group
            s = PicoSAT.inc_max_var(pico)
            v_p = vars[p]
            for i in excluded
                PicoSAT.add(pico, -s)
                PicoSAT.add(pico, -(v_p + i))
                PicoSAT.add(pico, 0)
            end
            push!(Hs, (s, Hold{P,V}(p, ver)))
        end

        # B(p, q, excluded) — third-party compat bound, coalesced by identical
        # excluded-index-set of q:  (¬s, ¬p@i, ¬q@j) for i ∈ p_idxs, j ∈ excluded
        bspecs = Tuple{P,P,Vector{Int},Vector{Int}}[]
        for p in names
            data_p = data[p]
            byq = Dict{P,Dict{Vector{Int},Vector{Int}}}()  # q => (excluded => p_idxs)
            for (i, v) in enumerate(data_p.versions)
                haskey(data_p.compat, v) || continue
                comp_pv = data_p.compat[v]
                for q in sort!(collect(keys(comp_pv)))
                    (q == p || q ∉ keys(data)) && continue
                    set = comp_pv[q]
                    excluded = Int[j for (j, w) in enumerate(data[q].versions) if w ∉ set]
                    isempty(excluded) && continue   # excludes nothing → no group
                    d = get!(() -> Dict{Vector{Int},Vector{Int}}(), byq, q)
                    push!(get!(() -> Int[], d, excluded), i)
                end
            end
            for q in sort!(collect(keys(byq)))
                d = byq[q]
                for excluded in sort!(collect(keys(d)))
                    push!(bspecs, (p, q, sort!(d[excluded]), excluded))
                end
            end
        end
        Bs = Tuple{Int,Bound{P,V}}[]
        for (p, q, p_idxs, excluded) in bspecs
            s = PicoSAT.inc_max_var(pico)
            v_p = vars[p]
            v_q = vars[q]
            for i in p_idxs, j in excluded
                PicoSAT.add(pico, -s)
                PicoSAT.add(pico, -(v_p + i))
                PicoSAT.add(pico, -(v_q + j))
                PicoSAT.add(pico, 0)
            end
            pvers = V[data[p].versions[i] for i in p_idxs]
            exset = Set(excluded)
            qallowed = V[w for (j, w) in enumerate(data[q].versions) if j ∉ exset]
            push!(Bs, (s, Bound{P,V}(p, pvers, q, qallowed)))
        end

        # selector => fact lookup
        facts = Dict{Int,Fact}()
        for groups in (Rs, Cs, Hs, Bs), (s, f) in groups
            facts[s] = f
        end

        D = typeof(data)
        return finalizer(finalize, DiagSAT{P,V,D}(
            data, pico, vars, reqs_sorted, Rs, Cs, Hs, Bs, facts))
    catch # on error free picosat solver
        PicoSAT.reset(pico)
        rethrow()
    end
end

function finalize(diag::DiagSAT)
    pico = diag.pico
    pico == C_NULL && return
    diag.pico = C_NULL
    PicoSAT.reset(pico)
end

# selector accessors

req_sels(diag::DiagSAT)    = Int[s for (s, _) in diag.Rs]
compat_sels(diag::DiagSAT) = Int[s for (s, _) in diag.Cs]
hold_sels(diag::DiagSAT)   = Int[s for (s, _) in diag.Hs]
bound_sels(diag::DiagSAT)  = Int[s for (s, _) in diag.Bs]

all_selectors(diag::DiagSAT) =
    Int[s for groups in (diag.Rs, diag.Cs, diag.Hs, diag.Bs) for (s, _) in groups]

# assume the given selectors on and test satisfiability
function sat_diag(diag::DiagSAT, assume)
    for s in assume
        PicoSAT.assume(diag.pico, s)
    end
    PicoSAT.sat(diag.pico) == PicoSAT.SATISFIABLE
end

# extract the current model's installed versions (call right after a SAT result)
function diag_solution(diag::DiagSAT{P,V}) where {P,V}
    sol = Dict{P,V}()
    for (p, v_p) in diag.vars
        PicoSAT.deref(diag.pico, v_p) < 0 && continue
        versions_p = diag.data[p].versions
        for i = 1:length(versions_p)
            if PicoSAT.deref(diag.pico, v_p + i) > 0
                sol[p] = versions_p[i]
                break
            end
        end
    end
    return sol
end

# SAT models may set package variables that nothing requires; prune a solution
# to the dependency-reachable closure of its root packages so reported
# solutions contain only packages the satisfied requirements actually need.
# Pruned solutions stay valid: every remaining version's dependencies are in
# the closure by construction, and compat constraints on absent packages are
# vacuous.
function prune_solution!(
    diag  :: DiagSAT{P,V},
    sol   :: Dict{P,V},
    roots,
) where {P,V}
    keep = Set{P}()
    queue = P[p for p in roots if haskey(sol, p)]
    while !isempty(queue)
        p = pop!(queue)
        p in keep && continue
        push!(keep, p)
        depends_p = diag.data[p].depends
        haskey(depends_p, sol[p]) || continue
        for q in depends_p[sol[p]]
            q in keys(diag.data) && haskey(sol, q) && push!(queue, q)
        end
    end
    filter!(kv -> kv.first in keep, sol)
    return sol
end

facts_of(diag::DiagSAT, sels) = Fact[diag.facts[s] for s in sels]

# shrinkable selectors in actionability-bias kind order: B → H → C → R, with the
# given requirement selectors `rsels` last (a cluster's, or all of them)
full_shrink(diag::DiagSAT, rsels = req_sels(diag)) =
    Int[bound_sels(diag); hold_sels(diag); compat_sels(diag); collect(rsels)]

# group MUS with actionability-biased deletion order.
#
# `sels_on` are all the selectors assumed on; `shrink` is the (ordered) subset we
# may delete, given in kind order B → H → C → R. Selectors in `sels_on ∖ shrink`
# are context: always kept on. Returns a minimal (w.r.t. `shrink`) subset of
# `sels_on` that is still unsatisfiable.
#
# Deleting B before C/H means a bound is dropped whenever the same conflict can
# be told through a user-owned constraint, so the surviving MUS prefers the
# actionable explanation. We start from the full assumption set rather than the
# failed-assumption core: that core is nondeterministic and can exclude the
# actionable explanation before the biased order gets to prefer it, whereas
# deleting from the full set is deterministic and always honors the bias.
#
# A single ordered pass suffices — no restart after deletions. Satisfiability
# is monotone under removing assumptions: if deleting `s` once made the
# instance SAT (so `s` was kept), any later deletion of `s` would test a
# subset of that SAT assumption set, which is necessarily still SAT. Kept
# elements provably stay kept, and the final set is minimal for the same
# reason: the final set minus any member is a subset of a tested SAT set.
function group_mus(
    diag    :: DiagSAT,
    sels_on :: AbstractVector{Int},
    shrink  :: AbstractVector{Int},
)
    sat_diag(diag, sels_on) && return Int[]
    mus = Set{Int}(sels_on)
    for s in shrink
        delete!(mus, s)
        sat_diag(diag, collect(mus)) && push!(mus, s) # needed: keep it
    end
    return sort!(collect(mus)) # sorted so downstream iteration is stable
end

# partition the unsatisfiable requirements into independent conflict clusters
# (requirement-level MICE): repeatedly extract a requirement-level MUS with all
# C/H/B held on as context, remove its members, repeat until satisfiable. Each
# disjoint MUS is one cluster = one story. Returns a vector of R-selector groups.
function cluster_reqs(diag::DiagSAT)
    rset = Set(req_sels(diag))
    context = Int[compat_sels(diag); hold_sels(diag); bound_sels(diag)]
    remaining = req_sels(diag)
    clusters = Vector{Int}[]
    while true
        sels_on = Int[remaining; context]
        # shrink only requirement selectors; C/H/B stay on as fixed context
        mus = group_mus(diag, sels_on, remaining)
        rmus = filter(in(rset), mus)
        isempty(rmus) && break
        push!(clusters, rmus)
        done = Set(rmus)
        remaining = filter(!in(done), remaining)
    end
    return clusters
end

# deterministic action ordering within a fix
kind_rank(::Requirement) = 1
kind_rank(::Hold)        = 2
kind_rank(::UserCompat)  = 3
kind_rank(::Bound)       = 4
order_actions(fs::Vector{Fact}) =
    sort(fs; by = f -> (kind_rank(f), fact_pkg(f)))

# enumerate globally-verified fixes as minimal correction sets. Over user-owned
# selectors U = R ∪ H ∪ C (all bounds B held on), a greedy keep-order pass
# (R, then H, then C) yields a maximal satisfiable subset; its complement is a
# minimal correction set. A permanent blocking clause ⋁_{s∈MCS} s forces each
# later fix to keep something a previous one dropped (undominated enumeration).
function enumerate_fixes(diag::DiagSAT{P,V}; max_fixes::Integer = 8) where {P,V}
    bsel = bound_sels(diag)
    rset = Set(req_sels(diag))
    U = Int[req_sels(diag); hold_sels(diag); compat_sels(diag)]  # keep order
    fixes = Fix{P,V}[]
    while length(fixes) < max_fixes
        # any model with all bounds on and the blocking clauses satisfied?
        sat_diag(diag, bsel) || break
        # greedy maximal satisfiable subset in keep order
        kept = Int[]
        for s in U
            push!(kept, s)
            sat_diag(diag, Int[bsel; kept]) || pop!(kept)
        end
        mcs = setdiff(U, kept)
        isempty(mcs) && break   # unreachable when the whole problem is UNSAT
        # the fix's concrete resulting versions, pruned to the closure of
        # the kept requirements
        sat_diag(diag, Int[bsel; kept])
        sol = diag_solution(diag)
        roots = P[diag.facts[s].pkg for s in kept if s in rset]
        prune_solution!(diag, sol, roots)
        push!(fixes, Fix{P,V}(order_actions(facts_of(diag, mcs)), sol))
        # block: the next fix must keep something this one dropped
        for s in mcs
            PicoSAT.add(diag.pico, s)
        end
        PicoSAT.add(diag.pico, 0)
    end
    return filter_dominated(fixes)
end

# a fix dominates another if it removes a proper subset of the restrictions the
# other removes; dominated fixes carry no information and are suppressed. The
# blocking clauses prevent duplicate fixes but can, in corner cases, force a
# later fix to drop a strict superset of an earlier fix's restrictions.
function filter_dominated(fixes::Vector{Fix{P,V}}) where {P,V}
    sets = [Set{Fact}(fix.actions) for fix in fixes]
    keep = [!any(sets[j] ⊊ sets[i] for j in eachindex(sets))
            for i in eachindex(sets)]
    return fixes[keep]
end

# probe each bound in a cluster's MUS: with the cluster's requirements on, all
# user compat/holds on, and every bound on except the probed one, is the
# cluster satisfiable? If so, relaxing that upstream compat declaration would
# resolve this conflict. Probing per-cluster (not globally) lets a suggestion
# fire even while other, independent conflicts remain unfixed.
function upstream_probes(
    diag    :: DiagSAT{P,V},
    cluster :: Vector{Int},
    mus     :: Vector{Int},
) where {P,V}
    roots = P[diag.facts[s].pkg for s in cluster]
    ups = UpstreamFix{P,V}[]
    for s in mus
        f = diag.facts[s]
        f isa Bound || continue
        sels = Int[cluster; compat_sels(diag); hold_sels(diag);
                   filter(!=(s), bound_sels(diag))]
        if sat_diag(diag, sels)
            sol = prune_solution!(diag, diag_solution(diag), roots)
            push!(ups, UpstreamFix{P,V}(f, sol))
        end
    end
    return sort!(ups; by = u -> (u.bound.pkg, u.bound.dep))
end

# order a cluster's MUS facts into story order: requirements first (sorted),
# then a BFS from the required packages along package links — introducing a
# package emits its C/H facts, then its B facts (each B introduces its dep
# target). Deterministic: the queue is processed in sorted order.
function order_chain(diag::DiagSAT{P,V}, facts::Vector{Fact}) where {P,V}
    reqfacts = sort!(Fact[f for f in facts if f isa Requirement]; by = fact_pkg)
    Cby = Dict{P,UserCompat{P,V}}()
    Hby = Dict{P,Hold{P,V}}()
    Bby = Dict{P,Vector{Bound{P,V}}}()
    for f in facts
        f isa UserCompat && (Cby[f.pkg] = f)
        f isa Hold       && (Hby[f.pkg] = f)
        f isa Bound      && push!(get!(() -> Bound{P,V}[], Bby, f.pkg), f)
    end
    foreach(v -> sort!(v; by = b -> b.dep), values(Bby))
    chain = copy(reqfacts)
    visited = Set{P}()
    queue = P[fact_pkg(f) for f in reqfacts]
    while !isempty(queue)
        sort!(queue)
        p = popfirst!(queue)
        p in visited && continue
        push!(visited, p)
        haskey(Cby, p) && push!(chain, Cby[p])
        haskey(Hby, p) && push!(chain, Hby[p])
        for b in get(Bby, p, Bound{P,V}[])
            push!(chain, b)
            push!(queue, b.dep)
        end
    end
    # emit any facts not reached by the BFS (deterministically)
    for f in sort!(Fact[f for f in facts if f ∉ chain];
                   by = f -> (kind_rank(f), fact_pkg(f)))
        push!(chain, f)
    end
    return chain
end

# rendering

# compress a subset of a package's versions into runs against the full list,
# e.g. [1,2,3,5] against 1:5 → "1–3, 5"; the whole list → "all versions"
function format_versions(whole::AbstractVector, subset::AbstractVector)
    isempty(subset) && return "no versions"
    idx = sort!(Int[findfirst(==(v), whole) for v in subset])
    length(idx) == length(whole) && return "all versions"
    # compress consecutive-in-`whole` versions into runs, then order the runs
    # and each run's endpoints by version value -- independent of whether
    # `whole` is sorted ascending or descending (the registry provider sorts
    # versions descending, which otherwise prints ranges backwards, e.g.
    # "1.11.0–1.0.0")
    runs = Tuple{eltype(whole),eltype(whole)}[]
    i = 1
    while i ≤ length(idx)
        j = i
        while j < length(idx) && idx[j+1] == idx[j] + 1
            j += 1
        end
        a, b = whole[idx[i]], whole[idx[j]]
        push!(runs, (min(a, b), max(a, b)))
        i = j + 1
    end
    sort!(runs; by = first)
    return join((lo == hi ? string(lo) : string(lo, "–", hi)
                 for (lo, hi) in runs), ", ")
end

render_fact(io::IO, f::Requirement, d::Diagnosis) =
    print(io, "you require ", f.pkg)
render_fact(io::IO, f::UserCompat, d::Diagnosis) =
    print(io, "your compat restricts ", f.pkg, " to ",
          format_versions(d.versions[f.pkg], f.allowed))
render_fact(io::IO, f::Hold, d::Diagnosis) =
    print(io, f.pkg, " is pinned at ", f.version)
function render_fact(io::IO, f::Bound, d::Diagnosis)
    all_p = d.versions[f.pkg]
    src = length(f.versions) == length(all_p) ? "(all versions)" :
          format_versions(all_p, f.versions)
    print(io, f.pkg, " ", src, " requires ", f.dep, " at ",
          format_versions(d.versions[f.dep], f.allowed))
end

render_action(f::Requirement) = string("drop requirement ", f.pkg)
render_action(f::UserCompat)  = string("relax your compat on ", f.pkg)
render_action(f::Hold)        = string("unpin ", f.pkg)

# packages worth naming in a fix's resulting solution: the requirements plus
# the dependencies actually implicated in some conflict. A verified solution
# lists the entire transitive closure at whatever versions it happened to pick,
# which is mostly noise -- the contested packages are what a fix meaningfully
# changes. The complete assignment is still available on `Fix.solution` /
# `UpstreamFix.solution` for programmatic use (e.g. emitting a manifest).
function relevant_pkgs(d::Diagnosis{P,V}) where {P,V}
    rel = Set{P}(d.reqs)
    for c in d.conflicts, f in c.chain
        push!(rel, fact_pkg(f))
        f isa Bound && push!(rel, f.dep)
    end
    return rel
end

function render_solution(d::Diagnosis{P,V}, sol::AbstractDict{P,V}) where {P,V}
    rel = relevant_pkgs(d)
    pkgs = collect(p for p in keys(sol) if p in rel)
    isempty(pkgs) && return "nothing"
    reqset = Set(d.reqs)
    sort!(pkgs)
    sort!(pkgs; by = p -> p ∉ reqset)   # requirements first (stable)
    join((string(p, " ", sol[p]) for p in pkgs), ", ")
end

# compact summary
function Base.show(io::IO, d::Diagnosis)
    nc = length(d.conflicts)
    nf = length(d.fixes)
    print(io, "Diagnosis(", nc, nc == 1 ? " conflict, " : " conflicts, ",
          nf, nf == 1 ? " fix)" : " fixes)")
end

# full report
function Base.show(io::IO, ::MIME"text/plain", d::Diagnosis)
    for (n, c) in enumerate(d.conflicts)
        n > 1 && println(io)
        rs = join(c.reqs, ", ")
        println(io, "Conflict ", n, ": ", rs,
                length(c.reqs) == 1 ? " cannot be satisfied." :
                " cannot be satisfied together.")
        for f in c.chain
            print(io, "  • ")
            render_fact(io, f, d)
            println(io)
        end
    end
    if !isempty(d.fixes)
        println(io)
        println(io, "Verified fixes:")
        for (n, fix) in enumerate(d.fixes)
            println(io, "  ", n, ". ",
                    join(map(render_action, fix.actions), " and "))
            println(io, "     → allows: ", render_solution(d, fix.solution))
        end
    end
    ups = [u for c in d.conflicts for u in c.upstream]
    if !isempty(ups)
        println(io)
        println(io, "Upstream fixes:")
        for u in ups
            println(io, "  • a release of ", u.bound.pkg,
                    " relaxing its compat on ", u.bound.dep)
            println(io, "    → allows: ", render_solution(d, u.solution))
        end
    end
end

# entry point

function diagnose(
    data   :: AbstractDict{P,<:PkgData{P,V}},
    reqs   :: SetOrVec{P};
    compat :: AbstractDict{P} = Dict{P,Vector{V}}(),
    holds  :: AbstractDict{P,V} = Dict{P,V}(),
    max_fixes :: Integer = 8,
) where {P,V}
    diag = DiagSAT(data, reqs; compat, holds)
    try
        if sat_diag(diag, all_selectors(diag))
            throw(ArgumentError(
                "resolvable with the given constraints: nothing to diagnose"))
        end
        context = Int[compat_sels(diag); hold_sels(diag); bound_sels(diag)]
        conflicts = Conflict{P,V}[]
        for cluster in cluster_reqs(diag)
            sels_on = Int[cluster; context]
            mus = group_mus(diag, sels_on, full_shrink(diag, cluster))
            chain = order_chain(diag, facts_of(diag, mus))
            upstream = upstream_probes(diag, cluster, mus)
            creqs = sort!(P[diag.facts[s].pkg for s in cluster])
            push!(conflicts, Conflict{P,V}(creqs, chain, upstream))
        end
        sort!(conflicts; by = c -> c.reqs)
        # fixes last: enumerate_fixes adds permanent blocking clauses
        fixes = enumerate_fixes(diag; max_fixes)
        versions = Dict{P,Vector{V}}(p => collect(data[p].versions) for p in keys(data))
        return Diagnosis{P,V}(copy(diag.reqs), conflicts, fixes, versions)
    finally
        finalize(diag)
    end
end
