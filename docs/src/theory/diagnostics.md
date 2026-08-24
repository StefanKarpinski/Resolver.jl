# Diagnosing an unsatisfiable resolve

When the requirements cannot be satisfied there is more to say than that
they cannot: which of them cannot hold together, why, and what a user
could change so that they can. Every one of those answers comes from
putting further questions to the universe the resolve was run against and
the instance it failed on, and this page states what licenses them.

Everything here is short, and two earlier results are why. Class space (see
[the layered page](layered.md#Preprocessing-for-caching)) makes a user
constraint's whole effect on the universe the set of classes it empties, so
there is exactly one kind of thing to relax and it is already a unit fact of
the instance. Relaxation-stable filtering (see the
[previous page](relaxation.md)) makes the filtered universe answer for every
relaxation, so there is exactly one universe to ask.

## Relaxations

Write ``Q`` for the failed query, ``\mathrm{reqs}_Q`` for its requirements,
and ``E_Q`` for the set of classes its constraints empty — a class ``Q``
admits no member of cannot be selected, and that is the whole of what a
constraint does to the universe.

A **withdrawal** from ``Q`` is a pair ``(A, B)`` with
``A \subseteq \mathrm{reqs}_Q`` and ``B \subseteq E_Q``: it requires only
``A`` and leaves empty only ``B``, so the requirements outside ``A`` are
dropped and the classes in ``E_Q \setminus B`` are given back. A *solution*
of it is one covering ``A`` and using no class in ``B``.

!!! note "Observation (a withdrawal is an assumption subset)"
    In the instance built for ``Q``, "``p`` is installed" is the unit
    assumption ``p`` and "class ``c`` of ``p`` is empty" is the unit
    assumption ``\lnot p@c``. Every other clause — dependency edges,
    registry conflicts, the structural at-most-one — is the registry's and
    is the same for every withdrawal. So ``(A, B)`` is satisfiable exactly
    when the instance is satisfiable assuming
    ``\{p : p \in A\} \cup \{\lnot p@c : p@c \in B\}``, and asking is one
    solve.

Every relaxation ``R \le Q`` in the previous page's sense is a withdrawal,
but not conversely: nothing says some set of constraints picks out an
arbitrary ``B \subseteq E_Q``. The questions below stay on the relaxations,
and the reason is the *package* granularity they are asked at.

## A package at a time

A report is written a package at a time — what has become of a package's
versions is one thing to say, and which of its constraints to relax is one
thing to do — so the questions are asked a package at a time too, over one
variable per package implying all of that package's assumptions
(`with_emptied_packages`). Defining such a variable is a conservative
extension: it occurs positively only in what is assumed, so every model of
the instance extends to one of the instance with the definition.

That granularity is also what keeps the questions on relaxations rather than
on withdrawals at large, and it is worth saying why, because a per-class
question would not.

!!! note "Observation (a package's worth of ``E_Q`` is a constraint's worth)"
    Let ``S`` be a set of packages and ``E_S`` the classes of ``E_Q``
    belonging to them. Then the withdrawal ``(A, E_Q \setminus E_S)`` is the
    relaxation ``R \le Q`` that requires ``A`` and lifts, for every
    ``p \in S``, every kind of ``Q`` that excludes any version of ``p``.

    *Proof.* After lifting those kinds nothing of ``Q`` forbids any version
    of a ``p \in S``, so every class of every such ``p`` has an admitted
    member and none of them is empty. The kinds are lifted for ``S`` alone,
    so every other package's classes are emptied exactly as ``Q`` emptied
    them. Hence ``R``'s empty set is ``E_Q \setminus E_S``, and its
    requirements are ``A``. ∎

The diagnosis only ever assumes a subset of the requirement literals and a
subset of the package literals, so every question it puts to the solver is
of that form — a genuine point of the previous page's lattice. A question
about *one* class of a package would not be: no constraint need pick out one
of a package's empty classes and leave its siblings empty.

## The filtered universe answers for every relaxation

The instance a resolve fails on is built over a universe that has been
filtered — reachability and redundancy elimination, run with ``Q``'s
constraints in force. The questions above are asked of that instance, so
what has to hold is that filtering for ``Q`` did not delete anything a
relaxation of ``Q`` would need.

That is Theorem C of the previous page: ``\mathrm{Lay}(F_Q(D), R) =
\mathrm{Lay}(D, R)`` for every ``R \le Q``, which is stronger than the
satisfiability statement the questions here need and covers it. The filter
is run under the rules that make it true, so nothing further is required of
a diagnosis than to stay below ``Q`` in the order — which, by the
observation above, everything it asks does.

What the proof turns on is a design decision worth naming here too: a class
the query empties is **deactivated, not deleted**. It keeps its row, its
column in every partner, and its dependency columns, which is exactly what
a question about giving it back needs to find still there.

## Conflicts are the coordinates of the repairs

A report hands out one menu per conflict and invites the reader to choose
from each of them independently. What that claims is that *every*
combination of choices repairs the query, and it is not something checked
after the fact: it is how the conflicts are found in the first place.

Write ``F`` for the minimal correction sets of the assumptions above — every
set of them whose withdrawal is satisfiable and no proper subset of which is
(`sat_mcses`) — and ``F_{\min} \subseteq F`` for the ones of smallest size.
A member of ``F_{\min}`` gives up as few of the user's own facts as any
repair of the query can.

Only ``F_{\min}`` is ever enumerated, and the way to enumerate just it is to
refuse to look at anything bigger. A repair is the set of assumptions a model
leaves unsatisfied, so "no repair bigger than ``k``" is an at-most-``k``
constraint over their negations — a few clauses over counting variables — and
the bound turns the enumeration into the one that is wanted.

!!! note "Observation (a bound cuts the family off rather than changing it)"
    Let ``B_k`` be the instance together with a constraint that at most ``k``
    of the assumptions go unsatisfied. Then ``B_k`` is satisfiable exactly
    when ``F`` has a member of size ``k`` or less, and the minimal correction
    sets of ``B_k`` are exactly those members.

    *Proof.* A model of ``B_k`` is a model of the instance leaving at most
    ``k`` assumptions unsatisfied, and the set it leaves unsatisfied is a
    correction set, which contains a minimal one; conversely a member of
    ``F`` of size ``k`` or less is witnessed by a model that leaves exactly
    it unsatisfied, and that model satisfies ``B_k``. So the two are
    satisfiable together.

    Let ``C`` be a minimal correction set of ``B_k``, witnessed by a model
    ``\mu`` of ``B_k``. The assumptions ``\mu`` leaves unsatisfied are a
    subset of ``C`` and a correction set of ``B_k``, so by minimality they
    are ``C`` itself, and ``|C| \le k``. And ``C`` is minimal for the
    instance: anything strictly inside it is smaller than ``k`` too, so it is
    a correction set of ``B_k`` if it is one at all, which minimality
    forbids. The converse is the same statement read backwards — a correction
    set of the instance of size ``k`` or less is one of ``B_k``, and
    minimality for the instance is minimality for ``B_k`` since ``B_k`` has
    the fewer models. ∎

So raising ``k`` from 1 until ``B_k`` becomes satisfiable finds the size of a
smallest repair — ``k`` is 2 to 6 on the registry queries this was measured on,
so it is a handful of solves, and each counter is as small as the answer — and
enumerating at that bound enumerates ``F_{\min}`` and nothing else.

!!! note "Observation (a family that counts right is a product)"
    Let ``C_1, \dots, C_k`` partition the literals ``F_{\min}`` uses, and
    suppose every ``M \in F_{\min}`` contains exactly one literal of each
    ``C_i``, and that ``|F_{\min}| = \prod_i |C_i|``. Then
    ``F_{\min} = C_1 \times \dots \times C_k``.

    *Proof.* Sending each ``M`` to the tuple of its literals is well defined
    by the first hypothesis and injective because a set is determined by its
    elements. An injection between finite sets of equal size is a bijection,
    so every tuple is some ``M``. ∎

Both hypotheses are cheap to check, and the partition to check them against
is forced: literals of the same ``C_i`` never share a repair, so the
candidate partition is the connected components of "never share a repair",
and if any partition works that one does. So the report offers ``C_i`` as
its ``i``-th conflict and the literals of ``C_i`` as that conflict's menu.
Every path through the menus is a member of ``F_{\min}`` — a repair, and one
of the cheapest. Nothing is left to compose and nothing is left to warn
about.

It also follows that there are exactly as many conflicts as a cheapest
repair has literals. The conflict count is not an artifact of how the search
went: it is the size of the smallest fix.

### The story of a coordinate

A menu says what to choose between without saying why. Fix any
``M \in F_{\min}``, and for the ``i``-th coordinate let ``c_i`` be its
literal of ``M``. Withdraw ``M \setminus \{c_i\}`` — every *other*
coordinate settled the way ``M`` settles it — and ask what is left.

It is unsatisfiable: ``M`` is minimal, so no proper subset of it repairs
anything. A MUS of it therefore exists, and it must contain ``c_i``, since
withdrawing ``c_i`` as well is withdrawing all of ``M``, which does repair.
And it contains no other literal of ``M``, because those were withdrawn
before the question was asked. So the ``i``-th conflict gets a story about
its own choice and nobody else's.

A MUS is not the whole story, though, and the gap is one a report cannot
afford. Several requirements can fail for one and the same reason; one of
them is already enough to make that reason unsatisfiable, so a MUS names one
of them — arbitrarily — and drops the rest as redundant. They are not
redundant to the user: they are required, they are broken, and this
conflict's fix is what rescues them. A report that omitted them would let
someone relax the bound and never learn what else had been failing.

So each requirement the MUS left out is asked one further question: is it
satisfiable on its own, against the packages as this coordinate finds them
(that is, with every *other* coordinate's choice already made)? If it is
not, then it is one this coordinate rescues — because applying ``c_i`` on
top is applying all of ``M``, and ``M`` repairs the query, so every
requirement of the query is satisfiable afterwards. Each such requirement
brings its own MUS, and the chain is the union.

The chain is therefore a union of minimal conflicts rather than a single
one: unsatisfiable, with every fact in it belonging to at least one of them,
but not with every fact load-bearing for the whole. A conflict whose
requirements fail *together* is the special case where there is only one.

### What the menus leave out

Two things can be true beyond what the menus say, and they are independent:

* ``F \ne F_{\min}``: repairs exist that give up more than these do. With
  ``F_{\min}`` in hand this is decided exactly, by one solve.
* ``F_{\min}`` is not a product: then the best that can be offered is the
  largest *rectangle* in it — a set of literals every repair in the
  rectangle shares, together with the alternatives that complete it — and
  the repairs outside the rectangle are cheapest ones the menus never reach.

!!! note "Observation (one solve decides whether anything larger exists)"
    Let ``F_{\min} = \{M_1, \dots, M_t\}`` be *all* the repairs of smallest
    size, and add to the instance the clause ``\bigvee_{c \in M_i} c`` for
    each ``i``. The result is satisfiable exactly when ``F \ne F_{\min}``.

    *Proof.* Write ``U(\mu)`` for the assumptions a model ``\mu`` leaves
    unsatisfied; it is a correction set, and the clause for ``M_i`` says
    exactly that ``M_i \not\subseteq U(\mu)``. So a model of the augmented
    instance has a ``U(\mu)`` containing no ``M_i``; the minimal correction
    set inside it is therefore none of them, and ``F`` has a member outside
    ``F_{\min}``. Conversely let ``N \in F \setminus F_{\min}``. A model
    witnessing ``N`` leaves exactly ``N`` unsatisfied, and no ``M_i`` is
    inside ``N``, since two distinct minimal correction sets are
    incomparable — so that model satisfies every added clause. ∎

Note which way round the clause goes: an assumption is stated as the user
would want it to hold, so the assumptions a model *satisfies* are the ones the
repair it witnesses leaves alone. Ruling ``M_i`` out is therefore asking for
one of its literals, not against it.

A report says nothing when neither holds, "Larger solutions also exist."
when only the first does, and "Other solutions also exist." whenever the second
does. Exploring those other solutions is a question this page does not
answer.

## Versions come from a resolve, not from the instance

Theorem C is about the answer, but the questions above only ask about
satisfiability. A report also says what the user would *get*, and that is
the layered answer, which depends on the version ordering.

The ordering a universe is laid out in belongs to the query. A class ranks
where its best admitted member ranks (`class_ranking`), and a relaxation
admits more members — so under a relaxation a class can compete at a better
version than it competes at here. A version read off the failed instance is
therefore a version of the query's universe, not of the relaxation's.

So no version in a report is ever read off the failed instance. A fix is
turned back into the withdrawal it is — the requirements it names dropped,
the constraint kinds it names lifted for the packages they name — and
[`relax`](@ref Resolver.relax) derives the relaxed problem from the query
and computes what the relaxation makes of the classes the query laid out:
the member each class stands for, and the order those members rank the
classes in. Resolving that on the failed resolve's own instance runs the
descent again under the relaxation's ranking, and the solution that comes
back is what the report shows. It is the solution the user will get,
because Theorem C says it is what a resolve of the relaxed problem from
scratch would give.

That a fix resolves at all is that theorem plus monotonicity: relaxing a
constraint admits at least the version the repair asked for, and admitting
more versions only adds models.

## Attribution needs no license

Why a class is empty is not a question about the instance. A constraint
forbids versions, a class is empty when every member of it is forbidden,
and which kinds forbid which version is read straight off the `Problem`
(`exclusion_kinds`, `class_exclusions`). No solver is involved, so nothing
has to be preserved by anything:

```math
\mathrm{reps}[p][c] = 0
\quad\Longleftrightarrow\quad
\text{every member of class } c \text{ of } p \text{ has an excluding kind.}
```

This is the one place a report can name a *user's* constraint at all —
registry compatibility is in the conflict matrices, where it is
indistinguishable from any other conflict — and it is why a story can say
which bound emptied a class without asking anything.
