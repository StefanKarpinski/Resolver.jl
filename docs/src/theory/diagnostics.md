# Diagnosing an unsatisfiable resolve

When the requirements cannot be satisfied there is more to say than that
they cannot: which of them cannot hold together, why, and what a user
could change so that they can. Every one of those answers comes from
putting further satisfiability questions to the instance the resolve
failed on, and this page states what licenses them.

Everything here is short, and class space (see the
[previous page](layered.md#Preprocessing-for-caching)) is why: a user
constraint's whole effect on the universe is the set of classes it
empties, so there is exactly one kind of thing to relax and it is already
a unit fact of the instance.

## Relaxations

A query is a set of requirements ``R`` together with constraints, and what
the constraints do to the universe is empty a set ``E`` of classes — a
class the query admits no member of cannot be selected, and that is the
whole of it.

A **relaxation** of the query is a pair ``(R', E')`` with
``R' \subseteq R`` and ``E' \subseteq E``: it requires only ``R'`` and
leaves empty only ``E'``, so the requirements outside ``R'`` are dropped
and the classes in ``E \setminus E'`` are given back. The query itself is
``(R, E)``, and a *solution* of the relaxation is one covering ``R'`` and
using no class in ``E'``.

!!! note "Observation (relaxations are assumption subsets)"
    In the instance built for the query, "``p`` is installed" is the unit
    assumption ``p`` and "class ``c`` of ``p`` is empty" is the unit
    assumption ``\lnot p@c``. Every other clause — dependency edges,
    registry conflicts, the structural at-most-one — is the registry's and
    is the same for every relaxation. So the relaxation ``(R', E')`` is
    satisfiable exactly when the instance is satisfiable assuming
    ``\{p : p \in R'\} \cup \{\lnot p@c : p@c \in E'\}``, and asking is one
    solve.

A report is written a package at a time — what has become of a package's
versions is one thing to say, and which of its constraints to relax is one
thing to do — so the questions are asked a package at a time too, over one
variable per package implying all of that package's assumptions. Defining
such a variable is a conservative extension: it occurs positively only in
what is assumed, so every model of the instance extends to one of the
instance with the definition, and the questions stay inside the class this
page covers.

## The filtered universe answers for every relaxation

The instance a resolve fails on is built over a universe that has been
filtered — reachability and redundancy elimination, run with the query's
constraints in force (`filter_pkg_info!`). The questions above are asked
of that instance, so what has to hold is that filtering for the query did
not delete anything a relaxation of the query would need.

!!! note "Proposition (filtering answers for every relaxation)"
    Let ``F`` be the universe obtained by filtering the raw universe for
    the query ``(R, E)``. For every relaxation ``(R', E')``, ``F`` has a
    solution of that relaxation iff the raw universe does.

    *Proof.* One direction is immediate: ``F``'s classes are a subset of
    the raw universe's, so its solutions are raw solutions.

    For the other, take a raw solution ``S`` of the relaxation.
    Reachability continues a package's prefix past a class the
    query emptied, since such a class is not somewhere the package can be
    installed. So the last kept class of an unsaturated package is not in
    ``E`` — had it been, the prefix would have continued past it — and is
    therefore admitted by the query, hence by the relaxation. Theorem 3's
    construction now applies verbatim: take last kept classes for
    unsaturated packages and ``S``'s own classes for saturated ones (which
    keep all of theirs), and Theorem 2's invariants make the splice
    conflict-free and dependency-closed. It covers ``R'``, and it uses no
    class in ``E'``: the unsaturated half uses none in ``E`` at all, and
    the saturated half takes its classes from ``S``.

    Redundancy elimination is handled by substitution, as in Theorem 6. It
    deletes a class only when a better class has a subset of its
    *registry* constraints — its dependency and conflict columns — and that
    comparison never mentions the query's constraints, so it holds for the
    relaxation unchanged. The dominating class is a candidate of the pass,
    and candidates exclude the classes the query emptied, so the
    relaxation admits it. Domination chains terminate at a class that
    survives, because reachability keeps prefixes. ∎

Both halves turn on the same design decision: a class the query empties is
**deactivated, not deleted**. It keeps its row, its column in every
partner, and its dependency columns, which is exactly what a question
about giving it back needs to find still there.

## Versions come from a resolve, not from the instance

The proposition is about satisfiability, and satisfiability is all the
questions above ask. A report also says what the user would *get*, and
that is a different quantity: the layered answer, which depends on the
version ordering.

The ordering a universe is laid out in belongs to the query. A class ranks
where its best admitted member ranks (`class_ranking`), and a relaxation
admits more members — so under a relaxation a class can compete at a
better version than it competes at here, and both filtering passes, which
are stated in terms of the rank order, would be reasoning about a
differently ordered universe. A version read off the failed instance is
therefore a version of the query's universe, not of the relaxation's.

So no version in a report is ever read off the failed instance. A fix is
carried out on the `Problem` — the requirements it names dropped, the
constraint sources it names taken off the packages they name — and the
modified problem is resolved against the T1 artifact
([`pkg_info`](@ref Resolver.pkg_info)), where the relaxation's own
representatives and its own class layout are computed from scratch. The
solution that comes back is what the report shows, and it is the solution
the user will get, because it is the same call the user's next resolve
makes.

That a fix resolves at all is the proposition above plus monotonicity:
relaxing a constraint source admits at least the version the repair asked
for, and admitting more versions only adds models.

## Attribution needs no license

Why a class is empty is not a question about the instance. A constraint
forbids versions, a class is empty when every member of it is forbidden,
and which sources forbid which version is read straight off the `Problem`
(`exclusion_sources`, `class_exclusions`). No solver is involved, so
nothing has to be preserved by anything:

```math
\mathrm{reps}[p][c] = 0
\quad\Longleftrightarrow\quad
\text{every member of class } c \text{ of } p \text{ has an excluding source.}
```

This is the one place a report can name a *user's* constraint at all —
registry compatibility is in the conflict matrices, where it is
indistinguishable from any other conflict — and it is why a story can say
which bound emptied a class without asking anything.
