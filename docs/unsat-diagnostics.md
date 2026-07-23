# Unsatisfiable diagnostics

When a set of requirements cannot be resolved, Resolver.jl can explain **why**
and list **verified fixes**. Rather than a bare "Unsatisfiable", it reports:

- **Conflicts** — minimal sets of requirements that cannot hold together (the
  *minimal unsatisfiable subsets*, or MUSes), rendered as a short chain of
  facts that walks from what you asked for to the incompatibility.
- **Verified fixes** — minimal changes that make the problem resolvable (the
  *minimal correction sets*, or MCSes), each checked by actually re-resolving,
  and shown with the versions it unlocks.
- **Upstream fixes** — where the resolution would be unblocked by a *new
  release* of a dependency that relaxed one of its own compat bounds.

Each fix is built from one or more of these actions:

| action | meaning |
|---|---|
| `drop requirement X` | stop requiring `X` |
| `relax your compat on X` | loosen a `[compat]` entry *you* imposed on `X` |
| `unpin X` | lift a pin holding `X` at one version |
| *upstream* | a dependency publishes a release with looser compat |

A resolve can have several **independent** conflicts at once; a single fix then
combines one action per conflict (joined with "and").

There are two ways to use it: the `diagnose` function (a library API over
abstract package/version data), and the `bin/resolve.jl` command line (which
resolves a real project against the registries and diagnoses automatically when
it fails).

## Synthetic examples

These use small hand-written data — packages `A`, `B`, … with integer versions
— so they reproduce with no registry. `PkgData(versions, depends, compat)`
describes one package: its versions, each version's dependencies, and each
version's compat bounds (allowed versions of a dependency).

```julia
using Resolver: PkgData, diagnose
D(v, deps, comp) = PkgData(v, deps, comp)
ND = Dict{Int,Vector{String}}            # empty "depends"
NC = Dict{Int,Dict{String,Vector{Int}}}  # empty "compat"
```

### drop a requirement

`A` needs `C` version 1, `B` needs `C` version 2, and `C` can only be one thing:

```julia
data = Dict(
    "A" => D([1], Dict(1 => ["C"]), Dict(1 => Dict("C" => [1]))),
    "B" => D([1], Dict(1 => ["C"]), Dict(1 => Dict("C" => [2]))),
    "C" => D([1, 2], ND(), NC()),
)
diagnose(data, ["A", "B"])
```

```
Conflict 1: A, B cannot be satisfied together.
  • you require A
  • you require B
  • A (all versions) requires C at 1
  • B (all versions) requires C at 2

Verified fixes:
  1. drop requirement B
     → allows: A 1, C 1
  2. drop requirement A
     → allows: B 1, C 2

Upstream fixes:
  • a release of A relaxing its compat on C
    → allows: A 1, B 1, C 2
  • a release of B relaxing its compat on C
    → allows: A 1, B 1, C 1
```

The `allows:` line shows only the versions the fix changes that are relevant to
the conflict (the requirements and the contested dependency), not the whole
resolved manifest.

### relax your compat

Here `A` needs `C` version 2, but you've restricted `C` to version 1 via
`compat`:

```julia
data = Dict(
    "A" => D([1], Dict(1 => ["C"]), Dict(1 => Dict("C" => [2]))),
    "C" => D([1, 2], ND(), NC()),
)
diagnose(data, ["A"]; compat = Dict("C" => [1]))
```

```
Conflict 1: A cannot be satisfied.
  • you require A
  • A (all versions) requires C at 2
  • your compat restricts C to 1

Verified fixes:
  1. relax your compat on C
     → allows: A 1, C 2
  2. drop requirement A
     → allows: nothing
```

### unpin a package

`C` is pinned at version 1, but `A` needs version 2:

```julia
data = Dict(
    "A" => D([1], Dict(1 => ["C"]), Dict(1 => Dict("C" => [2]))),
    "C" => D([1, 2], ND(), NC()),
)
diagnose(data, ["A"]; holds = Dict("C" => 1))
```

```
Conflict 1: A cannot be satisfied.
  • you require A
  • A (all versions) requires C at 2
  • C is pinned at 1

Verified fixes:
  1. unpin C
     → allows: A 1, C 2
  2. drop requirement A
     → allows: nothing
```

### several independent conflicts

`A`/`B` disagree about `C`, and independently `X`/`Y` disagree about `Z`. The
two conflicts are reported separately, and each verified fix combines one action
per conflict:

```julia
data = Dict(
    "A" => D([1], Dict(1 => ["C"]), Dict(1 => Dict("C" => [1]))),
    "B" => D([1], Dict(1 => ["C"]), Dict(1 => Dict("C" => [2]))),
    "C" => D([1, 2], ND(), NC()),
    "X" => D([1], Dict(1 => ["Z"]), Dict(1 => Dict("Z" => [1]))),
    "Y" => D([1], Dict(1 => ["Z"]), Dict(1 => Dict("Z" => [2]))),
    "Z" => D([1, 2], ND(), NC()),
)
diagnose(data, ["A", "B", "X", "Y"])
```

```
Conflict 1: A, B cannot be satisfied together.
  • you require A
  • you require B
  • A (all versions) requires C at 1
  • B (all versions) requires C at 2

Conflict 2: X, Y cannot be satisfied together.
  • you require X
  • you require Y
  • X (all versions) requires Z at 1
  • Y (all versions) requires Z at 2

Verified fixes:
  1. drop requirement B and drop requirement Y
     → allows: A 1, X 1, C 1, Z 1
  2. drop requirement B and drop requirement X
     → allows: A 1, Y 1, C 1, Z 2
  3. drop requirement A and drop requirement Y
     → allows: B 1, X 1, C 2, Z 1
  4. drop requirement A and drop requirement X
     → allows: B 1, Y 1, C 2, Z 2

Upstream fixes:
  • a release of A relaxing its compat on C  → allows: A 1, B 1, C 2
  • a release of B relaxing its compat on C  → allows: A 1, B 1, C 1
  • a release of X relaxing its compat on Z  → allows: X 1, Y 1, Z 2
  • a release of Y relaxing its compat on Z  → allows: X 1, Y 1, Z 1
```

## Real-registry examples (command line)

`bin/resolve.jl` resolves a project against the installed registries and, when
resolution fails, prints the same diagnosis automatically before exiting
nonzero.

The examples below are **reproducible**: they were produced against the General
registry at commit
[`6bb1cf9538b1db31c2f4e4517ede70bc4c44227e`](https://github.com/JuliaRegistries/General/commit/6bb1cf9538b1db31c2f4e4517ede70bc4c44227e).
To reproduce exactly, check that commit out and point `--registry` at it:

```console
$ git clone https://github.com/JuliaRegistries/General
$ git -C General checkout 6bb1cf9538b1db31c2f4e4517ede70bc4c44227e
$ julia --project=bin -e 'using Pkg; Pkg.instantiate()'
```

(Omit `--registry=./General` to resolve against your currently installed
registries instead; the diagnosis is the same feature, but versions may differ
as the registry evolves.)

The conflicts below all come from the same real split: many packages still
require `HTTP` 0.x while current ones require `HTTP` 1.x, and only one `HTTP`
version can be chosen.

### drop a requirement

`BulkSMS` (an `HTTP` 0.x package) and `AnthropicClient` (an `HTTP` 1.x package)
can't coexist:

```toml
# Project.toml
[deps]
BulkSMS = "5fa4ca3b-1c55-51a6-87f0-d5811ecf545c"
AnthropicClient = "e82deec3-d3b5-4b85-aba1-b8e1f470db1f"
```

```console
$ julia --project=bin bin/resolve.jl myproject --julia=1.11 --registry=./General
Unresolvable. Diagnosis:

Conflict 1: BulkSMS, AnthropicClient cannot be satisfied together.
  • you require BulkSMS
  • you require AnthropicClient
  • BulkSMS (all versions) requires HTTP at 0.6.14–0.8.19
  • AnthropicClient (all versions) requires HTTP at 1.0.0–1.11.0

Verified fixes:
  1. drop requirement AnthropicClient
     → allows: BulkSMS 0.0.1, HTTP 0.8.12
  2. drop requirement BulkSMS
     → allows: AnthropicClient 0.1.0, HTTP 1.0.0

Upstream fixes:
  • a release of BulkSMS relaxing its compat on HTTP
    → allows: BulkSMS 0.0.1, AnthropicClient 0.1.0, HTTP 1.0.0
  • a release of AnthropicClient relaxing its compat on HTTP
    → allows: BulkSMS 0.0.1, AnthropicClient 0.1.0, HTTP 0.8.12
```

### relax your compat

A project that needs `AnthropicClient` (which requires `HTTP` 1.x) but carries a
stale `[compat] HTTP = "0"`:

```toml
# Project.toml
[deps]
AnthropicClient = "e82deec3-d3b5-4b85-aba1-b8e1f470db1f"
HTTP = "cd3eb016-35fb-5094-929b-558a96fad6f3"

[compat]
HTTP = "0"
```

```console
$ julia --project=bin bin/resolve.jl myproject --julia=1.11 --registry=./General
Unresolvable. Diagnosis:

Conflict 1: AnthropicClient cannot be satisfied.
  • you require AnthropicClient
  • AnthropicClient (all versions) requires HTTP at 1.0.0–1.11.0
  • your compat restricts HTTP to 0.6.10–0.9.17

Verified fixes:
  1. relax your compat on HTTP
     → allows: HTTP 1.0.0, AnthropicClient 0.1.0
  2. drop requirement AnthropicClient
     → allows: HTTP 0.8.12

Upstream fixes:
  • a release of AnthropicClient relaxing its compat on HTTP
    → allows: HTTP 0.8.12, AnthropicClient 0.1.0
```

### unpin a package

Pin `HTTP` to an old version with `--pin` (repeatable, comma-separated; the
leading `v` is optional), then require a package that needs a newer one:

```toml
# Project.toml
[deps]
AnthropicClient = "e82deec3-d3b5-4b85-aba1-b8e1f470db1f"
```

```console
$ julia --project=bin bin/resolve.jl myproject --julia=1.11 --registry=./General --pin=HTTP@0.6.10
Unresolvable. Diagnosis:

Conflict 1: AnthropicClient cannot be satisfied.
  • you require AnthropicClient
  • AnthropicClient (all versions) requires HTTP at 1.0.0–1.11.0
  • HTTP is pinned at 0.6.10

Verified fixes:
  1. unpin HTTP
     → allows: AnthropicClient 0.1.0, HTTP 1.0.0
  2. drop requirement AnthropicClient
     → allows: nothing
```
