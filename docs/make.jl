using Documenter
using Resolver

makedocs(
    sitename = "Resolver.jl",
    modules = [Resolver],
    # the UnsatCores submodule is internal machinery for the diagnostics: its
    # docstrings document its own API for the code that consumes it, and belong
    # no more in this manual than the rest of src/ does. without this, checkdocs
    # demands an @docs block for every one of them
    checkdocs_ignored_modules = [Resolver.UnsatCores],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://karpinski.org/Resolver.jl",
    ),
    pages = [
        "Home" => "index.md",
        "API" => "api.md",
        "Theory" => [
            "Setting" => "theory/setting.md",
            "The layered solution & filtering" => "theory/layered.md",
            "Pareto front optimality" => "theory/front.md",
            "Relaxation-stable filtering" => "theory/relaxation.md",
            "Explaining an unsatisfiable resolve" => "theory/unsat-explanation.md",
        ],
    ],
)

deploydocs(
    repo = "github.com/StefanKarpinski/Resolver.jl",
    devbranch = "main",
    push_preview = true,
)
