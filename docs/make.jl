using Documenter
using Resolver

makedocs(
    sitename = "Resolver.jl",
    modules = [Resolver],
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
        ],
    ],
)

deploydocs(
    repo = "github.com/StefanKarpinski/Resolver.jl",
    devbranch = "main",
    push_preview = true,
)
