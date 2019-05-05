module Resolver

export resolve, PkgVer

struct PkgVer
    package::String
    version::VersionNumber
end

function resolve(
    required::Vector{String},
    versions::Dict{String,Vector{VersionNumber}},
    dependencies::Dict{PkgVer,Vector{String}},
    conflicts::Set{Tuple{PkgVer,PkgVer}},
)
    no_version = typemin(VersionNumber)
    packages = sort!(collect(keys(versions)), by = p->(p ∉ required, p))
    choices = [copy(versions[p]) for p in packages]

    for (i, pkg) in enumerate(packages)
        (pkg in required ? push! : pushfirst!)(choices[i], no_version)
    end

    function conflict(a::PkgVer, b::PkgVer)
        depends(v::PkgVer, d::String) =
            haskey(dependencies, v) && d in dependencies[v]
        (a, b) in conflicts || (b, a) in conflicts ||
        b.version == no_version && depends(a, b.package) ||
        a.version == no_version && depends(b, a.package)
    end

    function conflict(i::Int, j::Int)
        i > 0 && j > 0 || return false
        a = PkgVer(packages[i], choices[i][assigned[i]])
        b = PkgVer(packages[j], choices[j][assigned[j]])
        return conflict(a, b)
    end

    # working data structures
    n = length(packages)
    order = collect(0:n)     # indices into packages, zero is skipped
    assigned = zeros(Int, n) # indices into choices
    conflicted = [1]         # previously conflicted packages
    prioritized = 0          # already-prioritied packages, index into conflicted

    # output data
    solutions = Vector{Pair{String,VersionNumber}}[]

    while prioritized < length(conflicted)
        p₀ = conflicted[prioritized += 1]
        order[1], order[p₀+1] = order[p₀+1], order[1]
        feasible = true
        for (i, p) in enumerate(order)
            p == 0 && continue
            k₁ = 1 + (i == 1 && choices[p][1] == no_version)
            found = false
            for k = k₁:length(choices[p])
                assigned[p] = k
                if !any(conflict(p, order[j]) for j=1:i-1)
                    found = true
                    break
                end
                # for k > 2 this will already have been done
                if k ≤ 2 && choices[p][k] ≠ no_version && p ∉ conflicted
                    push!(conflicted, p)
                end
            end
            (feasible &= found) || break
        end
        order[1], order[p₀+1] = order[p₀+1], order[1]
        feasible || continue
        solution = Pair{String,VersionNumber}[]
        for i = 1:n
            v = choices[i][assigned[i]]
            v ≠ no_version && push!(solution, packages[i] => v)
        end
        push!(solutions, solution)
    end

    return solutions
end

end # module
