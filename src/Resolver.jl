module Resolver

export resolve

const PkgVer = Pair{String,VersionNumber}

function resolve(
    required::Vector{String},
    versions::Dict{String,Vector{VersionNumber}},
    dependencies::Dict{PkgVer,Vector{String}},
    conflicts::Set{Tuple{PkgVer,PkgVer}} = Set(NTuple{2,PkgVer}[]),
)
    # output data
    solutions = Vector{Pair{String,VersionNumber}}[]
    isempty(required) && return solutions

    no_version = typemin(VersionNumber)
    packages = sort!(collect(keys(versions)), by = p->(p ∉ required, p))
    choices = [copy(versions[p]) for p in packages]

    for (i, pkg) in enumerate(packages)
        (pkg in required ? push! : pushfirst!)(choices[i], no_version)
    end

    function conflict(a::PkgVer, b::PkgVer)
        depends(v::PkgVer, d::String) =
            haskey(dependencies, v) && d in dependencies[v]
        c = (a, b) in conflicts || (b, a) in conflicts ||
        b[2] == no_version && depends(a, b[1]) ||
        a[2] == no_version && depends(b, a[1])
        # @show a, b, c
        return c
    end

    n = length(packages)
    package = fill(0, n)        # i    : package assigned at index i
    version = fill(1, n, n)     # i, p : best version of p still possible
    reached = fill(false, n, n) # i, p : optimal[i,p] already reached

    function search!(i::Int)
        for p = 1:n
            # skip already assigned packages
            any(p == package[j] for j = 1:i-1) && continue
            # record ith package choice
            package[i] = p
            if i == n
                done = true
            else
                # optimal still-compatible version of p
                v = version[i, p]
                a = PkgVer(packages[p], choices[p][v])
                # optimal still-compatbile next versions
                done = false
                for q = 1:n
                    version[i+1, q] = 0
                    for w = version[i, q]:length(choices[q])
                        b = PkgVer(packages[q], choices[q][w])
                        conflict(a, b) && continue
                        version[i+1, q] = w
                        break
                    end
                    (done |= version[i+1, q] == 0) && break
                end
            end
            if done
                if i == n
                    @show package, version
                    # record solution
                    solution = Pair{String,VersionNumber}[]
                    for p = 1:n
                        v = choices[p][version[i, p]]
                        v ≠ no_version && push!(solution, packages[p] => v)
                    end
                    @show solution
                    push!(solutions, solution)
                    return # last p is unique
                end
            else
                search!(i+1)
            end
        end
    end
    search!(1)

    return solutions
end

end # module
