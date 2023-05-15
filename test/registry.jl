using Pkg: depots1
using Pkg.Registry: RegistryInstance, init_package_info!
using Pkg.Types: stdlibs
using Pkg.Versions: VersionSpec

const reg_path = joinpath(depots1(), "registries", "General")
const reg_inst = RegistryInstance(expanduser(reg_path))
const reg_dict = Dict(p.name => p for p in values(reg_inst.pkgs))
const excludes = push!(Set(first.(values(stdlibs()))), "julia")

dp = DepsProvider{String, VersionNumber, VersionSpec}() do pkg::String
    info = init_package_info!(reg_dict[pkg])
    vers = sort!(collect(keys(info.version_info)), rev=true)
    deps = Dict(v => String[] for v in vers)
    comp = Dict(v => Dict{String,VersionSpec}() for v in vers)
    # scan versions and populate deps & compat data
    for v in vers
        for (r, d) in info.deps
            v in r && union!(deps[v], keys(d))
        end
        for (r, c) in info.compat
            v in r && mergewith!(intersect, comp[v], c)
        end
    end
    foreach(sort!, values(deps))
    # scrub out excluded deps (stdlibs, julia itself)
    for d in values(deps)
        setdiff!(d, excludes)
    end
    for c in values(comp), x in excludes
        delete!(c, x)
    end
    # deduplicate data structures to save memory
    for i = 1:length(vers)-1, j = i+1:length(vers)
        v, w = vers[i], vers[j]
        deps[v] == deps[w] && (deps[v] = deps[w])
        comp[v] == comp[w] && (comp[v] = comp[w])
    end
    # return resolver PkgInfo data structure
    PkgInfo{String, VersionNumber, VersionSpec}(vers, deps, comp)
end

#=
all_names = sort!(collect(keys(reg_dict)))
filter!(!endswith("_jll"), all_names)
filter!(!in(excludes), all_names)
all_pkgs = find_packages(dp, all_names)
filter_reachable!(all_pkgs, all_names)
filter_redundant!(all_pkgs)
all_ix = find_interacts(all_pkgs)

const pairs = Tuple{String,String}[]
for p in all_names, q in get(all_ix, p, String[])
    p < q || continue
    reqs = [p, q]
    pkgs = find_packages(dp, reqs)
    filter_reachable!(pkgs, reqs)
    filter_redundant!(pkgs)
    all(length(pkgs[p].versions) > 1 for p in reqs) || continue
    push!(pairs, (p, q))
    @show reqs
end
=#

const picosat = expanduser("~/dev/picosat/picosat")
const picomus = expanduser("~/dev/picosat/picomus")

function decode_line(vars::Vector{String}, line::AbstractString)
    words = String.(split(line))
    if !isempty(words) && words[1] ∉ ("c", "p")
        words = map(words) do word
            i = tryparse(Int, word)
            i === nothing && return word
            i > 0 && return vars[i]
            i < 0 && return "!$(vars[-i])"
            return ""
        end
    end
    return strip(join(words, " "))
end

using TimerOutputs

function solve(reqs_str::AbstractString)
    to = TimerOutput()
    @timeit to "solve" begin
    reqs = String.(split(reqs_str, ','))
    @timeit to "collect packages" pkgs = find_packages(dp, reqs)
    @timeit to "filter reachable" filter_reachable!(pkgs, reqs)
    @timeit to "filter redundant" filter_redundant!(pkgs)

    # generate & solve problem
    @timeit to "generate SAT" problem = gen_sat("tmp/problem.cnf", pkgs, reqs)
    output = "tmp/output.txt"
    @timeit to "solve SAT" open(output, write=true) do io
        run(pipeline(ignorestatus(`$picomus $problem`), stdout=io))
    end

    # decode output
    @timeit to "decode" begin
    vars = String[]
    for p in sort!(collect(keys(pkgs)))
        push!(vars, "$p")
        for v in pkgs[p].versions
            push!(vars, "$p=$v")
        end
    end
    sat = nothing
    open(output) do io
        while !eof(io)
            line = readline(io)
            if line == "s SATISFIABLE"
                println("SATISFIABLE")
                sat = true
                seen = Set{String}()
                while !eof(io)
                    line = readline(io)
                    startswith(line, "v ") || continue
                    line = chop(line, head=2, tail=0)
                    line = decode_line(vars, line)
                    startswith(line, "!") && continue
                    contains(line, "=") || continue
                    pkg = String(split(line, "=")[1])
                    pkg in seen && continue
                    println(line)
                    push!(seen, pkg)
                end
            elseif line == "s UNSATISFIABLE"
                println("UNSATISFIABLE")
                sat = false
                core = Int[]
                while !eof(io)
                    line = readline(io)
                    startswith(line, "v ") || continue
                    line = chop(line, head=2, tail=0)
                    i = parse(Int, line)
                    i == 0 && break
                    push!(core, i)
                end
                lines = readlines(problem)
                popfirst!(lines)
                filter!(!isempty, lines)
                for line in lines[core]
                    println(decode_line(vars, line))
                end
            end
        end
    end
    end # @timeit "decode"
    end # @timeit "solve"
    show(to)
end
