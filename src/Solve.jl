# variables for each package:
#  - one for the package
#  - one for each version
#
# clauses for each package:
#  - one/zero if the package is required or not
#  - one specifying package versions
#  - one for each conflict
#  - optimality clauses, one per version

function gen_sat(
    out  :: Union{ArgWrite, Nothing},
    data :: Dict{P, PkgEntry{P}},
    reqs :: Vector{P},
) where {P}
    arg_write(out) do out
        # compute & output header
        names = sort!(collect(keys(data)))
        var = Dict{P, Int}() # variable indices
        v = 0 # number of variables
        x = 0 # number of conflicts
        for p in names
            var[p] = v + 1
            v += length(data[p].versions) + 1
            d = length(data[p].depends)
            x += sum(data[p].conflicts)
            x += sum(@view(data[p].conflicts[:, 1:d])) # count deps twice
        end
        # conflicts are double-counted
        @assert iseven(x)
        x >>= 1
        # number of clauses
        c = length(reqs) + v + x
        println(out, "p cnf $v $c")
        println(out)
        # output requirements clauses
        for p in sort(reqs)
            println(out, "$(var[p]) 0")
        end
        println(out)
        # output package version clauses
        for p in names
            print(out, "-$(var[p]) ")
            for i = 1:length(data[p].versions)
                print(out, "$(var[p]+i) ")
            end
            println(out, "0")
        end
        println(out)
        # output dependency clauses
        for p in names
            for i = 1:length(data[p].versions)
                for (j, q) in enumerate(data[p].depends)
                    data[p].conflicts[i, j] || continue
                    println(out, "-$(var[p]+i) $(var[q]) 0")
                end
            end
        end
        println(out)
        # output incompatibility clauses
        for p in names
            for q in sort!(collect(keys(data[p].interacts)))
                q < p || break
                b = data[p].interacts[q]
                for i = 1:length(data[p].versions),
                    j = 1:length(data[q].versions)
                    data[p].conflicts[i, b+j] || continue
                    println(out, "-$(var[p]+i) -$(var[q]+j) 0")
                end
            end
        end
        println(out)
        # output optimality clauses
        for p in names
            for i = 1:length(data[p].versions)
                print(out, "-$(var[p]) $(var[p]+i) ")
                # can be "excused" if a better version is chosen
                for j = 1:i-1
                    print(out, "$(var[p]+j) ")
                end
                # or if it conflicts with something else chosen
                for q in sort!(collect(keys(data[p].interacts)))
                    b = data[p].interacts[q]
                    for j = 1:length(data[q].versions)
                        data[p].conflicts[i, b+j] || continue
                        print(out, "$(var[q]+j) ")
                    end
                end
                println(out, "0")
            end
        end
    end
end

function gen_sat(
    data :: Dict{P, PkgEntry{P}},
    reqs :: Vector{P},
) where {P}
    gen_sat(nothing, data, reqs)
end

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

solve(reqs_str::AbstractString) = solve(String.(split(reqs_str, r"\s*,\s*")))

function solve(reqs::AbstractVector{<:AbstractString})
    to = TimerOutput()
    @timeit to "total" begin
        @timeit to "collect packages" pkgs = find_packages(dp, reqs)
        @timeit to "filter reachable" filter_reachable!(pkgs, reqs)
        @timeit to "filter redundant" filter_redundant!(pkgs)
        @timeit to "resolve versions" sat = solve(pkgs, reqs)
    end
    show(to)
    return sat
end

const picosat = expanduser("~/dev/picosat/picosat")
const picomus = expanduser("~/dev/picosat/picomus")

function solve(
    pkgs :: Dict{P, PkgInfo{P,V,S}},
    reqs :: AbstractVector{<:AbstractString},
) where {P,V,S}
    # generate & solve problem
    problem = gen_sat("tmp/problem.cnf", pkgs, reqs)
    output = "tmp/output.txt"
    open(output, write=true) do io
        run(pipeline(ignorestatus(`$picomus $problem`), stdout=io))
    end

    # decode output
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
    return sat
end

using LinearAlgebra: checksquare

function gen_max_clique_sat(
    out   :: Union{ArgWrite, Nothing},
    G     :: AbstractMatrix{Bool},
    parts :: Vector{Vector{Int}},
    force :: Vector{Int} = Int[],
)
    arg_write(out) do out
        # output the header
        n = checksquare(G)
        c = nnz(G)÷2 + length(parts) + length(force)
        println(out, "p cnf $n $c")
        println(out)
        # output vertex partitions (packages)
        for part in parts
            for i in part
                print(out, "$i ")
            end
            println(out, "0")
        end
        println(out)
        # output graph structure
        for (i, j, v) in zip(findnz(G)...)
            i < j && v || continue
            println(out, "-$i -$j 0")
        end
        isempty(force) || println(out)
        # output forced nodes
        for i in force
            println(out, "$i 0")
        end
    end
end

function gen_max_clique_sat(
    G     :: AbstractMatrix{Bool},
    parts :: Vector{Vector{Int}},
    force :: Vector{Int} = Int[],
)
    gen_max_clique_sat(nothing, G, parts, force)
end

function parse_picosat_output(io::ArgRead)
    arg_read(io) do io
        sat = nothing
        vals = Int[]
        while !eof(io)
            line = readline(io)
            if line == "s SATISFIABLE"
                sat = true
            elseif line == "s UNSATISFIABLE"
                sat = false
            elseif startswith(line, "v ")
                val = parse(Int, chop(line, head=2, tail=0))
                val != 0 && push!(vals, val)
            end
        end
        sat === nothing && error("error parsing picosat output")
        return sat, vals
    end
end

function solve_max_clique(
    G     :: AbstractMatrix{Bool},
    parts :: Vector{Vector{Int}},
    force :: Vector{Int} = Int[],
)
    problem = "tmp/problem.cnf"
    output = "tmp/output.txt"

    problem = gen_max_clique_sat(problem, G, parts, force)
    open(output, write=true) do io
        run(pipeline(ignorestatus(`$picomus $problem`), stdout=io))
    end
    sat, sol = parse_picosat_output(output)
    filter!(>(0), sol)
    return sat, sol
end

function complete_graph!(
    G::AbstractMatrix{Bool},
    parts::Vector{Vector{Int}},
)
    n = checksquare(G)
    C = Set{Tuple{Int,Int}}() # compat set
    for i = 1:n-1, j = i+1:n
        G[i, j] && continue # explicit incompatibility
        (i, j) in C && continue # known compatible
        sat, sol = solve_max_clique(G, parts, [i, j])
        if !sat
            # implicit incompatibility
            G[i, j] = G[j, i] = true
            continue
        end
        # solution including i & j
        # println(repr(sol))
        for k = 1:n
            sum(G[sol, k]) > 1 && continue
            # 0 or 1 incompatibilities with solution
            # 0 => k is already in the solution
            # 1 => k is swappable into solution
            for l in sol
                l < k && push!(C, (l, k))
                k < l && push!(C, (k, l))
            end
        end
    end
    return G
end

function complete_graph(
    G::AbstractMatrix{Bool},
    parts::Vector{Vector{Int}},
)
    complete_graph!(copy(G), parts)
end
