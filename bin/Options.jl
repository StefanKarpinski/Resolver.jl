module Options

export USAGE, usage, parse_opts!, handle_opts

const USAGE = Base.RefValue{String}()

function usage()
    println(stdout, USAGE[])
    exit(0)
end

function usage(msg)
    println(stderr, "[ERROR] $msg\n\n$(USAGE[])")
    exit(1)
end

const OPTS = Vector{Pair{Symbol,Union{String,Nothing}}}()

"""
Process `args`, checking for and extracting options. The `opts` vector is a list
of valid option names, which must contain only ASCII letters, numbers and dashes.
"""
function parse_opts!(
    args :: AbstractVector{<:AbstractString},
    opts :: AbstractVector{<:AbstractString},
)
    # check option names
    for opt in opts
        contains(opt, r"^[a-z0-9-]+$"i) ||
            error("Internal error: bad option name: " * repr(opt))
    end
    # process arguments
    i = 1
    while i ≤ length(args)
        arg = args[i]
        if arg[1] == '-'
            if arg == "--"
                deleteat!(args, i)
                break
            end
            arg in ("-h", "--help") && usage()
            m = match(r"^--([a-z-]+)(?:=(.+))?$", arg)
            !isnothing(m) && m[1] in opts || usage("Invalid option: $arg")
            opt = Symbol(replace(m[1], '-' => '_'))
            push!(OPTS, opt => m[2])
            deleteat!(args, i)
        else
            i += 1
        end
    end
end

function _handle_opts(body::Function, pred::Function, value::Any)
    for (opt, val) in OPTS
        pred(opt) || continue
        # check for unexpected arguments
        if val isa String
            hasmethod(body, Tuple{Symbol,String}) ||
            hasmethod(body, Tuple{String}) ||
            usage("Option $opt does not take an argument (got $(repr(val)))")
        end
        # check for required arguments
        if val isa Nothing
            hasmethod(body, Tuple{Symbol,Nothing}) ||
            hasmethod(body, Tuple{Nothing}) ||
            hasmethod(body, Tuple{Symbol}) ||
            hasmethod(body, Tuple{}) ||
                usage("Option $opt requires an argument")
        end
        # call the first applicable method of body
        value =
            applicable(body, opt, val) ? body(opt, val) :
            applicable(body, val) ? body(val) :
            applicable(body, opt) ? body(opt) :
            applicable(body) ? body() :
                error("Internal error: invalid handle_opts body")
    end
    return value
end

"""
    handle_opts(body::Function, opt::Symbol, value::Any=nothing)
    handle_opts(body::Function, re::Regex, value::Any=nothing)

For any option matching `name` or `re`, based on what methods `body` has, call
the first of the following methods that is applicable:
```
body(opt::Symbol, val::Union{String,Nothing})
body(val::Union{String,Nothing})
body(opt::Symbol)
body()
```
If `val` is a string and no method accepting a string value exists, an error
indicating an unexpected argument is raised. If `val` is `nothing` and all
methods require a string value, an error indicating a missing required argument
is raised.
"""
handle_opts(body::Function, opt::Symbol, value::Any=nothing) =
    _handle_opts(body, ==(opt), value)
handle_opts(body::Function, re::Regex, value::Any=nothing) =
    _handle_opts(body, opt->contains(string(opt), re), value)

end # module

using .Options
