module PicoSAT

using PicoSAT_jll

const UNKNOWN = 0
const SATISFIABLE = 10
const UNSATISFIABLE = 20

init() =
    ccall((:picosat_init, libpicosat), Ptr{Cvoid}, ())
reset(p::Ptr{Cvoid}) =
    ccall((:picosat_reset, libpicosat), Cvoid, (Ptr{Cvoid},), p)
adjust(p::Ptr{Cvoid}, N::Integer) =
    ccall((:picosat_adjust, libpicosat), Cvoid, (Ptr{Cvoid}, Cint), p, N)
add(p::Ptr{Cvoid}, lit::Integer) =
    ccall((:picosat_add, libpicosat), Cint, (Ptr{Cvoid}, Cint), p, lit)
sat(p::Ptr{Cvoid}, limit::Integer = -1) =
    ccall((:picosat_sat, libpicosat), Cint, (Ptr{Cvoid}, Cint), p, limit)
deref(p::Ptr{Cvoid}, lit::Integer) =
    ccall((:picosat_deref, libpicosat), Cint, (Ptr{Cvoid}, Cint), p, lit)
push(p::Ptr{Cvoid}) =
    ccall((:picosat_push, libpicosat), Cint, (Ptr{Cvoid},), p)
pop(p::Ptr{Cvoid}) =
    ccall((:picosat_pop, libpicosat), Cint, (Ptr{Cvoid},), p)

# for debugging by exporting CNF files
function print(p::Ptr{Cvoid}, path::AbstractString)
    f = ccall(:fopen, Ptr{Cvoid}, (Cstring, Cstring), path, "w")
    ccall((:picosat_print, libpicosat), Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}), p, f)
    @assert ccall(:fclose, Cint, (Ptr{Cvoid},), f) == 0
end

end