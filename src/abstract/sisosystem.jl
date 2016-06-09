abstract SisoSystem{T<:Real} <: LtiSystem

# Default behaviour wanted/implemented for `SisoSystem`s for catch-all purposes

# Dimention information
ndims(s::SisoSystem)          = 1
size(s::SisoSystem)           = 1
size(s::SisoSystem, dim::Int) = 1

# Iteration interface
start(::SisoSystem)               = 1
next(s::SisoSystem, state::Int)   = (s[state], state+1)
done(s::SisoSystem, state::Int)   = state > 1
eltype{T<:SisoSystem}(::Type{T})  = T
length(s::SisoSystem)             = 1
eachindex(s::SisoSystem)          = 1
endof(s::SisoSystem)              = 1

# Printing functions
function summary(s::SisoSystem)
  string("siso(nx=", numstates(s), (isdiscrete(s) ? string(",Ts=",
    samplingtime(s)) : ""), ")")
end

showcompact(io::IO, s::SisoSystem)  = print(io, summary(s))
show(io::IO, s::SisoSystem)         = print(io, summary(s))
showall(io::IO, s::SisoSystem)      = show(io, s)

# I/O mapping
numstates(s::SisoSystem)  = degree(numpoly(s)) + degree(denpoly(s))
numinputs(s::SisoSystem)  = 1
numoutputs(s::SisoSystem) = 1
