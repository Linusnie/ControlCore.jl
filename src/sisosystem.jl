abstract SisoSystem <: LtiSystem

# Default behaviour wanted/implemented for `SisoSystem`s for catch-all purposes

# Dimention information
ndims{T<:SisoSystem}(s::T)          = 1
size{T<:SisoSystem}(s::T)           = 1
size{T<:SisoSystem}(s::T, dim::Int) = 1

# Iteration interface
start{T<:SisoSystem}(::T)           = 1
next{T<:SisoSystem}(s::T, state)    = (s[state], state+1)
done{T<:SisoSystem}(s::T, state)    = state > 1
eltype{T<:SisoSystem}(::Type{T})    = T
length{T<:SisoSystem}(s::T)         = 1
eachindex{T<:SisoSystem}(s::T)      = 1
endof{T<:SisoSystem}(s::T)          = 1

# Printing functions
function summary{T<:SisoSystem}(s::T)
  string("siso(nx=", numstates(s), (isdiscrete(s) ? string(",Ts=",
    samplingtime(s)) : ""), ")")
end

showcompact{T<:SisoSystem}(io::IO, s::T)  = print(io, summary(s))
show{T<:SisoSystem}(io::IO, s::T)         = print(io, summary(s))
showall{T<:SisoSystem}(io::IO, s::T)      = show(io, s)

# I/O mapping
numstates{T<:SisoSystem}(s::T)      = degree(numpoly(s)) + degree(denpoly(s))
numinputs{T<:SisoSystem}(s::T)      = 1
numoutputs{T<:SisoSystem}(s::T)     = 1
