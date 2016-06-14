abstract SisoTf{T<:Real} <: SisoSystem{T}

# Printing functions
summary(s::SisoTf) = string("tf(nx=", numstates(s), (isdiscrete(s) ?
  string(",Ts=", samplingtime(s)) : ""), ")")

showcompact(io::IO, s::SisoTf)  = print(io, summary(s))
show(io::IO, s::SisoTf)         = print(io, summary(s))
showall(io::IO, s::SisoTf)      = print(io, summary(s))
