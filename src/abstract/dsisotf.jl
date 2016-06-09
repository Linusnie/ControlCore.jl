abstract DSisoTf{T} <: SisoTf{T}

# Sampling information
isdiscrete(s::DSisoTf)          = true
samplingtime{T}(s::DSisoTf{T})  = s.Ts

# Printing functions
summary(s::DSisoTf) = string("tf(nx=", numstates(s), ",Ts=", samplingtime(s),
  ")")

showcompact(io::IO, s::DSisoTf) = print(io, summary(s))
show(io::IO, s::DSisoTf)        = print(io, summary(s))
showall(io::IO, s::DSisoTf)     = print(io, summary(s))
