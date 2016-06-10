abstract CSisoTf{T} <: SisoTf{T}

# Sampling information
isdiscrete(s::CSisoTf)          = false
samplingtime{T}(s::CSisoTf{T})  = zero(T)

# Printing functions
summary(s::CSisoTf)             = string("tf(nx=", numstates(s), ")")

showcompact(io::IO, s::CSisoTf) = print(io, summary(s))
show(io::IO, s::CSisoTf)        = print(io, summary(s))
showall(io::IO, s::CSisoTf)     = print(io, summary(s))
