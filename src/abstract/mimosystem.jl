abstract MimoSystem{T<:SisoSystem} <: LtiSystem

# I/O mapping
numinputs(s::MimoSystem)  = size(getmatrix(s)::AbstractArray, 2)
numoutputs(s::MimoSystem) = size(getmatrix(s)::AbstractArray, 1)

# Dimension information
ndims(s::MimoSystem)              = ndims(getmatrix(s)::AbstractArray)
size(s::MimoSystem)               = size(getmatrix(s)::AbstractArray)
size(s::MimoSystem, dim::Int)     = size(getmatrix(s)::AbstractArray, dim)
size(s::MimoSystem, dims::Int...) = size(getmatrix(s)::AbstractArray, dims...)

# Iteration interface
start(::MimoSystem)               = 1
next(s::MimoSystem, state::Int)   = (s[state], state+1)
done(s::MimoSystem, state::Int)   = done(getmatrix(s)::AbstractArray, state)
eltype{T<:SisoSystem}(::Type{MimoSystem{T}})  = T
length(s::MimoSystem)             = length(getmatrix(s)::AbstractArray)
eachindex(s::MimoSystem)          = eachindex(getmatrix(s)::AbstractArray)
endof(s::MimoSystem)              = endof(getmatrix(s)::AbstractArray)

# Printing functions
summary(s::MimoSystem) = string("mimo(nx=", numstates(s), ",ny=", numoutputs(s),
  ",nu=", numinputs(s), (isdiscrete(s) ? ",Ts=", samplingtime(s), : ""), ")")

showcompact(io::IO, s::MimoSystem)  = print(io, summary(s))
show(io::IO, s::MimoSystem)         = print(io, summary(s))
showall(io::IO, s::MimoSystem)      = print(io, summary(s))
