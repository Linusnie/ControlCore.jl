abstract MimoSystem{T<:SisoSystem} <: LtiSystem

# common constructor

function mimo{M<:AbstractArray}(m::M)
  @assert eltype(m) <: SisoSystem
  if eltype(m) <: DSiso
    DMimo(m)
  elseif eltype(m) <: CSiso
    CMimo(m)
  else
    warn("Do not support mixed continuous and discrete MIMO")
    throw(DomainError)
  end
end

# I/O mapping
numstates(s::MimoSystem)            = map(numstates, getmatrix(s)::AbstractArray)
numinputs(s::MimoSystem)            = size(getmatrix(s)::AbstractArray, 2)
numoutputs(s::MimoSystem)           = size(getmatrix(s)::AbstractArray, 1)

# Dimension information
ndims(s::MimoSystem)                = ndims(getmatrix(s)::AbstractArray)
size(s::MimoSystem)                 = size(getmatrix(s)::AbstractArray)
size(s::MimoSystem, dim::Int)       = size(getmatrix(s)::AbstractArray, dim)
size(s::MimoSystem, dims::Int...)   = map(x->size(s, x), dims)

# Iteration interface
start(::MimoSystem)                 = 1
next(s::MimoSystem, state::Int)     = (s[state], state+1)
done(s::MimoSystem, state::Int)     = state > length(s)
eltype{T<:SisoSystem}(::Type{MimoSystem{T}})  = T
length(s::MimoSystem)               = length(getmatrix(s)::AbstractArray)
eachindex(s::MimoSystem)            = eachindex(getmatrix(s)::AbstractArray)
endof(s::MimoSystem)                = endof(getmatrix(s)::AbstractArray)

# indexing
Base.getindex(s::MimoSystem, i) = s.M[i]
Base.getindex(s::MimoSystem, i, j) = s.M[i,j]
Base.getindex(s::MimoSystem, inds...) = mimo(getindex(s.M, inds...))

Base.setindex!(s1::MimoSystem, s2::SisoSystem, i) = setindex!(s1.M, s2, i)
Base.setindex!(s1::MimoSystem, s2::SisoSystem, i, j) = setindex!(s1.M, s2, i, j)
Base.setindex!(s1::MimoSystem, s2::MimoSystem, inds...) = setindex!(s1.M, s2.M, inds...)

# Concatenation
f(x) = isa(x, MimoSystem) ? getfield(x, :m) : x
Base.vcat(H::LtiSystem...) = mimo(cat(1, map(f, H)...))
Base.hcat(H::LtiSystem...) = mimo(cat(2, map(f, H)...))



# Common type interface
zeros(s::MimoSystem)   = map(poles, getmatrix(s)::AbstractArray)
poles(s::MimoSystem)   = map(poles, getmatrix(s)::AbstractArray)
numvec(s::MimoSystem)  = map(numvec, getmatrix(s)::AbstractArray)
denvec(s::MimoSystem)  = map(denvec, getmatrix(s)::AbstractArray)
numpoly(s::MimoSystem) = map(numpoly, getmatrix(s)::AbstractArray)
denpoly(s::MimoSystem) = map(denpoly, getmatrix(s)::AbstractArray)
zpkdata(s::MimoSystem) = map(zpkdata, getmatrix(s)::AbstractArray)

# Printing functions
summary(s::MimoSystem) = string("mimo(nx=", numstates(s), ",ny=", numoutputs(s),
  ",nu=", numinputs(s), (isdiscrete(s) ? string(",Ts=", samplingtime(s)) : ""),
  ")")

showcompact(io::IO, s::MimoSystem)  = print(io, summary(s))
show(io::IO, s::MimoSystem)         = print(io, summary(s))
showall(io::IO, s::MimoSystem)      = print(io, summary(s))
