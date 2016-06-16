# definition of continuous mimo transfer function

CSiso = Union{CSisoTf}

immutable CMimo{T<:CSiso, M<:AbstractArray{T}} <: MimoSystem{T}
  m::M
  ny::Int
  nu::Int
  function call{M<:AbstractArray}(::Type{CMimo}, m::M)
    @assert 0 < ndims(m)
    @assert eltype(m) <: CSiso

    T = eltype(M)
    return new{T,M}(m, size(m, 1), size(m, 2))
  end
end

# creation of continuous mimo transfer functions

"""
    tf(num, den[, Ts])

Constructs a continuous transfer function with numerator `num`,
denominator `den` of type `Vector`.

A continuous transfer function can also be constructed from `Polynomials.Poly` objects.

A continuous MIMO transfer function can be constructed from a `Matrix` of `Vector` or
`Matrix` of `Polynomials.Poly` objects.

Discrete transfer functions are constructed using additional argument sampling
time `Ts`.

# Examples
```julia
julia> tf([1,0,3],[1,1,2])
  s^2 + 3
-----------
s^2 + s + 2
```

```julia
julia> tf(Polynomials.Poly([3,0,1]),Polynomials.Poly([2,1,1]))
  s^2 + 3
-----------
s^2 + s + 2
```

```julia
julia> tf([1,0,3],[1,1,2],1)
   z^-2 + 3
---------------
z^-2 + z^-1 + 2
Sampling time = 1
```
"""
function tf{M1<:AbstractArray, M2<:AbstractArray}(num::M1, den::M2)
  @assert eltype(eltype(num)) <: CSiso
  @assert eltype(eltype(den)) <: CSiso

  T = promote_type(eltype(num), eltype(den))
  if size(num, 1, 2) != size(den, 1, 2)
    warn("num and den dimensions must match")
    throw(DomainError())
  end
  m = similar(num, CSisoRational{T})
  for idx in eachindex(m)
    m[idx] = tf(num[idx], den[idx])
  end
  CMimo(m)
end

function tf{M<:AbstractArray}(gain::M)
  @assert eltype(eltype(gain)) <: Real

  m = similar(gain, CSisoRational{eltype(gain)})
  for idx in eachindex(m)
    m[idx] = tf(gain[idx])
  end
  CMimo(m)
end

function zpk{M1<:AbstractArray, M2<:AbstractArray, M3<:AbstractArray}(
  z::M1, p::M2, k::M3)
  @assert eltype(eltype(z)) <: Number
  @assert eltype(eltype(p)) <: Number
  @assert eltype(k)         <: Real

  T = promote_type(eltype(z), eltype(p))
  # Validate input and output dimensions match
  ny, nu = size(z, 1, 2)
  if (ny, nu) != size(p, 1, 2) || (ny, nu) != size(k, 1, 2)
    warn("num and den dimensions must match")
    throw(DomainError())
  end
  m = similar(z, CSisoZpk{T})
  for idx in eachindex(m)
    m[idx] = zpk(z_[idx], p_[idx], k[idx])
  end
  CMimo(m)
end

function zpk{M1<:AbstractArray}(gain::M1)
  @assert eltype(gain) <: Real

  m = similar(gain, CSisoZpk{eltype(gain)})
  for idx in eachindex(m)
    m[idx] = zpk(gain[idx])
  end
  CMimo(m)
end

# conversion and promotion

promote_rule{T<:Real,T1,T2}(::Type{CMimo{T1,T2}}, ::Type{T}) = CMimo
promote_rule{T<:Real}(::Type{CMimo}, ::Type{T}) = CMimo
convert{T<:Real}(::Type{CMimo}, x::T)           =  zpk(sparse([1],[1],[x]))

promote_rule{M1<:AbstractArray,T1,T2}(::Type{CMimo{T1,T2}}, ::Type{M1}) = CMimo
promote_rule{M1<:AbstractArray}(::Type{CMimo}, ::Type{M1}) = CMimo
convert{M1<:AbstractArray}(::Type{CMimo}, x::M1)           = zpk(x)

# overloading identities

one(s::CMimo)        = zpk(sparse([1],[1],[Int8(1)]))
one(::Type{CMimo})   = zpk(sparse([1],[1],[Int8(1)]))
zero(s::CMimo)       = zpk(sparse([1],[1],[Int8(0)]))
zero(::Type{CMimo})  = zpk(sparse([1],[1],[Int8(0)]))

# interface implementation

samplingtime(s::CMimo) = map(samplingtime, s.m)

getmatrix(s::CMimo) = s.m

isdiscrete(s::CMimo) = false

# overload printing functions

function show(io::IO, s::CMimo)
  println(io, "Continuous time transfer function model")
  println(io, "\ty = Gu")
  println(io, "with nu=", s.nu, ", ny=", s.ny, ".")
end

function showall{T1}(io::IO, s::CMimo{T1})
  show(io, s)
  println(io, "")
  for subs in s
    printtransferfunction(io::IO, subs)
    println(io, "")
  end
end

# overload mathematical operations

function +(s1::CMimo, s2::CMimo)
  if size(s1.m) != size(s2.m)
    warn("Systems have different shapes")
    throw(DomainError)
  end
  m = s1.m + s2.m
  CMimo(m)
end

+{T<:Real}(s::CMimo, n::T)          = CMimo(s.m+n)
+{T<:Real}(n::T, s::CMimo)          = +(s, n)

function +{T<:Real}(s::CMimo, n::Matrix{T})
  if size(s.m) != size(n)
    warn("Systems have different shapes")
    throw(DomainError)
  end
  CMimo(s.m + n)
end
+{T<:Real}(n::Matrix{T}, s::CMimo)  = (s, n)

.+{T<:Real}(n::T, s::CMimo)         = +(n, s)
.+{T<:Real}(s::CMimo, n::T)         = +(s, n)

.+{T<:Real}(n::Matrix{T}, s::CMimo) = +(n, s)
.+{T<:Real}(s::CMimo, n::Matrix{T}) = +(s, n)
.+(s1::CMimo, s2::CMimo)            = +(s1, s2)

-(s::CMimo)                         = CMimo(-s.m)

-{T<:Real}(s::CMimo, n::Matrix{T})  = +(s, -n)
-{T<:Real}(n::Matrix{T}, s::CMimo)  = +(-n, s)
-(s1::CMimo,s2::CMimo)              = +(s1,-s2)
-{T<:Real}(s::CMimo, n::T)          = +(s,-n)
-{T<:Real}(n::T, s::CMimo)          = +(-n, s)

.-{T<:Real}(n::T, s::CMimo)         = +(-n, s)
.-{T<:Real}(s::CMimo, n::T)         = +(s,-n)

.-{T<:Real}(n::Matrix{T}, s::CMimo) = -(n, s)
.-{T<:Real}(s::CMimo, n::Matrix{T}) = -(s, n)

function *(s1::CMimo, s2::CMimo)
  if size(s1.m, 2) != size(s2.m, 1)
    warn("s1*s2: s1 must have same number of inputs as s2 has outputs")
    throw(DomainError())
  end
  CMimo(s1.m*s2.m)
end

*{T<:Real}(s::CMimo, n::T)          = CMimo(s.m*n)
*{T<:Real}(n::T, s::CMimo)          = *(s, n)

*{T<:Real}(s::CMimo, n::Matrix{T})  = CMimo(s.m*n)
*{T<:Real}(n::Matrix{T}, s::CMimo)  = CMimo(n*s.m)

.*{T<:Real}(n::T, s::CMimo)         = *(n, s)
.*{T<:Real}(s::CMimo, n::T)         = *(s, n)
.*{T<:Real}(n::Matrix{T}, s::CMimo) = *(n, s)
.*{T<:Real}(s::CMimo, n::Matrix{T}) = *(s, n)
.*(s1::CMimo, s2::CMimo)            = *(s1,s2)

function /{T<:Real}(n::T, s::CMimo)
  warn("MIMO Transfer function inversion isn't implemented yet")
  throw(DomainError())
end
/(s1::CMimo, s2::CMimo)             = s1*(1/s2)
/{T<:Real}(s::CMimo, n::T)          = s*(1/n)
function /{T<:Real}(n::Matrix{T}, t::CMimo)
  warn("MIMO Transfer function inversion isn't implemented yet")
  throw(DomainError())
end
/{T<:Real}(s::CMimo, n::Matrix{T})  = CMimo(s.m/n)

./{T<:Real}(s::CMimo, n::T)         = s/n
./{T<:Real}(n::T, s::CMimo)         = n/s
./{T<:Real}(s::CMimo, n::Matrix{T}) = s/n
./{T<:Real}(n::Matrix{T}, s::CMimo) = n/s

==(s1::CMimo, s2::CMimo)            = s1.m != s2.m

!=(s1::CMimo, s2::CMimo)            = !(s1.m == s2.m)

function isapprox(
  s1::CMimo, s2::CMimo)
  # TODO: Implement
end
