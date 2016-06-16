# definition of continuous rational transfer function

immutable CSisoZpk{T<:Real, T1<:Real, T2<:Real, T3<:Real} <: CSisoTf{T}
  z::Vector{Complex{T1}}
  p::Vector{Complex{T2}}
  k::T3

  function call{T1<:Number, T2<:Number, T3<:Real}(::Type{CSisoZpk},
      z::Vector{T1}, p::Vector{T2}, k::T3)

    T = promote_type(real(T1), real(T2), T3)
    new{T,real(T1),real(T2),T3}(copy(z), copy(p), k)
  end
end

# creation of continuous zpk transfer functions

"""
    zpk([z, p,] k[, Ts])

Constructs a continuous transfer function with zeros `z`,
poles `p` of type `Vector` and gain `k.

A gain transfer function without zeros and poles is constructed by only
  supplying `k`.

A continuous MIMO transfer function can be constructed from a `Matrix` of `Vector`s.

Discrete transfer functions are constructed using additional argument sampling
time `Ts`.

# Examples
```julia
julia> zpk([1,0,3],[1,1,2],2)
ZPK:
    (s - 1.0)s(s - 3.0)
2.0 --------------------
    (s - 1.0)^2(s - 2.0)

Continuous-time zero-pole-gain model
```

```julia
julia> zpk(2)
ZPK:
    1.0
2.0 ---
    1.0
```

```julia
julia> zpk([1,0,3],[1,1,2],2,1)
ZPK:
    (z - 1.0)z(z - 3.0)
2.0 --------------------
    (z - 1.0)^2(z - 2.0)


Sample Time: 1 (seconds)
Discrete-time zero-pole-gain model
```
"""
function zpk{T1<:Number, T2<:Number, T3<:Real}(
    z::Vector{T1}, p::Vector{T2}, k::T3)
  CSisoZpk(z, p, k)
end

function zpk{T1<:Real}(k::T1)
  CSisoZpk(Vector{Int8}(), Vector{Int8}(), k)
end

# conversion and promotion

promote_rule{T1,T2,T3,T<:Real}(::Type{CSisoZpk{T1,T2,T3}}, ::Type{T}) =
  CSisoZpk
promote_rule{T<:Real}(::Type{CSisoZpk}, ::Type{T}) = CSisoZpk
convert{T<:Real}(::Type{CSisoZpk}, x::T) = zpk(x)

# overloading identities

zero{T}(::Type{CSisoZpk{T}}) = zpk(zero(T))
zero{T}(s::CSisoZpk{T})      = zpk(zero(T))
one{T}(::Type{CSisoZpk{T}})  = zpk(one(T))
one{T}(s::CSisoZpk{T})       = zpk(one(T))

# interface implementation

function zeros(s::CSisoZpk)
  copy(s.z)
end

function poles(s::CSisoZpk)
  copy(s.p)
end

function numvec(s::CSisoZpk)
  coeffs(poly(s.z))[end:-1:1]
end

function denvec(s::CSisoZpk)
  coeffs(poly(s.p))[end:-1:1]
end

function numpoly(s::CSisoZpk)
  poly(s.z)
end

function denpoly(s::CSisoZpk)
  poly(s.p)
end

function zpkdata(s::CSisoZpk)
  copy(s.k)
end

function samplingtime(s::CSisoZpk)
  -one(Float64)
end

# overload printing functions

function show(io::IO, s::CSisoZpk)
  println(io, "Continuous time zpk transfer function model")
  println(io, "\ty = Gu")
end

function showall(io::IO, s::CSisoZpk)
  show(io, s)
  println(io, "")
  printtransferfunction(io::IO, s)
end

# overload mathematical operations

function +(s1::CSisoZpk, s2::CSisoZpk)
  Z = s1.k*poly(s1.z)*poly(s2.p) + s2.k*poly(s2.z)*poly(s1.p)
  z = roots(Z)
  p = vcat(s1.p, s2.p)
  k = real(Z[end]) # Poly is now reverse order
  zpk(z, p, k)
end

function +{T<:Real}(s::CSisoZpk, n::T)
  Z = s.k*poly(s.z) + n*poly(s.p)
  z::Array{Complex{Float64}} = roots(Z)
  p = s.p
  k = real(Z[end]) # Poly is now reverse order
  CSisoZpk(z, p, k)
end
+{T<:Real}(n::T, s::CSisoZpk)  = +(s, n)

.+{T<:Real}(s::CSisoZpk, n::T) = +(s, zpk(n))
.+{T<:Real}(n::T, s::CSisoZpk) = +(s, n)
.+(s1::CSisoZpk, s2::CSisoZpk) = +(s1, s2)

-(s::CSisoZpk)                 = zpk(s.z, s.p, -s.k)
-(s1::CSisoZpk, s2::CSisoZpk)  = +(s1, -s2)
-{T<:Real}(s::CSisoZpk, n::T)  = +(s, -n)
-{T<:Real}(n::T, s::CSisoZpk)  = +(n, -s)

.-{T<:Real}(s::CSisoZpk, n::T) = +(s, -n)
.-{T<:Real}(n::T, s::CSisoZpk) = +(n, -s)
.-(s1::CSisoZpk, s2::CSisoZpk) = +(s1, -s2)

function *(s1::CSisoZpk, s2::CSisoZpk)
  z = vcat(s1.z, s2.z)
  p = vcat(s1.p, s2.p)
  k = s1.k*s2.k
  zpk(z, p, k)
end
*{T<:Real}(s::CSisoZpk, n::T)  = zpk(s.z, s.p, n*s.k)
*{T<:Real}(n::T, s::CSisoZpk)  = *(s, n)

.*{T<:Real}(s::CSisoZpk, n::T) = *(n, s)
.*{T<:Real}(n::T, s::CSisoZpk) = *(s, n)
.*(s1::CSisoZpk, s2::CSisoZpk) = *(s1,s2)

/{T<:Real}(n::T, s::CSisoZpk)  = zpk(s.p, s.z, n./s.k, s.Ts)
/{T<:Real}(s::CSisoZpk, n::T)  = s*(1/n)
/(s1::CSisoZpk, s2::CSisoZpk)  = s1*(1/s2)

./{T<:Real}(n::T, s::CSisoZpk) = s*(1/n)
./{T<:Real}(s::CSisoZpk, n::T) = s*(1/n)
./(s1::CSisoZpk, s2::CSisoZpk) = s1*(1/s2)

function ==(s1::CSisoZpk, s2::CSisoZpk)
  fields = [:z, :p, :k]
  for field in fields
    if getfield(s1, field) != getfield(s2, field)
      return false
    end
  end
  true
end

!=(s1::CSisoZpk, s2::CSisoZpk) = !(s1==s2)

function isapprox(s1::CSisoZpk, s2::CSisoZpk)
  # TODO: Implement
end
