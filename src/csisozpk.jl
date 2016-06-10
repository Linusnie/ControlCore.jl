immutable CSisoZpk{T<:AbstractFloat} <: CSisoTf{T}
  z::Vector{Complex{T}}
  p::Vector{Complex{T}}
  k::T
  function call{T}(::Type{CSisoZpk}, z::Vector{Complex{T}},
    p::Vector{Complex{T}}, k::T)
    new{T}(copy(z), copy(p), k)
  end
end

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
function zpk{T1<:AbstractFloat}(z::Vector{Complex{T1}},  p::Vector{Complex{T1}}, k::T1)
  CSisoZpk(z, p, k)
end

function zpk{T1<:Number, T2<:Number, T3<:Real}(z::Vector{T1}, p::Vector{T2}, k::T3)
  T   = promote_type(real(T1), real(T2), T3, Float16) # ensure AbstractFloat
  z_  = convert(Vector{Complex{T}},z)
  p_  = convert(Vector{Complex{T}},p)
  k_  = convert(T,k)
  CSisoZpk(z_, p_, k_)
end

function zpk{T1<:Real}(k::T1)
  zpk(Vector{T1}(), Vector{T1}(), k)
end

Base.promote_rule{T1<:AbstractFloat, T2<:AbstractFloat}(::CSisoZpk{T1}, ::CSisoZpk{T2}) = CSisoZpk{promote_type(T1, T2)}
Base.convert{T1<:AbstractFloat, T2<:AbstractFloat}(::Type{CSisoZpk{T1}}, s::CSisoZpk{T2}) = zpk(convert(Poly{T1}, s.z), convert(Poly{T1}, s.p),  convert(T1, s.k))

Base.promote_rule{T1<:AbstractFloat, T2<:Real}(::CSisoZpk{T1}, ::T2) = CSisoZpk{promote_type(T1, T2)}
Base.convert{T1<:AbstractFloat, T2<:Real}(::Type{CSisoZpk{T1}}, x::T2) = zpk([one(T1)], [one(T1)], convert(T1, x))

Base.one{T1<:AbstractFloat}(::Type{CSisoZpk{T1}}) = zpk(Vector{T1}(), Vector{T1}(), one(T1))
Base.one{T1<:AbstractFloat}(s::Type{CSisoZpk{T1}}) = zpk(Vector{T1}(), Vector{T1}(), one(T1))
Base.zero{T1<:AbstractFloat}(::Type{CSisoZpk{T1}}) = zpk(Vector{T1}(), Vector{T1}(), zero(T1))
Base.zero{T1<:AbstractFloat}(s::Type{CSisoZpk{T1}}) = zpk(Vector{T1}(), Vector{T1}(), zero(T1))

function zeros{T1<:AbstractFloat}(s::CSisoZpk{T1})
  return copy(s.z)
end

function poles{T1<:AbstractFloat}(s::CSisoZpk{T1})
  return copy(s.p)
end

function numvec{T1<:AbstractFloat}(s::CSisoZpk{T1})
  return coeffs(poly(s.z))[end:-1:1]
end

function denvec{T1<:AbstractFloat}(s::CSisoZpk{T1})
  return coeffs(poly(s.p))[end:-1:1]
end

function numpoly{T1<:AbstractFloat}(s::CSisoZpk{T1})
  return poly(s.z)
end

function denpoly{T1<:AbstractFloat}(s::CSisoZpk{T1})
  return poly(s.p)
end

function zpkdata{T1<:AbstractFloat}(s::CSisoZpk{T1})
  return copy(s.k)
end

function samplingtime{T1<:AbstractFloat}(s::CSisoZpk{T1})
  return -one(Float64)
end

function show(io::IO, s::CSisoZpk)
  println(io, "Discrete time zpk transfer function model")
  println(io, "\ty = Gu")
  if s.Ts > 0
    println(io, "with Ts=", s.Ts, ".")
  elseif s.Ts == 0
    println(io, "with Ts=unspecified.")
  end
end

function showall(io::IO, s::CSisoZpk)
  show(io, s)
  println(io, "")
  printtransferfunction(io::IO, s)
end

function +{T1<:AbstractFloat, T2<:AbstractFloat}(
  s1::CSisoZpk{T1}, s2::CSisoZpk{T2})
  T = promote_type(T1, T2)
  z::Array{T}
  p::Array{T}
  k::T

  Z = s1.k*poly(s1.z)*poly(s2.p) + s2.k*poly(s2.z)*poly(s1.p)
  z = roots(Z)
  p = vcat(s1.p, s2.p)
  k = Z[end] # Poly is now reverse order
  return zpk(z, p, k)
end

function +{T1<:AbstractFloat, T2<:Real}(
    s::CSisoZpk{T1}, n::T2)
  Z = s.k*poly(s.z) + n*poly(s.p)
  z = roots(Z)
  p = s.p
  k = real(Z[end]) # Poly is now reverse order
  return zpk(z_, p_, k)
end
#+{T1<:AbstractFloat, T2<:Real}(s::CSisoZpk{T1}, n::T2)  = +(s, zpk(n))
+{T1<:AbstractFloat, T2<:Real}(n::T2, s::CSisoZpk{T1})  = +(s, n)

.+{T1<:AbstractFloat, T2<:Real}(s::CSisoZpk{T1}, n::T2) = +(s, zpk(n))
.+{T1<:AbstractFloat, T2<:Real}(n::T2, s::CSisoZpk{T1}) = +(s, n)
.+{T1<:AbstractFloat, T2<:AbstractFloat}(s1::CSisoZpk{T1}, s2::CSisoZpk{T2}) = +(s1, s2)

-{T1<:AbstractFloat}(s::CSisoZpk{T1}) = zpk(s.z, s.p, -s.k)
function -{T1<:AbstractFloat, T2<:AbstractFloat}(
  s1::CSisoZpk{T1}, s2::CSisoZpk{T2})
  return +(s1, -s2)
end
-{T1<:AbstractFloat, T2<:Real}(s::CSisoZpk{T1}, n::T2)  = +(s, -n)
-{T1<:AbstractFloat, T2<:Real}(n::T2, s::CSisoZpk{T1})  = +(n, -t)

.-{T1<:AbstractFloat, T2<:Real}(s::CSisoZpk{T1}, n::T2) = +(s, -n)
.-{T1<:AbstractFloat, T2<:Real}(n::T2, s::CSisoZpk{T1}) = +(n, -t)
.-{T1<:AbstractFloat, T2<:AbstractFloat}(s1::CSisoZpk{T1}, s2::CSisoZpk{T2}) = +(s1, -s2)

function *{T1<:AbstractFloat, T2<:AbstractFloat}(
  s1::CSisoZpk{T1}, s2::CSisoZpk{T2})
  z = vcat(s1.z, s2.z)
  p = vcat(s1.p, s2.p)
  k = s1.k*s2.k
  return zpk(z, p, k)
end
*{T1<:AbstractFloat, T2<:Real}(s::CSisoZpk{T1}, n::T2)  = zpk(s.z, s.p, n*s.k)
*{T1<:AbstractFloat, T2<:Real}(n::T2, s::CSisoZpk{T1})  = *(s, n)

.*{T1<:AbstractFloat, T2<:Real}(s::CSisoZpk{T1}, n::T2) = *(n, t)
.*{T1<:AbstractFloat, T2<:Real}(n::T2, s::CSisoZpk{T1}) = *(s, n)
.*{T1<:AbstractFloat, T2<:AbstractFloat}(s1::CSisoZpk{T1}, s2::CSisoZpk{T2}) = *(s1,s2)

function /{T1<:AbstractFloat, T2<:Real}(
  n::T2, s::CSisoZpk{T1})
  return zpk(s.p, s.z, n./s.k, s.Ts)
end
/{T1<:AbstractFloat, T2<:Real}(s::CSisoZpk{T1}, n::T2)  = t*(1/n)
/{T1<:AbstractFloat, T2<:AbstractFloat}(s1::CSisoZpk{T1}, s2::CSisoZpk{T2}) = s1*(1/s2)
./{T1<:AbstractFloat, T2<:Real}(n::T2, s::CSisoZpk{T1}) = t*(1/n)
./{T1<:AbstractFloat, T2<:Real}(s::CSisoZpk{T1}, n::T2) = t*(1/n)
./{T1<:AbstractFloat, T2<:AbstractFloat}(s1::CSisoZpk{T1}, s2::CSisoZpk{T2}) = s1*(1/s2)

function =={T1<:AbstractFloat, T2<:AbstractFloat}(
  s1::CSisoZpk{T1}, s2::CSisoZpk{T2})
  fields = [:z, :p, :k]
  for field in fields
    if getfield(s1, field) != getfield(s2, field)
      return false
    end
  end
  true
end

!={T1<:AbstractFloat, T2<:AbstractFloat}(
  s1::CSisoZpk{T1}, s2::CSisoZpk{T2}) = !(s1==s2)

function isapprox{T1<:AbstractFloat, T2<:AbstractFloat}(
    s1::CSisoZpk{T1}, s2::CSisoZpk{T2})
  # TODO: Implement
end
