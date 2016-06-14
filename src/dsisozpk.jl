# definition of discrete rational transfer function

immutable DSisoZpk{T<:Real, T1<:Real, T2<:Real, T3<:Real} <: DSisoTf{T}
  z::Vector{Complex{T1}}
  p::Vector{Complex{T2}}
  k::T3
  Ts::Float64

  function call{T1<:Number, T2<:Number, T3<:Real}(::Type{DSisoZpk},
    z::Vector{T1}, p::Vector{T2}, k::T3, Ts::Float64)

    tmp = max(Ts, zero(Float64))
    Ts_ = tmp == zero(Float64) ? NaN : tmp
    T = promote_type(real(T1), real(T2), T3)
    new{T,real(T1),real(T2),T3}(copy(z), copy(p), k, Ts_)
  end
end

# creation of discrete zpk transfer functions

function zpk{T1<:Number, T2<:Number, T3<:Real, T4<:Real}(z::Vector{T1},
    p::Vector{T2}, k::T3, Ts::T4)
  DSisoZpk(z, p, k, Float64(Ts))
end

function zpk{T1<:Real, T2<:Real}(k::T1, Ts::T2)
  DSisoZpk(Vector{Int8}(), Vector{Int8}(), k, Float64(Ts))
end

# conversion and promotion

promote_rule{T1,T2,T3,T4,T5<:Real}(::Type{DSisoZpk{T1,T2,T3,T4}}, ::Type{T5}) =
  DSisoZpk
promote_rule{T1,T<:Real}(::Type{DSisoZpk{T1}}, ::Type{T}) = DSisoZpk
convert{T<:Real}(::Type{DSisoZpk}, x::T) = zpk(x, zero(Float64))

# overloading identities

zero{T}(::Type{DSisoZpk{T}}) = zpk(zero(T), zero(Float64))
zero{T}(s::DSisoZpk{T})      = zpk(zero(T), zero(Float64))
one{T}(::Type{DSisoZpk{T}})  = zpk(one(T), zero(Float64))
one{T}(s::DSisoZpk{T})       = zpk(one(T), zero(Float64))

# interface implementation

function zeros(s::DSisoZpk)
  copy(s.z)
end

function poles(s::DSisoZpk)
  copy(s.p)
end

function numvec(s::DSisoZpk)
  coeffs(poly(s.z))[end:-1:1]
end

function denvec(s::DSisoZpk)
  coeffs(poly(s.p))[end:-1:1]
end

function numpoly(s::DSisoZpk)
  poly(s.z)
end

function denpoly(s::DSisoZpk)
  poly(s.p)
end

function zpkdata(s::DSisoZpk)
  copy(s.k)
end

function samplingtime(s::DSisoZpk)
  copy(s.Ts)
end

# overload printing functions

function show(io::IO, s::DSisoZpk)
  println(io, "Discrete time zpk transfer function model")
  println(io, "\ty = Gu")
  if s.Ts > 0
    println(io, "with Ts=", s.Ts, ".")
  elseif s.Ts == 0
    println(io, "with Ts=unspecified.")
  end
end

function showall(io::IO, s::DSisoZpk)
  show(io, s)
  println(io, "")
  printtransferfunction(io::IO, s)
end

# overload mathematical operations

function +(s1::DSisoZpk, s2::DSisoZpk)
  if s1.Ts == s2.Ts || s2.Ts == NaN
    Z = s1.k*poly(s1.z)*poly(s2.p) + s2.k*poly(s2.z)*poly(s1.p)
    z::Array{Complex{Float64}} = roots(Z)
    p = vcat(s1.p, s2.p)
    k = real(Z[end]) # Poly is now reverse order
    return zpk(z, p, k, s1.Ts)
  elseif s1.Ts == NaN
    Z = s1.k*poly(s1.z)*poly(s2.p) + s2.k*poly(s2.z)*poly(s1.p)
    z = roots(Z)
    p = vcat(s1.p, s2.p)
    k = real(Z[end]) # Poly is now reverse order
    return zpk(z, p, k, s2.Ts)
  else
    warn("Sampling time mismatch")
    throw(DomainError())
  end
end

function +{T<:Real}(s::DSisoZpk, n::T)
  Z = s.k*poly(s.z) + n*poly(s.p)
  z::Array{Complex{Float64}} = roots(Z)
  p = s.p
  k = real(Z[end]) # Poly is now reverse order
  zpk(z, p, k, s.Ts)
end
+{T<:Real}(n::T, s::DSisoZpk)  = +(s, n)

.+{T<:Real}(s::DSisoZpk, n::T) = +(n, s)
.+{T<:Real}(n::T, s::DSisoZpk) = +(s, n)
.+(s1::DSisoZpk, s2::DSisoZpk) = +(s1, s2)

-(s::DSisoZpk)                 = zpk(s.z, s.p, -s.k, s.Ts)

-(s1::DSisoZpk, s2::DSisoZpk)  = +(s1, -s2)
-{T<:Real}(s::DSisoZpk, n::T)  = +(s, -n)
-{T<:Real}(n::T, s::DSisoZpk)  = +(n, -s)

.-{T<:Real}(s::DSisoZpk, n::T) = +(s, -n)
.-{T<:Real}(n::T, s::DSisoZpk) = +(n, -s)
.-(s1::DSisoZpk, s2::DSisoZpk) = +(s1, -s2)

function *(s1::DSisoZpk, s2::DSisoZpk)
  if s1.Ts == s2.Ts || s2.Ts == NaN
    z = vcat(s1.z, s2.z)
    p = vcat(s1.p, s2.p)
    k = s1.k*s2.k
    return zpk(z, p, k, s1.Ts)
  elseif s1.Ts == NaN
    z = vcat(s1.z, s2.z)
    p = vcat(s1.p, s2.p)
    k = s1.k*s2.k
    return zpk(z, p, k, s2.Ts)
  else
    warn("Sampling time mismatch")
    throw(DomainError())
  end
end
*{T<:Real}(s::DSisoZpk, n::T)  = zpk(s.z, s.p, n*s.k, s.Ts)
*{T<:Real}(n::T, s::DSisoZpk)  = *(s, n)

.*{T<:Real}(s::DSisoZpk, n::T) = *(n, s)
.*{T<:Real}(n::T, s::DSisoZpk) = *(s, n)
.*(s1::DSisoZpk, s2::DSisoZpk) = *(s1, s2)

/{T<:Real}(n::T, s::DSisoZpk)  = zpk(s.p, s.z, n./s.k, s.Ts)
/(s1::DSisoZpk, s2::DSisoZpk)  = s1*(1/s2)
/{T<:Real}(s::DSisoZpk, n::T)  = s*(1/n)

./{T<:Real}(n::T, s::DSisoZpk) = n*(1/s)
./{T<:Real}(s::DSisoZpk, n::T) = s*(1/n)
./(s1::DSisoZpk, s2::DSisoZpk) = s1*(1/s2)

function ==(s1::DSisoZpk, s2::DSisoZpk)
  fields = [:Ts, :z, :p, :k]
  for field in fields
    if getfield(s1, field) != getfield(s2, field)
      return false
    end
  end
  true
end

!=(s1::DSisoZpk, s2::DSisoZpk) = !(s1==s2)

function isapprox(s1::DSisoZpk, s2::DSisoZpk)
  # TODO: Implement
end
