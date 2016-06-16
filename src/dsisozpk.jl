# definition of discrete rational transfer function

immutable DSisoZpk{T<:Real, T1<:Real, T2<:Real, T3<:Real} <: DSisoTf{T}
  z::Vector{Complex{T1}}
  p::Vector{Complex{T2}}
  k::T3
  Ts::Float64

  function call{T1<:Number, T2<:Number, T3<:Real}(::Type{DSisoZpk},
    z::Vector{T1}, p::Vector{T2}, k::T3, Ts::Float64)

    Ts_ = Ts > zero(Float64) ? Ts : NaN
    T = promote_type(real(T1), real(T2), T3)
    new{T,real(T1),real(T2),T3}(z, p, k, Ts_)
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

zero{T1,T2,T3,T4}(::Type{DSisoZpk{T1,T2,T3,T4}}) = zpk(zero(Int8), zero(Float64))
zero{T}(::Type{DSisoZpk{T}})                     = zpk(zero(Int8), zero(Float64))
zero(::Type{DSisoZpk})                           = zpk(zero(Int8), zero(Float64))
zero{T}(s::DSisoZpk{T})                          = zpk(zero(Int8), zero(Float64))
one{T1,T2,T3,T4}(::Type{DSisoZpk{T1,T2,T3,T4}})  = zpk(one(Int8), zero(Float64))
one{T}(::Type{DSisoZpk{T}})                      = zpk(one(Int8), zero(Float64))
one(::Type{DSisoZpk})                            = zpk(one(Int8), zero(Float64))
one{T}(s::DSisoZpk{T})                           = zpk(one(Int8), zero(Float64))

# interface implementation

function zeros(s::DSisoZpk)
  copy(s.z)
end

function poles(s::DSisoZpk)
  copy(s.p)
end

function numvec(s::DSisoZpk)
  real(coeffs(poly(s.z))[end:-1:1])
end

function denvec(s::DSisoZpk)
  real(coeffs(poly(s.p))[end:-1:1])
end

function numpoly(s::DSisoZpk)
  poly(s.z)
end

function denpoly(s::DSisoZpk)
  poly(s.p)
end

function zpkdata(s::DSisoZpk)
  (s.z, s.p, s.k)
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

function +{T1,T2,T3,T4,T5,T6,T7,T8}(s1::DSisoZpk{T1,T2,T3,T4}, s2::DSisoZpk{T5,T6,T7,T8})
  Ts::Float64
  if s1.Ts == s2.Ts || isnan(s2.Ts)
    Ts = s1.Ts
  elseif isnan(s1.Ts)
    Ts = s2.Ts
  else
    warn("Sampling time mismatch")
    throw(DomainError())
  end
  #Tz = float(promote_type(T1, T5))
  #Tp = promote_type(T3, T7)
  #Tk = promote_type(T1, T5)
  p1,p2,pcommon = rmcommon(copy(s1.p), copy(s2.p))
  z1,z2,zcommon = rmcommon(copy(s1.z), copy(s2.z))
  Z = s1.k*poly(z1)*poly(p2) + s2.k*poly(z2)*poly(p1)
  z = vcat(convert(Vector{Complex{Float64}},roots(Z)), zcommon)
  p = vcat(p1, p2, pcommon)
  k = real(Z[end]) # Poly is now reverse order
  zpk(z, p, k, Ts) #::DSisoZpk{Tz,Tz,Tp,Tk}
end

function +{T1,T2,T3,T4,T<:Real}(s::DSisoZpk{T1,T2,T3,T4}, n::T)
  Tk = promote_type(T1, T)
  Tz = float(promote_type(T1, T))
  Z = s.k*poly(s.z) + n*poly(s.p)
  z = roots(Z)
  p = copy(s.p)
  k = real(Z[end])
  zpk(z, p, k, s.Ts)::DSisoZpk{Tz,Tz,T3,Tk}
end
+{T<:Real}(n::T, s::DSisoZpk)  = +(s, n)

.+{T<:Real}(s::DSisoZpk, n::T) = +(s, n)
.+{T<:Real}(n::T, s::DSisoZpk) = +(n, s)
.+(s1::DSisoZpk, s2::DSisoZpk) = +(s1, s2)

-(s::DSisoZpk)                 = zpk(copy(s.z), copy(s.p), -s.k, s.Ts)

-(s1::DSisoZpk, s2::DSisoZpk)  = +(s1, -s2)
-{T<:Real}(s::DSisoZpk, n::T)  = +(s, -n)
-{T<:Real}(n::T, s::DSisoZpk)  = +(n, -s)

.-{T<:Real}(s::DSisoZpk, n::T) = +(s, -n)
.-{T<:Real}(n::T, s::DSisoZpk) = +(n, -s)
.-(s1::DSisoZpk, s2::DSisoZpk) = +(s1, -s2)

function *(s1::DSisoZpk, s2::DSisoZpk)
  Ts::Float64
  if s1.Ts == s2.Ts || isnan(s2.Ts)
    Ts = s1.Ts
  elseif isnan(s1.Ts)
    s2.Ts
  else
    warn("Sampling time mismatch")
    throw(DomainError())
  end
  z1,p2,pcommon = rmcommon(s1.z, s2.p)
  p1,z2,zcommon = rmcommon(s1.p, s2.z)
  z = vcat(z1, z2)
  p = vcat(p1, p2)
  k = s1.k*s2.k
  zpk(z, p, k, Ts)
end
*{T<:Real}(s::DSisoZpk, n::T)  = zpk(copy(s.z), copy(s.p), n*s.k, s.Ts)
*{T<:Real}(n::T, s::DSisoZpk)  = *(s, n)

.*{T<:Real}(s::DSisoZpk, n::T) = *(n, s)
.*{T<:Real}(n::T, s::DSisoZpk) = *(s, n)
.*(s1::DSisoZpk, s2::DSisoZpk) = *(s1, s2)

/{T<:Real}(n::T, s::DSisoZpk)  = zpk(copy(s.p), copy(s.z), n./s.k, s.Ts)
/{T<:Real}(s::DSisoZpk, n::T)  = s*(1/n)
/(s1::DSisoZpk, s2::DSisoZpk)  = s1*(1/s2)

./{T<:Real}(n::T, s::DSisoZpk) = n*(1/s)
./{T<:Real}(s::DSisoZpk, n::T) = s*(1/n)
./(s1::DSisoZpk, s2::DSisoZpk) = s1*(1/s2)

function ==(s1::DSisoZpk, s2::DSisoZpk)
  s1.Ts == s2.Ts && s1.z == s2.z &&
    s1.p == s2.p && s1.k == s2.k
end

!=(s1::DSisoZpk, s2::DSisoZpk) = !(s1==s2)

function isapprox(s1::DSisoZpk, s2::DSisoZpk,
    rtol::Real=sqrt(eps()), atol::Real=0)
  sdiff = s2-s1
  norm(sdiff.k) < rtol
end
