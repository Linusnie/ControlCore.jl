# definition of discrete rational transfer function

immutable DSisoRational{T<:Real, T1<:Real, T2<:Real} <: DSisoTf{T}
  num::Poly{T1}
  den::Poly{T2}
  Ts::Float64

  function call{T1<:Real, T2<:Real}(::Type{DSisoRational}, num::Poly{T1},
      den::Poly{T2}, Ts::Float64)
    T = promote_type(eltype(num), eltype(den))
    Ts_ = Ts > zero(Float64) ? Ts : NaN
    new{T,eltype(num),eltype(den)}(num, den, Ts_)
  end
end

# creation of discrete rational transfer functions

function tf{V1<:AbstractVector, V2<:AbstractVector, T3<:Real}(num::V1, den::V2,
    Ts::T3)
  @assert eltype(num) <: Number string("num must be vector of T<:Number elements")
  @assert eltype(den) <: Number string("den must be vector of T<:Number elements")

  pnum = Poly(num[end:-1:1])
  pden = Poly(den[end:-1:1])
  DSisoRational(pnum, pden, Float64(Ts))
end

function tf{T1<:Real, T2<:Real, T3<:Real}(num::Poly{T1}, den::Poly{T2}, Ts::T3)
  DSisoRational(num, den, Float64(Ts))
end

function tf{T1<:Real, T2<:Real}(gain::T1, Ts::T2)
  DSisoRational(Poly([gain]), Poly([one(T1)]), Float64(Ts))
end

# conversion and promotion

promote_rule{T<:Real,T1,T2,T3}(::Type{DSisoRational{T1,T2,T3}}, ::Type{T}) =
  DSisoRational
promote_rule{T1,T<:Real}(::Type{DSisoRational{T1}}, ::Type{T}) = DSisoRational
convert{T<:Real}(::Type{DSisoRational}, x::T)   =
  tf([x], [one(T)], zero(Float64))

# overloading identities

zero{T1,T2,T3}(::Type{DSisoRational{T1,T2,T3}}) = tf([zero(Int8)], [one(Int8)], zero(Float64))
zero{T}(::Type{DSisoRational{T}})               = tf([zero(Int8)], [one(Int8)], zero(Float64))
zero(::Type{DSisoRational})                     = tf([zero(Int8)], [one(Int8)], zero(Float64))
zero{T}(s::DSisoRational{T})                    = tf([zero(Int8)], [one(Int8)], zero(Float64))
one{T1,T2,T3}(::Type{DSisoRational{T1,T2,T3}})  = tf([one(Int8)], [one(Int8)], zero(Float64))
one{T}(::Type{DSisoRational{T}})                = tf([one(Int8)], [one(Int8)], zero(Float64))
one(::Type{DSisoRational})                      = tf([one(Int8)], [one(Int8)], zero(Float64))
one{T}(s::DSisoRational{T})                     = tf([one(Int8)], [one(Int8)], zero(Float64))

# interface implementation

function zeros(s::DSisoRational)
  z::Vector{Complex{Float64}} = roots(s.num)
end

function poles(s::DSisoRational)
  p::Vector{Complex{Float64}} = roots(s.den)
end

function numvec(s::DSisoRational)
  coeffs(s.num)[end:-1:1]
end

function denvec(s::DSisoRational)
  coeffs(s.den)[end:-1:1]
end

function numpoly(s::DSisoRational)
  copy(s.num)
end

function denpoly(s::DSisoRational)
  copy(s.den)
end

function zpkdata(s::DSisoRational)
  (zeros(s), poles(s), s.num[1]/s.den[1])
end

function samplingtime(s::DSisoRational)
  s.Ts
end

# overload printing functions

function show(io::IO, s::DSisoRational)
  println(io, "Discrete time rational transfer function model")
  println(io, "\ty = Gu")
  if s.Ts > 0
    println(io, "with Ts=", s.Ts, ".")
  elseif s.Ts == 0
    println(io, "with Ts=unspecified.")
  end
end

function showall(io::IO, s::DSisoRational)
  show(io, s)
  println(io, "")
  printtransferfunction(io::IO, s)
end

# overload mathematical operations

function +(s1::DSisoRational, s2::DSisoRational)
  Ts::Float64
  if s1.Ts == s2.Ts || isnan(s2.Ts)
    Ts = s1.Ts
  elseif isnan(s1.Ts)
    Ts = s2.Ts
  else
    warn("Sampling time mismatch")
    throw(DomainError())
  end
  den1,den2,dengcd   = rmgcd(s1.den, s2.den)
  tf(s1.num*den2 + s2.num*den1, den1*den2*dengcd, Ts)
end
+{T<:Real}(s::DSisoRational, n::T)       = tf(s.num + n*s.den, copy(s.den), s.Ts)
+{T<:Real}(n::T, s::DSisoRational)       = s + n

.+{T<:Real}(s::DSisoRational, n::T)      = s + n
.+{T<:Real}(n::T, s::DSisoRational)      = s + n
.+(s1::DSisoRational, s2::DSisoRational) = +(s1,-s2)

-{T}(s::DSisoRational{T})                = tf(-s.num, copy(s.den), s.Ts)

-(s1::DSisoRational, s2::DSisoRational)  = +(s1,-s2)
-{T<:Real}(n::T, s::DSisoRational)       = +(n, -s)
-{T<:Real}(s::DSisoRational, n::T)       = +(s, -n)

.-{T<:Real}(s::DSisoRational, n::T)      = +(s, -n)
.-{T<:Real}(n::T, s::DSisoRational)      = +(-n, s)
.-(s1::DSisoRational, s2::DSisoRational) = +(s1,-s2)

function *(s1::DSisoRational, s2::DSisoRational)
  Ts::Float64
  if s1.Ts == s2.Ts || isnan(s2.Ts)
    Ts = s1.Ts
  elseif isnan(s1.Ts)
    Ts = s2.Ts
  else
    warn("Sampling time mismatch")
    throw(DomainError())
  end
  num1,den2,gcd1 = rmgcd(s1.num, s2.den)
  den1,num2,gcd2 = rmgcd(s1.den, s2.num)
  tf(num1*num2, den1*den2, Ts)
end
*{T<:Real}(s::DSisoRational, n::T)       = tf(s.num*n, copy(s.den), s.Ts)
*{T<:Real}(n::T, s::DSisoRational)       = *(s, n)

.*{T<:Real}(s::DSisoRational, n::T)      = *(s, n)
.*{T<:Real}(n::T, s::DSisoRational)      = *(n, s)
.*(s1::DSisoRational, s2::DSisoRational) = *(s1, s2)

/(s1::DSisoRational, s2::DSisoRational)  = s1*(1/s2)
/{T<:Real}(n::T, s::DSisoRational)       = tf(n*s.den, copy(s.num), s.Ts)
/{T<:Real}(s::DSisoRational, n::T)       = s*(1/n)

./{T<:Real}(n::T, s::DSisoRational)      = /(n, s)
./{T<:Real}(s::DSisoRational, n::T)      = /(s, n)
./(s1::DSisoRational, s2::DSisoRational) = /(s1, s2)

function ==(s1::DSisoRational, s2::DSisoRational)
  s1.Ts == s2.Ts && s1.num == s2.num &&
    (s1.den == s2.den || s1.num == zero(s1.num))
    # TODO scaling of num and den
end

!=(s1::DSisoRational, s2::DSisoRational) = !(s1 == s2)

function isapprox(s1::DSisoRational, s2::DSisoRational,
    rtol::Real=sqrt(eps()), atol::Real=0)
  sdiff = s2-s1
  return norm(sdiff.num) < rtol
end
