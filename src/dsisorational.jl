# definition of discrete rational transfer function

immutable DSisoRational{T<:Real, T1<:Real, T2<:Real} <: DSisoTf{T}
  num::Poly{T1}
  den::Poly{T2}
  Ts::Float64

  function call(::Type{DSisoRational},
      num::AbstractVector, den::AbstractVector, Ts::Float64)
    @assert eltype(num) <: Number string("num must be vector of T<:Number elements")
    @assert eltype(den) <: Number string("den must be vector of T<:Number elements")

    tmp = max(Ts, zero(Float64))
    Ts_ = tmp == zero(Float64) ? NaN : tmp
    T = promote_type(eltype(num), eltype(den))
    pnum = Poly(num)
    pden = Poly(den)
    new{T,eltype(num),eltype(den)}(pnum, pden, Ts_)
  end
end

# creation of discrete rational transfer functions

function tf{T1<:Real, T2<:Real, T3<:Real}(num::AbstractVector{T1},
    den::AbstractVector{T2}, Ts::T3)
  DSisoRational(num[end:-1:1], den[end:-1:1], Float64(Ts))
end

function tf{T1<:Real, T2<:Real, T3<:Real}(num::Poly{T1}, den::Poly{T2}, Ts::T3)
  tf(coeffs(num)[end:-1:1], coeffs(den)[end:-1:1], Float64(Ts))
end

function tf{T1<:Real, T2<:Real}(gain::T1, Ts::T2)
  DSisoRational([gain], [one(T1)], Float64(Ts))
end

# conversion and promotion

promote_rule{T<:Real,T1,T2,T3}(::Type{DSisoRational{T1,T2,T3}}, ::Type{T}) =
  DSisoRational
promote_rule{T<:Real}(::Type{DSisoRational}, ::Type{T}) = DSisoRational
convert{T<:Real}(::Type{DSisoRational}, x::T)   =
  tf([x], [one(T)], zero(Float64))

# overloading identities

zero{T}(::Type{DSisoRational{T}}) = tf([zero(T)],[one(T)], zero(Float64))
zero{T}(s::DSisoRational{T})      = tf([zero(T)],[one(T)], zero(Float64))
one{T}(::Type{DSisoRational{T}})  = tf([one(T)], [one(T)], zero(Float64))
one{T}(s::DSisoRational{T})       = tf([one(T)], [one(T)], zero(Float64))

# interface implementation

function zeros(s::DSisoRational)
  z::Vector{Complex{Float64}} = roots(s.num)
end

function poles(s::DSisoRational)
  p::Vector{Complex{Float64}} = roots(s.den)
end

function numvec(s::DSisoRational)
  copy(coeffs(s.num)[end:-1:1])
end

function denvec(s::DSisoRational)
  copy(coeffs(s.den)[end:-1:1])
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
  copy(s.Ts)
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
  if s1.Ts == s2.Ts || s2.Ts == NaN
    return tf(s1.num*s2.den + s2.num*s1.den, s1.den*s2.den, s1.Ts)
  elseif s1.Ts == NaN
    
    return tf(s1.num*s2.den + s2.num*s1.den, s1.den*s2.den, s2.Ts)
  else
    warn("Sampling time mismatch")
    throw(DomainError())
  end
end
+{T<:Real}(s::DSisoRational, n::T)       = tf(s.num + n*s.den, s.den, s.Ts)
+{T<:Real}(n::T, s::DSisoRational)       = s + n

.+{T<:Real}(s::DSisoRational, n::T)      = s + n
.+{T<:Real}(n::T, s::DSisoRational)      = s + n
.+(s1::DSisoRational, s2::DSisoRational) = +(s1,-s2)

-{T}(s::DSisoRational{T})                = tf(-s.num, s.den, s.Ts)

-(s1::DSisoRational, s2::DSisoRational)  = +(s1,-s2)
-{T<:Real}(n::T, s::DSisoRational)       = +(-n, s)
-{T<:Real}(s::DSisoRational, n::T)       = +(s, -n)

.-{T<:Real}(s::DSisoRational, n::T)      = +(s, -n)
.-{T<:Real}(n::T, s::DSisoRational)      = +(-n, s)
.-(s1::DSisoRational, s2::DSisoRational) = +(s1,-s2)

function *(s1::DSisoRational, s2::DSisoRational)
  if s1.Ts == s2.Ts || s2.Ts == NaN
    return tf(s1.num*s2.num, s1.den*s2.den, s1.Ts)
  elseif s1.Ts == NaN
    return tf(s1.num*s2.num, s1.den*s2.den, s2.Ts)
  else
    warn("Sampling time mismatch")
    throw(DomainError())
  end
end
*{T<:Real}(s::DSisoRational, n::T)       = tf(s.num*n, s.den, s.Ts)
*{T<:Real}(n::T, s::DSisoRational)       = *(s, n)

.*{T<:Real}(s::DSisoRational, n::T)      = *(s, n)
.*{T<:Real}(n::T, s::DSisoRational)      = *(n, s)
.*(s1::DSisoRational, s2::DSisoRational) = *(s1, s2)

/(s1::DSisoRational, s2::DSisoRational)  = s1*(1/s2)
/{T<:Real}(n::T, s::DSisoRational)       = tf(n*s.den, s.num, s.Ts)
/{T<:Real}(s::DSisoRational, n::T)       = s*(1/n)

./{T<:Real}(n::T, s::DSisoRational)      = n*(1/s)
./{T<:Real}(s::DSisoRational, n::T)      = s*(1/n)
./(s1::DSisoRational, s2::DSisoRational) = s1*(1/s2)

function ==(s1::DSisoRational, s2::DSisoRational)
  fields = [:num, :den, :Ts]
  for field in fields
      if getfield(s1, field) != getfield(s2, field)
          return false
      end
  end
  true
end

!=(s1::DSisoRational, s2::DSisoRational) = !(s1 == s2)

function isapprox(s1::DSisoRational, s2::DSisoRational,
    rtol::Real=sqrt(eps), atol::Real=0)
  sdiff = s2-s1
  return norm(sdiff.num) < rtol*max(norm(s1.num), norm(s2.num))
end
