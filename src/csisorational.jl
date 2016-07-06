# definition of continuous rational transfer function

immutable CSisoRational{T<:Real, T1<:Real, T2<:Real} <: CSisoTf{T}
  num::Poly{T1}
  den::Poly{T2}

  function call{T1<:Real, T2<:Real}(::Type{CSisoRational}, num::Poly{T1},
      den::Poly{T2})
    T = promote_type(eltype(num), eltype(den))
    new{T,eltype(num),eltype(den)}(num, den)
  end
end

# creation of continuous rational transfer functions

function tf{V1<:AbstractVector, V2<:AbstractVector}(num::V1, den::V2)
  @assert eltype(num) <: Number string("num must be vector of T<:Number elements")
  @assert eltype(den) <: Number string("den must be vector of T<:Number elements")
  pnum = Poly(num[end:-1:1])
  pden = Poly(den[end:-1:1])
  CSisoRational(pnum, pden)
end

function tf{T1<:Real, T2<:Real}(num::Poly{T1}, den::Poly{T2})
  CSisoRational(num, den)
end

function tf{T1<:Real}(gain::T1)
  CSisoRational(Poly([gain]), Poly([one(T1)]))
end

# conversion and promotion

promote_rule{T<:Real,T1,T2,T3}(::Type{CSisoRational{T1,T2,T3}}, ::Type{T}) =
  CSisoRational
promote_rule{T1,T<:Real}(::Type{CSisoRational{T1}}, ::Type{T}) = CSisoRational
convert{T<:Real}(::Type{CSisoRational}, x::T) = tf([x], [one(T)])

# overloading identities

zero{T1,T2,T3}(::Type{CSisoRational{T1,T2,T3}}) = tf([zero(Int8)], [one(Int8)])
zero{T}(::Type{CSisoRational{T}})               = tf([zero(Int8)], [one(Int8)])
zero(::Type{CSisoRational})                     = tf([zero(Int8)], [one(Int8)])
zero{T}(s::CSisoRational{T})                    = tf([zero(Int8)], [one(Int8)])
one{T1,T2,T3}(::Type{CSisoRational{T1,T2,T3}})  = tf([one(Int8)], [one(Int8)])
one{T}(::Type{CSisoRational{T}})                = tf([one(Int8)], [one(Int8)])
one(::Type{CSisoRational})                      = tf([one(Int8)], [one(Int8)])
one{T}(s::CSisoRational{T})                     = tf([one(Int8)], [one(Int8)])

# interface implementation

function zeros(s::CSisoRational)
  roots(s.num)
end

function poles(s::CSisoRational)
  roots(s.den)
end

function numvec(s::CSisoRational)
  coeffs(s.num)[end:-1:1]
end

function denvec(s::CSisoRational)
  coeffs(s.den)[end:-1:1]
end

function numpoly(s::CSisoRational)
  copy(s.num)
end

function denpoly(s::CSisoRational)
  copy(s.den)
end

function zpkdata(s::CSisoRational)
  (zeros(s), poles(s), s.num[1]/s.den[1])
end

function samplingtime(s::CSisoRational)
  -one(Float64)
end

# overload printing functions

function show(io::IO, s::CSisoRational)
  println(io, "Continuous time rational transfer function model")
  println(io, "\ty = Gu")
end

function showall(io::IO, s::CSisoRational)
  show(io, s)
  println(io, "")
  printtransferfunction(io::IO, s)
end

# overload mathematical operations

function +(s1::CSisoRational, s2::CSisoRational)
  den1,den2,dengcd   = rmgcd(s1.den, s2.den)
  tf(s1.num*den2 + s2.num*den1, den1*den2*dengcd)
end
+{T<:Real}(s::CSisoRational, n::T)       = tf(s.num + n*s.den, copy(s.den))
+{T<:Real}(n::T, s::CSisoRational)       = s + n

.+{T<:Real}(s::CSisoRational, n::T)      = s + n
.+{T<:Real}(n::T, s::CSisoRational)      = s + n
.+(s1::CSisoRational, s2::CSisoRational) = +(s1, s2)

-(s::CSisoRational)                      = tf(-s.num, copy(s.den))

-(s1::CSisoRational, s2::CSisoRational)  = +(s1,-s2)
-{T<:Real}(n::T, s::CSisoRational)       = +(n, -s)
-{T<:Real}(s::CSisoRational, n::T)       = +(s, -n)

.-{T<:Real}(s::CSisoRational, n::T)      = -(s, n)
.-{T<:Real}(n::T, s::CSisoRational)      = -(n, s)
.-(s1::CSisoRational, s2::CSisoRational) = +(s1, -s2)

function *(s1::CSisoRational, s2::CSisoRational)
  num1,den2,gcd1   = rmgcd(s1.num, s2.den)
  den1,num2,gcd2   = rmgcd(s1.den, s2.num)
  tf(num1*num2, den1*den2)
end

*{T<:Real}(s::CSisoRational, n::T)       = tf(s.num*n, copy(s.den))
*{T<:Real}(n::T, s::CSisoRational)       = *(s, n)

.*{T<:Real}(s::CSisoRational, n::T)      = *(s, n)
.*{T<:Real}(n::T, s::CSisoRational)      = *(n, s)
.*(s1::CSisoRational, s2::CSisoRational) = *(s1, s2)

/(s1::CSisoRational, s2::CSisoRational)  = s1*(1/s2)
/{T<:Real}(n::T, s::CSisoRational)       = tf(n*s.den, copy(s.num))
/{T<:Real}(s::CSisoRational, n::T)       = s*(1/n)

./{T<:Real}(n::T, s::CSisoRational)      = /(n, s)
./{T<:Real}(s::CSisoRational, n::T)      = /(s, n)
./(s1::CSisoRational, s2::CSisoRational) = /(s1, s2)

function ==(s1::CSisoRational, s2::CSisoRational)
  s1.num == s2.num && (s1.den == s2.den || s1.num == zero(s1.num))
end

!=(s1::CSisoRational, s2::CSisoRational) = !(s1 == s2)

function isapprox(s1::CSisoRational, s2::CSisoRational,
  rtol::Real=sqrt(eps()), atol::Real=0)
  sdiff = s2-s1
  norm(sdiff.num) < rtol
end
