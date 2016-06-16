# definition of continuous rational transfer function

immutable CSisoRational{T<:Real, T1<:Real, T2<:Real} <: CSisoTf{T}
  num::Poly{T1}
  den::Poly{T2}

  function call(::Type{CSisoRational}, num::AbstractVector,
      den::AbstractVector)
    @assert eltype(num) <: Number string("num must be vector of T<:Number elements")
    @assert eltype(den) <: Number string("den must be vector of T<:Number elements")
    pnum = Poly(num)
    pden = Poly(den)

    T = promote_type(eltype(z), eltype(p))
    new{T,eltype(num),eltype(den)}(pnum, pden)
  end
end

# creation of continuous rational transfer functions

function tf{T1<:Real, T2<:Real}(num::Vector{T1}, den::Vector{T2})
  CSisoRational(num_[end:-1:1], den_[end:-1:1])
end

function tf{T1<:Real, T2<:Real}(num::Vector{T1}, den::Vector{T2})
  CSisoRational(num_[end:-1:1], den_[end:-1:1])
end

function tf{T1<:Real}(gain::T1)
  CSisoRational([gain], [one(T1)])
end

# conversion and promotion

promote_rule{T<:Real,T1,T2,T3}(::Type{CSisoRational{T1,T2,T3}}, ::Type{T}) =
  CSisoRational
promote_rule{T<:Real}(::Type{CSisoRational}, ::Type{T}) = CSisoRational
convert{T<:Real}(::Type{CSisoRational}, x::T) = tf([x], [one(T)])

# overloading identities

zero{T}(::Type{CSisoRational{T}}) = tf([zero(T)], [one(T)])
zero{T}(s::CSisoRational{T})      = zero(CSisoRational{T})
one{T}(::Type{CSisoRational{T}})  = tf([one(T)], [one(T)])
one{T}(s::CSisoRational{T})       = one(CSisoRational{T})

# interface implementation

function zeros(s::CSisoRational)
  copy(roots(s.num))
end

function poles(s::CSisoRational)
  copy(roots(s.den))
end

function numvec(s::CSisoRational)
  copy(coeffs(s.num)[end:-1:1])
end

function denvec(s::CSisoRational)
  copy(coeffs(s.den)[end:-1:1])
end

function numpoly(s::CSisoRational)
  copy(s.num)
end

function denpoly(s::CSisoRational)
  copy(s.num)
end

function zpkdata(s::CSisoRational)
  (zeros(s), poles(s), num[1]/den[1])
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

+(s1::CSisoRational, s2::CSisoRational)  =
  tf(s1.num*s2.den + s2.num*s1.den, s1.den*s2.den)
+{T<:Real}(s::CSisoRational, n::T)       = tf(s.num + n*s.den, s.den)
+{T<:Real}(n::T, s::CSisoRational)       = s + n

.+{T<:Real}(s::CSisoRational, n::T)      = s + n
.+{T<:Real}(n::T, s::CSisoRational)      = s + n
.+(s1::CSisoRational, s2::CSisoRational) = +(s1, s2)

-(s::CSisoRational)                      = tf(-s.num, s.den)
-(s1::CSisoRational, s2::CSisoRational)  = +(s1,-s2)
-{T<:Real}(n::T, s::CSisoRational)       = +(-n, s)
-{T<:Real}(s::CSisoRational, n::T)       = +(s, -n)

.-{T<:Real}(s::CSisoRational, n::T)      = -(s, n)
.-{T<:Real}(n::T, s::CSisoRational)      = -(n, s)
.-(s1::CSisoRational, s2::CSisoRational) = +(s1, -s2)

*(s1::CSisoRational, s2::CSisoRational)  = tf(s1.num*s2.num, s1.den*s2.den)
*{T<:Real}(s::CSisoRational, n::T)       = tf(s.num*n, s.den)
*{T<:Real}(n::T, s::CSisoRational)       = *(s, n)

.*{T<:Real}(s::CSisoRational, n::T)      = *(s, n)
.*{T<:Real}(n::T, s::CSisoRational)      = *(n, s)
.*(s1::CSisoRational, s2::CSisoRational) = *(s1, s2)

/(s1::CSisoRational, s2::CSisoRational)  = s1*(1/s2)
/{T<:Real}(n::T, s::CSisoRational)       = tf(s.num/n, s.den)
/{T<:Real}(s::CSisoRational, n::T)       = s*(1/n)

./{T<:Real}(n::T, s::CSisoRational)      = tf(s.num/n, s.den)
./{T<:Real}(s::CSisoRational, n::T)      = s*(1/n)
./(s1::CSisoRational, s2::CSisoRational) = /(s1, s2)

function ==(s1::CSisoRational, s2::CSisoRational)
  fields = [:num, :den]
  for field in fields
      if getfield(s1, field) != getfield(s2, field)
          return false
      end
  end
  true
end

!=(s1::CSisoRational, s2::CSisoRational) = !(s1 == s2)

function isapprox(s1::CSisoRational, s2::CSisoRational,
  rtol::Real=sqrt(eps), atol::Real=0)
  sdiff = s2-s1
  norm(sdiff.num) < rtol*max(norm(s1.num), norm(s2.num))
end
