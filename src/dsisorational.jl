immutable DSisoRational{T<:AbstractFloat} <: DSisoTf{T}
  num::Poly{T}
  den::Poly{T}
  Ts::T
  function call{T}(::Type{DSisoRational}, num::Vector{T}, den::Vector{T}, Ts::T)
    Ts_ = max(Ts, zero(T))
    pnum = Poly(num,"z^-1")
    pden = Poly(den,"z^-1")
    new{T}(pnum, pden, Ts_)
  end
end

function tf{T<:AbstractFloat}(num::Vector{T}, den::Vector{T}, Ts::T)
  DSisoRational(num_[end:-1:1], den_[end:-1:1], Ts)
end

function tf{T<:AbstractFloat}(num::Poly{T}, den::Poly{T}, Ts::T)
  tf(coeffs(num)[end:-1:1], coeffs(den)[end:-1:1], Ts)
end

function tf{T1<:Real}(gain::T1, Ts::T1)
  DSisoRational([gain], [one(T1)], Ts)
end

function tf{T1<:Real, T2<:Real, T3<:Real}(num::Vector{T1}, den::Vector{T2}, Ts::T3)
  T = promote_type(T1, T2, T3, Float16)
  num_ = convert(Vector{T}, num)
  den_ = convert(Vector{T}, den)
  Ts_ = convert(T, Ts)
  DSisoRational(num_[end:-1:1], den_[end:-1:1], Ts_)
end

function tf{T1<:Real, T2<:Real, T3<:Real}(num::Poly{T1}, den::Poly{T2}, Ts::T3)
  tf(coeffs(num)[end:-1:1], coeffs(den)[end:-1:1], Ts)
end

function tf{T1<:Real, T2<:Real}(gain::T1, Ts::T2)
  tf([gain], [one(T1)], Ts)
end

Base.promote_rule{T1, T2}(::Type{DSisoRational{T1}}, ::Type{DSisoRational{T2}}) = DSisoRational{promote_type(T1, T2)}
Base.convert{T1, T2}(::Type{DSisoRational{T1}}, sys::DSisoRational{T2}) = tf(Poly{T1}(sys.num), Poly{T1}(sys.den), sys.Ts)

Base.promote_rule{T1, T2<:Real}(::Type{DSisoRational{T1}}, ::Type{T2}) = DSisoRational{promote_type(T1,T2)}
Base.convert{T1, T2<:Real}(::Type{DSisoRational{T1}}, x::T2) = tf([convert(T1,x)], [one(T1)], zero(Float64))

Base.zero{T1}(::Type{DSisoRational{T1}}) = tf(zero(Poly{T1}), one(Poly{T1}),zero(Float64))
Base.zero{T1}(s::DSisoRational{T1}) = tf(zero(Poly{T1}), one(Poly{T1}), s.Ts)
Base.one{T1}(::Type{DSisoRational{T1}}) = tf(one(Poly{T1}), one(Poly{T1}),zero(Float64))
Base.one{T1}(s::DSisoRational{T1}) = tf(one(Poly{T1}), one(Poly{T1}), s.Ts)

function zeros{T1}(s::DSisoRational{T1})
  return roots(s.num)::Vector{Complex{T1}}
end

function poles{T1}(s::DSisoRational{T1})
  return roots(s.den)::Vector{Complex{T1}}
end

function numvec{T1}(s::DSisoRational{T1})
  return copy(coeffs(s.num)[end:-1:1])
end

function denvec{T1}(s::DSisoRational{T1})
  return copy(coeffs(s.den)[end:-1:1])
end

function numpoly{T1}(s::DSisoRational{T1})
  return copy(s.num)
end

function denpoly{T1}(s::DSisoRational{T1})
  return copy(s.num)
end

function zpkdata{T1}(s::DSisoRational{T1})
  return (zeros(s), poles(s), num[1]/den[1])
end

function samplingtime{T1}(s::DSisoRational{T1})
  return copy(s.ts)
end

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

function +{T1, T2}(
    s1::DSisoRational{T1},
    s2::DSisoRational{T2})
  if s1.Ts == s2.Ts || s1.Ts == zero(T1) || s2.Ts == zero(T2)
    Ts, Ts2 = promote(s1.Ts, s2.Ts)
    return tf(s1.num*s2.den + s2.num*s1.den, s1.den*s2.den, Ts)
  else
    warn("Sampling time mismatch")
    throw(DomainError())
  end
end
function +{T1, T2<:Real}(s::DSisoRational{T1}, n::T2)
  return tf(s.num + n*s.den, s.den, s.Ts)
end
+{T1, T2<:Real}(n::T2, s::DSisoRational{T1})  = s + n

.+{T1, T2<:Real}(s::DSisoRational{T1}, n::T2) = s + n
.+{T1, T2<:Real}(n::T2, s::DSisoRational{T1}) = s + n
.+{T1, T2}(s1::DSisoRational{T1}, s2::DSisoRational{T2}) = +(s1,-s2)

-{T1}(s::DSisoRational{T1}) = tf(-s.num, s.den, s.Ts)
function -{T1, T2}(
  s1::DSisoRational{T1},
  s2::DSisoRational{T2})
  return +(s1,-s2)
end
-{T1, T2<:Real}(n::T2, s::DSisoRational{T1})  = +(-n, s)
-{T1, T2<:Real}(s::DSisoRational{T1}, n::T2)  = +(s, -n)

.-{T1, T2<:Real}(s::DSisoRational{T1}, n::T2) = +(s, -n)
.-{T1, T2<:Real}(n::T2, s::DSisoRational{T1}) = +(-n, s)
.-{T1, T2}(s1::DSisoRational{T1}, s2::DSisoRational{T2}) = +(s1,-s2)

function *{T1, T2}(
    s1::DSisoRational{T1},
    s2::DSisoRational{T2})
  if s1.Ts == s2.Ts || s1.Ts == zero(T1) || s2.Ts == zero(T2)
    Ts, Ts2 = promote(s1.Ts, s2.Ts)
    tf(s1.num*s2.num, s1.den*s2.den, Ts)
  else
    warn("Sampling time mismatch")
    throw(DomainError())
  end
end
*{T1, T2<:Real}(s::DSisoRational{T1}, n::T2)  = tf(s.num*n, s.den, s.Ts)
*{T1, T2<:Real}(n::T2, s::DSisoRational{T1})  = *(s, n)

.*{T1, T2<:Real}(s::DSisoRational{T1}, n::T2) = *(s, n)
.*{T1, T2<:Real}(n::T2, s::DSisoRational{T1}) = *(n, s)
.*{T1, T2}(s1::DSisoRational{T1}, s2::DSisoRational{T2}) = *(s1, s2)

function /{T1, T2}(
    s1::DSisoRational{T1},
    s2::DSisoRational{T2})
    return s1*(1/s2)
end
/{T1, T2<:Real}(n::T2, s::DSisoRational{T1})  = tf(n*s.den, s.num, s.Ts)
/{T1, T2<:Real}(s::DSisoRational{T1}, n::T2)  = s*(1/n)

./{T1, T2<:Real}(n::T2, s::DSisoRational{T1}) = n*(1/s)
./{T1, T2<:Real}(s::DSisoRational{T1}, n::T2) = s*(1/n)
./{T1, T2}(s1::DSisoRational{T1}, s2::DSisoRational{T2}) = s1*(1/s2)

function =={T1, T2}(
    s1::DSisoRational{T1},
    s2::DSisoRational{T2})
  fields = [:num, :den, :Ts]
  for field in fields
      if getfield(s1, field) != getfield(s2, field)
          return false
      end
  end
  true
end

!={T1, T2}(s1::DSisoRational{T1}, s2::DSisoRational{T2}) = !(s1 == s2)

function isapprox{T1, T2}(s1::DSisoRational{T1}, s2::DSisoRational{T2})
  # TODO: Implement
end
