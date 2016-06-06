immutable DiscreteSisoRational{T1<:Real} <: DiscreteSisoTf
  num::Poly{T1}
  den::Poly{T1}
  Ts::Float64
  function call{T1}(::Type{DiscreteSisoRational}, num::Vector{T1}, den::Vector{T1}, Ts::Float64)
    Ts_ = max(Ts, zero(Float64))
    pnum = Poly(num,"z^-1")
    pden = Poly(den,"z^-1")
    new{T1}(pnum, pden, Ts_)
  end
end

function tf{T1<:Real, T2<:Real}(num::Vector{T1}, den::Vector{T2}, Ts::Float64)
  num_, den_ = promote(num, den)
  DiscreteSisoRational(num_[end:-1:1], den_[end:-1:1], Ts)
end

function tf{T1<:Real, T2<:Real}(num::Poly{T1}, den::Poly{T2}, Ts::Float64)
  tf(coeffs(num)[end:-1:1], coeffs(den)[end:-1:1], Ts)
end

function tf{T1<:Real}(gain::T1, Ts::Float64)
  tf([gain], [one(T1)], Ts)
end

Base.promote_rule{T1<:Real, T2<:Real}(::Type{DiscreteSisoRational{T1}}, ::Type{DiscreteSisoRational{T2}}) = DiscreteSisoRational{promote_type(T1, T2)}
Base.convert{T1<:Real, T2<:Real}(::Type{DiscreteSisoRational{T1}}, sys::DiscreteSisoRational{T2}) = tf(Poly{T1}(sys.num), Poly{T1}(sys.den), sys.Ts)

Base.promote_rule{T1<:Real, T2<:Real}(::Type{DiscreteSisoRational{T1}}, ::Type{T2}) = DiscreteSisoRational{promote_type(T1,T2)}
Base.convert{T1<:Real, T2<:Real}(::Type{DiscreteSisoRational{T1}}, x::T2) = tf([convert(T1,x)], [one(T1)], zero(Float64))

## IDENTITIES ##
Base.zero{T1<:Real}(::Type{DiscreteSisoRational{T1}}) = tf(zero(Poly{T1}), one(Poly{T1}),zero(Float64))
Base.zero{T1<:Real}(s::DiscreteSisoRational{T1}) = tf(zero(Poly{T1}), one(Poly{T1}), s.Ts)
Base.one{T1<:Real}(::Type{DiscreteSisoRational{T1}}) = tf(one(Poly{T1}), one(Poly{T1}),zero(Float64))
Base.one{T1<:Real}(s::DiscreteSisoRational{T1}) = tf(one(Poly{T1}), one(Poly{T1}), s.Ts)

function zeros(s::DiscreteSisoRational)
  return roots(s.num)
end

function poles(s::DiscreteSisoRational)
  return roots(s.den)
end

showcompact(io::IO, s::DiscreteSisoRational) = print(io, summary(s))

function show(io::IO, s::DiscreteSisoRational)
  println(io, "Discrete time rational transfer function model")
  println(io, "\ty = Gu")
  if s.Ts > 0
    println(io, "with Ts=", s.Ts, ".")
  elseif s.Ts == 0
    println(io, "with Ts=unspecified.")
  end
end

function showall(io::IO, s::DiscreteSisoRational)
  show(io, s)
  println(io, "")
  printtransferfunction(io::IO, s)
end

Base.print(io::IO, s::DiscreteSisoRational) = show(io, s)

function printtransferfunction(io::IO, s::DiscreteSisoRational)
  numstr = print_poly_reverse(s.num)
  denstr = print_poly_reverse(s.den)

  # Figure out the length of the separating line
  len_num = length(numstr)
  len_den = length(denstr)
  dashcount = max(len_num, len_den)

  # Center the numerator or denominator
  if len_num < dashcount
    numstr = "$(repeat(" ", div(dashcount - len_num, 2)))$numstr"
  else
    denstr = "$(repeat(" ", div(dashcount - len_den, 2)))$denstr"
  end
  println(io, numstr)
  println(io, repeat("-", dashcount))
  println(io, denstr)
end

function summary(io::IO, s::DiscreteSisoRational)
  if s.Ts > 0
    println(io, string("tf(nu=1, ny=1, Ts=", s.Ts, "."))
  elseif s.Ts == 0
    println(io, string("tf(nu=1, ny=1, Ts=unspecified."))
  end
end


function +{T1<:Real, T2<:Real}(
    s1::DiscreteSisoRational{T1},
    s2::DiscreteSisoRational{T2})
  if s1.Ts == s2.Ts
    Ts, Ts2 = promote(s1.Ts, s2.Ts)
    return tf(s1.num*s2.den + s2.num*s1.den, s1.den*s2.den, Ts)
  else
    warn("Sampling time mismatch")
    throw(DomainError())
  end
end
function +{T1<:Real, T2<:Real}(s::DiscreteSisoRational{T1}, n::T2)
  return tf(s.num + n*s.den, s.den, s.Ts)
end
+{T1<:Real, T2<:Real}(n::T2, s::DiscreteSisoRational{T1})  = s + n

.+{T1<:Real, T2<:Real}(s::DiscreteSisoRational{T1}, n::T2) = s + n
.+{T1<:Real, T2<:Real}(n::T2, s::DiscreteSisoRational{T1}) = s + n
.+{T1<:Real, T2<:Real}(s1::DiscreteSisoRational{T1}, s2::DiscreteSisoRational{T2}) = +(s1,-s2)

-{T1<:Real}(s::DiscreteSisoRational{T1}) = tf(-s.num, s.den, s.Ts)
function -{T1<:Real, T2<:Real}(
  s1::DiscreteSisoRational{T1},
  s2::DiscreteSisoRational{T2})
  return +(s1,-s2)
end
-{T1<:Real, T2<:Real}(n::T2, s::DiscreteSisoRational{T1})  = +(-n, s)
-{T1<:Real, T2<:Real}(s::DiscreteSisoRational{T1}, n::T2)  = +(s, -n)

.-{T1<:Real, T2<:Real}(s::DiscreteSisoRational{T1}, n::T2) = +(s, -n)
.-{T1<:Real, T2<:Real}(n::T2, s::DiscreteSisoRational{T1}) = +(-n, s)
.-{T1<:Real, T2<:Real}(s1::DiscreteSisoRational{T1}, s2::DiscreteSisoRational{T2}) = +(s1,-s2)

function *{T1<:Real, T2<:Real}(
    s1::DiscreteSisoRational{T1},
    s2::DiscreteSisoRational{T2})
  if s1.Ts == s2.Ts
    Ts, Ts2 = promote(s1.Ts, s2.Ts)
    tf(s1.num*s2.num, s1.den*s2.den, Ts)
  else
    warn("Sampling time mismatch")
    throw(DomainError())
  end
end
*{T1<:Real, T2<:Real}(s::DiscreteSisoRational{T1}, n::T2)  = tf(s.num*n, s.den, s.Ts)
*{T1<:Real, T2<:Real}(n::T2, s::DiscreteSisoRational{T1})  = *(s, n)

.*{T1<:Real, T2<:Real}(s::DiscreteSisoRational{T1}, n::T2) = *(s, n)
.*{T1<:Real, T2<:Real}(n::T2, s::DiscreteSisoRational{T1}) = *(n, s)
.*{T1<:Real, T2<:Real}(s1::DiscreteSisoRational{T1}, s2::DiscreteSisoRational{T2}) = *(s1, s2)

function /{T1<:Real, T2<:Real}(
    s1::DiscreteSisoRational{T1},
    s2::DiscreteSisoRational{T2})
    return s1*(1/s2)
end
/{T1<:Real, T2<:Real}(n::T2, s::DiscreteSisoRational{T1})  = tf(n*s.den, s.num, s.Ts)
/{T1<:Real, T2<:Real}(s::DiscreteSisoRational{T1}, n::T2)  = s*(1/n)

./{T1<:Real, T2<:Real}(n::T2, s::DiscreteSisoRational{T1}) = n*(1/s)
./{T1<:Real, T2<:Real}(s::DiscreteSisoRational{T1}, n::T2) = s*(1/n)
./{T1<:Real, T2<:Real}(s1::DiscreteSisoRational{T1}, s2::DiscreteSisoRational{T2}) = s1*(1/s2)

function =={T1<:Real, T2<:Real}(
    s1::DiscreteSisoRational{T1},
    s2::DiscreteSisoRational{T2})
  fields = [:num, :den, :Ts]
  for field in fields
      if getfield(s1, field) != getfield(s2, field)
          return false
      end
  end
  true
end

!={T1<:Real, T2<:Real}(s1::DiscreteSisoRational{T1}, s2::DiscreteSisoRational{T2}) = !(s1 == s2)

function isapprox{T1<:Real, T2<:Real}(s1::DiscreteSisoRational{T1}, s2::DiscreteSisoRational{T2})
  # TODO: Implement
end
