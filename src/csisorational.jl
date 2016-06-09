immutable CSisoRational{T} <: CSisoTf{T}
  num::Poly{T1}
  den::Poly{T1}
  function call{T1}(::Type{DSisoRational}, num::Vector{T1}, den::Vector{T1})
    pnum = Poly(num,"s")
    pden = Poly(den,"s")
    new{T1}(pnum, pden)
  end
end

function tf{T1<:AbstractFloat}(num::Vector{T1}, den::Vector{T1})
  DSisoRational(num_[end:-1:1], den_[end:-1:1])
end

function tf{T1<:AbstractFloat}(num::Poly{T1}, den::Poly{T1})
  tf(coeffs(num)[end:-1:1], coeffs(den)[end:-1:1])
end

function tf{T1<:AbstractFloat}(gain::T1)
  CSisoRational([gain], [one(T1)])
end

function tf{T1<:Real, T2<:Real}(num::Vector{T1}, den::Vector{T2})
  num_, den_ = promote(num, den, Float16)
  DSisoRational(num_[end:-1:1], den_[end:-1:1])
end

function tf{T1<:Real, T2<:Real}(num::Poly{T1}, den::Poly{T2})
  tf(coeffs(num)[end:-1:1], coeffs(den)[end:-1:1])
end

function tf{T1<:Real}(gain::T1)
  T = promote_type(gain, Float16)
  tf([convert(T, gain)], [one(T)])
end

Base.promote_rule{T1, T2}(::Type{DSisoRational{T1}}, ::Type{DSisoRational{T2}}) = DSisoRational{promote_type(T1, T2)}
Base.convert{T1, T2}(::Type{DSisoRational{T1}}, sys::DSisoRational{T2}) = tf(Poly{T1}(sys.num), Poly{T1}(sys.den))

Base.promote_rule{T1, T2<:Real}(::Type{DSisoRational{T1}}, ::Type{T2}) = DSisoRational{promote_type(T1,T2)}
Base.convert{T1, T2<:Real}(::Type{DSisoRational{T1}}, x::T2) = tf([convert(T1,x)], [one(T1)])

Base.zero{T1}(::Type{DSisoRational{T1}}) = tf(zero(Poly{T1}), one(Poly{T1}))
Base.zero{T1}(t::DSisoRational{T1})      = Base.zero(DSisoRational{T1})
Base.one{T1}(::Type{DSisoRational{T1}})  = tf(one(Poly{T1}), one(Poly{T1}))
Base.one{T1}(t::DSisoRational{T1})       = Base.one(DSisoRational{T1})

function zeros{T1}(s::DSisoRational{T1})
  return copy(roots(s.num))
end

function poles{T1}(s::DSisoRational{T1})
  return copy(roots(s.den))
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
  return -one(Float64)
end

ndims(s::DSisoRational)  = 1
size(s::DSisoRational)   = 1

function getindex(s::DSisoRational, idx::Int)
  if idx != 1
    warn("A SISO transfer function only has one element")
    throw(DomainError())
  end
  s
end

function getindex(s::DSisoRational, rows, cols)
  if rows != 1 || cols != 1
    warn("A SISO transfer function only has one element")
    throw(DomainError())
  end
  s
end

getindex(s::DSisoRational, ::Colon, ::Colon) = s

function getindex(s::DSisoRational, ::Colon, cols)
  if cols != 1
    warn("A SISO transfer function only has one element")
    throw(DomainError())
  end
  s
end

function getindex(s::DSisoRational, rows, ::Colon)
  if rows != 1
    warn("A SISO transfer function only has one element")
    throw(DomainError())
  end
  s
end

start(s::DSisoRational)       = 1
next(s::DSisoRational, state) = (s.m[state], state+1)
done(s::DSisoRational, state) = state > length(s)
eltype{T1}(::Type{DSisoRational{T1}}) = DSisoRational{T1}
length(s::DSisoRational) = 1
eachindex(s::DSisoRational) = 1:length(s)
endof(s::DSisoRational) = length(s)

showcompact(io::IO, s::DSisoRational) = print(io, summary(s))

function show(io::IO, s::DSisoRational)
  println(io, "Continuous time rational transfer function model")
  println(io, "\ty = Gu")
end

function showall(io::IO, s::DSisoRational)
  show(io, s)
  println(io, "")
  printtransferfunction(io::IO, s)
end

Base.print(io::IO, s::DSisoRational) = show(io, s)

function printtransferfunction(io::IO, s::DSisoRational)
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

function summary(io::IO, s::DSisoRational)
  println(io, string("tf(nu=1, ny=1)."))
end

function +{T1, T2}(
    t1::DSisoRational{T1},
    t2::DSisoRational{T2})
    return tf(t1.num*t2.den + t2.num*t1.den, t1.den*t2.den)
end
+{T1, T2<:Real}(t::DSisoRational{T1}, n::T2)  =  tf(t.num + n*t.den, t.den)
+{T1, T2<:Real}(n::T2, t::DSisoRational{T1})  = t + n

.+{T1, T2<:Real}(t::DSisoRational{T1}, n::T2) = t + n
.+{T1, T2<:Real}(n::T2, t::DSisoRational{T1}) = t + n
.+{T1, T2}(t1::DSisoRational{T1}, t2::DSisoRational{T2}) = +(t1, t2)

-{T1}(t::DSisoRational{T1}) = tf(-t.num, t.den)
-{T1, T2}(t1::DSisoRational{T1}, t2::DSisoRational{T2}) = +(t1,-t2)
-{T1, T2<:Real}(n::T2, t::DSisoRational{T1}) = +(-n,t)
-{T1, T2<:Real}(t::DSisoRational{T1}, n::T2) = +(t, -n)

.-{T1, T2<:Real}(t::DSisoRational{T1}, n::T2) = -(t, n)
.-{T1, T2<:Real}(n::T2, t::DSisoRational{T1}) = -(n, t)
.-{T1, T2}(t1::DSisoRational{T1}, t2::DSisoRational{T2}) = +(t1, -t2)

*{T1, T2}(t1::DSisoRational{T1}, t2::DSisoRational{T2}) = tf(t1.num*t2.num, t1.den*t2.den)
*{T1, T2<:Real}(t::DSisoRational{T1}, n::T2) = tf(t.num*n, t.den)
*{T1, T2<:Real}(n::T2, t::DSisoRational{T1}) = *(t, n)

.*{T1, T2<:Real}(t::DSisoRational{T1}, n::T2) = *(t, n)
.*{T1, T2<:Real}(n::T2, t::DSisoRational{T1}) = *(n, t)
.*{T1, T2}(t1::DSisoRational{T1}, t2::DSisoRational{T2}) = *(t1, t2)

function /{T1, T2}(
    t1::DSisoRational{T1},
    t2::DSisoRational{T2})
    return t1*(1/t2)
end
/{T1, T2<:Real}(n::T2, t::DSisoRational{T1})  = tf(t.num/n, t.den)
/{T1, T2<:Real}(t::DSisoRational{T1}, n::T2)  = t*(1/n)
./{T1, T2<:Real}(n::T2, t::DSisoRational{T1}) = tf(t.num/n, t.den)
./{T1, T2<:Real}(t::DSisoRational{T1}, n::T2) = t*(1/n)
./{T1, T2}(t1::DSisoRational{T1}, t2::DSisoRational{T2}) = /(t1, t2)

function =={T1, T2}(
    t1::DSisoRational{T1},
    t2::DSisoRational{T2})
  fields = [:num, :den]
  for field in fields
      if getfield(t1, field) != getfield(t2, field)
          return false
      end
  end
  return true
end

!={T1, T2}(s1::DSisoRational{T1}, s2::DSisoRational{T2}) = !(s1 == s2)

function isapprox{T1, T2}(s1::DSisoRational{T1}, s2::DSisoRational{T2})
  # TODO: Implement
end
