immutable DSisoZpk{T} <: DSisoTf{T}
  z::Vector{Complex{T}}
  p::Vector{Complex{T}}
  k::T
  Ts::T
  function call{T}(::Type{DSisoZpk}, z::Vector{Complex{T}},
    p::Vector{Complex{T}}, k::T, Ts::T)
    Ts_ = max(Ts, zero(T))
    new{T}(copy(z), copy(p), k, Ts_)
  end
end

function zpk{T1<:AbstractFloat}(z::Vector{Complex{T1}},
  p::Vector{Complex{T1}}, k::T1, Ts::T1))
  DSisoZpk(z, p, k, Ts)
end

function zpk{T1<:Number, T2<:Number, T2<:Real, T4<:Real}(z::Vector{T1}, p::Vector{T2}, k::T2, Ts::T4)
  T   = promote_type(T1, T2, T2, T4, Float16) # ensure AbstractFloat
  z_  = convert(Vector{Complex{T}},z)
  p_  = convert(Vector{Complex{T}},p)
  k_  = convert(T,k)
  Ts_ = convert(T,Ts)
  DSisoZpk(z_, p_, k_, Ts_)
end

function zpk{T1<:Number, T2<:Real}(k::T1, Ts::T2)
  zpk(Vector{Complex{T1}}(), Vector{Complex{T1}}(), k, Float64(Ts))
end

Base.promote_rule{T1<:AbstractFloat, T2<:AbstractFloat}(::Type{DSisoZpk{T1}}, ::Type{DSisoZpk{T2}}) = DSisoZpk{promote_type(T1, T2), promote_type(T1, T2)}
Base.convert{T1<:AbstractFloat, T2<:AbstractFloat}(::Type{DSisoZpk{T1}}, s::DSisoZpk{T2}) =
  DSisoZpk(convert(Vector{Complex{T1}}, s.z), convert(Vector{Complex{T1}}, s.p),  convert(T1, s.k), convert(T1, s.Ts))

Base.promote_rule{T1<:AbstractFloat, T2<:Real}(::Type{DSisoZpk{T1}}, ::Type{T2}) = DSisoZpk{promote_type(T1, T2), promote_type(T1, T2)}
Base.convert{T1<:AbstractFloat, T2<:Real}(::Type{DSisoZpk{T1}}, x::T2) = zpk([one(T1)], [one(T1)], convert(T1, x), zero(T1))

Base.one{T1<:AbstractFloat, T2<:AbstractFloat}(::Type{DSisoZpk{T1}}) = zpk(Vector{Complex{T1}}(), Vector{Complex{T1}}(), one(T1), zero(T1))
Base.one{T1<:AbstractFloat, T2<:AbstractFloat}(s::DSisoZpk{T1}) = zpk(Vector{Complex{T1}}(), Vector{Complex{T1}}(), one(T1), s.Ts)
Base.zero{T1<:AbstractFloat, T2<:AbstractFloat}(::Type{DSisoZpk{T1}}) = zpk(Vector{Complex{T1}}(), Vector{Complex{T1}}(), zero(T1), zero(T1))
Base.zero{T1<:AbstractFloat, T2<:AbstractFloat}(s::DSisoZpk{T1}) = zpk(Vector{Complex{T1}}(), Vector{Complex{T1}}(), zero(T1), s.Ts)

function zeros{T1<:AbstractFloat}(s::DSisoZpk{T1})
  return copy(s.z)
end

function poles{T1<:AbstractFloat}(s::DSisoZpk{T1})
  return copy(s.p)
end

function numvec{T1<:AbstractFloat}(s::DSisoZpk{T1})
  return coeffs(poly(s.z))[end:-1:1]
end

function denvec{T1<:AbstractFloat}(s::DSisoZpk{T1})
  return coeffs(poly(s.p))[end:-1:1]
end

function numpoly{T1<:AbstractFloat}(s::DSisoZpk{T1})
  return poly(s.z)
end

function denpoly{T1<:AbstractFloat}(s::DSisoZpk{T1})
  return poly(s.p)
end

function zpkdata{T1<:AbstractFloat}(s::DSisoZpk{T1})
  return copy(s.k)
end

function samplingtime{T1<:AbstractFloat}(s::DSisoZpk{T1})
  return copy(s.ts)
end

ndims(s::DSisoZpk)  = 1
size(s::DSisoZpk)   = 1

function getindex(s::DSisoZpk, idx::Int)
  if idx != 1
    warn("A SISO transfer function only has one element")
    throw(DomainError())
  end
  s
end

function getindex(s::DSisoZpk, rows, cols)
  if rows != 1 || cols != 1
    warn("A SISO transfer function only has one element")
    throw(DomainError())
  end
  s
end

getindex(s::DSisoZpk, ::Colon, ::Colon) = s

function getindex(s::DSisoZpk, ::Colon, cols)
  if cols != 1
    warn("A SISO transfer function only has one element")
    throw(DomainError())
  end
  s
end

function getindex(s::DSisoZpk, rows, ::Colon)
  if rows != 1
    warn("A SISO transfer function only has one element")
    throw(DomainError())
  end
  s
end

start(s::DSisoZpk)           = 1
next(s::DSisoZpk, state)     = (s.m[state], state+1)
done(s::DSisoZpk, state)     = state > length(s)
eltype(::Type{DSisoZpk{T1}}) = DSisoZpk{T1}
length(s::DSisoZpk)          = 1
eachindex(s::DSisoZpk)       = 1:length(s)
endof(s::DSisoZpk)           = length(s)

showcompact(io::IO, s::DSisoZpk) = print(io, summary(s))

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

Base.print(io::IO, s::DSisoZpk) = show(io, t)

function printtransferfunction(io::IO, s::DSisoZpk)
  numstr = sprint(print_polyroots, s.z, "z")
  denstr = sprint(print_polyroots, s.p, "z")
  gainstr = s.k[1]==1.0 ? "" : "$(round(s.k[1], 6))"

  # Figure out the length of the separating line
  len_num = length(numstr)
  len_den = length(denstr)
  len_gain = length(gainstr)
  dashcount = max(len_num, len_den)

  # Center the numerator or denominator
  if len_num < dashcount
    numstr = "$(repeat(" ", div(dashcount - len_num, 2)))$numstr"
  else
    denstr = "$(repeat(" ", div(dashcount - len_den, 2)))$denstr"
  end
  println(io, repeat(" ", len_gain+1), numstr)
  println(io, gainstr, " ", repeat("-", dashcount))
  println(io, repeat(" ", len_gain+1), denstr)
end

function summary(io::IO, s::DSisoZpk)
  if s.Ts > 0
    println(io, string("zpk(nu=1, ny=1, Ts=", s.Ts, "."))
  elseif s.Ts == 0
    println(io, string("zpk(nu=1, ny=1, Ts=unspecified."))
  end
end

function +{T1<:AbstractFloat, T2<:AbstractFloat}(
  s1::DSisoZpk{T1}, s2::DSisoZpk{T2})
  if s1.Ts == s2.Ts
    Z = s1.k*poly(s1.z)*poly(s2.p) + s2.k*poly(s2.z)*poly(s1.p)
    z = roots(Z)
    p = vcat(s1.p, s2.p)
    k = Z[end] # Poly is now reverse order
    return zpk(z, p, k, s1.Ts)
  else
    warn("Sampling time mismatch")
    throw(DomainError())
  end
end

function +{T1<:AbstractFloat, T2<:Real}(
    s::DSisoZpk{T1}, n::T2)
  Z = s.k*poly(s.z) + n*poly(s.p)
  z = roots(Z)
  p = s.p
  k = Z[end] # Poly is now reverse order
  return zpk(z_, p_, k, s.Ts)
end
# +{T1<:AbstractFloat, T2<:Real}(s::DSisoZpk{T1}, n::T2)  = +(s, zpk(n,s.Ts))
+{T1<:AbstractFloat, T2<:Real}(n::T2, s::DSisoZpk{T1})  = +(s, n)

.+{T1<:AbstractFloat, T2<:Real}(s::DSisoZpk{T1}, n::T2) = +(n, s)
.+{T1<:AbstractFloat, T2<:Real}(n::T2, s::DSisoZpk{T1}) = +(s, n)
.+{T1<:AbstractFloat, T2<:AbstractFloat}(s1::DSisoZpk{T1}, s2::DSisoZpk{T2}) = +(s1, s2)

-{T1<:AbstractFloat}(s::DSisoZpk{T1}) = zpk(s.z, s.p, -s.k, s.Ts)
function -{T1<:AbstractFloat, T2<:AbstractFloat}(
  s1::DSisoZpk{T1}, s2::DSisoZpk{T2})
  return +(s1, -s2)
end
-{T1<:AbstractFloat, T2<:Real}(s::DSisoZpk{T1}, n::T2)  = +(s, -n)
-{T1<:AbstractFloat, T2<:Real}(n::T2, s::DSisoZpk{T1})  = +(n, -s)

.-{T1<:AbstractFloat, T2<:Real}(s::DSisoZpk{T1}, n::T2) = +(s, -n)
.-{T1<:AbstractFloat, T2<:Real}(n::T2, s::DSisoZpk{T1}) = +(n, -s)
.-{T1<:AbstractFloat, T2<:AbstractFloat}(s1::DSisoZpk{T1}, s2::DSisoZpk{T2}) = +(s1, -s2)

function *{T1<:AbstractFloat, T2<:AbstractFloat}(
  s1::DSisoZpk{T1}, s2::DSisoZpk{T2})
  if s1.Ts == s2.Ts
    Ts, Ts2 = promote(s1.Ts, s2.Ts)
    z = vcat(s1.z, s2.z)
    p = vcat(s1.p, s2.p)
    k = s1.k*s2.k
    return zpk(z, p, k, Ts)
  else
    warn("Sampling time mismatch")
    throw(DomainError())
  end
end
*{T1<:AbstractFloat, T2<:Real}(s::DSisoZpk{T1}, n::T2)  = zpk(s.z, s.p, n*s.k, s.Ts)
*{T1<:AbstractFloat, T2<:Real}(n::T2, s::DSisoZpk{T1})  = *(s, n)

.*{T1<:AbstractFloat, T2<:Real}(s::DSisoZpk{T1}, n::T2) = *(n, s)
.*{T1<:AbstractFloat, T2<:Real}(n::T2, s::DSisoZpk{T1}) = *(s, n)
.*{T1<:AbstractFloat, T2<:AbstractFloat}(s1::DSisoZpk{T1}, s2::DSisoZpk{T2}) = *(s1, s2)

function /{T1<:AbstractFloat, T2<:Real}(
  n::T2, s::DSisoZpk{T1})
  return zpk(s.p, s.z, n./s.k, s.Ts)
end
/{T1<:AbstractFloat, T2<:AbstractFloat}(s1::DSisoZpk{T1}, s2::DSisoZpk{T2}) = s1*(1/s2)
/{T1<:AbstractFloat, T2<:Real}(s::DSisoZpk{T1}, n::T2)  = s*(1/n)

./{T1<:AbstractFloat, T2<:Real}(n::T2, s::DSisoZpk{T1}) = n*(1/s)
./{T1<:AbstractFloat, T2<:Real}(s::DSisoZpk{T1}, n::T2) = s*(1/n)
./{T1<:AbstractFloat, T2<:AbstractFloat}(s1::DSisoZpk{T1}, s2::DSisoZpk{T2}) = s1*(1/s2)

function =={T1<:AbstractFloat, T2<:AbstractFloat}(
  s1::DSisoZpk{T1}, s2::DSisoZpk{T2})
  fields = [:Ts, :z, :p, :k]
  for field in fields
    if getfield(s1, field) != getfield(s2, field)
      return false
    end
  end
  true
end

!={T1<:AbstractFloat, T2<:AbstractFloat}(
  s1::DSisoZpk{T1}, s2::DSisoZpk{T2}) = !(s1==s2)

function isapprox{T1<:AbstractFloat, T2<:AbstractFloat}(
    s1::DSisoZpk{T1}, s2::DSisoZpk{T2})
  # TODO: Implement
end
