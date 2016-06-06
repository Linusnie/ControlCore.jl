immutable DiscreteSisoZpk{T1<:AbstractFloat, T2<:AbstractFloat} <: DiscreteSisoTf
  z::Vector{Complex{T1}}
  p::Vector{Complex{T1}}
  k::T2
  Ts::Float64
  function call{T1, T2}(::Type{DiscreteSisoZpk}, z::Vector{Complex{T1}},
    p::Vector{Complex{T1}}, k::T2, Ts::Float64)
    Ts_ = max(Ts, zero(Float64))
    new{T1, T2}(z, p, k, Ts_)
  end
end

function zpk{T1<:Number, T2<:Number, T3<:Real, T4<:Real}(z::Vector{T1}, p::Vector{T2}, k::T3, Ts::T4)
  T  = promote_type(real(T1), real(T2))
  z_ = convert(Vector{Complex{Float64}},z)
  p_ = convert(Vector{Complex{Float64}},p)
  DiscreteSisoZpk(z_, p_, Float64(k), Float64(Ts))
end

function zpk{T1<:Real, T2<:Real}(k::T1, Ts::T2)
  zpk(Vector{Complex{T1}}(), Vector{Complex{T1}}(), k, Float64(Ts))
end

Base.promote_rule{T1<:Real, T2<:Real, T3<:Real, T4<:Real}(::Type{DiscreteSisoZpk{T1, T2}}, ::Type{DiscreteSisoZpk{T3, T4}}) = DiscreteSisoZpk{promote_type(T1, T3), promote_type(T2, T4)}
Base.convert{T1<:Real, T2<:Real, T3<:Real, T4<:Real}(::Type{DiscreteSisoZpk{T1, T2}}, s::DiscreteSisoZpk{T3, T4}) = DiscreteSisoZpk(s.z, s.p,  convert(T2, s.k), s.Ts)

Base.promote_rule{T1<:Real, T2<:Real, T3<:Real}(::Type{DiscreteSisoZpk{T1, T2}}, ::Type{T3}) = DiscreteSisoZpk{T1, promote_type(T2, T3)}
Base.convert{T1<:Real, T2<:Real, T3<:Real}(::Type{DiscreteSisoZpk{T1, T2}}, x::T3) = zpk([one(T1)], [one(T1)], convert(T2, x), zero(Float64))

Base.one{T1<:Real, T2<:Real}(::Type{DiscreteSisoZpk{T1, T2}}) = zpk(Vector{Complex{T1}}(), Vector{Complex{T1}}(), one(T2), zero(Float64))
Base.one{T1<:Real, T2<:Real}(s::DiscreteSisoZpk{T1, T2}) = zpk(Vector{Complex{T1}}(), Vector{Complex{T1}}(), one(T2), s.Ts)
Base.zero{T1<:Real, T2<:Real}(::Type{DiscreteSisoZpk{T1, T2}}) = zpk(Vector{Complex{T1}}(), Vector{Complex{T1}}(), zero(T2), zero(Float64))
Base.zero{T1<:Real, T2<:Real}(s::DiscreteSisoZpk{T1, T2}) = zpk(Vector{Complex{T1}}(), Vector{Complex{T1}}(), zero(T2), s.Ts)

function zeros(s::DiscreteSisoZpk)
  return s.z
end

function poles(s::DiscreteSisoZpk)
  return s.p
end

showcompact(io::IO, s::DiscreteSisoZpk) = print(io, summary(s))

function show(io::IO, s::DiscreteSisoZpk)
  println(io, "Discrete time zpk transfer function model")
  println(io, "\ty = Gu")
  if s.Ts > 0
    println(io, "with Ts=", s.Ts, ".")
  elseif s.Ts == 0
    println(io, "with Ts=unspecified.")
  end
end

function showall(io::IO, s::DiscreteSisoZpk)
  show(io, s)
  println(io, "")
  printtransferfunction(io::IO, s)
end

Base.print(io::IO, s::DiscreteSisoZpk) = show(io, t)

function printtransferfunction(io::IO, s::DiscreteSisoZpk)
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

function summary(io::IO, s::DiscreteSisoZpk)
  if s.Ts > 0
    println(io, string("zpk(nu=1, ny=1, Ts=", s.Ts, "."))
  elseif s.Ts == 0
    println(io, string("zpk(nu=1, ny=1, Ts=unspecified."))
  end
end

function +{T1<:Real, T2<:Real, T3<:Real, T4<:Real}(
  s1::DiscreteSisoZpk{T1, T2}, s2::DiscreteSisoZpk{T3, T4})
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

function +{T1<:Real, T2<:Real, T3<:Real}(
    s::DiscreteSisoZpk{T1, T2}, n::T3)
  Z = s.k*poly(s.z) + n*poly(s.p)
  z = roots(Z)
  p = s.p
  k = Z[end] # Poly is now reverse order
  return zpk(z_, p_, k, s.Ts)
end
# +{T1<:Real, T2<:Real, T3<:Real}(s::DiscreteSisoZpk{T1, T2}, n::T3)  = +(s, zpk(n,s.Ts))
+{T1<:Real, T2<:Real, T3<:Real}(n::T3, s::DiscreteSisoZpk{T1, T2})  = +(s, n)

.+{T1<:Real, T2<:Real, T3<:Real}(s::DiscreteSisoZpk{T1, T2}, n::T3) = +(n, s)
.+{T1<:Real, T2<:Real, T3<:Real}(n::T3, s::DiscreteSisoZpk{T1, T2}) = +(s, n)
.+{T1<:Real, T2<:Real, T3<:Real, T4<:Real}(s1::DiscreteSisoZpk{T1, T2}, s2::DiscreteSisoZpk{T3, T4}) = +(s1, s2)

-{T1<:Real, T2<:Real}(s::DiscreteSisoZpk{T1, T2}) = zpk(s.z, s.p, -s.k, s.Ts)
function -{T1<:Real, T2<:Real, T3<:Real, T4<:Real}(
  s1::DiscreteSisoZpk{T1, T2}, s2::DiscreteSisoZpk{T3, T4})
  return +(s1, -s2)
end
-{T1<:Real, T2<:Real, T3<:Real}(s::DiscreteSisoZpk{T1, T2}, n::T3)  = +(s, -n)
-{T1<:Real, T2<:Real, T3<:Real}(n::T3, s::DiscreteSisoZpk{T1, T2})  = +(n, -s)

.-{T1<:Real, T2<:Real, T3<:Real}(s::DiscreteSisoZpk{T1, T2}, n::T3) = +(s, -n)
.-{T1<:Real, T2<:Real, T3<:Real}(n::T3, s::DiscreteSisoZpk{T1, T2}) = +(n, -s)
.-{T1<:Real, T2<:Real, T3<:Real, T4<:Real}(s1::DiscreteSisoZpk{T1, T2}, s2::DiscreteSisoZpk{T3, T4}) = +(s1, -s2)

function *{T1<:Real, T2<:Real, T3<:Real, T4<:Real}(
  s1::DiscreteSisoZpk{T1, T2}, s2::DiscreteSisoZpk{T3, T4})
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
*{T1<:Real, T2<:Real, T3<:Real}(s::DiscreteSisoZpk{T1, T2}, n::T3)  = zpk(s.z, s.p, n*s.k, s.Ts)
*{T1<:Real, T2<:Real, T3<:Real}(n::T3, s::DiscreteSisoZpk{T1, T2})  = *(s, n)

.*{T1<:Real, T2<:Real, T3<:Real}(s::DiscreteSisoZpk{T1, T2}, n::T3) = *(n, s)
.*{T1<:Real, T2<:Real, T3<:Real}(n::T3, s::DiscreteSisoZpk{T1, T2}) = *(s, n)
.*{T1<:Real, T2<:Real, T3<:Real, T4<:Real}(s1::DiscreteSisoZpk{T1, T2}, s2::DiscreteSisoZpk{T3, T4}) = *(s1, s2)

function /{T1<:Real, T2<:Real, T3<:Real}(
  n::T3, s::DiscreteSisoZpk{T1, T2})
  return zpk(s.p, s.z, n./s.k, s.Ts)
end
/{T1<:Real, T2<:Real, T3<:Real, T4<:Real}(s1::DiscreteSisoZpk{T1, T2}, s2::DiscreteSisoZpk{T3, T4}) = s1*(1/s2)
/{T1<:Real, T2<:Real, T3<:Real}(s::DiscreteSisoZpk{T1, T2}, n::T3)  = s*(1/n)

./{T1<:Real, T2<:Real, T3<:Real}(n::T3, s::DiscreteSisoZpk{T1, T2}) = n*(1/s)
./{T1<:Real, T2<:Real, T3<:Real}(s::DiscreteSisoZpk{T1, T2}, n::T3) = s*(1/n)
./{T1<:Real, T2<:Real, T3<:Real, T4<:Real}(s1::DiscreteSisoZpk{T1, T2}, s2::DiscreteSisoZpk{T3, T4}) = s1*(1/s2)

function =={T1<:Real, T2<:Real, T3<:Real, T4<:Real}(
  s1::DiscreteSisoZpk{T1, T2}, s2::DiscreteSisoZpk{T3, T4})
  fields = [:Ts, :z, :p, :k]
  for field in fields
    if getfield(s1, field) != getfield(s2, field)
      return false
    end
  end
  true
end

!={T1<:Real, T2<:Real, T3<:Real, T4<:Real}(
  s1::DiscreteSisoZpk{T1, T2}, s2::DiscreteSisoZpk{T3, T4}) = !(s1==s2)

function isapprox{T1<:Real, T2<:Real, T3<:Real, T4<:Real}(
    s1::DiscreteSisoZpk{T1, T2}, s2::DiscreteSisoZpk{T3, T4})
  # TODO: Implement
end
