immutable ContinuousSisoZpk{T1<:AbstractFloat, T2<:AbstractFloat} <: ContinuousSisoTf
  z::Vector{Complex{T1}}
  p::Vector{Complex{T1}}
  k::T2
  function call{T1, T2}(::Type{ContinuousSisoZpk}, z::Vector{Complex{T1}},
    p::Vector{Complex{T1}}, k::T2)
    new{T1, T2}(copy(z), copy(p), k)
  end
end

"""
    zpk([z, p,] k[, Ts])

Constructs a continuous transfer function with zeros `z`,
poles `p` of type `Vector` and gain `k.

A gain transfer function without zeros and poles is constructed by only
  supplying `k`.

A continuous MIMO transfer function can be constructed from a `Matrix` of `Vector`s.

Discrete transfer functions are constructed using additional argument sampling
time `Ts`.

# Examples
```julia
julia> zpk([1,0,3],[1,1,2],2)
ZPK:
    (s - 1.0)s(s - 3.0)
2.0 --------------------
    (s - 1.0)^2(s - 2.0)

Continuous-time zero-pole-gain model
```

```julia
julia> zpk(2)
ZPK:
    1.0
2.0 ---
    1.0
```

```julia
julia> zpk([1,0,3],[1,1,2],2,1)
ZPK:
    (z - 1.0)z(z - 3.0)
2.0 --------------------
    (z - 1.0)^2(z - 2.0)


Sample Time: 1 (seconds)
Discrete-time zero-pole-gain model
```
"""
function zpk{T1<:Number, T2<:Number, T3<:Real}(z::Vector{T1}, p::Vector{T2}, k::T3)
  z_, p_ = promote(z, p)
  ContinuousSisoZpk(z_, p_, k)
end

function zpk{T1<:Real}(k::T1)
  zpk(Vector{T1}(), Vector{T1}(), k)
end

Base.promote_rule{T1<:AbstractFloat, T2<:AbstractFloat, T3<:AbstractFloat, T4<:AbstractFloat}(::ContinuousSisoZpk{T1, T2}, ::ContinuousSisoZpk{T3, T4}) = ContinuousSisoZpk{promote_type(T1, T3), promote_type(T2, T4)}
Base.convert{T1<:AbstractFloat, T2<:AbstractFloat, T3<:AbstractFloat, T4<:AbstractFloat,}(::Type{ContinuousSisoZpk{T1, T2}}, t::ContinuousSisoZpk{T3, T4}) = zpk(Poly{T1}(t.z), Poly{T1}(t.p),  convert(T2, t.k))

Base.promote_rule{T1<:AbstractFloat, T2<:AbstractFloat, T3<:Real}(::ContinuousSisoZpk{T1, T2}, ::T3) = ContinuousSisoZpk{T1,T1, promote_type(T2, T3)}
Base.convert{T1<:AbstractFloat, T2<:AbstractFloat, T3<:Real}(::Type{ContinuousSisoZpk{T1, T2}}, x::T3) = zpk([one(T1)], [one(T1)], convert(T2, x))

Base.one{T1<:AbstractFloat, T2<:AbstractFloat}(::Type{ContinuousSisoZpk{T1, T2}}) = zpk(Vector{T1}(), Vector{T1}(), one(T2))
Base.one{T1<:AbstractFloat, T2<:AbstractFloat}(t::Type{ContinuousSisoZpk{T1, T2}}) = zpk(Vector{T1}(), Vector{T1}(), one(T2))
Base.zero{T1<:AbstractFloat, T2<:AbstractFloat}(::Type{ContinuousSisoZpk{T1, T2}}) = zpk(Vector{T1}(), Vector{T1}(), zero(T2))
Base.zero{T1<:AbstractFloat, T2<:AbstractFloat}(t::Type{ContinuousSisoZpk{T1, T2}}) = zpk(Vector{T1}(), Vector{T1}(), zero(T2))

function zeros{T1<:AbstractFloat, T2<:AbstractFloat}(s::ContinuousSisoZpk{T1, T2})
  return copy(s.z)
end

function poles{T1<:AbstractFloat, T2<:AbstractFloat}(s::ContinuousSisoZpk{T1, T2})
  return copy(s.p)
end

function numvec{T1<:AbstractFloat, T2<:AbstractFloat}(s::ContinuousSisoZpk{T1, T2})
  return coeffs(poly(s.z))[end:-1:1]
end

function denvec{T1<:AbstractFloat, T2<:AbstractFloat}(s::ContinuousSisoZpk{T1, T2})
  return coeffs(poly(s.p))[end:-1:1]
end

function numpoly{T1<:AbstractFloat, T2<:AbstractFloat}(s::ContinuousSisoZpk{T1, T2})
  return poly(s.z)
end

function denpoly{T1<:AbstractFloat, T2<:AbstractFloat}(s::ContinuousSisoZpk{T1, T2})
  return poly(s.p)
end

function zpkdata{T1<:AbstractFloat, T2<:AbstractFloat}(s::ContinuousSisoZpk{T1, T2})
  return copy(s.k)
end

function samplingtime{T1<:AbstractFloat, T2<:AbstractFloat}(s::ContinuousSisoZpk{T1, T2})
  return -one(Float64)
end

ndims(s::ContinuousSisoZpk)  = 1
size(s::ContinuousSisoZpk)   = 1

function getindex(s::ContinuousSisoZpk, idx::Int)
  if idx != 1
    warn("A SISO transfer function only has one element")
    throw(DomainError())
  end
  s
end

function getindex(s::ContinuousSisoZpk, rows, cols)
  if rows != 1 || cols != 1
    warn("A SISO transfer function only has one element")
    throw(DomainError())
  end
  s
end

getindex{T1<:AbstractFloat, T2<:AbstractFloat}(s::ContinuousSisoZpk{T1, T2}, ::Colon, ::Colon) = s

function getindex(s::ContinuousSisoZpk, ::Colon, cols)
  if cols != 1
    warn("A SISO transfer function only has one element")
    throw(DomainError())
  end
  s
end

function getindex(s::ContinuousSisoZpk, rows, ::Colon)
  if rows != 1
    warn("A SISO transfer function only has one element")
    throw(DomainError())
  end
  s
end

start(s::ContinuousSisoZpk)       = 1
next(s::ContinuousSisoZpk, state) = (s.m[state], state+1)
done(s::ContinuousSisoZpk, state) = state > length(s)
eltype{T1<:AbstractFloat, T2<:AbstractFloat}(::Type{ContinuousSisoZpk{T1,T2}}) = ContinuousSisoZpk{T1,T2}
length(s::ContinuousSisoZpk) = 1
eachindex(s::ContinuousSisoZpk) = 1:length(s)
endof(s::ContinuousSisoZpk) = length(s)

showcompact(io::IO, s::ContinuousSisoZpk) = print(io, summary(s))

function show(io::IO, s::ContinuousSisoZpk)
  println(io, "Discrete time zpk transfer function model")
  println(io, "\ty = Gu")
  if s.Ts > 0
    println(io, "with Ts=", s.Ts, ".")
  elseif s.Ts == 0
    println(io, "with Ts=unspecified.")
  end
end

function showall(io::IO, s::ContinuousSisoZpk)
  show(io, s)
  println(io, "")
  printtransferfunction(io::IO, s)
end

Base.print(io::IO, s::ContinuousSisoZpk) = show(io, t)

function printtransferfunction(io::IO, s::ContinuousSisoZpk)
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

function summary(io::IO, s::ContinuousSisoZpk)
  println(io, string("zpk(nu=1, ny=1)."))
end

function +{T1<:AbstractFloat, T2<:AbstractFloat, T3<:AbstractFloat, T4<:AbstractFloat}(
  s1::ContinuousSisoZpk{T1, T2}, s2::ContinuousSisoZpk{T3, T4})
  z = Array{T1}
  p = copy(z)
  k = Float64

  Z = s1.k*poly(s1.z)*poly(s2.p) + s2.k*poly(s2.z)*poly(s1.p)
  z = roots(Z)
  p = vcat(s1.p, s2.p)
  k = Z[end] # Poly is now reverse order
  return zpk(z, p, k)
end

function +{T1<:AbstractFloat, T2<:AbstractFloat, T3<:Real}(
    s::ContinuousSisoZpk{T1, T2}, n::T3)
  Z = s.k*poly(s.z) + n*poly(s.p)
  z = roots(Z)
  p = s.p
  k = Z[end] # Poly is now reverse order
  return zpk(z_, p_, k)
end
#+{T1<:AbstractFloat, T2<:AbstractFloat, T3<:Real}(t::ContinuousSisoZpk{T1, T2}, n::T3)  = +(t, zpk(n))
+{T1<:AbstractFloat, T2<:AbstractFloat, T3<:Real}(n::T3, t::ContinuousSisoZpk{T1, T2})  = +(t, n)

.+{T1<:AbstractFloat, T2<:AbstractFloat, T3<:Real}(t::ContinuousSisoZpk{T1, T2}, n::T3) = +(t, zpk(n))
.+{T1<:AbstractFloat, T2<:AbstractFloat, T3<:Real}(n::T3, t::ContinuousSisoZpk{T1, T2}) = +(t, n)
.+{T1<:AbstractFloat, T2<:AbstractFloat, T3<:AbstractFloat, T4<:AbstractFloat}(s1::ContinuousSisoZpk{T1, T2}, s2::ContinuousSisoZpk{T3, T4}) = +(s1, s2)

-{T1<:AbstractFloat, T2<:AbstractFloat}(t::ContinuousSisoZpk{T1, T2}) = zpk(t.z, t.p, -t.k)
function -{T1<:AbstractFloat, T2<:AbstractFloat, T3<:AbstractFloat, T4<:AbstractFloat}(
  s1::ContinuousSisoZpk{T1, T2}, s2::ContinuousSisoZpk{T3, T4})
  return +(s1, -s2)
end
-{T1<:AbstractFloat, T2<:AbstractFloat, T3<:Real}(t::ContinuousSisoZpk{T1, T2}, n::T3)  = +(t, -n)
-{T1<:AbstractFloat, T2<:AbstractFloat, T3<:Real}(n::T3, t::ContinuousSisoZpk{T1, T2})  = +(n, -t)

.-{T1<:AbstractFloat, T2<:AbstractFloat, T3<:Real}(t::ContinuousSisoZpk{T1, T2}, n::T3) = +(t, -n)
.-{T1<:AbstractFloat, T2<:AbstractFloat, T3<:Real}(n::T3, t::ContinuousSisoZpk{T1, T2}) = +(n, -t)
.-{T1<:AbstractFloat, T2<:AbstractFloat, T3<:AbstractFloat, T4<:AbstractFloat}(s1::ContinuousSisoZpk{T1, T2}, s2::ContinuousSisoZpk{T3, T4}) = +(s1, -s2)

function *{T1<:AbstractFloat, T2<:AbstractFloat, T3<:AbstractFloat, T4<:AbstractFloat}(
  s1::ContinuousSisoZpk{T1, T2}, s2::ContinuousSisoZpk{T3, T4})
  z = vcat(s1.z, s2.z)
  p = vcat(s1.p, s2.p)
  k = s1.k*s2.k
  return zpk(z, p, k)
end
*{T1<:AbstractFloat, T2<:AbstractFloat, T3<:Real}(t::ContinuousSisoZpk{T1, T2}, n::T3)  = zpk(t.z, t.p, n*t.k)
*{T1<:AbstractFloat, T2<:AbstractFloat, T3<:Real}(n::T3, t::ContinuousSisoZpk{T1, T2})  = *(t, n)

.*{T1<:AbstractFloat, T2<:AbstractFloat, T3<:Real}(t::ContinuousSisoZpk{T1, T2}, n::T3) = *(n, t)
.*{T1<:AbstractFloat, T2<:AbstractFloat, T3<:Real}(n::T3, t::ContinuousSisoZpk{T1, T2}) = *(t, n)
.*{T1<:AbstractFloat, T2<:AbstractFloat, T3<:AbstractFloat, T4<:AbstractFloat}(s1::ContinuousSisoZpk{T1, T2}, s2::ContinuousSisoZpk{T3, T4}) = *(s1,s2)

function /{T1<:AbstractFloat, T2<:AbstractFloat, T3<:Real}(
  n::T3, t::ContinuousSisoZpk{T1, T2})
  return zpk(t.p, t.z, n./t.k, t.Ts)
end
/{T1<:AbstractFloat, T2<:AbstractFloat, T3<:Real}(t::ContinuousSisoZpk{T1, T2}, n::T3)  = t*(1/n)
/{T1<:AbstractFloat, T2<:AbstractFloat, T3<:AbstractFloat, T4<:AbstractFloat}(s1::ContinuousSisoZpk{T1, T2}, s2::ContinuousSisoZpk{T3, T4}) = s1*(1/s2)
./{T1<:AbstractFloat, T2<:AbstractFloat, T3<:Real}(n::T3, t::ContinuousSisoZpk{T1, T2}) = t*(1/n)
./{T1<:AbstractFloat, T2<:AbstractFloat, T3<:Real}(t::ContinuousSisoZpk{T1, T2}, n::T3) = t*(1/n)
./{T1<:AbstractFloat, T2<:AbstractFloat, T3<:AbstractFloat, T4<:AbstractFloat}(s1::ContinuousSisoZpk{T1, T2}, s2::ContinuousSisoZpk{T3, T4}) = s1*(1/s2)

function =={T1<:AbstractFloat, T2<:AbstractFloat, T3<:AbstractFloat, T4<:AbstractFloat}(
  s1::ContinuousSisoZpk{T1, T2}, s2::ContinuousSisoZpk{T3, T4})
  fields = [:z, :p, :k]
  for field in fields
    if getfield(s1, field) != getfield(s2, field)
      return false
    end
  end
  true
end

!={T1<:AbstractFloat, T2<:AbstractFloat, T3<:AbstractFloat, T4<:AbstractFloat}(
  s1::ContinuousSisoZpk{T1, T2}, s2::ContinuousSisoZpk{T3, T4}) = !(s1==s2)

function isapprox{T1<:AbstractFloat, T2<:AbstractFloat, T3<:AbstractFloat, T4<:AbstractFloat}(
    s1::ContinuousSisoZpk{T1, T2}, s2::ContinuousSisoZpk{T3, T4})
  # TODO: Implement
end
