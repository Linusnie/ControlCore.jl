# Defined collection of similar siso systems
CSiso = Union{CSisoTf}

immutable CMimo{T<:CSiso} <: MimoSystem
  m::Matrix{T}
  ny::Int
  nu::Int
  function call{T}(::Type{CMimo}, m::Matrix{T})
    return new{T}(m, size(m, 1), size(m, 2))
  end
end


"""
    tf(num, den[, Ts])

Constructs a continuous transfer function with numerator `num`,
denominator `den` of type `Vector`.

A continuous transfer function can also be constructed from `Polynomials.Poly` objects.

A continuous MIMO transfer function can be constructed from a `Matrix` of `Vector` or
`Matrix` of `Polynomials.Poly` objects.

Discrete transfer functions are constructed using additional argument sampling
time `Ts`.

# Examples
```julia
julia> tf([1,0,3],[1,1,2])
  s^2 + 3
-----------
s^2 + s + 2
```

```julia
julia> tf(Polynomials.Poly([3,0,1]),Polynomials.Poly([2,1,1]))
  s^2 + 3
-----------
s^2 + s + 2
```

```julia
julia> tf([1,0,3],[1,1,2],1)
   z^-2 + 3
---------------
z^-2 + z^-1 + 2
Sampling time = 1
```
"""
# ```julia
# julia> tf(Polynomials.Poly([3,0,1]),Polynomials.Poly([2,1,1]),1)
#    z^-2 + 3
# ---------------
# z^-2 + z^-1 + 2
# Sampling time = 1
# ```
function tf{T1<:AbstractFloat}(num::Matrix{Vector{T1}}, den::Matrix{Vector{T1}})
  ny, nu = size(num, 1, 2)
  if (ny, nu) != size(den, 1, 2)
    warn("num and den dimensions must match")
    throw(DomainError())
  end
  m = Array(CSisoRational{T1}, ny, nu)
  for idx in eachindex(m)
    m[idx] = tf(num_[idx], den_[idx])
  end
  DMimoTf(m)
end

function tf{T1<:AbstractFloat}(gain::Matrix{T1})
  ny, nu = size(gain, 1, 2)
  m = Matrix{CSisoRational{T1}}(ny, nu)
  for idx in eachindex(m)
    m[idx] = tf(gain[idx])
  end
  DMimoTf(m)
end

function zpk{T1<:AbstractFloat}(z::Matrix{Vector{Complex{T1}}}, p::Matrix{Vector{Complex{T1}}}, k::Matrix{T1})
  ny, nu = size(z, 1, 2)
  if (ny, nu) != size(p, 1, 2) || (ny, nu) != size(k, 1, 2)
    warn("num and den dimensions must match")
    throw(DomainError())
  end
  m = Array(CSisoZpk{T1}, ny, nu)
  for idx in eachindex(m)
    m[idx] = zpk(z_[idx], p_[idx], k[idx])
  end
  DMimoTf(m)
end

function zpk{T1<:AbstractFloat}(gain::Matrix{T1})
  ny, nu = size(gain, 1, 2)
  m = Matrix{CSisoZpk{T1}}(ny, nu)
  for idx in eachindex(m)
    m[idx] = zpk(gain[idx], Ts)
  end
  DMimoTf(m)
end

function tf{T1<:Real, T2<:Real}(num::Matrix{Vector{T1}}, den::Matrix{Vector{T2}})
  ny, nu = size(num, 1, 2)
  if (ny, nu) != size(den, 1, 2)
    warn("num and den dimensions must match")
    throw(DomainError())
  end
  T = promote_type(T1, T2, Float16)
  num_ = convert(Matrix{Vector{T}}, num)
  den_ = convert(Matrix{Vector{T}}, den)
  m = Array(CSisoRational{T}, ny, nu)
  for idx in eachindex(s.m)
      m[idx] = tf(num_[idx], den_[idx])
  end
  return CMimo(m)
end

function tf{T1<:Real}(gain::Matrix{T1})
  ny, nu = size(gain, 1, 2)
  T = promote_type(T1, Float16)
  gain_ = convert(Matrix{T}, gain)
  m = Array(CSisoZpk{T1, T1}, ny, nu)
  for idx in eachindex(s.m)
    m[idx] = zpk(gain_[idx])
  end
  return CMimo(m)
end

function zpk{T1<:Number, T2<:Number, T3<:Real}(z::Matrix{Vector{T1}}, p::Matrix{Vector{T2}}, k::Matrix{T3})
  ny, nu = size(z, 1, 2)
  if (ny, nu) != size(p, 1, 2) || (ny, nu) != size(k, 1, 2)
    warn("num and den dimensions must match")
    throw(DomainError())
  end
  T  = promote_type(real(T1), real(T2), k, Float16)
  z_ = convert(Matrix{Vector{Complex{T}}}, z)
  p_ = convert(Matrix{Vector{Complex{T}}}, p)
  k_ = convert(Matrix{T}, k_)
  m = Array(CSisoZpk{T}, ny, nu)
  for idx in eachindex(m)
    m[idx] = zpk(z_[idx], p_[idx], k_[idx])
  end
  CMimo(m)
end

function zpk{T1<:Real}(gain::Matrix{T1})
  ny, nu = size(gain, 1, 2)
  T  = promote_type(T1, Float16)
  gain_ = convert(Matrix{T}, gain)
  m = Matrix{CSisoZpk{T}}(ny, nu)
  for idx in eachindex(m)
    m[idx] = zpk(gain_[idx])
  end
  CMimo(m, Ts)
end

function tf{T1<:CSisoTf}(m::Matrix{T1})
  return CMimo(m)
end

one{T1}(s::CMimo{T1})       = CMimo(ones(T1,1,1))
one{T1}(::Type{CMimo{T1}})  = CMimo(ones(T1,1,1))
zero{T1}(s::CMimo{T1})      = CMimo(zeros(T1,1,1))
zero{T1}(::Type{CMimo{T1}}) = CMimo(zeros(T1,1,1))

function zeros{T1}(s::CMimo{T1})
  z = Matrix{Vector}(s.ny, s.nu)
  for idx in eachindex(z)
    z[idx] = poles(t[idx])
  end
  return z
end

function poles{T1}(s::CMimo{T1})
  p = Matrix{Vector}(s.ny, s.nu)
  for idx in eachindex(p)
    p[idx] = poles(s.m[idx])
  end
  return p
end

function numvec{T1}(s::CMimo{T1})
  num = Matrix{Vector}(s.ny, s.nu)
  for idx in eachindex(s)
    num[idx] = numvec(s.m[idx])
  end
  return num
end

function denvec{T1}(s::CMimo{T1})
  den = Matrix{Vector}(s.ny, s.nu)
  for idx in eachindex(s)
    den[idx] = denvec(s.m[idx])
  end
  return den
end

function numpoly{T1}(s::CMimo{T1})
  num = Matrix{Vector}(s.ny, s.nu)
  for idx in eachindex(s)
    num[idx] = numpoly(s.m[idx])
  end
  return num
end

function denpoly{T1}(s::CMimo{T1})
  den = Matrix{Vector}(s.ny, s.nu)
  for idx in eachindex(s)
    den[idx] = denpoly(s.m[idx])
  end
  return den
end

function zpkdata{T1}(s::CMimo{T1})
  zpkdata_ = Matrix{Vector}(s.ny, s.nu)
  for idx in eachindex(s)
    zpkdata_[idx] = zpkdata(s.m[idx])
  end
  return zpkdata_
end

samplingtime{T1}(s::CMimo{T1}) = -one(Float64)

ndims(s::CMimo)  = 2
size(s::CMimo)   = (s.ny, s.nu)

function getindex(s::CMimo, idx::Int)
  row, col = divrem(idx-1, s.ny)
  CMimo(fill(s.m[row+1, col+1],1,1), s.Ts)
end

function getindex(s::CMimo, rows, cols)
  s2 = try
      [s.m[row, col] for row in rows, col in cols]
    catch exception
      warn("s[,j]: Index out of bounds")
      throw(exception)
    end
  CMimo(s2, s.Ts)
end

getindex(s::CMimo, ::Colon, ::Colon) = CMimo(s.m, s.Ts) # returns a copy of s

function getindex(s::CMimo, ::Colon, cols)
  s2 = try
      [s.m[row, col] for row in 1:s.ny, col in cols]
    catch exception
      warn("s[,j]: Index out of bounds")
      throw(exception)
    end
  CMimo(s2, s.Ts)
end

function getindex(s::CMimo, rows, ::Colon)
  s2 = try
      [s.m[row, col] for row in rows, col in 1:s.nu]
    catch exception
      warn("s[,j]: Index out of bounds")
      throw(exception)
    end
  CMimo(s2, s.Ts)
end

start(s::CMimo)       = 1
next(s::CMimo, state::Int) = (s.m[state], state+1)
done(s::CMimo, state::Int) = state > length(s)
eltype{T<:CSiso}(s::CMimo{T}) = T
length(s::CMimo) = length(s.m)
eachindex(s::CMimo) = 1:length(s)
endof(s::CMimo) = length(s)

showcompact(io::IO, s::CMimo) = print(io, summary(s))

function show(io::IO, s::CMimo)
  println(io, "Continuous time transfer function model")
  println(io, "\ty = Gu")
  println(io, "with nu=", s.nu, ", ny=", s.ny, ".")
end

function showall{T1}(io::IO, s::CMimo{T1})
  show(io, s)
  println(io, "")
  for subs in s
    printtransferfunction(io::IO, subs)
    println(io, "")
  end
end

function summary{T1<:Real}(s::CMimo{CSisoRational{T1}})
  string("tf(nu=", s.nu, ", ny=", s.ny, ")")
end
function summary{T1<:Number}(s::CMimo{CSisoZpk{T1}})
  string("zpk(nu=", s.nu, ", ny=", s.ny, ")")
end

function +{T1, T2}(
    s1::CMimo{T1}, s2::CMimo{T2})
  if size(s1.m) != size(s2.m)
    warn("Systems have different shapes")
    throw(DomainError)
  elseif s1.Ts != s2.Ts
    warn("Sampling time mismatch")
    throw(DomainError)
  end
  s1, s2 = promote(s1, s2)
  m = s1.m + s2.m
  return tf(m)
end

+{T1, T2<:Real}(s::CMimo{T1}, n::T2) = tf(s.m+n)
+{T1, T2<:Real}(n::T2, s::CMimo{T1}) = +(s, n)

function +{T1, T2<:Real}(s::CMimo{T1}, n::Matrix{T2})
  if size(s.m) != size(n)
    warn("Systems have different shapes")
    throw(DomainError)
  end
  CMimo(s.m + n)
end
+{T1, T2<:Real}(n::Matrix{T2}, s::CMimo{T1})  = (s, n)

.+{T1, T2<:Real}(n::T2, s::CMimo{T1}) = +(s, n)
.+{T1, T2<:Real}(s::CMimo{T1}, n::T2) = +(s, n)

.+{T1, T2<:Real}(n::Matrix{T2}, s::CMimo{T1}) = +(n, s)
.+{T1, T2<:Real}(s::CMimo{T1}, n::Matrix{T2}) = +(s, n)
.+{T1, T2}(s1::CMimo{T1}, s2::CMimo{T2}) = +(s1, s2)

-{T1}(s::CMimo{T1}) = tf(-s.m)

-{T1, T2<:Real}(s::CMimo{T1}, n::Matrix{T2}) = +(s, -n)
-{T1, T2<:Real}(n::Matrix{T2}, s::CMimo{T1}) = +(-n, s)
-{T1, T2}(s1::CMimo{T1}, s2::CMimo{T2}) = +(s1,-s2)
-{T1, T2<:Real}(s::CMimo{T1}, n::T2) = +(s,-n)
-{T1, T2<:Real}(n::T2, s::CMimo{T1}) = +(-n, s)

.-{T1, T2<:Real}(s::CMimo{T1}, n::T2) = +(s,-n)
.-{T1, T2<:Real}(n::T2, s::CMimo{T1}) = +(-n, s)

.-{T1, T2<:Real}(n::Matrix{T2}, s::CMimo{T1}) = -(n, s)
.-{T1, T2<:Real}(s::CMimo{T1}, n::Matrix{T2}) = -(s, n)

function *{T1, T2}(
  s1::CMimo{T1}, s2::CMimo{T2})
  if size(s1.m, 2) != size(s2.m, 1)
    warn("s1*s2: s1 must have same number of inputs as s2 has outputs")
    throw(DomainError())
  end
  m = s1.m*s2.m
  return tf(m)
end

function *{T1, T2<:Real}(
  s::CMimo{T1}, n::T2)
  tf(s.m*n)
end
*{T1, T2<:Real}(n::T2, s::CMimo{T1})  = *(s, n)
*{T1, T2<:Real}(s::CMimo{T1}, n::Matrix{T2}) = tf(s.m*n)
*{T1, T2<:Real}(n::Matrix{T2}, s::CMimo{T1}) = tf(n*s.m)

.*{T1, T2<:Real}(n::T2, s::CMimo{T1})         = *(n, s)
.*{T1, T2<:Real}(s::CMimo{T1}, n::T2)         = *(s, n)
.*{T1, T2<:Real}(n::Matrix{T2}, s::CMimo{T1}) = *(n, s)
.*{T1, T2<:Real}(s::CMimo{T1}, n::Matrix{T2}) = *(s, n)
.*{T1, T2}(s1::CMimo{T1}, s2::CMimo{T2}) = *(s1, s2)

function /{T1, T2<:Real}(n::T2, s::CMimo{T1})
  warn("MIMO TransferFunction inversion isn't implemented yet")
  throw(DomainError())
end
/{T1, T2}(s1::CMimo{T1}, s2::CMimo{T2}) = s1*(1/s2)
/{T1, T2<:Real}(s::CMimo{T1}, n::T2) = tf(s.m/n)

function /{T1, T2<:Real}(n::Matrix{T2}, s::CMimo{T1})
  warn("MIMO TransferFunction inversion isn't implemented yet")
  throw(DomainError())
end
/{T1, T2<:Real}(s::CMimo{T1}, n::Matrix{T2}) = s*(1/n)

./{T1, T2<:Real}(s::CMimo{T1}, n::T2) = s/n
./{T1, T2<:Real}(n::T2, s::CMimo{T1}) = n/s
./{T1, T2<:Real}(s::CMimo{T1}, n::Matrix{T2}) = s/n
./{T1, T2<:Real}(n::Matrix{T2}, s::CMimo{T1}) = n/s

function =={T1, T2}(
  s1::CMimo{T1}, s2::CMimo{T2})
  return s1.m != s2.m
end

!={T1, T2}(
  s1::CMimo{T1}, s2::CMimo{T2}) = !(s1.m == s2.m)

function isapprox{T1, T2}(
  s1::CMimo{T1}, s2::CMimo{T2})
  # TODO: Implement
end
