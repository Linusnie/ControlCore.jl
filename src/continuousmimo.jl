immutable ContinuousMimo{T1<:ContinuousSisoTf} <: TransferFunction
  m::Matrix{T1}
  ny::Int
  nu::Int
  function call{T1}(::Type{ContinuousMimo}, m::Matrix{T1})
    return new{T1}(m, size(m, 1), size(m, 2))
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
function tf{T1<:Real, T2<:Real}(num::Matrix{Vector{T1}}, den::Matrix{Vector{T2}})
  # Validate input and output dimensions match
  ny, nu = size(num, 1, 2)
  if (ny, nu) != size(den, 1, 2)
    warn("num and den dimensions must match")
    throw(DomainError())
  end
  num_, den_ = promote(num, den)
  m = Array(ContinuousSisoRational{eltype(num_)}, ny, nu)
  for idx in eachindex(s.m)
      m[idx] = tf(num_[idx], den_[idx])
  end
  return ContinuousMimo(m)
end

function tf{T1<:Real}(gain::Matrix{T1})
  ny, nu = size(gain, 1, 2)
  m = Array(ContinuousSisoZpk{T1, T1}, ny, nu)
  for idx in eachindex(s.m)
    m[idx] = zpk(gain[idx])
  end
  return ContinuousMimo(m)
end

function zpk{T1<:Number, T2<:Number, T3<:Real}(z::Matrix{Vector{T1}}, p::Matrix{Vector{T2}}, k::Matrix{T3})
  # Validate input and output dimensions match
  ny, nu = size(z, 1, 2)
  if (ny, nu) != size(p, 1, 2) || (ny, nu) != size(k, 1, 2)
    warn("num and den dimensions must match")
    throw(DomainError())
  end
  T  = promote_type(real(T1), real(T2))
  z_ = convert(Matrix{Vector{Complex{Float64}}},z)
  p_ = convert(Matrix{Vector{Complex{Float64}}},p)
  m = Array(ContinuousSisoZpk{Float64, Float64}, ny, nu)
  for idx in eachindex(m)
    m[idx] = zpk(z_[idx], p_[idx], Float64(k[idx]))
  end
  ContinuousMimo(m, Ts)
end

function zpk{T1<:Real}(gain::Matrix{T1})
  ny, nu = size(gain, 1, 2)
  m = Matrix{ContinuousSisoZpk{Float64, Float64}}(ny, nu)
  for idx in eachindex(m)
    m[idx] = zpk(gain[idx])
  end
  ContinuousMimo(m, Ts)
end

function tf{T1<:ContinuousSisoTf}(m::Matrix{T1})
  return ContinuousMimo(m)
end

one{T1<:ContinuousSisoTf}(s::ContinuousMimo{T1}) = ContinuousMimo(fill(one(T1),1,1))

function one{T1<:ContinuousSisoTf}(::Type{ContinuousMimo{T1}})
  s = one(T1)
  return ContinuousMimo(fill(s,1,1))
end
function zero{T1<:ContinuousSisoTf}(s::ContinuousMimo{T1})
  ContinuousMimo(fill(zero(T1),1,1))
end
function zero{T1<:ContinuousSisoTf}(::Type{ContinuousMimo{T1}})
  s = zero(T1)
  return ContinuousMimo(fill(s,1,1))
end

function zeros{T1<:ContinuousSisoTf}(s::ContinuousMimo{T1})
  z = Matrix{Vector}(s.ny, s.nu)
  for idx in eachindex(z)
    z[idx] = poles(t[idx])
  end
  return z
end

function poles{T1<:ContinuousSisoTf}(s::ContinuousMimo{T1})
  p = Matrix{Vector}(s.ny, s.nu)
  for idx in eachindex(p)
    p[idx] = poles(s.m[idx])
  end
  return p
end

function numvec{T1<:ContinuousSisoTf}(s::ContinuousMimo{T1})
  num = Matrix{Vector}(s.ny, s.nu)
  for idx in eachindex(s)
    num[idx] = numvec(s.m[idx])
  end
  return num
end

function denvec{T1<:ContinuousSisoTf}(s::ContinuousMimo{T1})
  den = Matrix{Vector}(s.ny, s.nu)
  for idx in eachindex(s)
    den[idx] = denvec(s.m[idx])
  end
  return den
end

function numpoly{T1<:ContinuousSisoTf}(s::ContinuousMimo{T1})
  num = Matrix{Vector}(s.ny, s.nu)
  for idx in eachindex(s)
    num[idx] = numpoly(s.m[idx])
  end
  return num
end

function denpoly{T1<:ContinuousSisoTf}(s::ContinuousMimo{T1})
  den = Matrix{Vector}(s.ny, s.nu)
  for idx in eachindex(s)
    den[idx] = denpoly(s.m[idx])
  end
  return den
end

function zpkdata{T1<:ContinuousSisoTf}(s::ContinuousMimo{T1})
  zpkdata_ = Matrix{Vector}(s.ny, s.nu)
  for idx in eachindex(s)
    zpkdata_[idx] = zpkdata(s.m[idx])
  end
  return zpkdata_
end

samplingtime{T1<:ContinuousSisoTf}(s::ContinuousMimo{T1}) = -one(Float64)

ndims(s::ContinuousMimo)  = 2
size(s::ContinuousMimo)   = (s.ny, s.nu)

function getindex(s::ContinuousMimo, idx::Int)
  row, col = divrem(idx-1, s.ny)
  ContinuousMimo(fill(s.m[row+1, col+1],1,1), s.Ts)
end

function getindex(s::ContinuousMimo, rows, cols)
  s2 = try
      [s.m[row, col] for row in rows, col in cols]
    catch exception
      warn("s[,j]: Index out of bounds")
      throw(exception)
    end
  ContinuousMimo(s2, s.Ts)
end

getindex(s::ContinuousMimo, ::Colon, ::Colon) = ContinuousMimo(s.m, s.Ts) # returns a copy of s

function getindex(s::ContinuousMimo, ::Colon, cols)
  s2 = try
      [s.m[row, col] for row in 1:s.ny, col in cols]
    catch exception
      warn("s[,j]: Index out of bounds")
      throw(exception)
    end
  ContinuousMimo(s2, s.Ts)
end

function getindex(s::ContinuousMimo, rows, ::Colon)
  s2 = try
      [s.m[row, col] for row in rows, col in 1:s.nu]
    catch exception
      warn("s[,j]: Index out of bounds")
      throw(exception)
    end
  ContinuousMimo(s2, s.Ts)
end

start(s::ContinuousMimo)       = 1
next(s::ContinuousMimo, state) = (s.m[state], state+1)
done(s::ContinuousMimo, state) = state > length(s)
eltype{T1<:ContinuousSisoTf}(s::ContinuousMimo{T1}) = T1
length(s::ContinuousMimo) = length(s.m)
eachindex(s::ContinuousMimo) = 1:length(s)
endof(s::ContinuousMimo) = length(s)

showcompact(io::IO, s::ContinuousMimo) = print(io, summary(s))

function show(io::IO, s::ContinuousMimo)
  println(io, "Continuous time transfer function model")
  println(io, "\ty = Gu")
  println(io, "with nu=", s.nu, ", ny=", s.ny, ".")
end

function showall{T1<:ContinuousSisoTf}(io::IO, s::ContinuousMimo{T1})
  show(io, s)
  println(io, "")
  for subs in s
    printtransferfunction(io::IO, subs)
    println(io, "")
  end
end

function summary{T1<:Real}(s::ContinuousMimo{ContinuousSisoRational{T1}})
  string("tf(nu=", s.nu, ", ny=", s.ny, ")")
end
function summary{T1<:Number, T2<:Real}(s::ContinuousMimo{ContinuousSisoZpk{T1, T2}})
  string("zpk(nu=", s.nu, ", ny=", s.ny, ")")
end

function +{T1<:ContinuousSisoTf, T2<:ContinuousSisoTf}(
    s1::ContinuousMimo{T1}, s2::ContinuousMimo{T2})
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

+{T1<:ContinuousSisoTf} (s::ContinuousMimo{T1}) = tf(s.m+n)
+{T1<:ContinuousSisoTf, T2<:Real}(n::T2, s::ContinuousMimo{T1}) = +(s, n)

function +{T1<:ContinuousSisoTf, T2<:Real}(s::ContinuousMimo{T1}, n::Matrix{T2})
  if size(s.m) != size(n)
    warn("Systems have different shapes")
    throw(DomainError)
  end
  ContinuousMimo(s.m + n)
end
+{T1<:ContinuousSisoTf, T2<:Real}(n::Matrix{T2}, s::ContinuousMimo{T1})  = (s, n)

.+{T1<:ContinuousSisoTf, T2<:Real}(n::T2, s::ContinuousMimo{T1}) = +(s, n)
.+{T1<:ContinuousSisoTf, T2<:Real}(s::ContinuousMimo{T1}, n::T2) = +(s, n)

.+{T1<:ContinuousSisoTf, T2<:Real}(n::Matrix{T2}, s::ContinuousMimo{T1}) = +(n, s)
.+{T1<:ContinuousSisoTf, T2<:Real}(s::ContinuousMimo{T1}, n::Matrix{T2}) = +(s, n)
.+{T1<:ContinuousSisoTf, T2<:ContinuousSisoTf}(s1::ContinuousMimo{T1}, s2::ContinuousMimo{T2}) = +(s1, s2)

-{T1<:ContinuousSisoTf}(s::ContinuousMimo{T1}) = tf(-s.m)

-{T1<:ContinuousSisoTf, T2<:Real}(s::ContinuousMimo{T1}, n::Matrix{T2}) = +(s, -n)
-{T1<:ContinuousSisoTf, T2<:Real}(n::Matrix{T2}, s::ContinuousMimo{T1}) = +(-n, s)
-{T1<:ContinuousSisoTf, T2<:ContinuousSisoTf}(s1::ContinuousMimo{T1}, s2::ContinuousMimo{T2}) = +(s1,-s2)
-{T1<:ContinuousSisoTf, T2<:Real}(s::ContinuousMimo{T1}, n::T2) = +(s,-n)
-{T1<:ContinuousSisoTf, T2<:Real}(n::T2, s::ContinuousMimo{T1}) = +(-n, s)

.-{T1<:ContinuousSisoTf, T2<:Real}(s::ContinuousMimo{T1}, n::T2) = +(s,-n)
.-{T1<:ContinuousSisoTf, T2<:Real}(n::T2, s::ContinuousMimo{T1}) = +(-n, s)

.-{T1<:ContinuousSisoTf, T2<:Real}(n::Matrix{T2}, s::ContinuousMimo{T1}) = -(n, s)
.-{T1<:ContinuousSisoTf, T2<:Real}(s::ContinuousMimo{T1}, n::Matrix{T2}) = -(s, n)

function *{T1<:ContinuousSisoTf, T2<:ContinuousSisoTf}(
  s1::ContinuousMimo{T1}, s2::ContinuousMimo{T2})
  if size(s1.m, 2) != size(s2.m, 1)
    warn("s1*s2: s1 must have same number of inputs as s2 has outputs")
    throw(DomainError())
  end
  m = s1.m*s2.m
  return tf(m)
end

function *{T1<:ContinuousSisoTf, T2<:Real}(
  s::ContinuousMimo{T1}, n::T2)
  tf(s.m*n)
end
*{T1<:ContinuousSisoTf, T2<:Real}(n::T2, s::ContinuousMimo{T1})  = *(s, n)
*{T1<:ContinuousSisoTf, T2<:Real}(s::ContinuousMimo{T1}, n::Matrix{T2}) = tf(s.m*n)
*{T1<:ContinuousSisoTf, T2<:Real}(n::Matrix{T2}, s::ContinuousMimo{T1}) = tf(n*s.m)

.*{T1<:ContinuousSisoTf, T2<:Real}(n::T2, s::ContinuousMimo{T1})         = *(n, s)
.*{T1<:ContinuousSisoTf, T2<:Real}(s::ContinuousMimo{T1}, n::T2)         = *(s, n)
.*{T1<:ContinuousSisoTf, T2<:Real}(n::Matrix{T2}, s::ContinuousMimo{T1}) = *(n, s)
.*{T1<:ContinuousSisoTf, T2<:Real}(s::ContinuousMimo{T1}, n::Matrix{T2}) = *(s, n)
.*{T1<:ContinuousSisoTf, T2<:ContinuousSisoTf}(s1::ContinuousMimo{T1}, s2::ContinuousMimo{T2}) = *(s1, s2)

function /{T1<:ContinuousSisoTf, T2<:Real}(n::T2, s::ContinuousMimo{T1})
  warn("MIMO TransferFunction inversion isn't implemented yet")
  throw(DomainError())
end
/{T1<:ContinuousSisoTf, T2<:ContinuousSisoTf}(s1::ContinuousMimo{T1}, s2::ContinuousMimo{T2}) = s1*(1/s2)
/{T1<:ContinuousSisoTf, T2<:Real}(s::ContinuousMimo{T1}, n::T2) = tf(s.m/n)

function /{T1<:ContinuousSisoTf, T2<:Real}(n::Matrix{T2}, s::ContinuousMimo{T1})
  warn("MIMO TransferFunction inversion isn't implemented yet")
  throw(DomainError())
end
/{T1<:ContinuousSisoTf, T2<:Real}(s::ContinuousMimo{T1}, n::Matrix{T2}) = s*(1/n)

./{T1<:ContinuousSisoTf, T2<:Real}(s::ContinuousMimo{T1}, n::T2) = s/n
./{T1<:ContinuousSisoTf, T2<:Real}(n::T2, s::ContinuousMimo{T1}) = n/s
./{T1<:ContinuousSisoTf, T2<:Real}(s::ContinuousMimo{T1}, n::Matrix{T2}) = s/n
./{T1<:ContinuousSisoTf, T2<:Real}(n::Matrix{T2}, s::ContinuousMimo{T1}) = n/s

function =={T1<:ContinuousSisoTf, T2<:ContinuousSisoTf}(
  s1::ContinuousMimo{T1}, s2::ContinuousMimo{T2})
  return s1.m != s2.m
end

!={T1<:ContinuousSisoTf, T2<:ContinuousSisoTf}(
  s1::ContinuousMimo{T1}, s2::ContinuousMimo{T2}) = !(s1.m == s2.m)

function isapprox{T1<:ContinuousSisoTf, T2<:ContinuousSisoTf}(
  s1::ContinuousMimo{T1}, s2::ContinuousMimo{T2})
  # TODO: Implement
end
