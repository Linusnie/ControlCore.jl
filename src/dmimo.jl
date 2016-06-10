# Defined collection of similar siso systems
DSiso = Union{DSisoTf}

immutable DMimo{T<:DSiso} <: MimoSystem{T}
  m::Matrix{T}
  ny::Int
  nu::Int
  function call{T}(::Type{DMimo}, m::Matrix{T})
    return new{T}(m, size(m,1), size(m,2))
  end
end

function tf{T1<:AbstractFloat}(num::Matrix{Vector{T1}}, den::Matrix{Vector{T1}}, Ts::T1)
  ny, nu = size(num, 1, 2)
  if (ny, nu) != size(den, 1, 2)
    warn("num and den dimensions must match")
    throw(DomainError())
  end
  m = Array(DSisoRational{T1}, ny, nu)
  for idx in eachindex(m)
    m[idx] = tf(num_[idx], den_[idx], Ts)
  end
  DMimo(m)
end

function tf{T1<:AbstractFloat}(gain::Matrix{T1}, Ts::T1)
  ny, nu = size(gain, 1, 2)
  m = Matrix{DSisoRational{T1}}(ny, nu)
  for idx in eachindex(m)
    m[idx] = tf(gain[idx], Ts)
  end
  DMimo(m)
end

function zpk{T1<:AbstractFloat}(z::Matrix{Vector{Complex{T1}}}, p::Matrix{Vector{Complex{T1}}}, k::Matrix{T1}, Ts::T1)
  # Validate input and output dimensions match
  ny, nu = size(z, 1, 2)
  if (ny, nu) != size(p, 1, 2) || (ny, nu) != size(k, 1, 2)
    warn("num and den dimensions must match")
    throw(DomainError())
  end
  m = Array(DSisoZpk{T1}, ny, nu)
  for idx in eachindex(m)
    m[idx] = zpk(z_[idx], p_[idx], k[idx], Ts)
  end
  DMimo(m)
end

function zpk{T1<:AbstractFloat}(gain::Matrix{T1}, Ts::T1)
  ny, nu = size(gain, 1, 2)
  m = Matrix{DSisoZpk{T1}}(ny, nu)
  for idx in eachindex(m)
    m[idx] = zpk(gain[idx], Ts)
  end
  DMimo(m, Ts)
end

function tf{T1<:Real, T2<:Real, T3<:Real}(num::Matrix{Vector{T1}},
  den::Matrix{Vector{T2}}, Ts::T3)
  ny, nu = size(num, 1, 2)
  if (ny, nu) != size(den, 1, 2)
    warn("num and den dimensions must match")
    throw(DomainError())
  end
  T = promote_type(T1, T2, T3, Float16)
  num_ = convert(Matrix{Vector{T}}, num)
  den_ = convert(Matrix{Vector{T}}, den)
  Ts_ = convert(T, Ts)
  m = Array(DSisoRational{T}, ny, nu)
  for idx in eachindex(m)
    m[idx] = tf(num_[idx], den_[idx], Ts_)
  end
  DMimo(m)
end

function tf{T1<:Real, T2<:Real}(gain::Matrix{T1}, Ts::T2)
  ny, nu = size(gain, 1, 2)
  T = promote_type(T1, T2, Float16)
  gain_ = convert(Matrix{T}, gain)
  Ts_ = convert(T, Ts)
  m = Matrix{DSisoRational{T}}(ny, nu)
  for idx in eachindex(m)
    m[idx] = tf(gain_[idx], Ts_)
  end
  DMimo(m)
end

function zpk{T1<:Number, T2<:Number, T3<:Real, T4<:Real}(z::Matrix{Vector{T1}},
  p::Matrix{Vector{T2}}, k::Matrix{T3}, Ts::T4)
  ny, nu = size(z, 1, 2)
  if (ny, nu) != size(p, 1, 2) || (ny, nu) != size(k, 1, 2)
    warn("num and den dimensions must match")
    throw(DomainError())
  end
  T  = promote_type(real(T1), real(T2), T3, T4, Float16)
  z_ = convert(Matrix{Vector{Complex{T}}}, z)
  p_ = convert(Matrix{Vector{Complex{T}}}, p)
  k_ = convert(Matrix{T}, k)
  Ts_ = convert(T, Ts)
  m = Array(DSisoZpk{T}, ny, nu)
  for idx in eachindex(m)
    m[idx] = zpk(z_[idx], p_[idx], k_[idx], Ts_)
  end
  DMimo(m)
end

function zpk{T1<:Real, T2<:Real}(gain::Matrix{T1}, Ts::T2)
  ny, nu = size(gain, 1, 2)
  T = promote_type(T1, T2, Float16)
  gain_ = convert(Matrix{T}, gain)
  Ts_ = convert(T, Ts)
  m = Matrix{DSisoZpk{T}}(ny, nu)
  for idx in eachindex(m)
    m[idx] = zpk(gain_[idx], Ts_)
  end
  DMimo(m)
end

function tf{T1<:DSiso}(m::Matrix{T1})
  Ts = m[1,1].Ts
  for idx in eachindex(m)
    if Ts != m[idx].Ts
      warn("Sampling time mismatch")
      throw(DomainError())
    end
  end
  DMimo(m)
end

one{T1<:DSiso}(s::DMimo{T1})      = DMimo(fill(one(s), 1, 1))
one{T1<:DSiso}(::Type{DMimo{T1}}) = DMimo(fill(one(T1), 1, 1))
zero{T1<:DSiso}(s::DMimo{T1})     = DMimo(fill(zero(s),1,1))
zero{T1}(::Type{DMimo{T1}})       = DMimo(fill(zero(T1), 1, 1))

function zeros{T1}(s::DMimo{T1})
  z = Matrix{Vector}(s.ny, s.nu)
  for idx in eachindex(z)
    z[idx] = poles(t[idx])
  end
  return z
end

function poles{T1}(s::DMimo{T1})
  p = Matrix{Vector}(s.ny, s.nu)
  for idx in eachindex(p)
    p[idx] = poles(s.m[idx])
  end
  return p
end

function numvec{T1}(s::DMimo{T1})
  num = Matrix{Vector}(s.ny, s.nu)
  for idx in eachindex(s)
    num[idx] = numvec(s.m[idx])
  end
  return num
end

function denvec{T1}(s::DMimo{T1})
  den = Matrix{Vector}(s.ny, s.nu)
  for idx in eachindex(s)
    den[idx] = denvec(s.m[idx])
  end
  return den
end

function numpoly{T1}(s::DMimo{T1})
  num = Matrix{Vector}(s.ny, s.nu)
  for idx in eachindex(s)
    num[idx] = numpoly(s.m[idx])
  end
  return num
end

function denpoly{T1}(s::DMimo{T1})
  den = Matrix{Vector}(s.ny, s.nu)
  for idx in eachindex(s)
    den[idx] = denpoly(s.m[idx])
  end
  return den
end

function zpkdata{T1}(s::DMimo{T1})
  zpkdata_ = Matrix{Vector}(s.ny, s.nu)
  for idx in eachindex(s)
    zpkdata_[idx] = zpkdata(s.m[idx])
  end
  return zpkdata_
end

samplingtime{T1}(s::DMimo{T1}) = samplingtime(s)

ndims(s::DMimo)  = 2
size(s::DMimo)   = (s.ny, s.nu)

function getindex(s::DMimo, idx::Int)
  row, col = divrem(idx-1, s.ny)
  DMimo(fill(s.m[row+1, col+1],1,1), samplingtime(s))
end

function getindex(s::DMimo, rows, cols)
  s2 = try
      [s.m[row, col] for row in rows, col in cols]
    catch exception
      warn("s[,j]: Index out of bounds")
      throw(exception)
    end
  DMimo(s2, samplingtime(s))
end

getindex(s::DMimo, ::Colon, ::Colon) = DMimo(s.m, samplingtime(s)) # returns a copy of s

function getindex(s::DMimo, ::Colon, cols)
  s2 = try
      [s.m[row, col] for row in 1:s.ny, col in cols]
    catch exception
      warn("s[,j]: Index out of bounds")
      throw(exception)
    end
  DMimo(s2, samplingtime(s))
end

function getindex(s::DMimo, rows, ::Colon)
  s2 = try
      [s.m[row, col] for row in rows, col in 1:s.nu]
    catch exception
      warn("s[,j]: Index out of bounds")
      throw(exception)
    end
  DMimo(s2, samplingtime(s))
end

start(s::DMimo)       = 1
next(s::DMimo, state::Int) = (s.m[state], state+1)
done(s::DMimo, state::Int) = state > length(s)
eltype{T1}(s::DMimo{T1}) = T1
length(s::DMimo) = length(s.m)
eachindex(s::DMimo) = 1:length(s)
endof(s::DMimo) = length(s)

showcompact(io::IO, s::DMimo) = print(io, summary(s))

function show(io::IO, s::DMimo)
  println(io, "Discrete time transfer function model")
  println(io, "\ty = Gu")
  if samplingtime(s) > 0
    println(io, "with nu=", s.nu, ", ny=", s.ny, ", Ts=", samplingtime(s),".")
  elseif samplingtime(s) == 0
    println(io, "with nu=", s.nu, ", ny=", s.ny, ", Ts=unspecified", samplingtime(s),".")
  end
end

function showall{T1}(io::IO, s::DMimo{T1})
  show(io, s)
  println(io, "")
  for subs in s
    printtransferfunction(io::IO, subs)
    println(io, "")
  end
end

function summary{T1<:Real}(s::DMimo{DSisoRational{T1}})
  string("tf(nu=", s.nu, ", ny=", s.ny, ", Ts=", samplingtime(s), ")")
end
function summary{T1<:AbstractFloat}(s::DMimo{DSisoZpk{T1}})
  string("zpk(nu=", s.nu, ", ny=", s.ny, ", Ts=", samplingtime(s), ")")
end

function +{T1, T2}(
    s1::DMimo{T1}, s2::DMimo{T2})
  if size(s1.m) != size(s2.m)
    warn("Systems have different shapes")
    throw(DomainError)
  elseif s1.Ts != s2.Ts
    warn("Sampling time mismatch")
    throw(DomainError)
  end
  m = s1.m + s2.m
  return tf(m)
end

+{T1, T2<:Real}(s::DMimo{T1}, n::T2) =
  DMimo(s.m+n, samplingtime(s))
+{T1, T2<:Real}(n::T2, s::DMimo{T1})  = +(s, n)

function +{T1, T2<:Real}(s::DMimo{T1}, n::Matrix{T2})
  if size(s.m) != size(n)
    warn("Systems have different shapes")
    throw(DomainError)
  end
  DMimo(s.m + n, samplingtime(s))
end
+{T1, T2<:Real}(n::Matrix{T2}, s::DMimo{T1})  = (s, n)

.+{T1, T2<:Real}(n::T2, s::DMimo{T1})         = +(n, s)
.+{T1, T2<:Real}(s::DMimo{T1}, n::T2)         = +(s, n)

.+{T1, T2<:Real}(n::Matrix{T2}, s::DMimo{T1}) = +(n, s)
.+{T1, T2<:Real}(s::DMimo{T1}, n::Matrix{T2}) = +(s, n)
.+{T1, T2}(
      s1::DMimo{T1}, s2::DMimo{T2})                    = +(s1, s2)

-{T1}(s::DMimo{T1})                   = tf(-s.m)

-{T1, T2<:Real}(s::DMimo{T1}, n::Matrix{T2}) = +(s, -n)
-{T1, T2<:Real}(n::Matrix{T2}, s::DMimo{T1}) = +(-n, s)
-{T1, T2}(s1::DMimo{T1},s2::DMimo{T2}) = +(s1,-s2)
-{T1, T2<:Real}(s::DMimo{T1}, n::T2)  = +(s,-n)
-{T1, T2<:Real}(n::T2, s::DMimo{T1})  = +(-n, s)

.-{T1, T2<:Real}(n::T2, s::DMimo{T1}) = +(-n, s)
.-{T1, T2<:Real}(s::DMimo{T1}, n::T2) = +(s,-n)

.-{T1, T2<:Real}(n::Matrix{T2}, s::DMimo{T1}) = -(n, s)
.-{T1, T2<:Real}(s::DMimo{T1}, n::Matrix{T2}) = -(s, n)

function *{T1, T2}(
  s1::DMimo{T1}, s2::DMimo{T2})
  if size(s1.m, 2) != size(s2.m, 1)
    warn("s1*s2: s1 must have same number of inputs as s2 has outputs")
    throw(DomainError())
  elseif s1.Ts != s2.Ts
    warn("Sampling time mismatch")
    throw(DomainError())
  end
  tf(s1.m*s2.m)
end

function *{T1, T2<:Real}(
  s::DMimo{T1}, n::T2)
  tf(s.m*n)
end
*{T1, T2<:Real}(n::T2, s::DMimo{T1})  = *(s, n)

*{T1, T2<:Real}(s::DMimo{T1}, n::Matrix{T2}) =
      DMimo(s.m*n, samplingtime(s))
*{T1, T2<:Real}(n::Matrix{T2}, s::DMimo{T1}) =
      DMimo(n*s.m, samplingtime(s))

.*{T1, T2<:Real}(n::T2, s::DMimo{T1}) = *(n, s)
.*{T1, T2<:Real}(s::DMimo{T1}, n::T2) = *(s, n)
.*{T1, T2<:Real}(n::Matrix{T2}, s::DMimo{T1}) = *(n, s)
.*{T1, T2<:Real}(s::DMimo{T1}, n::Matrix{T2}) = *(s, n)
.*{T1, T2}(
  s1::DMimo{T1}, s2::DMimo{T2}) = *(s1,s2)

function /{T1, T2<:Real}(n::T2, s::DMimo{T1})
  warn("MIMO TransferFunction inversion isn't implemented yet")
  throw(DomainError())
end
/{T1, T2}(
  s1::DMimo{T1}, s2::DMimo{T2}) = s1*(1/s2)
/{T1, T2<:Real}(s::DMimo{T1}, n::T2) = s*(1/n)
function /{T1, T2<:Real}(n::Matrix{T2}, t::DMimo{T1})
  warn("MIMO TransferFunction inversion isn't implemented yet")
  throw(DomainError())
end
/{T1, T2<:Real}(s::DMimo{T1}, n::Matrix{T2}) =
  tf(s.m/n)

./{T1, T2<:Real}(s::DMimo{T1}, n::T2) = s/n
./{T1, T2<:Real}(n::T2, s::DMimo{T1}) = n/s
./{T1, T2<:Real}(s::DMimo{T1}, n::Matrix{T2}) = s/n
./{T1, T2<:Real}(n::Matrix{T2}, s::DMimo{T1}) = n/s

function =={T1, T2}(
  s1::DMimo{T1}, s2::DMimo{T2})
  s1.m != s2.m
end

!={T1, T2}(
  s1::DMimo{T1}, s2::DMimo{T2}) = !(s1.m == s2.m)

function isapprox{T1, T2}(
  s1::DMimo{T1}, s2::DMimo{T2})
  # TODO: Implement
end
