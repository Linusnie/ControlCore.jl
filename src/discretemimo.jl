immutable DiscreteMimo{T1<:DiscreteSisoTf} <: TransferFunction
  m::Matrix{T1}
  ny::Int
  nu::Int
  Ts::Float64
  function call{T1}(::Type{DiscreteMimo}, m::Matrix{T1}, Ts::Float64)
    return new{T1}(m, size(m,1), size(m,2), max(Ts,0.))
  end
end

function tf{T1<:Real, T2<:Real}(num::Matrix{Vector{T1}}, den::Matrix{Vector{T2}}, Ts::Float64)
  # Validate input and output dimensions match
  ny, nu = size(num, 1, 2)
  if (ny, nu) != size(den, 1, 2)
    warn("num and den dimensions must match")
    throw(DomainError())
  end
  num_, den_ = promote(num, den)
  m = Array(DiscreteSisoRational{promote_type(T1,T2)}, ny, nu)
  for idx in eachindex(m)
    m[idx] = tf(num_[idx], den_[idx], Ts)
  end
  DiscreteMimo(m, Ts)
end

function tf{T1<:Real}(gain::Matrix{T1}, Ts::Float64)
  ny, nu = size(gain, 1, 2)
  m = Matrix{DiscreteSisoRational{T1}}(ny, nu)
  for idx in eachindex(m)
    m[idx] = tf(gain[idx], Ts)
  end
  DiscreteMimo(m, Ts)
end

function zpk{T1<:Number, T2<:Number, T3<:Real}(z::Matrix{Vector{T1}}, p::Matrix{Vector{T2}}, k::Matrix{T3}, Ts::Float64)
  # Validate input and output dimensions match
  ny, nu = size(z, 1, 2)
  if (ny, nu) != size(p, 1, 2) || (ny, nu) != size(k, 1, 2)
    warn("num and den dimensions must match")
    throw(DomainError())
  end
  T  = promote_type(real(T1), real(T2))
  z_ = convert(Matrix{Vector{Complex{Float64}}},z)
  p_ = convert(Matrix{Vector{Complex{Float64}}},p)
  m = Array(DiscreteSisoZpk{Float64, Float64}, ny, nu)
  for idx in eachindex(m)
    m[idx] = zpk(z_[idx], p_[idx], Float64(k[idx]), Ts)
  end
  DiscreteMimo(m, Ts)
end

function zpk{T1<:Real}(gain::Matrix{T1}, Ts::Float64)
  ny, nu = size(gain, 1, 2)
  m = Matrix{DiscreteSisoZpk{Float64, Float64}}(ny, nu)
  for idx in eachindex(m)
    m[idx] = zpk(gain[idx], Ts)
  end
  DiscreteMimo(m, Ts)
end

function tf{T1<:DiscreteSisoTf}(m::Matrix{T1})
  Ts = m[1,1].Ts
  for idx in eachindex(m)
    if Ts != m[idx].Ts
      warn("Sampling time mismatch")
      throw(DomainError())
    end
  end
  DiscreteMimo(m, Ts)
end

one{T1<:DiscreteSisoTf}(s::DiscreteMimo{T1}) =
  length(s) > 0 ? DiscreteMimo(fill(one(s.m[1]),1,1), s.Ts) :
                  DiscreteMimo(fill(one(T1),1,1), s.Ts)

function one{T1<:DiscreteSisoTf}(::Type{DiscreteMimo{T1}})
  s = one(T1)
  return DiscreteMimo(fill(s,1,1), s.Ts)
end
function zero{T1<:DiscreteSisoTf}(s::DiscreteMimo{T1})
  if length(s) > 0
    return DiscreteMimo(fill(zero(s.m[1]),1,1), s.Ts)
  else
    return DiscreteMimo(fill(zero(T1),1,1), s.Ts)
  end
end
function zero{T1<:DiscreteSisoTf}(::Type{DiscreteMimo{T1}})
  s = zero(T1)
  return DiscreteMimo(fill(s,1,1), s.Ts)
end

function zeros(s::DiscreteMimo)
  z = Matrix{Vector}(s.ny, s.nu)
  for idx in eachindex(z)
    z[idx] = poles(t[idx])
  end
  return z
end

function poles(s::DiscreteMimo)
  p = Matrix{Vector}(s.ny, s.nu)
  for idx in eachindex(p)
    p[idx] = poles(s.m[idx])
  end
  return p
end

ndims(s::DiscreteMimo)  = 2
size(s::DiscreteMimo)   = (s.ny, s.nu)

function getindex(s::DiscreteMimo, idx::Int)
  row, col = divrem(idx-1, s.ny)
  DiscreteMimo(fill(s.m[row+1, col+1],1,1), s.Ts)
end

function getindex(s::DiscreteMimo, rows, cols)
  s2 = try
      [s.m[row, col] for row in rows, col in cols]
    catch exception
      warn("s[,j]: Index out of bounds")
      throw(exception)
    end
  DiscreteMimo(s2, s.Ts)
end

getindex(s::DiscreteMimo, ::Colon, ::Colon) = DiscreteMimo(s.m, s.Ts) # returns a copy of s

function getindex(s::DiscreteMimo, ::Colon, cols)
  s2 = try
      [s.m[row, col] for row in 1:s.ny, col in cols]
    catch exception
      warn("s[,j]: Index out of bounds")
      throw(exception)
    end
  DiscreteMimo(s2, s.Ts)
end

function getindex(s::DiscreteMimo, rows, ::Colon)
  s2 = try
      [s.m[row, col] for row in rows, col in 1:s.nu]
    catch exception
      warn("s[,j]: Index out of bounds")
      throw(exception)
    end
  DiscreteMimo(s2, s.Ts)
end

start(s::DiscreteMimo)       = 1
next(s::DiscreteMimo, state) = (s.m[state], state+1)
done(s::DiscreteMimo, state) = state > length(s)
eltype{T1<:Real}(::Type{DiscreteMimo{DiscreteSisoRational{T1}}}) = DiscreteMimo{DiscreteSisoRational{T1}}
eltype{T1<:Number, T2<:Real}(::Type{DiscreteMimo{DiscreteSisoZpk{T1,T2}}}) = DiscreteMimo{DiscreteSisoZpk{T1,T2}}
length(s::DiscreteMimo) = length(s.m)
eachindex(s::DiscreteMimo) = 1:length(s)
endof(s::DiscreteMimo) = length(s)

showcompact(io::IO, s::DiscreteMimo) = print(io, summary(s))

function show(io::IO, s::DiscreteMimo)
  println(io, "Discrete time transfer function model")
  println(io, "\ty = Gu")
  if s.Ts > 0
    println(io, "with nu=", s.nu, ", ny=", s.ny, ", Ts=", s.Ts,".")
  elseif s.Ts == 0
    println(io, "with nu=", s.nu, ", ny=", s.ny, ", Ts=unspecified", s.Ts,".")
  end
end

function showall{T1<:DiscreteSisoTf}(io::IO, s::DiscreteMimo{T1})
  show(io, s)
  println(io, "")
  for subs in s
    printtransferfunction(io::IO, subs)
    println(io, "")
  end
end

function summary{T1<:Real}(s::DiscreteMimo{DiscreteSisoRational{T1}})
  string("tf(nu=", s.nu, ", ny=", s.ny, ", Ts=", s.Ts, ")")
end
function summary{T1<:Number, T2<:Real}(s::DiscreteMimo{DiscreteSisoZpk{T1, T2}})
  string("zpk(nu=", s.nu, ", ny=", s.ny, ", Ts=", s.Ts, ")")
end

function +{T1<:DiscreteSisoTf, T2<:DiscreteSisoTf}(
    s1::DiscreteMimo{T1}, s2::DiscreteMimo{T2})
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

+{T1<:DiscreteSisoTf, T2<:Real}(s::DiscreteMimo{T1}, n::T2) =
  DiscreteMimo(s.m+n, s.Ts)
+{T1<:DiscreteSisoTf, T2<:Real}(n::T2, s::DiscreteMimo{T1})  = +(s, n)

function +{T1<:DiscreteSisoTf, T2<:Real}(s::DiscreteMimo{T1}, n::Matrix{T2})
  if size(s.m) != size(n)
    warn("Systems have different shapes")
    throw(DomainError)
  end
  DiscreteMimo(s.m + n, s.Ts)
end
+{T1<:DiscreteSisoTf, T2<:Real}(n::Matrix{T2}, s::DiscreteMimo{T1})  = (s, n)

.+{T1<:DiscreteSisoTf, T2<:Real}(n::T2, s::DiscreteMimo{T1})         = +(n, s)
.+{T1<:DiscreteSisoTf, T2<:Real}(s::DiscreteMimo{T1}, n::T2)         = +(s, n)

.+{T1<:DiscreteSisoTf, T2<:Real}(n::Matrix{T2}, s::DiscreteMimo{T1}) = +(n, s)
.+{T1<:DiscreteSisoTf, T2<:Real}(s::DiscreteMimo{T1}, n::Matrix{T2}) = +(s, n)
.+{T1<:DiscreteSisoTf, T2<:DiscreteSisoTf}(
      s1::DiscreteMimo{T1}, s2::DiscreteMimo{T2})                    = +(s1, s2)

-{T1<:DiscreteSisoTf}(s::DiscreteMimo{T1})                   = tf(-s.m)
-{T1<:DiscreteSisoTf, T2<:Real}(s::DiscreteMimo{T1}, n::Matrix{T2}) = +(s, -n)
-{T1<:DiscreteSisoTf, T2<:Real}(n::Matrix{T2}, s::DiscreteMimo{T1}) = +(-n, s)
-{T1<:DiscreteSisoTf, T2<:DiscreteSisoTf}(s1::DiscreteMimo{T1},s2::DiscreteMimo{T2}) = +(s1,-s2)
-{T1<:DiscreteSisoTf, T2<:Real}(s::DiscreteMimo{T1}, n::T2)  = +(s,-n)
-{T1<:DiscreteSisoTf, T2<:Real}(n::T2, s::DiscreteMimo{T1})  = +(-n, s)

.-{T1<:DiscreteSisoTf, T2<:Real}(n::T2, s::DiscreteMimo{T1}) = +(-n, s)
.-{T1<:DiscreteSisoTf, T2<:Real}(s::DiscreteMimo{T1}, n::T2) = +(s,-n)

function *{T1<:DiscreteSisoTf, T2<:DiscreteSisoTf}(
  s1::DiscreteMimo{T1}, s2::DiscreteMimo{T2})
  if size(s1.m, 2) != size(s2.m, 1)
    warn("s1*s2: s1 must have same number of inputs as s2 has outputs")
    throw(DomainError())
  elseif s1.Ts != s2.Ts
    warn("Sampling time mismatch")
    throw(DomainError())
  end
  tf(s1.m*s2.m)
end

function *{T1<:DiscreteSisoTf, T2<:Real}(
  s::DiscreteMimo{T1}, n::T2)
  tf(s.m*n)
end
*{T1<:DiscreteSisoTf, T2<:Real}(n::T2, s::DiscreteMimo{T1})  = *(s, n)

*{T1<:DiscreteSisoTf, T2<:Real}(s::DiscreteMimo{T1}, n::Matrix{T2}) =
      DiscreteMimo(s.m*n, s.Ts)
*{T1<:DiscreteSisoTf, T2<:Real}(n::Matrix{T2}, s::DiscreteMimo{T1}) =
      DiscreteMimo(n*s.m, s.Ts)

.*{T1<:DiscreteSisoTf, T2<:Real}(n::T2, s::DiscreteMimo{T1}) = *(n, s)
.*{T1<:DiscreteSisoTf, T2<:Real}(s::DiscreteMimo{T1}, n::T2) = *(s, n)
.*{T1<:DiscreteSisoTf, T2<:Real}(s::DiscreteMimo{T1}, n::Matrix{T2}) =
      DiscreteMimo(s.m*n, s.Ts)
.*{T1<:DiscreteSisoTf, T2<:Real}(n::Matrix{T2}, s::DiscreteMimo{T1}) =
      DiscreteMimo(n*s.m, s.Ts)

function /{T1<:DiscreteSisoTf, T2<:Real}(n::T2, s::DiscreteMimo{T1})
  warn("MIMO TransferFunction inversion isn't implemented yet")
  throw(DomainError())
end
/{T1<:DiscreteSisoTf, T2<:DiscreteSisoTf}(
  s1::DiscreteMimo{T1}, s2::DiscreteMimo{T2}) = s1*(1/s2)

./{T1<:DiscreteSisoTf, T2<:Real}(s::DiscreteMimo{T1}, n::T2) = s*(1/n)
./{T1<:DiscreteSisoTf, T2<:Real}(n::T2, s::DiscreteMimo{T1}) = DiscreteMimo(n./s.m, s.Ts)

function =={T1<:DiscreteSisoTf, T2<:DiscreteSisoTf}(
  s1::DiscreteMimo{T1}, s2::DiscreteMimo{T2})
  s1.m != s2.m
end

!={T1<:DiscreteSisoTf, T2<:DiscreteSisoTf}(
  s1::DiscreteMimo{T1}, s2::DiscreteMimo{T2}) = !(s1.m == s2.m)

function isapprox{T1<:DiscreteSisoTf, T2<:DiscreteSisoTf}(
  s1::DiscreteMimo{T1}, s2::DiscreteMimo{T2})
  # TODO: Implement
end
