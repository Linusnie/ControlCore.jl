immutable ContinuousSsMimo{T<:Real} <: ContinuousSs
  A::Matrix{T}
  B::Matrix{T}
  C::Matrix{T}
  D::Matrix{T}
  nx::Int
  nu::Int
  ny::Int

  function call{T}(::Type{ContinuousSsMimo}, A::Matrix{T}, B::Matrix{T},
    C::Matrix{T}, D::Matrix{T})
    na, ma = size(A)
    nb, mb = size(B)
    nc, mc = size(C)
    nd, md = size(D)

    if na != ma
      warn("A must be square")
      throw(DomainError())
    elseif nb != na
      warn("B must have the same row size as that of A")
      throw(DomainError())
    elseif mc != ma
      warn("C must have the same column size as that of A")
      throw(DomainError())
    elseif md != mb
      warn("D must have the same column size as that of B")
      throw(DomainError())
    elseif nd != nc
      warn("D must have the same row size as that of C")
      throw(DomainError())
    end

    new{T}(A, B, C, D, na, mb, nc)
  end
end

# creation of continuous state space types

function ss{T1<:Real, n1, T2<:Real, n2, T3<:Real, n3, T4<:Real, n4}(
  A::Array{T1,n1}, B::Array{T2,n2}, C::Array{T3,n3}, D::Array{T4,n4})
  @assert 0 < n1 < 3 "A can have at most 2 dimensions"
  @assert 0 < n2 < 3 "B can have at most 2 dimensions"
  @assert 0 < n3 < 3 "C can have at most 2 dimensions"
  @assert 0 < n4 < 3 "D can have at most 2 dimensions"

  # if A is a vector, it should be upgraded to a matrix
  na, ma = size(A, 1, 2)
  if na != ma
    warn("A must be square")
    throw(DomainError())
  end

  a = (n1 == one(n1)) ? reshape(A, na, ma) : A

  # if B is a vector, it should be upgraded to a matrix
  nb, mb = size(B, 1, 2)

  b = (n2 == one(n2)) ? ((n1 == one(n1)) ? reshape(B, mb, nb) :
    reshape(B, nb, mb)) : B

  # if C is a vector, it should be upgraded to a matrix
  nc, mc = size(C, 1, 2)

  c = (n3 == one(n3)) ? ((n1 == one(n1)) ? reshape(C, nc, mc) :
    reshape(C, mc, nc)) : C

  # if D is a vector, it should be upgraded to a matrix
  d = isempty(D) ? zeros(size(c, 1), size(b, 2)) :
    ((n4 == one(n4)) ? reshape(D, size(D, 1), size(D, 2)) : D)

  ContinuousSsMimo(promote(a,b,c,d)...)
end

ss{T<:Real}(g::T) = ss(g*ones(T, 1, 1))

function ss{T<:Real}(g::Vector{T})
  a = zeros(T, 0, 0)
  b = zeros(T, 0, 1)
  c = zeros(T, size(g, 1), 0)
  d = reshape(g, size(g, 1), 1)

  ContinuousSsMimo(a,b,c,d)
end


function ss{T<:Real}(g::Matrix{T})
  a = zeros(T, 0, 0)
  b = zeros(T, 0, size(g, 2))
  c = zeros(T, size(g, 1), 0)

  ContinuousSsMimo(a,b,c,g)
end

# overloading identities

one{T<:Real}(s::ContinuousSsMimo{T})        = ss(ones(T, 1, 1))
one{T<:Real}(::Type{ContinuousSsMimo{T}})   = ss(ones(T, 1, 1))
zero{T<:Real}(s::ContinuousSsMimo{T})       = ss(zeros(T, 1, 1))
zero{T<:Real}(::Type{ContinuousSsMimo{T}})  = ss(zeros(T, 1, 1))

# overload inv and zeros

function inv{T<:Real}(s::ContinuousSsMimo{T})
  # TODO: Implement
end

function zeros{T<:Real}(s::ContinuousSsMimo{T})
  # TODO: Implement
end

function poles{T<:Real}(s::ContinuousSsMimo{T})
  # TODO: Implement
end

# overload slicing functions

ndims{T<:Real}(s::ContinuousSsMimo{T})          = 2
size{T<:Real}(s::ContinuousSsMimo{T})           = (s.ny, s.nu)
size{T<:Real}(s::ContinuousSsMimo{T}, dim::Int) = dim > ndims(s) ? 1 : size(s)[dim]

function getindex{T<:Real}(s::ContinuousSsMimo{T}, idx::Int)
  if idx < 1 || idx > length(s.D)
    warn("s[idx]: Trying to access idx < 1 or idx > length(s.D)")
    throw(BoundsError(s.D, idx))
  end

  col, row = divrem(idx-1, s.ny)

  ContinuousSsMimo(s.A, hcat(s.B[:, col+1]), s.C[row+1, :],
    fill(s.D[row+1, col+1], 1, 1))
end

getindex{T<:Real}(s::ContinuousSsMimo{T}, ::Colon) = s

function getindex{T<:Real}(s::ContinuousSsMimo{T}, row::Int, col::Int)
  if row < 1 || row > s.ny
    warn("s[i,]: Trying to access non-existent outputs")
    throw(BoundsError(s.C, row))
  elseif col < 1 || col > s.nu
    warn("s[,j]: Trying to access non-existent inputs")
    throw(BoundsError(s.B, col))
  end

  ContinuousSsMimo(s.A, hcat(s.B[:, col]), s.C[row, :],
    fill(s.D[row, col], 1, 1))
end

function getindex{T<:Real}(s::ContinuousSsMimo{T}, rows, cols)
  b = try
    [s.B[row, col] for row in 1:s.nx, col in cols]
  catch exception
    warn("s[,j]: Trying to access non-existent inputs")
    throw(exception)
  end

  c = try
    [s.C[row, col] for row in rows, col in 1:s.nx]
  catch exception
    warn("s[i,]: Trying to access non-existent outputs")
    throw(exception)
  end

  d = [s.D[row, col] for row in rows, col in cols]

  ContinuousSsMimo(s.A, b, c, d)
end

getindex{T<:Real}(s::ContinuousSsMimo{T}, ::Colon, ::Colon) = s

function getindex{T<:Real}(s::ContinuousSsMimo{T}, ::Colon, cols)
  b = try
    [s.B[row, col] for row in 1:s.nx, col in cols]
  catch exception
    warn("s[,j]: Trying to access non-existent inputs")
    throw(exception)
  end

  d = [s.D[row, col] for row in 1:s.ny, col in cols]

  ContinuousSsMimo(s.A, b, s.C, d)
end

function getindex{T<:Real}(s::ContinuousSsMimo{T}, rows, ::Colon)
  c = try
    [s.C[row, col] for row in rows, col in 1:s.nx]
  catch exception
    warn("s[i,]: Trying to access non-existent outputs")
    throw(exception)
  end

  d = [s.D[row, col] for row in rows, col in 1:s.nu]

  ContinuousSsMimo(s.A, s.B, c, d)
end

# overload iteration interface

start{T<:Real}(::ContinuousSsMimo{T})         = 1
next{T<:Real}(s::ContinuousSsMimo{T}, state)  = (s[state], state+1)
done{T<:Real}(s::ContinuousSsMimo{T}, state)  = state > length(s.D)
eltype{T<:Real}(::Type{ContinuousSsMimo{T}})  = ContinuousSsSiso{T}
length{T<:Real}(s::ContinuousSsMimo{T})       = length(s.D)
eachindex{T<:Real}(s::ContinuousSsMimo{T})    = eachindex(s.D)
endof{T<:Real}(s::ContinuousSsMimo{T})        = endof(s.D)

# overload printing functions

showcompact{T<:Real}(io::IO, s::ContinuousSsMimo{T}) = print(io, summary(s))

function show{T<:Real}(io::IO, s::ContinuousSsMimo{T})
  println(io, "Continuous time state space model")
  println(io, "\tẋ = Ax + Bu")
  println(io, "\ty = Cx + Du")
  println(io, "with nx=", s.nx, ", nu=", s.nu, ", ny=", s.ny, ".")
end

function showall{T<:Real}(io::IO, s::ContinuousSsMimo{T})
  show(io, s)
  println(io, "System matrix (A):")
  println(io, s.A)
  println(io, "Input matrix (B):")
  println(io, s.B)
  println(io, "Output matrix (C):")
  println(io, s.C)
  println(io, "Feedforward matrix (D):")
  println(io, s.D)
end

function summary{T<:Real}(s::ContinuousSsMimo{T})
  string("ss(nx=", s.nx, ",nu=", s.nu, ",ny=", s.ny, ")")
end

# overload mathematical operations

-{T<:Real}(s::ContinuousSsMimo{T}) = ContinuousSsMimo(s.A, s.B, -s.C, -s.D)

function +{T1<:Real, T2<:Real}(s1::ContinuousSsMimo{T1}, s2::ContinuousSsMimo{T2})
  # Ensure systems have same shapes
  if size(s1) != size(s2)
    warn("s1+s2: size(s1) != size(s2)")
    throw(DomainError())
  end

  T3 = promote_type(T1, T2)

  a = vcat(hcat(s1.A, zeros(T3, s1.nx, s2.nx)),
          hcat(zeros(T3, s2.nx, s1.nx), s2.A))
  b = vcat(s1.B, s2.B)
  c = hcat(s1.C, s2.C)
  d = s1.D + s2.D

  ContinuousSsMimo(a,b,c,d)
end

+{T1<:Real, T2<:Real}(s::ContinuousSsMimo{T1}, g::T2) =
  ContinuousSsMimo(promote(copy(s.A), copy(s.B), copy(s.C), s.D + g)...)
+{T1<:Real, T2<:Real}(g::T2, s::ContinuousSsMimo{T1}) = +(s, g)

function +{T1<:Real, T2<:Real}(s::ContinuousSsMimo{T1}, g::Matrix{T2})
  if size(s.D) != size(g)
    warn("s+g: size(s.D) != size(g)")
    throw(DomainError())
  end

  ContinuousSsMimo(promote(copy(s.A), copy(s.B), copy(s.C), s.D + g)...)
end
+{T1<:Real, T2<:Real}(g::Matrix{T2}, s::ContinuousSsMimo{T1}) = +(s, g)

.+{T1<:Real, T2<:Real}(s1::ContinuousSsMimo{T1}, s2::ContinuousSsMimo{T2}) = +(s1, s2)
.+{T1<:Real, T2<:Real}(s::ContinuousSsMimo{T1}, g::T2)  = +(s, g)
.+{T1<:Real, T2<:Real}(g::T2, s::ContinuousSsMimo{T1})  = +(s, g)

function -{T1<:Real, T2<:Real}(s1::ContinuousSsMimo{T1}, s2::ContinuousSsMimo{T2})
  # Ensure systems have same shapes
  if size(s1) != size(s2)
    warn("s1-s2: size(s1) != size(s2)")
    throw(DomainError())
  end

  T3 = promote_type(T1, T2)

  a = vcat(hcat(s1.A, zeros(T3, s1.nx, s2.nx)),
          hcat(zeros(T3, s2.nx, s1.nx), s2.A))
  b = vcat(s1.B, s2.B)
  c = hcat(s1.C, -s2.C)
  d = s1.D - s2.D

  ContinuousSsMimo(a,b,c,d)
end

-{T1<:Real, T2<:Real}(s::ContinuousSsMimo{T1}, g::T2) =
  ContinuousSsMimo(promote(copy(s.A), copy(s.B), copy(s.C), s.D - g)...)
-{T1<:Real, T2<:Real}(g::T2, s::ContinuousSsMimo{T1}) = +(g, -s)

function -{T1<:Real, T2<:Real}(s::ContinuousSsMimo{T1}, g::Matrix{T2})
  if size(s.D) != size(g)
    warn("s+g: size(s.D) != size(g)")
    throw(DomainError())
  end

  ContinuousSsMimo(promote(copy(s.A), copy(s.B), copy(s.C), s.D - g)...)
end
-{T1<:Real, T2<:Real}(g::Matrix{T2}, s::ContinuousSsMimo{T1}) = +(g, -s)

.-{T1<:Real, T2<:Real}(s1::ContinuousSsMimo{T1}, s2::ContinuousSsMimo{T2}) = -(s1, s2)
.-{T1<:Real, T2<:Real}(s::ContinuousSsMimo{T1}, g::T2)  = -(s, g)
.-{T1<:Real, T2<:Real}(g::T2, s::ContinuousSsMimo{T1})  = +(g, -s)

function *{T1<:Real, T2<:Real}(s1::ContinuousSsMimo{T1}, s2::ContinuousSsMimo{T2})
  # Remark: s1*s2 implies u -> s2 -> s1 -> y
  if s1.nu != s2.ny
    warn("s1*s2: s1.nu != s2.ny")
    throw(DomainError())
  end

  T3 = promote_type(T1, T2)

  a = vcat(hcat(s1.A, s1.B*s2.C),
          hcat(zeros(T3, s2.nx, s1.nx), s2.A))
  b = vcat(s1.B*s2.D, s2.B)
  c = hcat(s1.C, s1.D*s2.C)
  d = s1.D * s2.D

  ContinuousSsMimo(a,b,c,d)
end

*{T1<:Real, T2<:Real}(s::ContinuousSsMimo{T1}, g::T2) =
  ContinuousSsMimo(promote(copy(s.A), s.B*g, copy(s.C), s.D*g)...)
*{T1<:Real, T2<:Real}(g::T2, s::ContinuousSsMimo{T1}) =
  ContinuousSsMimo(promote(copy(s.A), copy(s.B), g*s.C, g*s.D)...)

function *{T1<:Real, T2<:Real}(s::ContinuousSsMimo{T1}, g::Matrix{T2})
  if s.nu != size(g, 1)
    warn("s*g: s.nu != size(g, 1)")
    throw(DomainError())
  end

  ContinuousSsMimo(promote(copy(s.A), s.B*g, copy(s.C), s.D*g))
end

function *{T1<:Real, T2<:Real}(g::Matrix{T2}, s::ContinuousSsMimo{T1})
  if s.ny != size(g, 2)
    warn("g*s: s.ny != size(g, 2)")
    throw(DomainError())
  end

  ContinuousSsMimo(promote(copy(s.A), copy(s.B), g*s.C, g*s.D))
end

.*{T1<:Real, T2<:Real}(s1::ContinuousSsMimo{T1}, s2::ContinuousSsMimo{T2}) = *(s1, s2)
.*{T1<:Real, T2<:Real}(s::ContinuousSsMimo{T1}, g::T2)  = *(s, g)
.*{T1<:Real, T2<:Real}(g::T2, s::ContinuousSsMimo{T1})  = *(g, s)

/{T1<:Real, T2<:Real}(s1::ContinuousSsMimo{T2}, s2::ContinuousSsMimo{T1}) =
  *(s1, inv(s2))

/{T1<:Real, T2<:Real}(s::ContinuousSsMimo{T1}, g::T2) =
  ContinuousSsMimo(promote(copy(s.A), s.B/g, copy(s.C), s.D/g)...)
/{T1<:Real, T2<:Real}(g::T2, s::ContinuousSsMimo{T1}) = *(g, inv(s))

/{T1<:Real, T2<:Real}(s::ContinuousSsMimo{T1}, g::Matrix{T2}) =
  ContinuousSsMimo(promote(copy(s.A), s.B*inv(g), copy(s.C), s.D*inv(g))...)
/{T1<:Real, T2<:Real}(g::Matrix{T2}, s::ContinuousSsMimo{T1}) = *(g, inv(s))

./{T1<:Real, T2<:Real}(s1::ContinuousSsMimo{T1}, s2::ContinuousSsMimo{T2}) = /(s1, s2)
./{T1<:Real, T2<:Real}(s::ContinuousSsMimo{T1}, g::T2)  = /(s, g)
./{T1<:Real, T2<:Real}(g::T2, s::ContinuousSsMimo{T1})  = /(g, s)

function =={T1<:Real, T2<:Real}(s1::ContinuousSsMimo{T1}, s2::ContinuousSsMimo{T2})
  # TODO: Implement
end

!={T1<:Real, T2<:Real}(s1::ContinuousSsMimo{T1}, s2::ContinuousSsMimo{T2}) = !(s1 == s2)

function isapprox{T1<:Real, T2<:Real}(s1::ContinuousSsMimo{T1},
  s2::ContinuousSsMimo{T2})
  # TODO: Implement
end
