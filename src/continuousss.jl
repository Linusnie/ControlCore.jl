# continuous state space type definition

immutable ContinuousSs{T<:Real} <: StateSpace
  A::Matrix{T}
  B::Matrix{T}
  C::Matrix{T}
  D::Matrix{T}
  nx::Int
  nu::Int
  ny::Int

  function call{T}(::Type{ContinuousSs}, A::Matrix{T}, B::Matrix{T},
    C::Matrix{T}, D::Matrix{T})
    na, ma = size(A)
    nb, mb = size(B)
    nc, mc = size(C)
    nd, md = size(D)

    if na != ma
      warn("A must be square")
      throw(DomainError())
    elseif nb != na
      warn("B must have the same row size as A")
      throw(DomainError())
    elseif mc != ma
      warn("C must have the same column size as A")
      throw(DomainError())
    elseif md != mb
      warn("D must have the same column size as B")
      throw(DomainError())
    elseif nd != nc
      warn("D must have the same row size as C")
      throw(DomainError())
    end

    new{T}(copy(A), copy(B), copy(C), copy(D), na, mb, nc)
  end
end

# creation of continuous state space types

function ss{T1<:Real, n1, T2<:Real, n2, T3<:Real, n3, T4<:Real, n4}(
  A::Array{T1,n1}, B::Array{T2,n2}, C::Array{T3,n3},
  D::Array{T4,n4} = zeros(Int8,0,0))
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

  ContinuousSs(promote(a,b,c,d)...)
end

function ss{T<:Real, n}(D::Array{T,n})
  @assert 0 < n < 3 "Gain matrix can have at most 2 dimensions"

  d = (n == one(n)) ? reshape(D, size(D, 1), 1) : D

  a = zeros(T, 0, 0)
  b = zeros(T, 0, size(D, 2))
  c = zeros(T, size(D, 1), 0)

  ContinuousSs(a,b,c,d)
end

ss{T<:Real}(D::T) = ss(D*ones(T, 1, 1))

# overloading Base functions

one{T<:Real}(s::ContinuousSs{T})        = ss(ones(T, 1, 1))
one{T<:Real}(::Type{ContinuousSs{T}})   = ss(ones(T, 1, 1))
zero{T<:Real}(s::ContinuousSs{T})       = ss(zeros(T, 1, 1))
zero{T<:Real}(::Type{ContinuousSs{T}})  = ss(zeros(T, 1, 1))

function zeros{T<:Real}(s::ContinuousSs{T})
  # TODO: Implement
end

function poles{T<:Real}(s::ContinuousSs{T})
  # TODO: Implement
end

ndims{T<:Real}(s::ContinuousSs{T})          = 2
size{T<:Real}(s::ContinuousSs{T})           = (s.ny, s.nu)
size{T<:Real}(s::ContinuousSs{T}, dim::Int) = dim > ndims(s) ? 1 : size(s)[dim]

function getindex{T<:Real}(s::ContinuousSs{T}, outidx::Int, inidx::Int)
  ContinuousSs(s.A, hcat(s.B[:, inidx]), s.C[outidx, :],
    fill(s.D[outidx, inidx], 1, 1))
end

function show{T<:Real}(io::IO, s::ContinuousSs{T})
    println(io, "Continuous time state space model")
    println(io, "\tẋ = Ax + Bu")
    println(io, "\ty = Cx + Du")
    println(io, "with nx=", s.nx, ", nu=", s.nu, ", ny=", s.ny, ".")
end

showcompact{T<:Real}(io::IO, s::ContinuousSs{T}) = print(io, summary(s))

function summary{T<:Real}(s::ContinuousSs{T})
  string("ss(nx=", s.nx, ",nu=", s.nu, ",ny=", s.ny, ")")
end

-{T<:Real}(s::ContinuousSs{T}) = ContinuousSs(s.A, s.B, -s.C, -s.D)

function +{T1<:Real, T2<:Real}(s1::ContinuousSs{T1}, s2::ContinuousSs{T2})
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

  ContinuousSs(a,b,c,d)
end

function +{T1<:Real, T2<:Real}(s::ContinuousSs{T1}, g::T2)
  d = s.D + g
  ContinuousSs(promote(s.A, s.B, s.C, d)...)
end
+{T1<:Real, T2<:Real}(g::T2, s::ContinuousSs{T1}) = +(s, g)

function -{T1<:Real, T2<:Real}(s1::ContinuousSs{T1}, s2::ContinuousSs{T2})
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

  ContinuousSs(a,b,c,d)
end

function -{T1<:Real, T2<:Real}(s::ContinuousSs{T1}, g::T2)
  d = s.D - g
  ContinuousSs(promote(s.A, s.B, s.C, d)...)
end
-{T1<:Real, T2<:Real}(g::T2, s::ContinuousSs{T1}) = +(g, -s)

function *{T1<:Real, T2<:Real}(s1::ContinuousSs{T1}, s2::ContinuousSs{T2})
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

  ContinuousSs(a,b,c,d)
end

*{T1<:Real, T2<:Real}(s::ContinuousSs{T1}, g::T2) =
  ContinuousSs(promote(s.A, s.B, s.C*g, s.D*g)...)

*{T1<:Real, T2<:Real}(g::T2, s::ContinuousSs{T1}) = *(s, g)

function /{T1<:Real, T2<:Real}(n::T2, s::ContinuousSs{T1})
  # Ensure s.D is invertible
  Dinv = try
    inv(s.D)
  catch e
    warn("n/s: s.D not invertible")
    throw(e)
  end

  a = s.A - s.B*Dinv*s.C
  b = s.B*Dinv
  c = -n*Dinv*s.C
  d = n*Dinv

  ContinuousSs(promote(a,b,c,d)...)
end

/{T1<:Real, T2<:Real}(s::ContinuousSs{T1}, n::T2) =
  ContinuousSs(promote(s.A, s.B, s.C/n, s.D/n)...)

function =={T1<:Real, T2<:Real}(s1::ContinuousSs{T1}, s2::ContinuousSs{T2})
    for field in fieldnames(ContinuousSs)
        if getfield(s1, field) != getfield(s2, field)
            return false
        end
    end
    return true
end

!={T1<:Real, T2<:Real}(s1::ContinuousSs{T1}, s2::ContinuousSs{T2}) = !(s1 == s2)
