# continuous state space type definition

immutable ContinuousSsSiso{T<:Real} <: ContinuousSs
  A::Matrix{T}
  B::Vector{T}
  C::Vector{T}
  D::T
  nx::Int

  function call{T}(::Type{ContinuousSsSiso}, A::Matrix{T}, B::Vector{T},
    C::Vector{T}, D::T)
    na, ma  = size(A)
    nb      = length(B)
    nc      = length(C)

    if na != ma
      warn("A must be square")
      throw(DomainError())
    elseif nb != na
      warn("B must have a length of ", na)
      throw(DomainError())
    elseif nc != ma
      warn("C must have a length of ", na)
      throw(DomainError())
    end

    new{T}(A, B, C, D, na)
  end
end

# creation of continuous state space types

ss{T<:Real}(A::Matrix{T}, B::Vector{T}, C::Vector{T}, D::T = zero(T)) =
  ContinuousSsSiso(A, B, C, D)

ss{T<:Real}(A::T, B::Vector{T}, C::Vector{T}, D::T = zero(T)) =
  ContinuousSsSiso(fill(A, 1, 1), B, C, D)

function ss{T1<:Real, T2<:Real, T3<:Real, T4<:Real}(A::Matrix{T1}, B::Vector{T2},
  C::Vector{T3}, D::T4 = zero(T4))
  T = promote_type(T1, T2, T3, T4)
  ContinuousSsSiso(convert(Matrix{T}, A), convert(Vector{T}, B),
                  convert(Vector{T}, C), convert(T, D))
end

function ss{T1<:Real, T2<:Real, T3<:Real, T4<:Real}(A::T1, B::Vector{T2},
  C::Vector{T3}, D::T4 = zero(T4))
  T = promote_type(T1, T2, T3, T4)
  ContinuousSsSiso(fill(convert(T, A), 1, 1), convert(Vector{T}, B),
                  convert(Vector{T}, C), convert(T, D))
end

function ss{T1<:Real, n1, T2<:Real, n2, T3<:Real, n3, T4<:Real}(
  A::Array{T1,n1}, B::Array{T2,n2}, C::Array{T3,n3}, D::T4 = zero(T4))
  @assert 0 < n1 < 3 "A can have at most 2 dimensions"
  @assert 0 < n2 < 3 "B can have at most 2 dimensions"
  @assert 0 < n3 < 3 "C can have at most 2 dimensions"

  T = promote_type(T1, T2, T3, T4)

  # if A is a vector, it should be upgraded to a matrix
  na, ma = size(A, 1, 2)
  if na != ma
    warn("A must be square")
    throw(DomainError())
  end

  a = (n1 == one(n1)) ? reshape(A, na, ma) : A

  # if B is a matrix, it should be downgraded to a vector
  nb, mb = size(B, 1, 2)
  if mb > 1 && nb > 1
    warn("B must be a vector of size ", na)
    throw(DomainError())
  end

  b = (n2 == one(n2)) ? B : vec(B)

  # if C is a matrix, it should be downgraded to a vector
  nc, mc = size(C, 1, 2)
  if nc > 1 && mc > 1
    warn("C must be a vector of size ", na)
    throw(DomainError())
  end

  c = (n3 == one(n3)) ? C : vec(C)

  ContinuousSsSiso(convert(Matrix{T}, a), convert(Vector{T}, b),
                  convert(Vector{T}, c), convert(T, D))
end

ss{T<:Real}(g::T) = ContinuousSsSiso(zeros(T, 0, 0), zeros(T, 0), zeros(T, 0), g)

# conversion and promotion

promote_rule{T1<:Real, T2<:Real}(::Type{ContinuousSsSiso{T1}},
  ::Type{T2}) = ContinuousSsSiso{promote_type(T1, T2)}
convert{T1<:Real, T2<:Real}(::Type{ContinuousSsSiso{T1}}, g::T2) =
  ss(convert(T1, g))

# overloading identities

one{T<:Real}(s::ContinuousSsSiso{T})        = ss(one(T))
one{T<:Real}(::Type{ContinuousSsSiso{T}})   = ss(one(T))
zero{T<:Real}(s::ContinuousSsSiso{T})       = ss(zero(T))
zero{T<:Real}(::Type{ContinuousSsSiso{T}})  = ss(zero(T))

# overload inv and zeros

function inv{T<:Real}(s::ContinuousSsSiso{T})
  # TODO: Implement
  throw(ErrorException("inv(s) for s::ContinuousSsSiso is not implemented"))
end

function zeros{T<:Real}(s::ContinuousSsSiso{T})
  # TODO: Implement
  throw(ErrorException("zeros(s) for s::ContinuousSsSiso is not implemented"))
end

function poles{T<:Real}(s::ContinuousSsSiso{T})
  # TODO: Implement
  throw(ErrorException("poles(s) for s::ContinuousSsSiso is not implemented"))
end

# overload slicing functions

ndims{T<:Real}(s::ContinuousSsSiso{T})          = 1
size{T<:Real}(s::ContinuousSsSiso{T})           = 1
size{T<:Real}(s::ContinuousSsSiso{T}, dim::Int) = 1

function getindex{T<:Real}(s::ContinuousSsSiso{T}, idx::Int)
  if idx != 1
    warn("s[idx]: Trying to access idx != 1")
    throw(BoundsError(s.D, idx))
  end

  ContinuousSsSiso(s.A, copy(s.B), copy(s.C), s.D)
end

getindex{T<:Real}(s::ContinuousSsSiso{T}, ::Colon) = s

# overload iteration interface

start{T<:Real}(::ContinuousSsSiso{T})         = 1
next{T<:Real}(s::ContinuousSsSiso{T}, state)  = (s[state], state+1)
done{T<:Real}(s::ContinuousSsSiso{T}, state)  = state > 1
eltype{T<:Real}(::Type{ContinuousSsSiso{T}})  = ContinuousSsSiso{T}
length{T<:Real}(s::ContinuousSsSiso{T})       = 1
eachindex{T<:Real}(s::ContinuousSsSiso{T})    = 1
endof{T<:Real}(s::ContinuousSsSiso{T})        = endof(s.D)

# overload printing functions

showcompact{T<:Real}(io::IO, s::ContinuousSsSiso{T}) = print(io, summary(s))

function show{T<:Real}(io::IO, s::ContinuousSsSiso{T})
  println(io, "Continuous time state space model")
  println(io, "\tẋ = Ax + Bu")
  println(io, "\ty = Cx + Du")
  println(io, "with nx=", s.nx, ".")
end

function showall{T<:Real}(io::IO, s::ContinuousSsSiso{T})
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

function summary{T<:Real}(s::ContinuousSsSiso{T})
  string("ss(nx=", s.nx, ")")
end

# overload mathematical operations

-{T<:Real}(s::ContinuousSsSiso{T}) = ContinuousSsSiso(s.A, s.B, -s.C, -s.D)

function +{T1<:Real, T2<:Real}(s1::ContinuousSsSiso{T1}, s2::ContinuousSsSiso{T2})
  T = promote_type(T1, T2)

  a = vcat(hcat(s1.A, zeros(T, s1.nx, s2.nx)),
          hcat(zeros(T, s2.nx, s1.nx), s2.A))
  b = vcat(s1.B, s2.B)
  c = hcat(s1.C, s2.C)
  d = s1.D + s2.D

  ContinuousSsSiso(a,b,c,d)
end

function +{T1<:Real, T2<:Real}(s::ContinuousSsSiso{T1}, g::T2)
  T = promote_type(T1, T2)

  ContinuousSsSiso(convert(Matrix{T}, copy(s.A)), convert(Vector{T}, copy(s.B)),
                  convert(Vector{T}, copy(s.C)), s.D + g)
end
+{T1<:Real, T2<:Real}(g::T2, s::ContinuousSsSiso{T1}) = +(s, g)

.+{T1<:Real, T2<:Real}(s1::ContinuousSsSiso{T1}, s2::ContinuousSsSiso{T2}) = +(s1, s2)
.+{T1<:Real, T2<:Real}(s::ContinuousSsSiso{T1}, g::T2)  = +(s, g)
.+{T1<:Real, T2<:Real}(g::T2, s::ContinuousSsSiso{T1})  = +(s, g)

function -{T1<:Real, T2<:Real}(s1::ContinuousSsSiso{T1}, s2::ContinuousSsSiso{T2})
  T = promote_type(T1, T2)

  a = vcat(hcat(s1.A, zeros(T, s1.nx, s2.nx)),
          hcat(zeros(T, s2.nx, s1.nx), s2.A))
  b = vcat(s1.B, s2.B)
  c = hcat(s1.C, -s2.C)
  d = s1.D - s2.D

  ContinuousSsSiso(a,b,c,d)
end

function -{T1<:Real, T2<:Real}(s::ContinuousSsSiso{T1}, g::T2)
  T = promote_type(T1, T2)

  ContinuousSsSiso(convert(Matrix{T}, copy(s.A)), convert(Vector{T}, copy(s.B)),
                  convert(Vector{T}, copy(s.C)), s.D - g)
end
-{T1<:Real, T2<:Real}(g::T2, s::ContinuousSsSiso{T1}) = +(g, -s)

.-{T1<:Real, T2<:Real}(s1::ContinuousSsSiso{T1}, s2::ContinuousSsSiso{T2}) = -(s1, s2)
.-{T1<:Real, T2<:Real}(s::ContinuousSsSiso{T1}, g::T2)  = -(s, g)
.-{T1<:Real, T2<:Real}(g::T2, s::ContinuousSsSiso{T1})  = +(g, -s)

function *{T1<:Real, T2<:Real}(s1::ContinuousSsSiso{T1}, s2::ContinuousSsSiso{T2})
  # Remark: y = (s1*s2) u implies u -> s2 -> s1 -> y

  T = promote_type(T1, T2)

  a = vcat(hcat(s1.A, s1.B*s2.C'),
          hcat(zeros(T, s2.nx, s1.nx), s2.A))
  b = vcat(s1.B*s2.D, s2.B)
  c = hcat(s1.C, s1.D*s2.C)
  d = s1.D * s2.D

  ContinuousSsSiso(a,b,c,d)
end

function *{T1<:Real, T2<:Real}(s::ContinuousSsSiso{T1}, g::T2)
  T = promote_type(T1, T2)

  ContinuousSsSiso(convert(Matrix{T}, copy(s.A)), s.B*g,
                  convert(Vector{T}, copy(s.C)), s.D*g)
end
function *{T1<:Real, T2<:Real}(g::T2, s::ContinuousSsSiso{T1})
  T = promote_type(T1, T2)

  ContinuousSsSiso(convert(Matrix{T}, copy(s.A)), convert(Vector{T}, copy(s.B)),
                  g*s.C, g*s.D)
end

function *{T1<:Real, T2<:Real}(s::ContinuousSsSiso{T1}, g::Matrix{T2})
  if s.nu != size(g, 1)
    warn("s*g: s.nu != size(g, 1)")
    throw(DomainError())
  end

  ContinuousSsSiso(promote(copy(s.A), s.B*g, copy(s.C), s.D*g)...)
end

function *{T1<:Real, T2<:Real}(g::Matrix{T2}, s::ContinuousSsSiso{T1})
  if s.ny != size(g, 2)
    warn("g*s: s.ny != size(g, 2)")
    throw(DomainError())
  end

  ContinuousSsSiso(promote(copy(s.A), copy(s.B), g*s.C, g*s.D)...)
end

.*{T1<:Real, T2<:Real}(s1::ContinuousSsSiso{T1}, s2::ContinuousSsSiso{T2}) = *(s1, s2)
.*{T1<:Real, T2<:Real}(s::ContinuousSsSiso{T1}, g::T2)  = *(s, g)
.*{T1<:Real, T2<:Real}(g::T2, s::ContinuousSsSiso{T1})  = *(g, s)

/{T1<:Real, T2<:Real}(s1::ContinuousSsSiso{T2}, s2::ContinuousSsSiso{T1}) =
  *(s1, inv(s2))

/{T1<:Real, T2<:Real}(s::ContinuousSsSiso{T1}, g::T2) =
  ContinuousSsSiso(promote(copy(s.A), s.B/g, copy(s.C), s.D/g)...)
/{T1<:Real, T2<:Real}(g::T2, s::ContinuousSsSiso{T1}) = *(g, inv(s))

function /{T1<:Real, T2<:Real}(s::ContinuousSsSiso{T1}, g::Matrix{T2})
  ginv = inv(g)

  ContinuousSsSiso(promote(copy(s.A), s.B*ginv, copy(s.C), s.D*ginv)...)
end
/{T1<:Real, T2<:Real}(g::Matrix{T2}, s::ContinuousSsSiso{T1}) = *(g, inv(s))

./{T1<:Real, T2<:Real}(s1::ContinuousSsSiso{T1}, s2::ContinuousSsSiso{T2}) = /(s1, s2)
./{T1<:Real, T2<:Real}(s::ContinuousSsSiso{T1}, g::T2)  = /(s, g)
./{T1<:Real, T2<:Real}(g::T2, s::ContinuousSsSiso{T1})  = /(g, s)

function =={T1<:Real, T2<:Real}(s1::ContinuousSsSiso{T1}, s2::ContinuousSsSiso{T2})
  # TODO: Implement
  throw(ErrorException("==(s1,s2) for ContinuousSsSiso is not implemented"))
end

!={T1<:Real, T2<:Real}(s1::ContinuousSsSiso{T1}, s2::ContinuousSsSiso{T2}) = !(s1 == s2)

function isapprox{T1<:Real, T2<:Real}(s1::ContinuousSsSiso{T1},
  s2::ContinuousSsSiso{T2})
  # TODO: Implement
  throw(ErrorException("isapprox(s1,s2) ContinuousSsSiso is not implemented"))
end
