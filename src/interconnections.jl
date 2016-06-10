# Series interconnection
series{T1, T2}(s1::T1, s2::T2)                            = *(promote(s1,s2)...)
series{T1<:LtiSystem, T2<:LtiSystem}(s1::T1, s2::T2)      = *(promote(s1,s2)...)
series{T1<:LtiSystem}(s1::T1, s2::T1)                     = *(s1, s2)

series{T1<:SisoSystem, T2<:Real}(s1::T1, n::T2)           = *(s1, n)
series{T1<:SisoSystem, T2<:Real}(n::T2, s1::T1)           = *(n ,s1)

series{T1<:SisoSystem, T2<:Real}(s1::T1, n::Matrix{T2})   = *(s1, n)
series{T1<:SisoSystem, T2<:Real}(n::Matrix{T2}, s1::T1)   = *(n ,s1)

series{T1<:MimoSystem, T2<:Real}(s1::T1, n::T2)           = *(s1, n)
series{T1<:MimoSystem, T2<:Real}(n::T2, s1::T1)           = *(n ,s1)

series{T1<:MimoSystem, T2<:Real}(s1::T1, n::Matrix{T2})   = *(s1, n)
series{T1<:MimoSystem, T2<:Real}(n::Matrix{T2}, s1::T1)   = *(n ,s1)

# Paralell interconnection
paralell{T1, T2}(s1::T1, s2::T2)                          = +(promote(s1,s2)...)
paralell{T1<:LtiSystem, T2<:LtiSystem}(s1::T1, s2::T2)    = +(promote(s1,s2)...)
paralell{T1<:LtiSystem}(s1::T1, s2::T1)                   = *(s1, s2)

paralell{T1<:SisoSystem, T2<:Real}(s1::T1, n::T2)         = *(s1, n)
paralell{T1<:SisoSystem, T2<:Real}(n::T2, s1::T1)         = *(n ,s1)

paralell{T1<:SisoSystem, T2<:Real}(s1::T1, n::Matrix{T2}) = *(s1, n)
paralell{T1<:SisoSystem, T2<:Real}(n::Matrix{T2}, s1::T1) = *(n ,s1)

paralell{T1<:MimoSystem, T2<:Real}(s1::T1, n::T2)         = *(s1, n)
paralell{T1<:MimoSystem, T2<:Real}(n::T2, s1::T1)         = *(n ,s1)

paralell{T1<:MimoSystem, T2<:Real}(s1::T1, n::Matrix{T2}) = *(s1, n)
paralell{T1<:MimoSystem, T2<:Real}(n::Matrix{T2}, s1::T1) = *(n ,s1)

# Feedback interconnection
function feedback{T1, T2}(s1::T1, s2::T2)
  s1_,s2_ = promote(s1,s2)
  /(s1_,series(one(T2),s2_))
end
function feedback{T1<:LtiSystem, T2<:LtiSystem}(s1::T1, s2::T2)
  s1_,s2_ = promote(s1,s2)
  /(s1,series(one(T2),s2))
end
function feedback{T1<:LtiSystem}(s1::T1, s2::T1)
  /(s1,series(one(T1),s2))
end

feedback{T1<:SisoSystem, T2<:Real}(s1::T1, n::T2)         =
  /(s1, series(one(T2), n))
feedback{T1<:SisoSystem, T2<:Real}(n::T2, s1::T1)         =
  /(n, series(one(T1), s1))

feedback{T1<:SisoSystem, T2<:Real}(s1::T1, n::Matrix{T2}) =
  /(s1, series(one(T2), n))
feedback{T1<:SisoSystem, T2<:Real}(n::Matrix{T2}, s1::T1) =
  /(n, series(one(T1), s1))

feedback{T1<:MimoSystem, T2<:Real}(s1::T1, n::T2)         =
  /(s1, series(one(T2), n))
feedback{T1<:MimoSystem, T2<:Real}(n::T2, s1::T1)         =
  /(n, series(one(T1), s1))

feedback{T1<:MimoSystem, T2<:Real}(s1::T1, n::Matrix{T2}) =
  /(s1, series(one(T2), n))
feedback{T1<:MimoSystem, T2<:Real}(n::Matrix{T2}, s1::T1) =
  /(n, series(one(T1), s1))
