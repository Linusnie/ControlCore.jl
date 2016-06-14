include("ControlCore.jl")
using ControlCore
using Polynomials

s1 = zpk(sparse([1],[1],[5.]), 1)

for sub in s1
  println("hej")
end

s2 = zpk(Array{Float64,2 }(0,0), 1)

denpoly(s1)
