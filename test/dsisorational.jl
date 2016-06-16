println("Starting discrete SISO rational transfer function type tests...")

using Polynomials

# Construction
s1 = tf([1, 2], [1, 2, 1], 1)
s2 = tf([1, 2], [1, 2, 1], Float64(1))
s3 = tf(Poly([Float64(2), 1.0]), Poly([1.0, Float64(2), 1]), 1)
s4 = tf(5.0, 1)
@test typeof(s1)     == ControlCore.DSisoRational{Int,Int,Int}
@test typeof(s2)     == ControlCore.DSisoRational{Int,Int,Int}
@test typeof(s3)     == ControlCore.DSisoRational{Float64,Float64,Float64}
@test typeof(s4)     == ControlCore.DSisoRational{Float64,Float64,Float64}

for s in [s1,s2,s3]
# I/O mapping
  @test numstates(s)    == 3
  @test numinputs(s)    == 1
  @test numoutputs(s)   == 1

# Dimension information
  @test ndims(s)        == 1
  @test size(s)         == 1
  @test size(s, 2)      == 1
  @test size(s, 1, 2)   == (1,1)

# Iteration interface
  @test s[1]            == s
  @test s[:]            == s
  @test_throws BoundsError s[2]

# poles and zeros
  @test zeros(s)         ≈ Array{Complex{eltype(s.num)},1}([-2+0im])
  @test poles(s)         ≈ Array{Complex{eltype(s.den)},1}([-1+0im, -1+0im])

# return vectors
  @test numvec(s)       == [1, 2]
  @test denvec(s)       == [1, 2, 1]
  @test numpoly(s)      == Poly([2, 1])
  @test denpoly(s)      == Poly([1, 2, 1])
  @test zpkdata(s)      == (zeros(s), poles(s), s.num[1]/s.den[1])
  @test samplingtime(s) == s.Ts
end

# Addition
