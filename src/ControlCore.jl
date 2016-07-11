module ControlCore

# Import conversion and promotion functions for overloading
import Base: convert, promote_rule

# Import identities for overloading
import Base: one, zero

# Import inv and zeros
import Base: inv, zeros

# Import slicing functions
import Base: ndims, size, getindex, setindex!

# Import iteration interface functions
import Base: start, next, done, eltype, length, eachindex, endof

# Import printing functions
import Base: showcompact, show, showall, summary

# Import mathematical operations for overloading
import Base: +, .+, -, .-, *, .*, /, ./, ==, !=, isapprox, transpose

# Export only the useful functions
export
  tf,
  zpk,
  ss,
  mimo,
  series,
  paralell,
  feedback,
  numstates,
  numinputs,
  numoutputs,
  getmatrix,
  degree,
  poles,
  numvec,
  denvec,
  numpoly,
  denpoly,
  zpkdata,
  samplingtime

# Polynomials package is needed
using Polynomials

include("abstract/ltisystem.jl")
include("abstract/sisosystem.jl")
include("abstract/sisotf.jl")
include("abstract/csisotf.jl")
include("abstract/dsisotf.jl")
include("abstract/sisoss.jl")
include("abstract/mimosystem.jl")
include("dsisorational.jl")
include("dsisozpk.jl")
include("dmimo.jl")
include("csisorational.jl")
include("csisozpk.jl")
include("cmimo.jl")
include("csisoss.jl")
include("dsisoss.jl")
include("cmimoss.jl")
include("dmimoss.jl")
include("display.jl")
include("interconnections.jl")
include("conversions.jl")

end # module
