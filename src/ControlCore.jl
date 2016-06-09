module ControlCore

# Import conversion and promotion functions for overloading
import Base: convert, promote_rule

# Import identities for overloading
import Base: one, zero

# Import inv and zeros
import Base: inv, zeros

# Import slicing functions
import Base: ndims, size, getindex

# Import iteration interface functions
import Base: start, next, done, eltype, length, eachindex, endof

# Import printing functions
import Base: showcompact, show, showall, summary

# Import mathematical operations for overloading
import Base: +, .+, -, .-, *, .*, /, ./, ==, !=, isapprox

# Export only the useful functions
export
  tf,
  zpk,
  ss,
  series,
  parallel,
  feedback,
  degree,
  poles,
  numvec,
  denvec,
  numpoly,
  denpoly

# Polynomials package is needed
using Polynomials

include("abstract.jl")
include("display.jl")
include("discretesisorational.jl")
include("discretesisozpk.jl")
include("discretemimo.jl")
include("continuoussisorational.jl")
include("continuoussisozpk.jl")
include("continuousmimo.jl")
include("continuoussssiso.jl")
include("continuousssmimo.jl")

end # module
