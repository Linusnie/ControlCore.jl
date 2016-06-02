module ControlCore

# Import mathematical operations for overloading
import Base: +, -, *, /, .+, .-, .*, ./, ==, !=

# Import conversion and promotion functions for overloading
import Base: convert, promote_rule

# Import identities for overloading
import Base: one, zero, zeros

# Import helper functions
import Base: ndims, size, getindex, showcompact, show, showall, summary

# Export only the useful functions
export
  tf,
  ss,
  series,
  parallel,
  feedback

# Polynomials package is needed
using Polynomials

include("abstract.jl")
include("continuousss.jl")

end # module
