# Some Thoughts on the Type Hierarchy

The type system will look like this:

- abstract LtiSystem
  - abstract SisoSystem
    - abstract SisoTf
      - abstract CSisoTf
        - ContinuousRational{T<:AbstractFloat}
        - ContinuousZpk{T<:AbstractFloat}
      - abstract DSisoTf
        - DiscreteRational{T<:AbstractFloat}
        - DiscreteZpk{T<:AbstractFloat}
    - abstract SisoSs
      - ContinuousSsSiso{T<:AbstractFloat}
      - DiscreteSsSiso{T<:AbstractFloat}
  - abstract MimoSystem
    - MimoTf{S<:SisoSystem}
    - MimoSs{T<:AbstractFloat}

# Interface Requirements

For **SisoSystem**s, the following functions are required:

- numpoly{S<:SisoSystem}(sys::S) should return Poly{T}
- denpoly{S<:SisoSystem}(sys::S) should return Poly{T}
- numvec{S<:SisoSystem}(sys::S) should return Vector{T}
- denvec{S<:SisoSystem}(sys::S) should return Vector{T}
- zpkdata{S<:SisoSystem}(sys::S) should return a tuple of (z, p, k)
- poles{S<:SisoSystem}(sys::S) should return roots(denpoly(sys))
- zeros{S<:SisoSystem}(sys::S) should return roots(denpoly(sys))

for doing, for example, basic mathematical operations among different
types of transfer functions.

For **MimoSystem**s, we need to have:

- getmatrix{S<:MimoSystem}(sys::S) should return the input-output mapping (mat)
- ...

for functions to be defined later on. Then, we can have mappings from the
corresponding Siso types directly to the Mimo versions.

We need to require the **LtiSystem** to have:

- numstates{S<:LtiSystem}(sys::S) should return nx in Int
- numinputs{S<:LtiSystem}(sys::S) should return nu in Int
- numoutputs{S<:LtiSystem}(sys::S) should return ny in Int
- isdiscrete(sys::LtiSystem) should return true/false for respective systems
- samplingtime{S<:LtiSystem}(sys::S) should return zero(T) for continuous and
  gain-only discrete time systems, and Ts::T for the rest of the discrete time
  systems.
- eltype
- start
- next
- done
- ...

for iteration/slicing/etc. functions to work properly.
