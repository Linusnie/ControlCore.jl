# ControlCore

**Build Status**

-  Unix/OSX: [![Build Status][travis-ci-img]][travis-ci-link]
-  Windows: [![Build status][appveyor-ci-img]][appveyor-ci-link]

Core control systems functionality for **analysis**, **design** and
**identification** tools to be implemented later on.

This repository is meant to provide a basic set of tools, *i.e.*, *transfer
function* and *state space* types as well as the basic mathematical operations
defined on them, in a way that other tools for **analyzing** or **identifying**
control systems would use the same set of functionality in a transparent,
coherent way.

The planned functionality to include in this toolbox is (ticked boxes are
implemented):

**Creation of types**

- [ ] `tf`,
- [ ] `zpk`,
- [ ] `ss`.

**Conversions and promotions**

- [ ] `convert`,
- [ ] `promote_rule`.

**Identity overloading**

- [ ] `one`,
- [ ] `zero`,
- [ ] `inv`.

**Slicing functions**

- [ ] `ndims`,
- [ ] `size`,
- [ ] `getindex`.

**Iteration interface**

- [ ] `start`,
- [ ] `next`,
- [ ] `done`,
- [ ] `eltype`,
- [ ] `length`,
- [ ] `eachindex`,
- [ ] `endof`.

**Printing functions**

- [ ] `showcompact`,
- [ ] `show`,
- [ ] `showall`,
- [ ] `summary`.

**Basic operations**

- [ ] `+`,
- [ ] `.+`,
- [ ] `-`,
- [ ] `.-`,
- [ ] `*`,
- [ ] `.*`,
- [ ] `/`,
- [ ] `./`,
- [ ] `==`,
- [ ] `!=`,
- [ ] `isapprox`.

**Basic functionality**

- [ ] `degree`
- [ ] `zeros`
- [ ] `poles`
- [ ] `numvec`
- [ ] `denvec`
- [ ] `numpoly`
- [ ] `denpoly`

**Interconnections**

- [ ] `series`,
- [ ] `parallel`,
- [ ] `feedback`.

# Some Thoughts on the Type Hierarchy

The type system will look like this:

- `abstract LtiSystem`
  - `abstract SisoSystem{T<:Real}`
    - `abstract SisoTf{T<:AbstractFloat} <: SisoSystem{T}`
      - `abstract CSisoTf{T} <: SisoTf{T}`
        - `CSisoRational{T} <: CSisoTf{T}`
        - `CSisoZpk{T} <: CSisoTf{T}`
      - `abstract DSisoTf{T} <: SisoTf{T}`
        - `DSisoRational{T} <: DSisoTf{T}`
        - `DSisoZpk{T} <: DSisoTf{T}`
    - `abstract SisoSs{T} <: SisoSystem{T}`
      - `CSisoSs{T} <: SisoSs{T}`
      - `DSisoSs{T} <: SisoSs{T}`
  - `abstract MimoSystem`
    - `CMimo{T<:CSiso}`
    - `CMimoSs{T<:Real}`
    - `DMimo{T<:DSiso}`
    - `DMimoSs{T<:Real}`
- `CSiso = Union{CSisoTf,CSisoSs}`
- `DSiso = Union{DSisoTf,DSisoSs}`

# Interface Requirements

For `SisoSystem`s, the following functions are required:

- `numpoly{S<:SisoSystem}(sys::S)` should return `Poly{T}`,
- `denpoly{S<:SisoSystem}(sys::S)` should return `Poly{T}`,
- `numvec{S<:SisoSystem}(sys::S)` should return `Vector{T}`,
- `denvec{S<:SisoSystem}(sys::S)` should return `Vector{T}`,
- `zpkdata{S<:SisoSystem}(sys::S)` should return a tuple of `(z, p, k)`,
- `poles{S<:SisoSystem}(sys::S)` should return `roots(denpoly(sys))`,
- `zeros{S<:SisoSystem}(sys::S)` should return `roots(denpoly(sys))`.

for doing, for example, basic mathematical operations among different
types of transfer functions.

For `MimoSystem`s, we need to have:

- `getmatrix{S<:MimoSystem}(sys::S)` should return the input-output mapping,
- ...

for functions to be defined later on. Then, we can have mappings from the
corresponding Siso types directly to the Mimo versions.

We need to require the `LtiSystem` to have:

- `numstates{S<:LtiSystem}(sys::S)` should return `nx::Int`,
- `numinputs{S<:LtiSystem}(sys::S)` should return `nu::Int`,
- `numoutputs{S<:LtiSystem}(sys::S)` should return `ny::Int`,
- `isdiscrete(sys::LtiSystem)` should return `true`/`false` for respective
  systems,
- `samplingtime{S<:LtiSystem}(sys::S)` should return `zero(T)` for continuous
  and gain-only discrete time systems, and `Ts::T` for the rest of the discrete
  time systems,
- `eltype`,
- `start`,
- `next`,
- `done`,
- ...

for iteration/slicing/etc. functions to work properly.

[travis-ci-img]: https://travis-ci.org/KTH-AC/ControlCore.jl.svg?branch=master
[travis-ci-link]: https://travis-ci.org/KTH-AC/ControlCore.jl
[appveyor-ci-img]: https://ci.appveyor.com/api/projects/status/geqrrlwve5ycjh0a/branch/master?svg=true
[appveyor-ci-link]: https://ci.appveyor.com/project/aytekinar/controlcore-jl/branch/master
