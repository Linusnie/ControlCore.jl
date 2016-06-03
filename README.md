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
- [ ] `evalfreq`

**Interconnections**

- [ ] `series`,
- [ ] `parallel`,
- [ ] `feedback`.

[travis-ci-img]: https://travis-ci.org/KTH-AC/ControlCore.jl.svg?branch=master
[travis-ci-link]: https://travis-ci.org/KTH-AC/ControlCore.jl
[appveyor-ci-img]: https://ci.appveyor.com/api/projects/status/geqrrlwve5ycjh0a/branch/master?svg=true
[appveyor-ci-link]: https://ci.appveyor.com/project/aytekinar/controlcore-jl/branch/master
