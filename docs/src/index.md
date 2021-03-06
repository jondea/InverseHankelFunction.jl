# Inverse Hankel Function

Provides a function which finds `z` in the equation `hankelh1(ν,z) = h` for a
given `h`.

[![Build Status](https://travis-ci.org/jondea/InverseHankelFunction.jl.svg?branch=master)](https://travis-ci.org/jondea/InverseHankelFunction.jl)
[![Coverage Status](https://coveralls.io/repos/jondea/InverseHankelFunction.jl/badge.svg?branch=master)](https://coveralls.io/r/jondea/InverseHankelFunction.jl?branch=master)
[![codecov.io](http://codecov.io/github/jondea/InverseHankelFunction.jl/coverage.svg?branch=master)](http://codecov.io/github/jondea/InverseHankelFunction.jl?branch=master)

## Get started
To install this package, call
```julia
import Pkg
Pkg.add("https://github.com/jondea/InverseHankelFunction.jl")
```
or alternatively type `] add https://github.com/jondea/InverseHankelFunction.jl`
in the REPL.

## Branches

For a given `h` there are typically many solutions to the equation, so to define
a single valued function, we take two approaches:
* Define a "normalised" Hankel function `hbar(z) = h(z)/h(z_0)`, and analytically
  continue our inverse from the point `z_0`.
  This is currently the best studied and most completely implemented approach,
  and we discuss it [here](hankelh1n).
  In a related way, we also define the inverse Hankel function which "passes
  through" `z_0` using the interface `invhankelh1(ν, h, PassingThrough(z_0))`.
* A more general approach is to define a branch index (which we denote as `b`)
  and find a way to ennumerate them.
  This approach is less well developed, and we discuss it [here](branches).
