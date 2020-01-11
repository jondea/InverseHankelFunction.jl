# Inverse Hankel Function

Provides a function which finds `z` in the equation `hankelh1(ν,z) = h` for a
given `h`.

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://jondea.github.io/InverseHankelFunction.jl/dev)
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
For a given `h` there are typically many solutions to this equation, so to define
a single valued function, we take two approaches:
* Define a "normalised" Hankel function hbar(z) = h(z)/h(z_0), and analytically
  continue our inverse from the point `z_0`
* Define a branch index (which we denote as `b`) which is 
