# Normalised Hankel Function of First Kind and its Inverse

One of the ways to define a single valued inverse of the Hankel function is to
consider the inverse around a specific known point (`z₀`).
To do this, we first define the normalised Hankel function, that is
```julia
hankelh1n(ν, z₀, z) = hankelh1(ν, z) / hankelh1(ν, z₀)

```

See below for an example of this function plotted below (argument of function across complex domain)
for `ν = 0, z₀ = 1`.
![](../plots/hankelh1n_nu_0_z_0_1_complex_arg.png)

In white, we have plotted the path which the inverse of the this function takes in complex space.
The same function is plotted below, with the real and imaginary parts plotted against the argument of the inverse function.

![](../plots/invhankelh1n_nu_0_z_0_1_real_arg.svg)

We repeat these two plots for `ν = 3, z₀ = 2`.
![](../plots/hankelh1n_nu_3_z_0_2_complex_arg.png)

![](../plots/invhankelh1n_nu_3_z_0_2_real_arg.svg)

From [Deakin (2020)](https://www.research.manchester.ac.uk/portal/en/theses/optimal-pml-transformations-for-the-helmholtz-equation(2617fdfb-06e9-4fbf-9bfc-934f6b361572).html) we know that the inverse of the normalised Hankel function is single valued for 
all real positive `z₀` and `ν` and ` 0 <=  z <= 1`.
Also we know that `hankelh1n(ν, z₀, z) → i∞ as z → 0`

