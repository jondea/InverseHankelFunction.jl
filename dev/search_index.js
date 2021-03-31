var documenterSearchIndex = {"docs":
[{"location":"functions/#Functions-1","page":"Functions","title":"Functions","text":"","category":"section"},{"location":"functions/#Normalised-Hankel-Function-1","page":"Functions","title":"Normalised Hankel Function","text":"","category":"section"},{"location":"functions/#","page":"Functions","title":"Functions","text":"hankelh1n\nHankelH1N","category":"page"},{"location":"functions/#InverseHankelFunction.hankelh1n","page":"Functions","title":"InverseHankelFunction.hankelh1n","text":"hankelh1n(ν, z₀ [,z])\n\nHankel function normalised at some z₀, H^{(1)}\\nu(z)/H^{(1)}\\nu(z₀), returns a function if z is not provided\n\n\n\n\n\n","category":"function"},{"location":"functions/#InverseHankelFunction.HankelH1N","page":"Functions","title":"InverseHankelFunction.HankelH1N","text":"Struct to cache the value of the Hankel function of the first kind at z₀, can be called like a function after construction\n\n\n\n\n\n","category":"type"},{"location":"functions/#Derivatives-of-Hankel-function-1","page":"Functions","title":"Derivatives of Hankel function","text":"","category":"section"},{"location":"functions/#","page":"Functions","title":"Functions","text":"diffhankelh1\ndiff2hankelh1\ndiff3hankelh1","category":"page"},{"location":"functions/#InverseHankelFunction.diffhankelh1","page":"Functions","title":"InverseHankelFunction.diffhankelh1","text":"diffhankelh1(ν, z [, h[, hm1]])\n\nReturns the first derivative of hankelh1 with respect to z. Can optionally pass in a precomputed H_ν(z) and H_{ν-1}(z).\n\n\n\n\n\n","category":"function"},{"location":"functions/#InverseHankelFunction.diff2hankelh1","page":"Functions","title":"InverseHankelFunction.diff2hankelh1","text":"diff2hankelh1(ν, z [, h[, hm1]])\n\nReturns the first and second derivatives of hankelh1 with respect to z as tuple Can optionally pass in a precomputed H_ν(z) and H_{ν-1}(z).\n\n\n\n\n\n","category":"function"},{"location":"functions/#InverseHankelFunction.diff3hankelh1","page":"Functions","title":"InverseHankelFunction.diff3hankelh1","text":"diff3hankelh1(ν, z [, h[, hm1]])\n\nReturns the first, second and third derivatives of hankelh1 with respect to z as tuple Can optionally pass in a precomputed H_ν(z) and H_{ν-1}(z).\n\n\n\n\n\n","category":"function"},{"location":"functions/#Inverse-of-Normalised-Hankel-Function-1","page":"Functions","title":"Inverse of Normalised Hankel Function","text":"","category":"section"},{"location":"functions/#","page":"Functions","title":"Functions","text":"invhankelh1n\ninvhankelh1n_sortedvec\ndiffinvhankelh1n_sortedvec\ninvhankelh1n_adaptive_solve","category":"page"},{"location":"functions/#InverseHankelFunction.invhankelh1n","page":"Functions","title":"InverseHankelFunction.invhankelh1n","text":"invhankelh1n(ν::Integer, z₀::Number, h̄::Number)\n\nFind z such that H^(1)_nu(z)H^(1)_nu(z₀) = barh for the branch continued from z₀ See also: hankelh1n\n\n\n\n\n\n","category":"function"},{"location":"functions/#Asymptotics-of-Hankel-function-1","page":"Functions","title":"Asymptotics of Hankel function","text":"","category":"section"},{"location":"functions/#","page":"Functions","title":"Functions","text":"hankel_arg_asymptotic_scale\nsmall_arg_hankelh1\ninv_small_arg_hankelh1\nlarge_arg_hankelh1\ninv_large_arg_hankelh1","category":"page"},{"location":"functions/#InverseHankelFunction.hankel_arg_asymptotic_scale","page":"Functions","title":"InverseHankelFunction.hankel_arg_asymptotic_scale","text":"hankel_arg_asymptotic_scale(ν, z)\n\nDetermines largeness or smallness of argument to determine the validity of the asymptotic forms of the hankel functions.\n\n\n\n\n\n","category":"function"},{"location":"functions/#InverseHankelFunction.small_arg_hankelh1","page":"Functions","title":"InverseHankelFunction.small_arg_hankelh1","text":"small_arg_hankelh1(ν::Integer, z::Number) -> Complex\n\nThe small argument (z) asymptotic of hankelh1\n\n\n\n\n\n","category":"function"},{"location":"functions/#InverseHankelFunction.inv_small_arg_hankelh1","page":"Functions","title":"InverseHankelFunction.inv_small_arg_hankelh1","text":"inv_small_arg_hankelh1(ν::Integer, h::Complex, b::Integer) -> Complex\n\nThe bth branch of the inverse of the small argument (z) asymptotic of hankelh1 Not valid for b > ν\n\n\n\n\n\n","category":"function"},{"location":"functions/#InverseHankelFunction.large_arg_hankelh1","page":"Functions","title":"InverseHankelFunction.large_arg_hankelh1","text":"large_arg_hankelh1(ν::Integer, z::Number) -> Complex\n\nThe large argument (z) asymptotic of hankelh1\n\n\n\n\n\n","category":"function"},{"location":"functions/#InverseHankelFunction.inv_large_arg_hankelh1","page":"Functions","title":"InverseHankelFunction.inv_large_arg_hankelh1","text":"inv_large_arg_hankelh1(ν::Integer, h::Complex, b::Integer) -> Complex\n\nThe bth branch of the inverse of the large argument (z) asymptotic of hankelh1\n\n\n\n\n\n","category":"function"},{"location":"branches/#Inverse-Hankel-function-branches-1","page":"Branches","title":"Inverse Hankel function branches","text":"","category":"section"},{"location":"branches/#","page":"Branches","title":"Branches","text":"To understand the branches in the inverse of the Hankel function, we look at the behaviour of several Hankel functions near the origin. The white contours are constant argument, and the black lines are constant magnitude. We can tell that there are multiple branches in all of these cases. For example, in the first plot, following a line of constant magnitude in the upper half, we can see multiple points with the same argument. These points are multiple solutions (z) for a particular value of the the Hankel function (h).","category":"page"},{"location":"branches/#","page":"Branches","title":"Branches","text":"As we increase the order of the Hankel function (which we will denote by m), we see the small argument (z≈0) (asymptotics)[asymptotics] become clearer, with a zero of multiplicity equal to m. Therefore around z=0 we see that there will be at least m branches of the inverse of the Hankel function. This motivates us to start our enumeration of our branches here, indexed by b","category":"page"},{"location":"branches/#Contour-plots-of-constant-argument-and-magnitude-for-Hankel-functions-against-complex-z-1","page":"Branches","title":"Contour plots of constant argument and magnitude for Hankel functions against complex z","text":"","category":"section"},{"location":"branches/#th-order-Hankel-function-1","page":"Branches","title":"0th order Hankel function","text":"","category":"section"},{"location":"branches/#","page":"Branches","title":"Branches","text":"(Image: )","category":"page"},{"location":"branches/#st-order-Hankel-function-1","page":"Branches","title":"1st order Hankel function","text":"","category":"section"},{"location":"branches/#","page":"Branches","title":"Branches","text":"(Image: )","category":"page"},{"location":"branches/#nd-order-Hankel-function-1","page":"Branches","title":"2nd order Hankel function","text":"","category":"section"},{"location":"branches/#","page":"Branches","title":"Branches","text":"(Image: )","category":"page"},{"location":"branches/#rd-order-Hankel-function-1","page":"Branches","title":"3rd order Hankel function","text":"","category":"section"},{"location":"branches/#","page":"Branches","title":"Branches","text":"(Image: )","category":"page"},{"location":"branches/#th-order-Hankel-function-2","page":"Branches","title":"8th order Hankel function","text":"","category":"section"},{"location":"branches/#","page":"Branches","title":"Branches","text":"(Image: )","category":"page"},{"location":"hankelh1n/#Normalised-Hankel-Function-of-First-Kind-and-its-Inverse-1","page":"Normalised Hankel Function","title":"Normalised Hankel Function of First Kind and its Inverse","text":"","category":"section"},{"location":"hankelh1n/#","page":"Normalised Hankel Function","title":"Normalised Hankel Function","text":"One of the ways to define a single valued inverse of the Hankel function is to consider the inverse around a specific known point (z₀). To do this, we first define the normalised Hankel function, that is","category":"page"},{"location":"hankelh1n/#","page":"Normalised Hankel Function","title":"Normalised Hankel Function","text":"hankelh1n(ν, z₀, z) = hankelh1(ν, z) / hankelh1(ν, z₀)\n","category":"page"},{"location":"hankelh1n/#","page":"Normalised Hankel Function","title":"Normalised Hankel Function","text":"See below for an example of this function plotted below (argument of function across complex domain) for ν = 0, z₀ = 1. (Image: )","category":"page"},{"location":"hankelh1n/#","page":"Normalised Hankel Function","title":"Normalised Hankel Function","text":"In white, we have plotted the path which the inverse of the this function takes in complex space. The same function is plotted below, with the real and imaginary parts plotted against the argument of the inverse function.","category":"page"},{"location":"hankelh1n/#","page":"Normalised Hankel Function","title":"Normalised Hankel Function","text":"(Image: )","category":"page"},{"location":"hankelh1n/#","page":"Normalised Hankel Function","title":"Normalised Hankel Function","text":"We repeat these two plots for ν = 3, z₀ = 2. (Image: )","category":"page"},{"location":"hankelh1n/#","page":"Normalised Hankel Function","title":"Normalised Hankel Function","text":"(Image: )","category":"page"},{"location":"hankelh1n/#","page":"Normalised Hankel Function","title":"Normalised Hankel Function","text":"From Deakin (2020) we know that the inverse of the normalised Hankel function is single valued for  all real positive z₀ and ν and 0 <=  z <= 1. Also we know that hankelh1n(ν, z₀, z) → i∞ as z → 0","category":"page"},{"location":"#Inverse-Hankel-Function-1","page":"Home","title":"Inverse Hankel Function","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"Provides a function which finds z in the equation hankelh1(ν,z) = h for a given h.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"(Image: Build Status) (Image: Coverage Status) (Image: codecov.io)","category":"page"},{"location":"#Get-started-1","page":"Home","title":"Get started","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"To install this package, call","category":"page"},{"location":"#","page":"Home","title":"Home","text":"import Pkg\nPkg.add(\"https://github.com/jondea/InverseHankelFunction.jl\")","category":"page"},{"location":"#","page":"Home","title":"Home","text":"or alternatively type ] add https://github.com/jondea/InverseHankelFunction.jl in the REPL.","category":"page"},{"location":"#Branches-1","page":"Home","title":"Branches","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"For a given h there are typically many solutions to the equation, so to define a single valued function, we take two approaches:","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Define a \"normalised\" Hankel function hbar(z) = h(z)/h(z_0), and analytically continue our inverse from the point z_0. This is currently the best studied and most completely implemented approach, and we discuss it (here)[hankelh1n].\nA more general approach is to define a branch index (which we denote as b) and find a way to ennumerate them. This approach is less well developed, and we discuss it (here)[branches].","category":"page"},{"location":"asymptotics/#Hankel-function-asymptotics-1","page":"Asymptotics","title":"Hankel function asymptotics","text":"","category":"section"},{"location":"asymptotics/#","page":"Asymptotics","title":"Asymptotics","text":"We make use of asymptotic forms (approximations) of the Hankel functions, which are valid for small and large argument.","category":"page"},{"location":"asymptotics/#","page":"Asymptotics","title":"Asymptotics","text":"[Maths here]","category":"page"},{"location":"asymptotics/#","page":"Asymptotics","title":"Asymptotics","text":"Reference: Abramowitz and Stegun/ Digital Library of Mathematical Functions.","category":"page"},{"location":"asymptotics/#","page":"Asymptotics","title":"Asymptotics","text":"(Image: )","category":"page"},{"location":"asymptotics/#","page":"Asymptotics","title":"Asymptotics","text":"(Image: )","category":"page"},{"location":"asymptotics/#","page":"Asymptotics","title":"Asymptotics","text":"(Image: )","category":"page"}]
}
