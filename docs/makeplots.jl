using Plots, SpecialFunctions, LaTeXStrings
pyplot()

include("src/branchesplots.jl")
include("src/asymptotics.jl")
include("src/hankelh1nplots.jl")

function makeplots()

    println("Generating plots")

    make_asymptotic_plots()

    make_branches_plots()

    make_hankelh1n_plots()

end
