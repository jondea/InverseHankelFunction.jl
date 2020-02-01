using Plots, SpecialFunctions, LaTeXStrings
pyplot()

include("src/branchesplots.jl")
include("src/asymptotics.jl")

function makeplots()

    println("Generating plots")

    make_asymptotic_plots()

    make_branches_plots()

end
