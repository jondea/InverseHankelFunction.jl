using Documenter, InverseHankelFunction

include("makeplots.jl")

makedocs(
    modules=[InverseHankelFunction],
    format=Documenter.HTML(),
    sitename="InverseHankelFunction.jl",
    pages=[
        "Home" => "index.md",
        "Normalised Hankel Function" => "hankelh1n.md",
        "Asymptotics" => "asymptotics.md",
        "Branches" => "branches.md",
        "Functions" => "functions.md",
    ]
)

if get(ENV, "TRAVIS", "") == ""
    makeplots()
end

# Only build plots in travis if we are deploying
# And dont install the dependencies unless we are deploying
function myDeps()
    if get(ENV, "TRAVIS", "") != ""
        println("Installing deploy dependencies")
        makeplots()
    end
end

deploydocs(
    repo = "github.com/jondea/InverseHankelFunction.jl.git",
    deps = myDeps
)
