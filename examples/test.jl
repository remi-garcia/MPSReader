using GLPK

include("../src/MPSReader.jl")

files = [
    "lpsolve.mps",
    "min4m3.mps",
    "min4m3i.mps",
    "min4m3b.mps",
    "min4m3mip.mps",
    "min4m3c0.mps",
    "max4m3.mps"
]
for file in files
    println("--- "*file*" ---")
    myModel, myVariables = mpstomodel("data/"*file, GLPK.Optimizer)
    @info myModel
    println()
end






#
