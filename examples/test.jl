#using GLPK

include("../src/MPSReader.jl")

fixed = true

files = Dict{String, Bool}([
    "lpsolve.mps" => fixed,
    "min4m3.mps" => fixed,
    "min4m3i.mps" => fixed,
    "min4m3b.mps" => fixed,
    "min4m3mip.mps" => fixed,
    "min4m3c0.mps" => fixed,
    "max4m3.mps" => !fixed
])
for file in keys(files)
    println("--- "*file*" ---")
    myModel, myVariables = mpstomodel("data/"*file, files[file])
    @info myModel
    #JuMP.optimize!(myModel, with_optimizer(GLPK.Optimizer))
    #println(JuMP.termination_status(myModel))
    #println(JuMP.objective_value(myModel))
    println()
end






#
