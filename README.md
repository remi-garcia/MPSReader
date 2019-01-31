# MPSReader

MPSReader is a Julia package that transform a .mps file (MPS format) into a model using JuMP.

The MPS format represents problem of the following form :
<p align="center">
optimize &nbsp; cᵀ x + c₀
&nbsp;&nbsp;
subject to &nbsp; L ≤ Ax ≤ U
&nbsp;&nbsp; l ≤ x ≤ u
</p>

optimize stands for minimize or maximize, if the objective sense is not defined MPSReader will use minimize as default. c₀ ∈ ℝ is a constant term and c ∈ ℝⁿ is the cost function.

L and U are vectors of constraints bounds. By default the lower bound of a constraint is -Inf and the upper bound is +Inf.
l and u are vectors of bounds of x. By default the lower bound of a variable is 0 and the upper bound is +Inf.
x can have simultaneously binary elements, integer elements and continuous elements. Semi-continuous variables are not supported and using MPSReader on a file containing it will throw an error.

This package exports 2 functions: `readmps()` and `mpstomodel()`

### Prerequisites

julia 1.0.0
MathOptInterface v0.8.1
JuMP v0.18.5+ [`~/.julia/dev/JuMP`]
LinearAlgebra

This package is using the development version of JuMP, it will be JuMP v0.19.0. If you are using julia 1.0.0 or newer you can install the needed package as follow:
```julia
julia> using Pkg
julia> Pkg.add("MathOptInterface")
julia> Pkg.add("LinearAlgebra")
julia> Pkg.develop("JuMP")
```

### Installing

TODO

## Running the tests

TODO

### Usage

```julia
julia> varTypes, bounds, objsense, c, c0, A, b, conTypes = readmps("mympsfile.mps")
julia> myModel, myVariables = mpstomodel("mympsfile.mps", mySolver) # mySolver = GLPK
```

`varTypes` is a vector of booleans, elements are true if the variable is continuous and at false if it is an integer or a boolean. `bounds` is a vector of pairs, the first is the lower bound of the variable and the second the upper bound. `objsense` is a Symbol for the optimizing sense, (`:min` by default). `conTypes` is a vector that gives the type of the constraint, 1 for ≤, 2 for = and 3 for >=. Each constraint having a lower and an upper bound is duplicated into a constraint with a lower bound and another with an upper bound.


### Problem Collections

* MIPLIB: https://miplib.zib.de/

### References

* MPS file format: http://lpsolve.sourceforge.net/5.5/mps-format.htm
