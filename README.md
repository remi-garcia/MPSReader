# MPSReader

MPSReader is a Julia package that transform a .mps file (MPS format) into a model using JuMP.

The MPS format represents problem of the following form :</br>
optimize cᵀ x + c₀</br>
subject to L ≤ Ax ≤ U</br>
l ≤ x ≤ u</br>

optimize stands for minimize or maximize, MPSReader will use minimize as default but it can be specified in Free MPS Format. c₀ ∈ ℝ is a constant term and c ∈ ℝⁿ is the cost function.

L and U are vectors of constraints bounds. By default the lower bound of a constraint is -Inf and the upper bound is +Inf.
l and u are vectors of bounds of x. By default the lower bound of a variable is 0 and the upper bound is +Inf.
x can have simultaneously binary elements, integer elements and continuous elements. Semi-continuous variables are not supported and using MPSReader on a file containing it will throw an error.

This package exports 2 functions: `readmps()` and `mpstomodel()`

### Prerequisites

julia 1.0.0</br>
MathOptInterface v0.8.1</br>
JuMP v0.18.5+ [`~/.julia/dev/JuMP`]</br>
LinearAlgebra</br>
GLPK (only to run examples</br>

This package is using the development version of JuMP, it will be JuMP v0.19.0. If you are using julia 1.0.0 or newer you can install the needed package as follow:
```julia
julia> using Pkg
julia> Pkg.add("MathOptInterface")
julia> Pkg.add("LinearAlgebra")
julia> Pkg.develop("JuMP")
julia> Pkg.add("GLPK") #That's optional
```

### Installing

TODO

## Running the tests

TODO

### Usage

```julia
julia> fixed = true #false
julia> varTypes, bounds, objsense, c, c0, A, b, conTypes = readmps("mympsfile.mps", fixed)
julia> myModel, myVariables = mpstomodel("mympsfile.mps", mySolver, fixed) # mySolver = GLPK.Optimizer
```

`varTypes` is a vector of booleans, elements are true if the variable is continuous and at false if it is an integer or a boolean. `bounds` is a vector of pairs, the first is the lower bound of the variable and the second the upper bound. `objsense` is a Symbol for the optimizing sense, (`:min` by default). `conTypes` is a vector that gives the type of the constraint, 1 for ≤, 2 for = and 3 for ≥. Each constraint having a lower and an upper bound is duplicated into a constraint with a lower bound and another with an upper bound.

`myModel` is the model in minimization with constraints of type ≥ and =. Therefore it should have this form:</br>
minimize cᵀ x + c₀</br>
subject to Ax ≥ b</br>
           Aeq x ≥ beq</br>
l ≤ x ≤ u</br>

By default fixed is true but if you are using a Free MPS Format you can set fixed variable to false. As lpsolve authors we suggest to try the fixed format and if it fails switch for the free format.

### Problem Collections

* MIPLIB: https://miplib.zib.de/

### References

* MPS file format: http://lpsolve.sourceforge.net/5.5/mps-format.htm
