#module MPSReader

using JuMP
using MathOptInterface
using LinearAlgebra

#export readmps, mpstomodel

function mpstomodel(filename::String, fixed::Bool = true)
    varBin, varInt, varFloat, bounds, objsense, c, c0, Aeq, beq, Aineq, bineq = mpstomatrices(filename, fixed)

    nVar = length(varBin) + length(varInt) + length(varFloat)

    myModel = Model()
    @variable(myModel, x[1:nVar])
    if !isapprox(c0, 0)
        if objsense == :min
            @objective(myModel, Min, dot(c, x) + c0)
        else
            @objective(myModel, Max, dot(c, x) + c0)
        end
    else
        if objsense == :min
            @objective(myModel, Min, dot(c, x))
        else
            @objective(myModel, Max, dot(c, x))
        end
    end
    # Get type to model and bounds/fixed
    for i in varBin
        JuMP.set_binary(x[i])
        if bounds[i][1] == bounds[i][2]
            JuMP.fix(x[i], bounds[i][1])
        end
    end
    for i in varInt
        JuMP.set_integer(x[i])
        if bounds[i][1] == bounds[i][2]
            JuMP.fix(x[i], bounds[i][1])
        else
            if isfinite(bounds[i][1])
                JuMP.set_lower_bound(x[i], bounds[i][1])
            end
            if isfinite(bounds[i][2])
                JuMP.set_upper_bound(x[i], bounds[i][2])
            end
        end
    end
    for i in varFloat
        if bounds[i][1] == bounds[i][2]
            JuMP.fix(x[i], bounds[i][1])
        else
            if isfinite(bounds[i][1])
                JuMP.set_lower_bound(x[i], bounds[i][1])
            end
            if isfinite(bounds[i][2])
                JuMP.set_upper_bound(x[i], bounds[i][2])
            end
        end
    end

    # Add constraints
    m = length(bineq)
    if objsense == :min
        @constraint(myModel, [i = 1:m], dot(Aineq[i,:], x) >= bineq[i])
    else
        @constraint(myModel, [i = 1:m], dot(Aineq[i,:], x) <= bineq[i])
    end
    m = length(beq)
    @constraint(myModel, [i = 1:m], dot(Aeq[i,:], x) == beq[i])

    return myModel, x
end


function mpstomatrices(filename::String, fixed::Bool = true)
    varTypes, bounds, objsense, c, c0, A, b, conTypes = readmps(filename, fixed)

    # Get precise var type
    varBin = Vector{Int64}()
    varInt = Vector{Int64}()
    varFloat = Vector{Int64}()

    nVar = length(varTypes)
    nCon = length(conTypes)

    for i in 1:nVar
        isFloat = varTypes[i]
        if isFloat
            push!(varFloat, i)
        else
            if isapprox(0, bounds[i][1]) && isapprox(1, bounds[i][2])
                push!(varBin, i)
            else
                push!(varInt, i)
            end
        end
    end

    # Aineq x >= bineq in min
    ismin = 1
    # Aineq x <= bineq in max
    if objsense == :max
        ismin = -1
    end

    # Split A and b in equality and inequality
    Aeq = Matrix{Float64}(undef, 0, nVar)
    beq = Vector{Float64}()
    Aineq = Matrix{Float64}(undef, 0, nVar)
    bineq = Vector{Float64}()

    for i in 1:nCon
        if conTypes[i] == 2
            Aeq = vcat(Aeq, A[i,:]')
            push!(beq, b[i])
        elseif conTypes[i] == 1
            Aineq = vcat(Aineq, -ismin * A[i,:]')
            push!(bineq, -ismin * b[i])
        elseif conTypes[i] == 3
            Aineq = vcat(Aineq, ismin * A[i,:]')
            push!(bineq, ismin * b[i])
        else
            error("Should never happen")
        end
    end

    return varBin, varInt, varFloat, bounds, objsense, c, c0, Aeq, beq, Aineq, bineq
end


function readmps(filename::String, fixed::Bool = true)
    sections = Dict{String, Bool}([
        "NAME" => false,
        "ROWS" => false,
        "COLUMNS" => false,
        "RHS" => false,
        "BOUNDS" => false,
        "RANGES" => false,
        "SOS" => false,
        "ENDATA" => false
    ])
    if !fixed
        sections["OBJSENSE"] = false
        sections["OBJNAME"] = false
    end

    varNames = Dict{String, Int}()
    conNames = Dict{String, Int}()
    varTypes = BitArray{1}() # true : float   false : int
    conTypes = Vector{Int}() # 1 : <=, 2 : ==, 3 : >=

    # Default values
    name = ""
    objectiveName = ""
    objsense = :min
    nVar = 0
    nCon = 0
    A = Matrix{Float64}(undef,0,0)
    b = Vector{Float64}()
    c = Vector{Float64}()
    c0 = 0.0
    bounds = Vector{Tuple{Float64, Float64}}()

    mps = open(filename, "r")
    line = ""
    while !eof(mps) && !sections["ENDATA"]
        line = readline(mps)
        # Read the section name and avoid empty lines and comments
        while strip(line) == "" || line[1] == '*'
            line = readline(mps)
        end
        words = split(line)
        section = words[1]

        # Errors in sections names and order
        if !haskey(sections, section)
            error("Section ", section, " is not a valide section")
        elseif sections[section]
            error("Section ", section, " appears twice")
        elseif section == "COLUMNS" && !sections["ROWS"]
            error("ROWS must come before COLUMNS")
        elseif section == "RHS" && (!sections["COLUMNS"] || !sections["NAME"])
            error("NAME and COLUMNS must come before RHS")
        elseif section == "BOUNDS" && !sections["COLUMNS"]
            error("COLUMNS must come before BOUNDS")
        elseif section == "RANGES" && !sections["RHS"]
            error("RHS must come before RANGES")
        end

        # Read info
        if section == "NAME"
            name = words[2]
        elseif section == "OBJSENSE"
            objsense = readObjsense(mps)
        elseif section == "OBJNAME"
            objectiveName = readObjname(mps)
        elseif section == "ROWS"
            objectiveName, conNames, conTypes = readRows(mps, objectiveName, fixed)
            nCon = length(conTypes)
        elseif section == "COLUMNS"
            varNames, c, A, varTypes = readColumns(mps, objectiveName, conNames, fixed)
            nVar = length(keys(varNames))
            bounds = fill((0.0, Inf), nVar)
        elseif section == "RHS"
            b, c0 = readRhs(mps, objectiveName, conNames, fixed)
        elseif section == "RANGES"
            A, b, conTypes = readRanges(mps, A, b, conNames, conTypes, fixed)
            nCon = length(conTypes)
        elseif section == "BOUNDS"
            bounds, varTypes = readBounds(mps, varNames, varTypes, fixed)
        elseif section == "SOS"
            error("SOS section reader not ready yet")
        end
        sections[section] = true
    end
    if !sections["ENDATA"]
        error("No ENDATA section in the file")
    end
    close(mps)

    return varTypes, bounds, objsense, c, c0, A, b, conTypes
end


function readBounds(mps, varNames, varTypesinit, fixed)
    varTypes = copy(varTypesinit)
    pos = position(mps)
    bounds = fill((0.0, Inf), length(keys(varNames)))

    line = readline(mps)
    while strip(line) == "" || line[1] == '*'
        line = readline(mps)
    end
    while line[1] == ' '
        pos = position(mps)
        words = Vector{String}()

        if fixed
            lenLine = length(line)
            words = strip.([
                line[2:min(lenLine, 3)],
                line[5:min(lenLine, 12)],
                line[15:min(lenLine, 22)],
                line[25:min(lenLine, 36)]
            ])
        else
            words = split(line)
        end

        boundType = words[1]
        if !fixed
            boundType = uppercase(boundType)
        end
        name = words[3]

        if boundType == "FR"
            # Free variable
            bounds[varNames[name]] = (-Inf, Inf)
        elseif boundType == "MI"
            # Lower bound -inf
            if isfinite(bounds[varNames[name]][2])
                bounds[varNames[name]] = (-Inf, bounds[varNames[name]][2])
            else
                bounds[varNames[name]] = (-Inf, 0.0)
            end
        elseif boundType == "PL"
            # Default value
            bounds[varNames[name]] = (bounds[varNames[name]][1], Inf)
        elseif boundType == "BV"
            # Binary variable
            bounds[varNames[name]] = (0.0, 1.0)
            varTypes[varNames[name]] = false
        elseif boundType == "SC"
            # Semi-continuous variable
            error("SC bound not ready yet")
        elseif boundType == "LO"
            # Lower bound
            boundVal = parse(Float64, words[4])
            bounds[varNames[name]] = (boundVal, bounds[varNames[name]][2])
        elseif boundType == "UP"
            # Upper bound
            boundVal = parse(Float64, words[4])
            bounds[varNames[name]] = (bounds[varNames[name]][1], boundVal)
        elseif boundType == "FX"
            # Fixed variable
            boundVal = parse(Float64, words[4])
            bounds[varNames[name]] = (boundVal, boundVal)
        elseif boundType == "LI"
            # Integer variable / Lower bound
            boundVal = parse(Float64, words[4])
            bounds[varNames[name]] = (boundVal, bounds[varNames[name]][2])
            varTypes[varNames[name]] = false
        elseif boundType == "UI"
            # Integer variable / Upper bound
            boundVal = parse(Float64, words[4])
            bounds[varNames[name]] = (bounds[varNames[name]][1], boundVal)
            varTypes[varNames[name]] = false
        else
            error("Unknown bound type ", boundType)
        end

        line = readline(mps)
        while strip(line) == "" || line[1] == '*'
            pos = position(mps)
            line = readline(mps)
        end
    end

    seek(mps, pos)
    return bounds, varTypes
end


function readRanges(mps, Ainit, binit, conNames, conTypesinit, fixed)
    A = copy(Ainit)
    b = copy(binit)
    conTypes = copy(conTypesinit)
    pos = position(mps)

    line = readline(mps)
    while strip(line) == "" || line[1] == '*'
        line = readline(mps)
    end
    while line[1] == ' '
        pos = position(mps)
        words = Vector{String}()

        if fixed
            lenLine = length(line)
            words = strip.([
                line[5:min(lenLine, 12)],
                line[15:min(lenLine, 22)],
                line[25:min(lenLine, 36)],
                line[40:min(lenLine, 47)],
                line[50:min(lenLine, 61)]
            ])
            while words[end] == ""
                pop!(words)
            end
        else
            words = split(line)
        end
        lenWords = length(words)

        for i in 2:2:lenWords
            name = words[i]
            if name != ""
                if conTypes[conNames[name]] == 2
                    error("Range on an equality constraint")
                elseif conTypes[conNames[name]] == 1
                    push!(conTypes, 3)
                    b0 = -abs(parse(Float64, words[i+1]))
                elseif conTypes[conNames[name]] == 3
                    push!(conTypes, 1)
                    b0 = abs(parse(Float64, words[i+1]))
                else
                    error("Should never happen")
                end
                A = vcat(A, A[conNames[name],:]')
                push!(b, b0)
            end
        end

        line = readline(mps)
        while strip(line) == "" || line[1] == '*'
            pos = position(mps)
            line = readline(mps)
        end
    end

    seek(mps, pos)
    return A, b, conTypes
end


function readRhs(mps, objectiveName, conNames, fixed)
    pos = position(mps)
    b = zeros(Float64, length(keys(conNames)))
    c0 = 0.0

    line = readline(mps)
    while strip(line) == "" || line[1] == '*'
        line = readline(mps)
    end
    while line[1] == ' '
        pos = position(mps)
        words = Vector{String}()

        if fixed
            lenLine = length(line)
            words = strip.([
                line[5:min(lenLine, 12)],
                line[15:min(lenLine, 22)],
                line[25:min(lenLine, 36)],
                line[40:min(lenLine, 47)],
                line[50:min(lenLine, 61)]
            ])
            while words[end] == ""
                pop!(words)
            end
        else
            words = split(line)
        end
        lenWords = length(words)

        for i in 2:2:lenWords
            name = words[i]
            if name == objectiveName
                if c0 != 0
                    error("Multiple presence of ", objectiveName, " in RHS")
                else
                    c0 = -parse(Float64, words[i+1])
                end
            elseif name != ""
                if b[conNames[name]] != 0
                    error("Multiple presence of ", name, " in RHS")
                else
                    b[conNames[name]] = parse(Float64, words[i+1])
                end
            end
        end

        line = readline(mps)
        while strip(line) == "" || line[1] == '*'
            pos = position(mps)
            line = readline(mps)
        end
    end

    seek(mps, pos)
    return b, c0
end


function readColumns(mps, objectiveName, conNames, fixed)
    pos = position(mps)
    varNames = Dict{String, Int}()
    varTypes = BitArray{1}() # true : float   false : int
    nCon = length(keys(conNames))
    A = zeros(Float64, nCon, 0)
    initVar = zeros(Float64, nCon, 1)
    c = Vector{Float64}()
    nVar = 1
    isFloat = true

    line = readline(mps)
    while strip(line) == "" || line[1] == '*'
        line = readline(mps)
    end
    while line[1] == ' '
        pos = position(mps)
        words = Vector{String}()

        if fixed
            lenLine = length(line)
            words = strip.([
                line[5:min(lenLine, 12)],
                line[15:min(lenLine, 22)],
                line[25:min(lenLine, 36)],
                line[40:min(lenLine, 47)],
                line[50:min(lenLine, 61)]
            ])
            while words[end] == ""
                pop!(words)
            end
        else
            words = split(line)
        end
        lenWords = length(words)

        # MARKER for INT
        if words[2] == "'MARKER'"
            ind = 3
            if fixed
                ind = 4
            end
            if words[ind] == "'INTORG'"
                isFloat = false
            elseif words[ind] == "'INTEND'"
                isFloat = true
            else
                error("Unknown MARKER in MPS file")
            end
        else
            name = words[1]
            if !haskey(varNames, name)
                varNames[name] = nVar
                nVar += 1
                A = [A initVar]
                push!(varTypes, isFloat)
                push!(c, 0.0)
            end

            for i in 2:2:lenWords
                if words[i] == objectiveName
                    # Case add to objective
                    c[varNames[name]] = parse(Float64, words[i+1])
                else
                    # Case add to A
                    A[conNames[words[i]], varNames[name]] = parse(Float64, words[i+1])
                end
            end
        end

        line = readline(mps)
        while strip(line) == "" || line[1] == '*'
            pos = position(mps)
            line = readline(mps)
        end
    end

    seek(mps, pos)
    return varNames, c, A, varTypes
end


function readRows(mps, objectiveName::String, fixed::Bool)
    pos = position(mps)
    objective = false
    conNames = Dict{String, Int}()
    conTypes = Vector{Int}() # 1 : <=, 2 : ==, 3 : >=
    nCon = 1

    line = readline(mps)
    while strip(line) == "" || line[1] == '*'
        line = readline(mps)
    end
    while line[1] == ' '
        pos = position(mps)
        if fixed
            typeTmp = strip(line[2:3])
            name = strip(line[5:min(length(line), 12)])
        else
            typeTmp, name = split(line)
            typeTmp = uppercase(typeTmp)
        end
        if typeTmp == "N"
            if objective
                error("MultiObjectives")
            else
                objective = true
                objectiveName = name
            end
        else
            if typeTmp == "G"
                push!(conTypes, 3)
            elseif typeTmp == "L"
                push!(conTypes, 1)
            elseif typeTmp == "E"
                push!(conTypes, 2)
            else
                error("Error in type of the row")
            end
            conNames[name] = nCon
            nCon += 1
        end

        line = readline(mps)
        while strip(line) == "" || line[1] == '*'
            pos = position(mps)
            line = readline(mps)
        end
    end

    seek(mps, pos)
    return objectiveName, conNames, conTypes
end


function readObjname(mps)
    line = readline(mps)
    while strip(line) == "" || line[1] == '*'
        line = readline(mps)
    end

    seek(mps, position(mps))
    return strip(line)
end


function readObjsense(mps)
    objectiveSense = ["MIN", "MAX"]

    line = readline(mps)
    while strip(line) == "" || line[1] == '*'
        line = readline(mps)
    end

    if !(strip(line) in objectiveSense)
        error("No valid information in OBJSENSE section")
    end

    seek(mps, position(mps))
    return strip(line) == "MIN" ? :min : :max
end

#end # module
