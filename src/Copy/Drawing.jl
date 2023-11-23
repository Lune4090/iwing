#= --- Comments--- =#
# In general, original structs or functions are minimal for easy maintanance
# Original structs are designed to be compatible with Makie & GeometryBasics

#= --- Structs --- =#

#= --- Functions --- =#

function draw_ViDARs_result(fig, AxisDict::Dict, AllObjDict::Dict, flags::Dict,
    center::Point3f, center_PDVec::Point3f,
    Llim::Float64, Rlim::Float64, scanres::Int, ScanningGrid::Vector{NamedTuple})
    # Visualize

    IsAnyNewAxCreated = false

    # Draw Stage
    if !flags["IsStageVisualised"] # Generate new Ax as Stage only if StageAx is not created
        AxisDict["StageAx"] = Axis(fig[1, 1], aspect=1)
        xlims!(AxisDict["StageAx"], -10, 10)
        ylims!(AxisDict["StageAx"], -10, 10)
        flags["IsStageVisualised"] = true
        IsAnyNewAxCreated = true
    else # Delete plots to initialize axis
        empty!(AxisDict["StageAx"])
    end

    # Draw Obj in StageAx
    # Draw polygons
    for key in keys(AllObjDict)
        lines!(AxisDict["StageAx"], AllObjDict[key].collisionmesh_world, color=:lightblue)
    end

    # Draw sight bar
    leftlimitsline = RotZ(Llim) * center_PDVec
    rightlimitsline = RotZ(-Rlim) * center_PDVec

    lines!(AxisDict["StageAx"],
        [center[1] + rightlimitsline[1], center[1], center[1] + leftlimitsline[1]],
        [center[2] + rightlimitsline[2], center[2], center[2] + leftlimitsline[2]],
        [center[3] + rightlimitsline[3], center[3], center[3] + leftlimitsline[3]],
        color=:lightgreen)

    # Draw ViDARs result
    # initial setting
    if !flags["IsViDARsVisualised"] # Generate new Ax as ViDARs only if ViDARsAx is not created
        AxisDict["ViDARsAxDict"] = Dict("ViDARsAx" => PolarAxis(fig[1, 2], rlimits=(0, 10), thetalimits=(-Llim, Rlim)), "ViDARsResult" => [])
        flags["IsViDARsVisualised"] = true
        IsAnyNewAxCreated = true
    else # Delete plots to initialize axis
        for plots in AxisDict["ViDARsAxDict"]["ViDARsResult"]
            delete!(AxisDict["ViDARsAxDict"]["ViDARsAx"], plots)
            AxisDict["ViDARsAxDict"]["ViDARsResult"] = []
        end
    end

    rarr = Vector{Float64}(undef, scanres)
    θarr = -Rlim:(Llim+Rlim)/(scanres-1):Llim
    for i in 1:scanres
        if ScanningGrid[i].dist != Inf
            rarr[i] = ScanningGrid[i].dist
        else
            rarr[i] = 0
        end
    end

    push!(
        AxisDict["ViDARsAxDict"]["ViDARsResult"],
        scatter!(AxisDict["ViDARsAxDict"]["ViDARsAx"], θarr, rarr,
            color=:lightblue, markersize=5)
    )

end