#= --- Comments--- =#
# In general, original structs or functions are minimal for easy maintanance
# Original structs are designed to be compatible with Makie & GeometryBasics

#= --- Structs --- =#

#= --- Functions --- =#

function draw_ViDARs_result(fig, AxisDict::Dict, AllObjDict::Dict, flags::Dict,
    player_position::Point3f, player_direction::Point3f,
    Llim::Float64, Rlim::Float64, scanres::Int, ScanningGrid::Vector{NamedTuple}, search_direction::Point3f)
    # Visualize

    IsAnyNewAxCreated = false
    lim_main = 7.5
    lim_vidars = 50

    # Draw Stage
    if !flags["IsStageVisualised"] # Generate new Ax as Stage only if StageAx is not created
        AxisDict["StageAx"] = Axis(fig[1, 1], aspect=1)
        xlims!(AxisDict["StageAx"], -lim_main, lim_main)
        ylims!(AxisDict["StageAx"], -lim_main, lim_main)
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
    leftlimitsline = RotZ(Llim) * search_direction
    rightlimitsline = RotZ(Rlim) * search_direction

    lines!(AxisDict["StageAx"],
        [player_position[1] + rightlimitsline[1], player_position[1], player_position[1] + leftlimitsline[1]],
        [player_position[2] + rightlimitsline[2], player_position[2], player_position[2] + leftlimitsline[2]],
        [player_position[3] + rightlimitsline[3], player_position[3], player_position[3] + leftlimitsline[3]],
        color=:lightgreen)

    lines!(AxisDict["StageAx"],
        [player_position[1], player_position[1] + player_direction[1]],
        [player_position[2], player_position[2] + player_direction[2]],
        [player_position[3], player_position[3] + player_direction[3]],
        color=:orange)

    lines!

    # Draw ViDARs result
    # initial setting
    if !flags["IsViDARsVisualised"] # Generate new Ax as ViDARs only if ViDARsAx is not created
        AxisDict["ViDARsAxDict"] = Dict("ViDARsAx" => PolarAxis(fig[1, 2], rlimits=(0, lim_vidars)), "ViDARsResult" => [])
        flags["IsViDARsVisualised"] = true
        IsAnyNewAxCreated = true
    else # Delete plots to initialize axis
        for plots in AxisDict["ViDARsAxDict"]["ViDARsResult"]
            delete!(AxisDict["ViDARsAxDict"]["ViDARsAx"], plots)
            AxisDict["ViDARsAxDict"]["ViDARsResult"] = []
        end
    end

    rarr = Vector{Float64}(undef, scanres)
    Llim += decide_arg2D(search_direction)
    Rlim += decide_arg2D(search_direction)
    θarr = LinRange(Rlim, Llim, scanres + 1)
    θarr = θarr[1:scanres]
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
            color=:blue, markersize=10)
    )

end