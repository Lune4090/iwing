#= --- Comments--- =#
# In general, original structs or functions are minimal for easy maintanance
# Original structs are designed to be compatible with Makie & GeometryBasics

#= --- Structs --- =#

#= --- Functions --- =#

function draw_main(fig, AxisDict::Dict, AllObjDict::Dict, flags::Dict,
    player_position::Point3f, player_direction::QuatRotation, lim_main::Float32,
    Llim::Float32, Rlim::Float32, search_direction::QuatRotation)

    # Draw Stage
    if !flags["IsStageVisualised"] # Generate new Ax as Stage only if StageAx is not created
        AxisDict["StageAx"] = Axis(fig[1, 1], aspect=1)
        xlims!(AxisDict["StageAx"], -lim_main, lim_main)
        ylims!(AxisDict["StageAx"], -lim_main, lim_main)
        flags["IsStageVisualised"] = true
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
        [player_position[1], player_position[1] + cos(rotation_angle(player_direction))],
        [player_position[2], player_position[2] + sin(rotation_angle(player_direction))],
        [player_position[3], player_position[3]],
        color=:orange)
end

function draw_ViDARs_result(fig, AxisDict::Dict, flags::Dict,
    Llim::Float32, Rlim::Float32, scanres::Int, lim_vidars::Float32,
    AnalyzedGrid::Vector{Vector}, search_direction::QuatRotation)

    # Draw ViDARs result
    # initial setting
    if !flags["IsViDARsVisualised"] # Generate new Ax as ViDARs only if ViDARsAx is not created
        AxisDict["ViDARsAxDict"] = Dict("ViDARsAx" => PolarAxis(fig[1, 2], rlimits=(0, lim_vidars)), "ViDARsResult" => [])
        flags["IsViDARsVisualised"] = true
    else # Delete plots to initialize axis
        for plots in AxisDict["ViDARsAxDict"]["ViDARsResult"]
            delete!(AxisDict["ViDARsAxDict"]["ViDARsAx"], plots)
            AxisDict["ViDARsAxDict"]["ViDARsResult"] = []
        end
    end

    rarr = Vector{Float32}(undef, 0)
    θarr = Vector{Float32}(undef, 0)
    Llim += rotation_angle(search_direction)
    Rlim += rotation_angle(search_direction)
    tmpθarr = LinRange(Rlim, Llim, scanres + 1)
    tmpθarr = tmpθarr[1:scanres]
    for i in 1:scanres
        for j in 1:length(AnalyzedGrid[i])
            if AnalyzedGrid[i][j].dist != Inf
                push!(rarr, AnalyzedGrid[i][j].dist)
                push!(θarr, tmpθarr[i])
            end
        end
    end

    push!(
        AxisDict["ViDARsAxDict"]["ViDARsResult"],
        scatter!(AxisDict["ViDARsAxDict"]["ViDARsAx"], θarr, rarr,
            color=:blue, markersize=10)
    )

end