#= --- Comments--- =#
# In general, original structs or functions are minimal for easy maintanance
# Original structs are designed to be compatible with Makie & GeometryBasics

#= --- Structs --- =#

#= --- Functions --- =#

function draw_ViDARs_result(fig, AxisDict::Dict, AllObjDict::Dict, flags::Dict,
    center::Vector{Float64}, center_PDVec::Vector{Float64},
    Llim::Float64, Rlim::Float64, scanres::Int, ScanningGrid::Vector{Dict{String,Any}})
    # Visualize

    IsAnyNewAxCreated = false

    # Draw Stage
    if !flags["IsStageVisualised"] # Generate new Ax as Stage only if StageAx is not created
        AxisDict["StageAxDict"] = Dict("StageAx" => Axis(fig[1, 1], aspect=1), "nonplayerplots" => [], "playerplots" => [])
        #maincam = maindisplay.scene.camera
        xlims!(AxisDict["StageAxDict"]["StageAx"], -10, 10)
        ylims!(AxisDict["StageAxDict"]["StageAx"], -10, 10)
        #zlims!(AxisDict["StageAxDict"]["StageAx"], -1, 1)
        flags["IsStageVisualised"] = true
        IsAnyNewAxCreated = true
    else # Delete plots to initialize axis
        for plots in AxisDict["StageAxDict"]["nonplayerplots"]
            delete!(AxisDict["StageAxDict"]["StageAx"], plots)
        end
        for plots in AxisDict["StageAxDict"]["playerplots"]
            delete!(AxisDict["StageAxDict"]["StageAx"], plots)
        end
        AxisDict["StageAxDict"]["nonplayerplots"] = []
        AxisDict["StageAxDict"]["playerplots"] = []
    end

    #=
        # Draw Obj in StageAx
        # Draw polygons
        facenum = AllObjDict[1].CollisionMesh.face_num
        for i in 1:facenum
            push!(AxisDict["StageAxDict"]["nonplayerplots"], lines!(
                AxisDict["StageAxDict"]["StageAx"],
                [AllObjDict[1].CollisionMesh.polygons[1, 1, i], AllObjDict[1].CollisionMesh.polygons[1, 2, i]],
                [AllObjDict[1].CollisionMesh.polygons[2, 1, i], AllObjDict[1].CollisionMesh.polygons[2, 2, i]],
                [AllObjDict[1].CollisionMesh.polygons[3, 1, i], AllObjDict[1].CollisionMesh.polygons[3, 2, i]],
                color=:lightblue
            ))
        end
    	=#

    # Draw player
    push!(AxisDict["StageAxDict"]["playerplots"], scatter!(
        AxisDict["StageAxDict"]["StageAx"],
        center[1],
        center[2],
        center[3],
        color=:orange
    ))
    leftlimitsline = quaternion2vector(rotate_vector_by_quaternion(MyRotationQuaternion([0, 0, 1], Llim), center_PDVec))
    rightlimitsline = quaternion2vector(rotate_vector_by_quaternion(MyRotationQuaternion([0, 0, 1], -Rlim), center_PDVec))

    push!(AxisDict["StageAxDict"]["playerplots"], lines!(
        AxisDict["StageAxDict"]["StageAx"],
        [center[1], center[1] + leftlimitsline[1]],
        [center[2], center[2] + leftlimitsline[2]],
        [center[3], center[3] + leftlimitsline[3]],
        color=:lightgreen
    ))
    push!(AxisDict["StageAxDict"]["playerplots"], lines!(
        AxisDict["StageAxDict"]["StageAx"],
        [center[1], center[1] + rightlimitsline[1]],
        [center[2], center[2] + rightlimitsline[2]],
        [center[3], center[3] + rightlimitsline[3]],
        color=:lightgreen
    ))

    # Draw ViDARs result
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

    # Draw ViDARs result
    rarr = Vector{Float64}(undef, scanres)
    θarr = Vector{Float64}(undef, scanres)
    for i in 1:scanres
        if ScanningGrid[i]["Dist"] != Inf
            rarr[i] = ScanningGrid[i]["Dist"]
        else
            rarr[i] = 0
        end
        θarr[i] = i * (Llim + Rlim) / scanres - Rlim
    end

    push!(
        AxisDict["ViDARsAxDict"]["ViDARsResult"],
        scatter!(AxisDict["ViDARsAxDict"]["ViDARsAx"], θarr, rarr,
            color=:lightblue, markersize=5)
    )

end