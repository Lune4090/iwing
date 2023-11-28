#= --- Comments--- =#
# In general, original structs or functions are minimal for easy maintanance
# Original structs are designed to be compatible with Makie & GeometryBasics
# This file mainly aims to pseudo physics calculation

#= --- Structs --- =#

#= --- Functions --- =#

function propertycalc(
    ScanningGrid::T;
    wavetype::String="Normal", detect_lim::Float32=0.1f0, standard_len::Float32=20.0f0,
    player_velocity::Point3f=Point3f(0)) where {T<:Vector{Vector}}
    AnalyzedGrid = similar(ScanningGrid)
    stepθ = 0
    for gridθ in ScanningGrid
        stepθ += 1
        gridθ, isallInf = erase_inf(gridθ)
        isallInf && sort!(gridθ, by=x -> x[1])
        transparency = 1.0f0
        analayzedresults = []
        isallInf && for eachresult in gridθ
            analyzedres = reflectioncalc(player_velocity,
                eachresult.dist, detect_lim, standard_len,
                wavetype, transparency,
                eachresult.normv, eachresult.scandir, eachresult.attr)
            transparency *= analyzedres.transparency
            push!(analayzedresults, analyzedres)
        end
        AnalyzedGrid[stepθ] = analayzedresults
    end
    return AnalyzedGrid
end

function erase_inf(gridθ)
    tmp_gridθ = similar(gridθ)
    for idx in eachindex(tmp_gridθ)
        tmp_gridθ[idx] = (dist=Inf, normv=nothing, scandir=nothing, num=nothing, objkey=nothing, attr=nothing)
    end
    num_nonInfres = 0
    for eachresult in gridθ
        if eachresult.dist != Inf && eachresult.dist != -Inf
            num_nonInfres += 1
            tmp_gridθ[num_nonInfres] = eachresult
        end
    end
    num_nonInfres == 0 && return tmp_gridθ[1:1], false
    return tmp_gridθ[1:num_nonInfres], true
end

function reflectioncalc(player_velocity::Point3f,
    dist::Float32, detect_lim::Float32, standard_len::Float32,
    wavetype::String, transparency::Float32,
    scandir::Point3f, normvec::Point3f, attr::Dict)
    reflection = 1.0f0 * transparency
    # distance decay
    reflection -= (1.0f0 - detect_lim) * (dist / standard_len)
    # size correction
    reflection *= attr["size"]
    # wavetype correction
    reflection *= attr["reflectionspectrum"][wavetype]
    # calc if reflection over detect_lim
    coefficient_reflectangle = normvec ⋅ -scandir # roughnessとか定義して使うかも
    analyzed_res = (
        dist=reflection > detect_lim ? dist : Inf,
        relative_velocity=dopplercalc(player_velocity, attr["velocity"]),
        transparency=attr["transparency"],
        reflection_amp=reflection)
    return analyzed_res
end

# まだ出来てない
function dopplercalc(player_velocity, obj_velocity)
    relative_velocity = player_velocity - obj_velocity
    velX = abs(relative_velocity[1])
    velY = abs(relative_velocity[2])
    velZ = abs(relative_velocity[3])
    return [velX, velY, velZ]
end