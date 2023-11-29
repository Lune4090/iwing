clear() = run(`cmd /c cls`)

using GLMakie
using LinearAlgebra
using GeometryBasics
using GeometryBasics: Mesh
using Rotations

include("Drawing.jl")
include("GameObject.jl")
include("UtlityFuncs.jl")
include("ViDARs.jl")
include("PseudoPhysics.jl")

#= --- Comments--- =#
# In general, original structs or functions are minimal for easy maintanance
# Original structs are designed to be compatible with Makie & GeometryBasics
# These codes control game progress by combining a lot of function&structs

#= --- Structs --- =#

#= --- Functions --- =#


# InitialSetting to start game
# Define variable and get some outside data
function InitialSetting(savedata::NamedTuple)
    # Generate Window
    fig = Figure()
    AxisDict = Dict{String,Any}()
    display(fig)

    # Object initialize
    AllObjDict = Dict{Int,MyGameObject}()

    # Flag initialize
    flags = Dict{String,Any}()
    flags["IsFieldVisualised"] = false
    flags["IsStageVisualised"] = false
    flags["IsViDARsVisualised"] = false

    # Set variable based on saved data
    AllObjDict = savedata.AllObjDict
    player_pos = savedata.player_pos
    player_direction = savedata.player_direction
    search_direction = savedata.search_direction
    chase_direction = savedata.chase_direction

    return fig, AxisDict,
    AllObjDict, flags,
    player_pos, player_direction,
    search_direction, chase_direction
end

# Call function which should be done in a frame
function EachFrame(fig, AxisDict::Dict, AllObjDict::Dict, flags::Dict,
    player_inputs)
    #println("----------------------------------------------")
    #println("EachFrameExecution Started")

    # Update all object
    for key in keys(AllObjDict)
        obj = AllObjDict[key]
        #println("----------------------------------------------")
        #println("Start processing object : $key")
        object_mover!(obj, Point3f(0, 0.05, 0), 60.0f0)
        update_globalizedmesh2D!(obj)
    end

    # Detect globalized(collision)mesh
    @time ScanningGrid_search = ViDARsLoop(
        AllObjDict, player_position, search_direction,
        Rlim, Llim, scanres)

    # Calculate detected things property
    @time AnalyzedGrid_search = propertycalc(
        ScanningGrid_search; wavetype="IR")

    # Draw display
    @time draw_main(fig, AxisDict, AllObjDict, flags,
        player_position, player_direction, 10.0f0,
        Llim, Rlim, search_direction)

    @time draw_ViDARs_result(fig, AxisDict, flags,
        Llim, Rlim, scanres, 50.0f0,
        AnalyzedGrid_search, search_direction)

    return search_direction
end


#=-------- Main Game loop ---------=#
function GameLoop(fig, AxisDict::Dict, AllObjDict::Dict, flags::Dict,
    player_position::Point3f, player_direction::QuatRotation,
    search_direction::QuatRotation, chase_direction::QuatRotation)

    Llim = π / 36.0f0
    Rlim = π / -36.0f0
    scanres = 5
    search_direction = RotZ((Llim - Rlim) / 2) * search_direction |> QuatRotation

    player_inputs = EachFrame(fig, AxisDict, AllObjDict, flags, player_inputs)

    return search_direction
end


#=------- End Game Execution --------=#
function EndGame()
    "Game end"
end

#=-------- GameMain --------=#

function GameMain(savedata::NamedTuple)
    args... = InitialSetting(savedata)

    #=-------- Begin GameLoop --------=#
    while to_value(events(fig).window_open)
        GameLoop(args...)
    end

    #=-------- End Game execution--------=#
    EndGame()
end

#GameMain()

#= --- Comment Zone ---=#

#= --- TestCodes --- =#

#=
現時点でのViDARSについての結論
・@threadsは効果なし(むしろ1/2位に落ち込んだ)→いやこれは立ち上げ時にスレッド数指定してないからか
・内部の計算が結構複雑だからGPUのレンダリングパイプラインを上手く使わない限り高速化はほぼ不可能
・取り敢えずCPUでいく
・目安として，キャラクター当たり判定のエッジの基準長を1[m]とした時，θres=30で走査角π/6(: dθ= 2°)だと50m先のオブジェクトがガタガタして映る感じ
・dθ=2°の場合，1000エッジくらいまでは60FPSを安定して出せる
=#