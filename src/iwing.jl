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
function InitialSetting()
    # defince new variable
    fig = Figure()
    AxisDict = Dict{String,Any}()
    display(fig)
    AllObjDict = Dict{Int,MyGameObject}()
    flags = Dict{String,Any}()
    # variables initial setting
    flags["IsFieldVisualised"] = false
    flags["IsStageVisualised"] = false
    flags["IsViDARsVisualised"] = false

    return fig, AxisDict, AllObjDict, flags
end

# Call function which should be done in a frame
function EachFrame(fig, AxisDict::Dict, AllObjDict::Dict, flags::Dict,
    player_position::Point3f, player_direction::QuatRotation,
    search_direction::QuatRotation, chase_direction::QuatRotation)
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

    Llim = π / 36.0f0
    Rlim = π / -36.0f0
    scanres = 5
    search_direction = RotZ((Llim - Rlim) / 2) * search_direction |> QuatRotation

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

# Main Game loop
function GameLoop(fig, AxisDict::Dict, AllObjDict::Dict, flags::Dict,
    player_position::Point3f, player_direction::QuatRotation,
    search_direction::QuatRotation, chase_direction::QuatRotation)
    #println("----------------------------------------------")
    #println("GameLoop Started")
    search_direction = EachFrame(fig, AxisDict, AllObjDict, flags, player_position, player_direction, search_direction, chase_direction)
    return search_direction
end

#= --- Comment Zone ---=#

#= --- TestCodes --- =#

GLMakie.activate!(framerate=120, render_on_demand=false)

function GameMain()
    fig, AxisDict, AllObjDict, flags = InitialSetting()
    # 全てのオブジェクトは作られたのちにこのDictに追記されることで初めて名前と存在をゲームから認められる

    p1 = Point3f(1.0, 0.0, 0.0)
    p2 = Point3f(2.0, 2.0, 0.0)
    tmpvertices2D = [p1, p2]
    tmpfaces2D = [TriangleFace([1, 2, 1]), TriangleFace([2, 1, 2])]

    # 辞書に登録
    tmpdict = Dict("Normal" => 1.0f0, "IR" => 1.0f0)
    for i in 1:10
        AllObjDict[i] = MyGameObject(
            Point3f(2.0, 0.0, 0.0),
            QuatRotation(RotZ(0)),
            Mesh(tmpvertices2D, tmpfaces2D),
            Mesh(tmpvertices2D, tmpfaces2D),
            Dict("size" => 1.0f0, "reflectionspectrum" => tmpdict, "transparency" => 1.0f0),
            Mesh(tmpvertices2D, tmpfaces2D),
            Mesh(tmpvertices2D, tmpfaces2D))
    end

    @show length(AllObjDict)

    # PlayerCharacter position (最終的にはObjとして引き渡す)
    player_pos = Point3f(0.0)
    player_direction = QuatRotation(RotZ(0))
    search_direction = QuatRotation(RotZ(0))
    chase_direction = QuatRotation(RotZ(0))

    # plotting object initialization

    while to_value(events(fig).window_open)
        player_direction = player_direction |> x -> RotZ(π / 180) * x |> QuatRotation
        println("------------")
        search_direction = GameLoop(fig, AxisDict, AllObjDict, flags, player_pos, player_direction, search_direction, chase_direction)
    end
end

GameMain()

#=
現時点でのViDARSについての結論
・@threadsは効果なし(むしろ1/2位に落ち込んだ)→いやこれは立ち上げ時にスレッド数指定してないからか
・内部の計算が結構複雑だからGPUのレンダリングパイプラインを上手く使わない限り高速化はほぼ不可能
・取り敢えずCPUでいく
・目安として，キャラクター当たり判定のエッジの基準長を1[m]とした時，θres=30で走査角π/6(: dθ= 2°)だと50m先のオブジェクトがガタガタして映る感じ
・dθ=2°の場合，1000エッジくらいまでは60FPSを安定して出せる
=#