clear() = run(`cmd /c cls`)

using GLMakie
using GeometryBasics
using GeometryBasics: Mesh
using Rotations

include("Drawing.jl")
include("GameObject.jl")
include("UtlityFuncs.jl")
include("ViDARs.jl")

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

# Main Game loop
function GameLoop(fig, AxisDict::Dict, AllObjDict::Dict, flags::Dict,
    center::Point3f, center_direction::Point3f)
    #println("----------------------------------------------")
    #println("GameLoop Started")
    EachFrame(fig, AxisDict, AllObjDict, flags, center, center_direction)
end

# Call function which should be done in a frame
function EachFrame(fig, AxisDict::Dict, AllObjDict::Dict, flags::Dict,
    center::Point3f, center_direction::Point3f)
    #println("----------------------------------------------")
    #println("EachFrameExecution Started")

    # Polygon projection
    for key in keys(AllObjDict)
        obj = AllObjDict[key]
        #println("----------------------------------------------")
        #println("Start processing object : $key")
        generate_globalizedmesh2D!(obj)
    end

    # put them into flags
    Llim = π / 4
    Rlim = π / 4
    scanres = 100

    # detect globalized(collision)mesh
    @time ScanningGrid = ViDARsLoop(
        AllObjDict, center, center_direction,
        Rlim, Llim, scanres)

    # Draw display
    @time draw_ViDARs_result(fig, AxisDict, AllObjDict, flags,
        center, center_direction,
        Llim, Rlim, scanres, ScanningGrid)

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

    Mygameobj = MyGameObject(
        Point3f(2.0, 2.0, 0.0),
        Point3f(1.0, 0.0, 0.0),
        Mesh(tmpvertices2D, tmpfaces2D),
        Mesh(tmpvertices2D, tmpfaces2D),
        Dict(),
        Mesh(tmpvertices2D, tmpfaces2D),
        Mesh(tmpvertices2D, tmpfaces2D))

    # 辞書に登録
    AllObjDict[1] = Mygameobj

    # PlayerCharacter position (最終的にはObjとして引き渡す)
    player_pos = Point3f(0.0)
    player_direction = Point3f(0.0, 1.0, 0.0)

    # plotting object initialization

    while to_value(events(fig).window_open)
        player_direction = player_direction |> x -> RotZ(π / 120) * x |> Point3f
        println("------------")
        GameLoop(fig, AxisDict, AllObjDict, flags, player_pos, player_direction)
    end
end

GameMain()