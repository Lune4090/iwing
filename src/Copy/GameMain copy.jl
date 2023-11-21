using GLMakie
using LinearAlgebra

include("AboutQuarternion.jl")
include("AboutDrawing.jl")
include("AboutMesh.jl")
include("AboutObject.jl")
include("UtlityFuncs.jl")
include("ViDARs.jl")

GLMakie.activate!(framerate=120, render_on_demand=false)

# Call function which should be done in a frame
function EachFrame(fig, AllObjDict::Dict, FlagBoolDict::Dict, FlagIntDict::Dict, FlagFloatDict::Dict, AxisDict::Dict,
    center::Vector{Float64}, center_PDVec::Vector{Float64})
    #println("----------------------------------------------")
    #println("EachFrameExecution Started")

    # Polygon projection
    for key in keys(AllObjDict)
        IGO = AllObjDict[key]
        #println("----------------------------------------------")
        #println("Start processing object : $key")
        TranslatePolygonLocal2Global!(IGO.CollisionMesh, IGO)
        #println("Finish translating CollisionMesh (object : $key)")
        #TranslatePolygonLocal2Global!(IGO.AppearanceMesh, IGO)
        #println("Finish translating AppearanceMesh (object : $key)")

    end

    Llim = π / 4
    Rlim = π / 4
    scanres = 1000

    # 上記関数でローカル座標情報から生成されたグローバル座標上のポリゴンを処理
    ScanningGrid = ViDARsLoop2D(
        AllObjDict::Dict, center::Vector{Float64},
        center_PDVec::Vector{Float64};
        θRlim=Rlim, θLlim=Llim, θres=scanres)

    # Draw display
    DrawObjects(fig, AllObjDict, FlagBoolDict, FlagIntDict, FlagFloatDict, AxisDict,
        center, center_PDVec, Llim, Rlim, scanres, ScanningGrid)

end

# Main Game loop
function GameLoop(fig, AllObjDict::Dict, FlagBoolDict::Dict, FlagIntDict::Dict, FlagFloatDict::Dict, AxisDict::Dict,
    center::Vector{Float64}, center_PDVec::Vector{Float64})
    #println("----------------------------------------------")
    #println("GameLoop Started")
    EachFrame(fig, AllObjDict, FlagBoolDict, FlagIntDict, FlagFloatDict, AxisDict, center, center_PDVec)
end

# InitialSetting to start game
function InitialSetting()
    fig = Figure()
    display(fig)
    AllObjDict = Dict{Int,MyGameObject}()
    FlagBoolDict = Dict{String,Bool}()
    FlagIntDict = Dict{String,Int}()
    FlagFloatDict = Dict{String,Float64}()
    AxisDict = Dict{String,Any}()
    FlagBoolDict["IsFieldVisualised"] = false
    FlagBoolDict["IsStageVisualised"] = false
    FlagBoolDict["IsViDARsVisualised"] = false
    println("InitialSetting is finished without any problem")
    return fig, AllObjDict, FlagBoolDict, FlagIntDict, FlagFloatDict, AxisDict
end

#= --- Comment Zone ---=#

#= --- TestCodes --- =#

function GameMain()
    fig, AllObjDict, FlagBoolDict, FlagIntDict, FlagFloatDict, AxisDict = InitialSetting()
    # 全てのオブジェクトは作られたのちにこのDictに追記されることで初めて名前と存在をゲームから認められる
    tmpvertices2D = [
        0.0 0.0 0.0;
        1.0 1.0 0.0
    ]

    tmpfaces2D = [
        1 2
    ]

    tmpvertices2D = permutedims(tmpvertices2D)
    tmpfaces2D = permutedims(tmpfaces2D)

    vertnum = 20
    scale = 10
    polynum = 10

    tmpvertices2D = rand(3, vertnum) * scale
    # convert to pseudo 2D
    for i in 1:vertnum
        tmpvertices2D[3, i] = 0.0
    end
    # generate face(∼polygon) randomly
    tmpfaces2D = rand(2, polynum) * scale
    # convert to 1-start idx
    tmpfaces2D = map(tmpfaces2D) do x
        floor(x) .+ 1.0
    end
    tmpfaces2D = map(tmpfaces2D) do x
        convert(Int, x)
    end

    #=
        tmpvertices2D = rand(3, 10) * 10
    	=#
    tmpMeshPDVec = [-1.0, 0.0, 0.0]
    Mymesh = CollisionMesh(tmpvertices2D, tmpfaces2D, tmpMeshPDVec, 2)
    tmpIGOpos = [1.0, 1.0, 0.0]
    tmpIGOPDVec = [1.0, 0.0, 0.0]
    tmpattrdict1 = Dict()
    tmpAttrArr = [tmpattrdict1]
    Mygameobj = MyGameObject(
        tmpIGOpos,
        tmpIGOPDVec,
        Mymesh,
        Nothing,
        tmpAttrArr)

    # 辞書に登録
    AllObjDict[1] = Mygameobj

    # PlayerCharacter position (最終的にはObjとして引き渡す)
    player_pos = [0.0, 0.0, 0.0]
    player_PDVec = [1.0, 0.0, 0.0]

    # plotting object initialization
    #GameLoop(fig, AllObjDict, FlagBoolDict, FlagIntDict, FlagFloatDict, AxisDict, player_pos, player_PDVec)


    while to_value(events(fig).window_open)
        player_PDVec = quaternion2vector(rotate_vector_by_quaternion(MyRotationQuaternion([0, 0, 1], π / 120), player_PDVec))
        println("------------")
        GameLoop(fig, AllObjDict, FlagBoolDict, FlagIntDict, FlagFloatDict, AxisDict, player_pos, player_PDVec)
        println("------------")
    end
end

GameMain()
