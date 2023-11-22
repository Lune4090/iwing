using GLMakie
using LinearAlgebra

GLMakie.activate!(framerate=120, render_on_demand=false)

#= --- Original Structures Definition --- =#
abstract type MyAbstractQuaternion end
abstract type MyAbstractNormalizedQuaternion <: MyAbstractQuaternion end
abstract type MyAbstractRotationQuaternion <: MyAbstractNormalizedQuaternion end

mutable struct MyQuaternion <: MyAbstractQuaternion
    x::Float64
    y::Float64
    z::Float64
    θ::Float64
end

MyQuaternion(x::Real, y::Real, z::Real, θ::Real) = MyQuaternion(convert(Float64, x), convert(Float64, y), convert(Float64, z), convert(Float64, θ))
# Vector -> Quarternion
MyQuaternion(vec::Vector) = MyQuaternion(vec[1], vec[2], vec[3], 0.0)
# Vector, θ -> Quarternion
MyQuaternion(vec::Vector, θ) = MyQuaternion(vec[1], vec[2], vec[3], θ)

mutable struct MyNormalizedQuaternion <: MyAbstractNormalizedQuaternion
    x::Float64
    y::Float64
    z::Float64
    θ::Float64
    MyNormalizedQuaternion(x, y, z, θ) = x^2 + y^2 + z^2 != 1 ? error("MyNormalizedQuaternion's Norm != 1!") : new(x, y, z, θ)
end

# Int->Float64 Conversion
MyNormalizedQuaternion(x::Int, y::Int, z::Int, θ::Int) = MyNormalizedQuaternion(convert(Float64, x), convert(Float64, y), convert(Float64, z), convert(Float64, θ))
# Vector -> Quarternion
MyNormalizedQuaternion(vec::Vector) = MyNormalizedQuaternion(vec[1], vec[2], vec[3], 0.0)
# Vector, θ -> Quarternion
MyNormalizedQuaternion(vec::Vector, θ) = MyNormalizedQuaternion(vec[1], vec[2], vec[3], θ)

mutable struct MyRotationQuaternion <: MyAbstractRotationQuaternion
    x::Float64
    y::Float64
    z::Float64
    θ::Float64
end
# Vector, θ -> Quarternion, this is special constructor for Rotation Quarternion
# https://qiita.com/drken/items/0639cf34cce14e8d58a5
MyRotationQuaternion(vec::Vector, θ) = MyRotationQuaternion(vec[1] * sin(θ / 2), vec[2] * sin(θ / 2), vec[3] * sin(θ / 2), cos(θ / 2))

#= --- Original Structures' functions --- =#

function RotateQuaternion(Rotator::MyAbstractQuaternion, Rotated::MyAbstractQuaternion)
    return MyQuaternion(
        Rotator.θ * Rotated.x - Rotator.z * Rotated.y + Rotator.y * Rotated.z + Rotator.x * Rotated.θ,
        Rotator.z * Rotated.x + Rotator.θ * Rotated.y - Rotator.x * Rotated.z + Rotator.y * Rotated.θ,
        -Rotator.y * Rotated.x + Rotator.x * Rotated.y + Rotator.θ * Rotated.z + Rotator.z * Rotated.θ,
        -Rotator.x * Rotated.x - Rotator.y * Rotated.y - Rotator.z * Rotated.z + Rotator.θ * Rotated.θ
    )
end

# https://qiita.com/kenjihiranabe/items/945232fbde58fab45681
function RotateVectorbyQuaternion(q::MyAbstractQuaternion, vec::Vector)
    tmp = MyQuaternion(
        q.θ * vec[1] - q.z * vec[2] + q.y * vec[3] + q.x * 1,
        q.z * vec[1] + q.θ * vec[2] - q.x * vec[3] + q.y * 1,
        -q.y * vec[1] + q.x * vec[2] + q.θ * vec[3] + q.z * 1,
        -q.x * vec[1] - q.y * vec[2] - q.z * vec[3] + q.θ * 1
    )
    return RotateQuaternion(tmp, MyQuaternion(-q.x, -q.y, -q.z, q.θ))

end

function Quaternion2Vector(q::MyQuaternion)
    v = Vector{Float64}(undef, 3)
    v[1] = q.x
    v[2] = q.y
    v[3] = q.z
    return v
end

function Vector2Quarternion(v::Vector{Real})
    return MyQuaternion(v[1], v[2], v[3], 1)
end

abstract type MyAbstractMesh end
abstract type MyAbstractCollisionMesh <: MyAbstractMesh end
abstract type MyAbstractAppearanceMesh <: MyAbstractMesh end

@kwdef mutable struct CollisionMesh <: MyAbstractCollisionMesh
    Vertices::Matrix{Float64}
    Faces::Matrix{Int}
    PostureDirectionVector::Vector{Float64}
    dim::Int
    polygons = Nothing
    face_num::Int64 = 0
    is_Decomposed = false
    function CollisionMesh(
        Vertices::Matrix{Float64},
        Faces::Matrix{Int},
        PostureDirectionVector::Vector{Float64},
        dim::Int,
        polygons=Nothing,
        face_num::Int=0,
        is_Decomposed=false)
        if size(Vertices, 1) != 3
            error("Vertice should be 3Dim!")
        elseif size(Faces, 1) != dim
            error("The number of vertices should be matched to defiend mesh dim($(dim))!")
        else
            new(
                Vertices,
                Faces,
                PostureDirectionVector,
                dim,
                polygons,
                face_num,
                is_Decomposed)
        end
    end
end

# VerticesからFacesに面毎に格納された頂点番号を引数にpolygonを生成
# linearindexingによるアンローリングよりforの二重ループの方が速かったのでこのまま
function DecomposeMesh!(mesh::MyAbstractMesh)
    num_faces = size(mesh.Faces, 2)
    dim = size(mesh.Faces, 1)
    polygons = Array{Float64}(undef, 3, dim, num_faces)
    for EachFace in 1:num_faces
        for EachVertex in 1:dim
            polygons[:, EachVertex, EachFace] .= mesh.Vertices[:, mesh.Faces[EachVertex, EachFace]]
        end
    end
    mesh.polygons = polygons
    mesh.face_num = num_faces
    mesh.is_Decomposed = true
end

@kwdef mutable struct InGameObj
    Pos::Vector{Float64} # これがグローバル座標
    PostureDirectionVector::Vector{Float64}
    CollisionMesh::CollisionMesh
    AppearanceMesh::Any
    # AppearanceMeshは外部ツールで製作した.stl形式のオブジェクトのloadを想定
    # AppearanceMeshは触らない(Makie.meshscatter!()にmarkerとして渡すだけ)
    Attributes::Vector{Dict{Any,Any}}
    # AttributesはIGOのCollisionMeshの各面毎に算出される
    # Attributes，Mesh側が持つべきという説もある
    # KinecticProperty
    # Velocity
    # ArgVelocity
    # PassiveMaterialProperty
    # Roughness
    # ReflectionSpectrum
    # ActiveMaterialProperty
    # RadiationSpectrum
    # Raise error when PDvec != 1
    function InGameObj(Pos::Vector{Float64}, PostureDirectionVector::Vector{Float64}, CollisionMesh::Any, AppearanceMesh::Any, Attributes::Vector{Dict{Any,Any}})
        norm_PDVec = PostureDirectionVector[1]^2 + PostureDirectionVector[2]^2 + PostureDirectionVector[3]^2
        if norm_PDVec != 1
            error("PostureDirectionVector should be equal to 1 !!!")
        else
            new(Pos, PostureDirectionVector, CollisionMesh, AppearanceMesh, Attributes)
        end
    end
end

function CalcNormalVector(polygon::Matrix{Float64})
    # ここ、2Dポリゴンを3Dポリゴンの3つ目の頂点を2つ目に合わせることで実装している影響で代えられない
    vec_a = polygon[:, 2] - polygon[:, 1]
    vec_b = polygon[:, 3] - polygon[:, 1]
    NormalVector = cross(vec_a, vec_b)
    return NormalVector / norm(NormalVector)
end

# 基本的に2D、3Dを問わずポリゴンは3次元座標系で保持される為流用可能
function TranslatePolygonLocal2Global!(mesh, IGO::InGameObj)
    DecomposeMesh!(mesh)
    # Local2Global Coordinate Translation
    # Rotation (by Quarternion[λx*sin(θ/2),λy*sin(θ/2),λz*sin(θ/2),cos(θ/2)])
    RotateAxis = cross(mesh.PostureDirectionVector, IGO.PostureDirectionVector)
    # 通常は回転軸を中心に回転
    if RotateAxis != [0.0, 0.0, 0.0]
        RotateAxis = RotateAxis / norm(RotateAxis)
        # CounterClockwiseを仮定
        RotateAngle = acos(mesh.PostureDirectionVector ⋅ IGO.PostureDirectionVector)
        for face in 1:mesh.face_num
            for vert in 1:mesh.dim
                #println("Axis vec : $RotateAxis")
                #println("Angle : $RotateAngle")
                mesh.polygons[:, vert, face] .= Quaternion2Vector(
                    RotateVectorbyQuaternion(
                        MyRotationQuaternion(RotateAxis, RotateAngle),
                        mesh.polygons[:, vert, face])
                )
            end
        end
        # 真反対を向いていて回転軸が定まらない時は各座標を反転
    elseif mesh.PostureDirectionVector != IGO.PostureDirectionVector
        for face in 1:mesh.face_num
            for vert in 1:mesh.dim
                mesh.polygons[:, vert, face] .= -mesh.polygons[:, vert, face]
            end
        end
        # そうでもない場合は一致しているので何もしないで
    end

    # ParallelTranslation
    for face in 1:mesh.face_num
        for vert in 1:mesh.dim
            mesh.polygons[:, vert, face] .= mesh.polygons[:, vert, face] .+ IGO.Pos
        end
    end

end

function DecideArg2D(Vec::Vector{Float64}; isAgainstX=true)
    tmpVec = Vec / norm(Vec) # Normを1に戻す
    # グローバルのx軸を基準(θ = 0)とし、反時計回りを正とする
    # -π<θ<=πで走査範囲を定義
    θ = acos(tmpVec[1])
    # y座標が負なら負になる
    if tmpVec[2] < 0
        θ = -θ
    end
    return θ
end

#= --- GameLoop & Inner caller functions--- =#

function ViDARsLoop2D(
    AllObjDict::Dict,
    center::Vector{Float64}, center_PDVec::Vector{Float64};
    θRlim=π / 4,
    θLlim=π / 4,
    θres=1000,
    CalcReflection=false)

    #println("----------------------------------------------")
    #println("Start ViDARsLoop")
    # 描画情報を格納するScanningGridをθres個のDictを持ったVectorとして生成
    ScanningGrid = Vector{Dict{String,Any}}(undef, θres)
    # 走査偏角の決定
    dθ = (θRlim + θLlim) / θres
    # 2D走査の準備
    center_PDVec[3] = 0 # Z座標を強制的にゼロにする
    # X軸正方向を基準に基準角を決定
    θstart = DecideArg2D(center_PDVec)
    #println("----------------------------------------------")
    #println("Start Scanning")
    # 反時計回りに走査を行う
    for step_θ in 1:θres
        ScanningGrid[step_θ] = Dict("Dist" => Inf)
        # 計算はベクトルを利用
        # https://risalc.info/src/line-plane-intersection-point.html
        # 走査線の方向ベクトルScanDirを導出，走査線方程式はcenter + t*ScanDir (0<t)
        ScanDir = [cos(step_θ * dθ + θstart - θRlim), sin(step_θ * dθ + θstart - θRlim), 0]
        #println("----------------------------------------------")
        #println("Start Scanning Step : $step_θ/$θres")
        for key in keys(AllObjDict)
            IGO = AllObjDict[key]
            #println("----------------------------------------------")
            #println("Start processing object : $key")
            # ここからポリゴン毎の処理
            for face_num in 1:IGO.CollisionMesh.face_num
                #println("----------------------------------------------")
                #println("Start Processing polygon No.$face_num")
                ReturnDict = Dict{String,Any}()
                polygon = IGO.CollisionMesh.polygons[:, :, face_num]
                # ポリゴンの頂点1から2へのベクトルをV_v1_v2，視点からポリゴン頂点1へのベクトルをV_c_v1とする
                V_v1_v2 = polygon[:, 2] - polygon[:, 1]
                V_c_v1 = polygon[:, 1] - center
                # 内積とノルムから角度ϕと距離hを導出，平面方程式n ⋅ x = hを導く
                θpoly = DecideArg2D(V_v1_v2)
                ϕ = DecideArg2D(-V_c_v1) - θpoly
                h = norm(V_c_v1) * sin(ϕ)
                n = [cos(θpoly - π / 2), sin(θpoly - π / 2), 0.0]
                # 交点座標は以下の通り
                tmp_dp = n ⋅ ScanDir
                is_parallel = tmp_dp == 0 ? true : false
                is_inner = false
                is_sameside = false
                if !is_parallel
                    t = (h - n ⋅ center) / tmp_dp
                    is_sameside = t >= 0 ? true : false
                    Xpt = center + ScanDir * t
                    # 交点の線内判定
                    # 両頂点からのベクトルとの内積の符号が異なる時のみ交点が線内となりその他の処理を行う
                    V_v1_Xpt = (Xpt - polygon[:, 1])
                    V_v2_Xpt = (Xpt - polygon[:, 2])
                    is_inner = dot(V_v1_Xpt, V_v2_Xpt) <= 0 ? true : false
                end
                if is_inner && is_sameside
                    ReturnDict["Dist"] = norm(Xpt - center)
                    if ScanningGrid[step_θ]["Dist"] > ReturnDict["Dist"]
                        ReturnDict["NormalVec"] = n
                        ReturnDict["FaceNum"] = face_num
                        ReturnDict["ObjectKey"] = key
                        ScanningGrid[step_θ] = ReturnDict
                    end
                end
            end
        end
    end
    return ScanningGrid
end

function DrawObjects(fig, AllObjDict::Dict, FlagBoolDict::Dict, FlagIntDict::Dict, FlagFloatDict::Dict, AxisDict::Dict,
    center::Vector{Float64}, center_PDVec::Vector{Float64}, Llim::Float64, Rlim::Float64, scanres::Int, ScanningGrid::Vector{Dict{String,Any}})
    # Visualize

    IsAnyNewAxCreated = false

    # Draw Stage
    if !FlagBoolDict["IsStageVisualised"] # Generate new Ax as Stage only if StageAx is not created
        AxisDict["StageAxDict"] = Dict("StageAx" => Axis(fig[1, 1], aspect=1), "nonplayerplots" => [], "playerplots" => [])
        #maincam = maindisplay.scene.camera
        xlims!(AxisDict["StageAxDict"]["StageAx"], -10, 10)
        ylims!(AxisDict["StageAxDict"]["StageAx"], -10, 10)
        #zlims!(AxisDict["StageAxDict"]["StageAx"], -1, 1)
        FlagBoolDict["IsStageVisualised"] = true
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

    # Draw player
    push!(AxisDict["StageAxDict"]["playerplots"], scatter!(
        AxisDict["StageAxDict"]["StageAx"],
        center[1],
        center[2],
        center[3],
        color=:orange
    ))
    leftlimitsline = Quaternion2Vector(RotateVectorbyQuaternion(MyRotationQuaternion([0, 0, 1], Llim), center_PDVec))
    rightlimitsline = Quaternion2Vector(RotateVectorbyQuaternion(MyRotationQuaternion([0, 0, 1], -Rlim), center_PDVec))

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
    if !FlagBoolDict["IsViDARsVisualised"] # Generate new Ax as ViDARs only if ViDARsAx is not created
        AxisDict["ViDARsAxDict"] = Dict("ViDARsAx" => PolarAxis(fig[1, 2], rlimits=(0, 10), thetalimits=(-Llim, Rlim)), "ViDARsResult" => [])
        FlagBoolDict["IsViDARsVisualised"] = true
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
    @time ScanningGrid = ViDARsLoop2D(
        AllObjDict::Dict, center::Vector{Float64},
        center_PDVec::Vector{Float64};
        θRlim=Rlim, θLlim=Llim, θres=scanres)

    # Draw display
    @time DrawObjects(fig, AllObjDict, FlagBoolDict, FlagIntDict, FlagFloatDict, AxisDict,
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
    AllObjDict = Dict{Int,InGameObj}()
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
    polynum = 30

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
    Mygameobj = InGameObj(
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
        player_PDVec = Quaternion2Vector(RotateVectorbyQuaternion(MyRotationQuaternion([0, 0, 1], π / 120), player_PDVec))
        println("------------")
        GameLoop(fig, AllObjDict, FlagBoolDict, FlagIntDict, FlagFloatDict, AxisDict, player_pos, player_PDVec)
        println("------------")
    end
end

GameMain()

#= ------ Test of 3D Translation ------ =#
#=
player_PDVec = [1.0, 0.0, 0.0]
MyRotationQuaternion([0, 0, 1], π / 6)
RotateVectorbyQuaternion(MyRotationQuaternion([0, 0, 1], π / 3), player_PDVec)

Quaternion2Vector(RotateVectorbyQuaternion(MyRotationQuaternion([0, 0, 1], π / 4), player_PDVec))
=#

#= ------ Speed Test ------=#
#=
tmpvertices2D = [
        0.0 0.0 0.0;
        1.0 0.0 0.0;
        0.0 1.0 0.0;
        1.0 1.0 0.0;
        -1.0 0.0 0.0
    ]

    tmpfaces2D = [
        1 2;
        1 3;
        2 4;
        5 1;
        3 2;
        4 1;
        5 4;
        4 3
    ]
tmpvertices2D = rand(3, 20)
tmpfaces2D = 10 .* (1 .+ rand(2, 10))
tmpfaces2D = map(tmpfaces2D) do x
    floor(x)
end
tmpfaces2D = map(tmpfaces2D) do x
    convert(Int, x)
end

tmpMeshPDVec = [0.0, 1.0, 0.0]
Mymesh = CollisionMesh(tmpvertices2D, tmpfaces2D, tmpMeshPDVec, 2)

function tmp1()
    for i in 1:100
        DecomposeMesh!(Mymesh)
    end
end
=#

#=
#Mygameobj.PostureDirectionVector = [1.0, 0.0, 0.0]
# いやこれ勝手にノルム1のベクトルを外から書き換えられるのまずいな


事前に配列サイズを指定するのとしないのとの差
for i in 1:100
    @time a = MakeScanningGrid(i * 100, π / 2, π / 2)
end
# 0.001725 seconds (θres = 10000, alloc ∼2.5kB)

for i in 1:100
    @time a = MakeScanningGrid_static(i * 100, π / 2, π / 2)
end
# 0.000457 seconds (θres = 10000, alloc ∼1kB)

何故かlinearindexの方が数割ほど遅いしメモリも食う
=#