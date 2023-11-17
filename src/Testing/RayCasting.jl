using GLMakie
using LinearAlgebra

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
MyRotationQuaternion(vec::Vector, θ) = MyRotationQuaternion(vec[1] * cos(θ / 2), vec[2] * cos(θ / 2), vec[3] * cos(θ / 2), sin(θ / 2))

#= --- Original Structures' functions --- =#
function RotateQuaternion(Rotator::MyAbstractQuaternion, Rotated::MyAbstractQuaternion)
    return MyQuaternion(
        Rotator.θ * Rotated.x - Rotator.z * Rotated.y + Rotator.y * Rotated.z + Rotator.x * Rotated.θ,
        Rotator.z * Rotated.x + Rotator.θ * Rotated.y - Rotator.x * Rotated.z + Rotator.y * Rotated.θ,
        -Rotator.y * Rotated.x + Rotator.x * Rotated.y + Rotator.θ * Rotated.z + Rotator.z * Rotated.θ,
        -Rotator.x * Rotated.x - Rotator.y * Rotated.y - Rotator.z * Rotated.z + Rotator.θ * Rotated.θ
    )
end

function RotateVectorbyQuaternion(q::MyAbstractQuaternion, vec::Vector)
    return MyQuaternion(
        q.θ * vec[1] - q.z * vec[2] + q.y * vec[3] + q.x * 0,
        q.z * vec[1] + q.θ * vec[2] - q.x * vec[3] + q.y * 0,
        -q.y * vec[1] + q.x * vec[2] + q.θ * vec[3] + q.z * 0,
        -q.x * vec[1] - q.y * vec[2] - q.z * vec[3] + q.θ * 0
    )
end

function Quaternion2Vector(q::MyQuaternion)
    v = Vector{Float64}(undef, 3)
    v[1] = q.x
    v[2] = q.y
    v[3] = q.z
    return v
end

function Vector2Quarternion(v::Vector{Real})
    return MyQuaternion(v[1], v[2], v[3], 0)
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

# NormalVectorの向きはZ軸正側からXY平面を見てv2-v1をCounterClockwiseに回した方向とする
# つまり、コリジョンメッシュを作るときは同方向から見て時計回りに頂点番号を振らないと
# 法線が内向きになってしまう
function CalcNormalVector2D(polygon2D::Matrix{Float64})
    normalilzed_edge = (polygon2D[:, 2] - polygon2D[:, 1]) / norm(polygon2D[:, 2] - polygon2D[:, 1])
    NormalVector = cross([0.0, 0.0, 1.0], normalilzed_edge)
    return NormalVector / norm(NormalVector)
end

function CalcDistAndNormalizedvec(eyepos::Vector, polygon::Matrix{Float64})
    v1 = polygon[:, 1] - eyepos
    v2 = polygon[:, 2] - eyepos
    v3 = polygon[:, 3] - eyepos
    norm_v1 = norm(v1)
    norm_v2 = norm(v2)
    norm_v3 = norm(v3)
    NormalizedV1 = v1 / norm_v1
    NormalizedV2 = v2 / norm_v2
    NormalizedV3 = v3 / norm_v3
    return [norm_v1, norm_v2, norm_v3, NormalizedV1, NormalizedV2, NormalizedV3]
end

# 基本的に2D、3Dを問わずポリゴンは3次元座標系で保持される為流用可能
function TranslatePolygonLocal2Global!(mesh, IGO::InGameObj)
    DecomposeMesh!(mesh)
    # Local2Global Coordinate Translation
    # ParallelTranslation
    for face in 1:mesh.face_num
        for vert in 1:mesh.dim
            mesh.polygons[:, vert, face] .= mesh.polygons[:, vert, face] .+ IGO.Pos
        end
    end
    # Roation (by Quarternion[λx*cos(θ/2),λy*cos(θ/2),λz*cos(θ/2),sin(θ/2)])
    RotateAxis = cross(mesh.PostureDirectionVector, IGO.PostureDirectionVector)
    # 通常は回転軸を中心に回転
    if RotateAxis != 0
        RotateAxis = RotateAxis / norm(RotateAxis)
        # CounterClockwiseを仮定
        RotateAngle = acos(mesh.PostureDirectionVector ⋅ IGO.PostureDirectionVector)
        for face in 1:mesh.face_num
            for vert in 1:mesh.dim
                println("Axis vec : $RotateAxis")
                println("Angle : $RotateAngle")
                mesh.polygons[:, vert, face] .= Quaternion2Vector(
                    RotateVectorbyQuaternion(
                        MyRotationQuaternion(RotateAxis, RotateAngle),
                        mesh.polygons[:, vert, face])
                )
            end
        end
        # 真反対を向いていて回転軸が定まらない時は各座標を反転
    elseif acos(mesh.PostureDirectionVector ⋅ IGO.PostureDirectionVector) == -1
        for face in 1:mesh.face_num
            for vert in 1:mesh.dim
                mesh.polygons[:, vert, face] .= -mesh.polygons[:, vert, face]
            end
        end
        # そうでもない場合は一致しているので何もしないで
    end
end

#= --- GameLoop & Inner caller functions--- =#

function ViDARsLoop2D(
    AllObjDict::Dict,
    center::Vector{Float64}, center_PDVec::Vector{Float64};
    θRlim=π / 4,
    θLlim=π / 4,
    θres=1000,
    CalcReflection=false)

    println("----------------------------------------------")
    println("Start ViDARsLoop")
    # 描画情報を格納するScanningGridをθres個のDictを持ったVectorとして生成
    ScanningGrid = Vector{Dict{String,Any}}(undef, θres)
    # 走査角の決定
    dθ = (θRlim + θLlim) / θres
    # 2D走査の準備
    center_PDVec[3] = 0 # Z座標を強制的にゼロにする
    center_PDVec = center_PDVec / norm(center_PDVec) # Normを1に戻す
    # グローバルのx軸を基準(θ = 0)とし、反時計回りを正とする
    # -π<θ<πで走査範囲を定義
    θstart = acos(center_PDVec[1])
    # y座標が負なら負になる
    if center_PDVec[2] < 0
        θstart = -θstart
        # θ=πとθ=0はここで区別する
    elseif center_PDVec[2] == 0 && center_PDVec[1] < 0
        θstart = π
    end
    println("----------------------------------------------")
    println("Start Scanning")
    # 反時計回りに走査を行う
    for step_θ in 1:θres
        ScanningGrid[step_θ] = Dict("Dist" => Inf)
        # 走査点のグローバル座標を導出
        scanning_pointX = center[1] + cos(step_θ * dθ + θstart - θRlim)
        scanning_pointY = center[2] + sin(step_θ * dθ + θstart - θRlim)
        # 走査線方程式の係数を導出
        # y = α_1 *X + β_1
        α_1 = (scanning_pointY - center[2]) / (scanning_pointX - center[1])
        β_1 = scanning_pointY - α_1 * scanning_pointX
        println("----------------------------------------------")
        println("Start Scanning Step : $step_θ/$θres")
        for key in keys(AllObjDict)
            IGO = AllObjDict[key]
            # ViDARsLoop内ではPolygonを直接変形する処理はしないとしてcopyを取っている
            polygons = IGO.CollisionMesh.polygons
            println("----------------------------------------------")
            println("Start processing object : $key")
            # ここからポリゴン毎の処理
            for face_num in 1:IGO.CollisionMesh.face_num
                println("----------------------------------------------")
                println("Start Processing polygon No.$face_num")
                ReturnDict = Dict{String,Any}()
                polygon = polygons[:, :, face_num]
                # y = α_2 *X + β_2
                α_2 = (polygon[2, 2] - polygon[2, 1]) / (polygon[1, 2] - polygon[1, 1])
                β_2 = polygon[2, 1] - α_2 * polygon[1, 1]
                # 解はX = (β_2 - β_1)/(α_1 - α_2)
                # 分母0によるDiv0で発生するNaN回避が必要
                if α_1 != α_2
                    CrossingPointX = (β_2 - β_1) / (α_1 - α_2)
                    CrossingPointY = CrossingPointX * α_1 + β_1
                    # 交点の線内判定
                    VecVert2Vert = polygon[:, 2] - polygon[:, 1]
                    VecVert2CrossPoint = [CrossingPointX, CrossingPointY, 0.0] - polygon[:, 1]
                    DotProduct = dot(VecVert2CrossPoint, VecVert2Vert)
                    # 内積が0<=DotProduct<=1の時のみ交点が線内となりその他の処理を行う
                    if 0 <= DotProduct <= 1
                        ReturnDict["Dist"] = norm([(CrossingPointX - scanning_pointX), (CrossingPointY - scanning_pointY)])
                        if ScanningGrid[step_θ]["Dist"] > ReturnDict["Dist"]
                            ReturnDict["NormalVec"] = CalcNormalVector2D(polygon)
                            ReturnDict["FaceNum"] = face_num
                            ReturnDict["ObjectKey"] = key
                            ScanningGrid[step_θ] = ReturnDict
                        end
                    end
                end
            end
        end
    end
    return ScanningGrid
end

# Call function which should be done in a frame
function EachFrame(fig, AllObjDict::Dict,
    center::Vector{Float64}, center_PDVec::Vector{Float64})
    println("----------------------------------------------")
    println("EachFrameExecution Started")

    # Polygon projection
    for key in keys(AllObjDict)
        IGO = AllObjDict[key]
        println("----------------------------------------------")
        println("Start processing object : $key")
        TranslatePolygonLocal2Global!(IGO.CollisionMesh, IGO)
        println("Finish translating CollisionMesh (object : $key)")
        #TranslatePolygonLocal2Global!(IGO.AppearanceMesh, IGO)
        println("Finish translating AppearanceMesh (object : $key)")

    end

    Llim = π / 4
    Rlim = π / 4
    scanres = 1000

    # 上記関数でローカル座標情報から生成されたグローバル座標上のポリゴンを処理
    ScanningGrid = ViDARsLoop2D(
        AllObjDict::Dict, center::Vector{Float64},
        center_PDVec::Vector{Float64};
        θRlim=Rlim, θLlim=Llim, θres=scanres)

    # Visualize

    # Draw Maindisplay
    maindisplay = Axis3(fig[1, 1], aspect=(1, 1, 1))
    maincam = maindisplay.scene.camera
    xlims!(maindisplay, -5, 5)
    ylims!(maindisplay, -5, 5)
    zlims!(maindisplay, -5, 5)

    facenum = AllObjDict[1].CollisionMesh.face_num

    for i in 1:facenum
        lines!(
            maindisplay,
            AllObjDict[1].CollisionMesh.polygons[1, :, i],
            AllObjDict[1].CollisionMesh.polygons[2, :, i],
            AllObjDict[1].CollisionMesh.polygons[3, :, i],
            color=:lightblue
        )
    end

    # point players
    scatter!(
        maindisplay,
        center[1],
        center[2],
        center[3],
        color=:orange)
    lines!(
        maindisplay,
        [center[1], center[1] + center_PDVec[1]],
        [center[2], center[2] + center_PDVec[2]],
        [center[3], center[3] + center_PDVec[3]],
        color=:lightgreen
    )
    # Draw PolarAxis
    polax = PolarAxis(fig[1, 2])
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
    scatobject = scatter!(polax, θarr, rarr)
    display(fig)
end

# Main Game loop
function GameLoop(fig, AllObjDict::Dict,
    center::Vector{Float64}, center_PDVec::Vector{Float64})
    println("----------------------------------------------")
    println("GameLoop Started")
    EachFrame(fig, AllObjDict, center, center_PDVec)
end

# InitialSetting to start game
function InitialSetting()
    fig = Figure()
    display(fig)
    println("InitialSetting is finished without any problem")
    return fig
end

#= --- Comment Zone ---=#

#= --- TestCodes --- =#

AllObjDict = Dict()

function Main_()
    fig = InitialSetting()
    # 全てのオブジェクトは作られたのちにこのDictに追記されることで初めて名前と存在をゲームから認められる
    AllObjDict = Dict{Int,InGameObj}()
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
    tmpvertices2D = permutedims(tmpvertices2D)
    tmpfaces2D = permutedims(tmpfaces2D)

    tmpMeshPDVec = [0.0, 1.0, 0.0]
    Mymesh = CollisionMesh(tmpvertices2D, tmpfaces2D, tmpMeshPDVec, 2)
    tmpIGOpos = [2.0, 2.0, 0.0]
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
    player_PDVec = [1.0, 0.0, 1.0]
    GameLoop(fig, AllObjDict, player_pos, player_PDVec)
end

Main_()


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
tmpvertices2D = rand(3, 2000)
tmpfaces2D = 1000 .* (1 .+ rand(2, 100000))
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