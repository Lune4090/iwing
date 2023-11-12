using GLMakie
using LinearAlgebra

#= --- Original Structures Definition --- =#
abstract type MyAbstractQuaternion end
abstract type MyAbstractNormalizedQuaternion <: MyAbstractQuaternion end

mutable struct MyQuaternion <: MyAbstractQuaternion
    x::Float64
    y::Float64
    z::Float64
    θ::Float64
end

MyQuaternion(x::Real, y::Real, z::Real, θ::Real) = MyQuaternion(convert(Float64, x), convert(Float64, y), convert(Float64, z), convert(Float64, θ))
MyQuaternion(vec::Vector) = MyQuaternion(vec[1], vec[2], vec[3], 0.0)

mutable struct MyNormalizedQuaternion <: MyAbstractNormalizedQuaternion
    x::Float64
    y::Float64
    z::Float64
    θ::Float64
    MyNormalizedQuaternion(x, y, z, θ) = x^2 + y^2 + z^2 != 1 ? error("RotationQuaternion's Norm != 1!") : new(x, y, z, θ)
end

# Int->Float64 Conversion
MyNormalizedQuaternion(x::Int, y::Int, z::Int, θ::Int) = MyNormalizedQuaternion(convert(Float64, x), convert(Float64, y), convert(Float64, z), convert(Float64, θ))

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

@kwdef mutable struct CollisionMesh
    Vertices::Matrix{Float64} # ここでローカル座標が決まる
    Faces::Matrix{Int}
end

@kwdef mutable struct InGameObj
    Pos::Vector{Float64} # これがグローバル座標
    PostureDirectionVector::Vector{Float64}
    RotationQuaternion::MyNormalizedQuaternion
    CollisionMesh::CollisionMesh
    AppearanceMesh::Any
    # AppearanceMeshは外部ツールで製作した.stl形式のオブジェクトのloadを想定
    # AppearanceMeshは触らない(Makie.meshscatter!()にmarkerとして渡すだけ)
    Attributes::Array{Dict,1}
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
    InGameObj(Pos, PostureDirectionVector, RotationQuaternion, CollisionMesh, AppearanceMesh, Attributes) = PostureDirectionVector[1]^2 + PostureDirectionVector[2]^2 + PostureDirectionVector[3]^2 != 1 ? error("PostureDirectionVector should be equal to 1!!!") : new(Pos, PostureDirectionVector, RotationQuaternion, CollisionMesh, AppearanceMesh, Attributes)
end

# IGO間の関係性がメインのコンテナ
@kwdef mutable struct StructuredInGameObj
    InGameObjArr::Vector{InGameObj}
    InGameObjsRelation
end

function CalcNormalVector(polygon::Matrix{Float64})
    # ここ、2Dポリゴンを3Dポリゴンの3つ目の頂点を2つ目に合わせることで実装している影響で代えられない
    vec_a = polygon[2, :] - polygon[1, :]
    vec_b = polygon[3, :] - polygon[1, :]
    NormalVector = cross(vec_a, vec_b)
    return NormalVector / norm(NormalVector)
end

function CalcDistPoint2Polygon(eyepos::Vec{3,Float32}, polygon::Matrix{Float64})
    return norm(transpose(eyepos) - reduce(+, polygon; dims=1))
end

function CalcDistNormalizedvec(eyepos::Vector, polygon::Matrix{Float64})
    v1 = polygon[1, :] - eyepos
    v2 = polygon[2, :] - eyepos
    v3 = polygon[3, :] - eyepos
    norm_v1 = norm(v1)
    norm_v2 = norm(v2)
    norm_v3 = norm(v3)
    NormalizedV1 = v1 / norm_v1
    NormalizedV2 = v2 / norm_v2
    NormalizedV3 = v3 / norm_v3
    return [norm_v1, norm_v2, norm_v3, NormalizedV1, NormalizedV2, NormalizedV3]
end

# VerticesからFacesに面毎に格納された頂点番号を引数にpolygonを生成
# Facesがn*3のArrayであることを仮定してるから注意
function DecomposeMesh(mesh::CollisionMesh)
    num_faces = size(mesh.Faces, 1)
    polygons = Array{Float64}(undef, 3, 3, num_faces)
    for EachFace in 1:num_faces
        for EachVertex in 1:3
            polygons[EachVertex, :, EachFace] .= mesh.Vertices[mesh.Faces[EachFace, EachVertex], :]
        end
    end
    return polygons, num_faces
end

# 現時点ではxy平面上の走査グリッド(つまり線)生成しか実装していない
function MakeScanningGrid(θres::Int, θRlim::Float64, θLlim::Float64; GridDim=1)
    dθ = (θRlim + θLlim) / θres
    if !(GridDim == 1 || GridDim == 2)
        error("dim should be 1 or 2")
    end
    if θRlim >= π || θRlim <= 0
        error("θRlim should be 0 < θ < π")
    end
    if θLlim >= π || θLlim <= 0
        error("θLlim should be 0 < θ < π")
    end

    if GridDim == 1
        ScanningGrid = []
        step_θ = 1
        while step_θ <= θres
            #println("scanningstep: $step_θ")
            # Grid is filled as Clockwise
            # θ,cos(nθ),sin(nθ),cos(n+1θ),cos(n+1θ)
            append!(ScanningGrid, [
                step_θ * dθ,
                cos(step_θ * dθ - θLlim), sin(step_θ * dθ - θLlim),
                cos((step_θ + 1) * dθ - θLlim), cos((step_θ + 1) * dθ - θLlim)
            ])
            step_θ += 1
        end
    end
    if GridDim == 2
        error("sorry, 2D grid is not ready...")
    end
    return ScanningGrid
end

# Grid should be filled as CounterClockwise
function MakeScanningGrid_static(θres::Int, θRlim::Float64, θLlim::Float64; GridDim=1)
    #dθ = (θRlim + θLlim) / θres
    if !(GridDim == 1 || GridDim == 2)
        error("dim should be 1 or 2")
    end
    if θRlim >= π || θRlim <= 0
        error("θRlim should be 0 < θ < π")
    end
    if θLlim >= π || θLlim <= 0
        error("θLlim should be 0 < θ < π")
    end

    if GridDim == 1

        ScanningGrid = zeros(θres)
        #=
                step_θ = 1
                while step_θ <= θres
                    #println("scanningstep: $step_θ")
                    # θ,cos(nθ),sin(nθ),cos(n+1θ),cos(n+1θ)
                    ScanningGrid[step_θ] = [
                        step_θ * dθ,
                        cos(step_θ * dθ - θLlim), sin(step_θ * dθ - θLlim),
                        cos((step_θ + 1) * dθ - θLlim), cos((step_θ + 1) * dθ - θLlim)
                    ]
                    step_θ += 1
                end
        		=#
    end
    if GridDim == 2
        error("sorry, 2D grid is not ready...")
    end
    return ScanningGrid
end

function ProjectPolygons()

end


# 以下のループをR-bufferに書き直す(まあほぼ変わらない)
function ViDARsLoop(center::Vector{Float64}, center_PDVec::Vector{Float64}, RayCastedObject::StructuredInGameObj;
    DetectionLengthUpperLim::Float64=300.0,
    DetectionAngleRightLim::Float64=π / 4,
    DetectionAngleLeftLim::Float64=π / 4,
    θres=1000,
    CalcReflection=false)

    println("----------------------------------------------")
    # θres個のVector型の要素を持ったVectorとしてScanningGridを受け取る
    dθ = (DetectionAngleRightLim + DetectionAngleLeftLim) / θres
    # 一旦愚直にforで実装
    for IGO in RayCastedObject.InGameObjArr
        println("----------------------------------------------")
        println("Start processing IGO")
        println("----------------------------------------------")
        polygons, num_faces = DecomposeMesh(IGO.CollisionMesh)
        # Local2Global Coordinate Translation
        polygons = polygons .+ IGO.Pos
        # Rotation (ここ違うから再度変更の必要あり)
        #=
                for i in 1:num_faces
                    for j in 1:3
                        polygons[:, j, i] = Quaternion2Vector(RotateVectorbyQuaternion(
                            IGO.RotationQuaternion, polygons[:, j, i]))
                    end
                end
        		=#
        NewAttributes = []

        # ここからポリゴン毎の処理
        for face_num in 1:num_faces
            println("----------------------------------------------")
            println("Start Processing polygon: $face_num")
            println("----------------------------------------------")
            ReturnDict = Dict()

            # θres個のVector型の要素を持ったVectorとしてScanningGridを受け取る
            ReturnDict["ScanningGrid"] = MakeScanningGrid_static(θres, DetectionAngleRightLim, DetectionAngleLeftLim)
            # Dist, norm算出
            ReturnDict["DisV1"], ReturnDict["DisV2"], ReturnDict["DisV3"], ReturnDict["V1"], ReturnDict["V2"], ReturnDict["V3"] = CalcDistNormalizedvec(center, polygons[:, :, face_num])

            ReturnDict["NormVec"] = [0.0, 0.0, 0.0]
            # 距離で捜査範囲内か確認
            if ReturnDict["DisV1"] > DetectionLengthUpperLim &&
               ReturnDict["DisV2"] > DetectionLengthUpperLim &&
               ReturnDict["DisV3"] > DetectionLengthUpperLim
                println("Too far to search")
            else
                # Polygon毎に法線ベクトルを算出
                ReturnDict["NormVec"] = CalcNormalVector(polygons[:, :, face_num])

                #= --- ここから2Dの走査処理に特化した内容に変わる ---=#

                center_PDVec[3] = 0 # Z座標を強制的にゼロにする
                center_PDVec = center_PDVec / norm(center_PDVec) # Normを1に戻す

                # 各頂点についての処理対象判定を開始
                is_vert1_inner = false
                is_vert1_outer_left = false
                is_vert1_outer_right = false
                is_vert2_inner = false
                is_vert2_outer_left = false
                is_vert2_outer_right = false
                θ1 = 0
                θ2 = 0
                is_calc_valid = true # NaN発生時の計算回避用
                println("polygons: $(polygons[:, :, face_num])")
                for vertex_num in 1:2
                    # 2Dであることを利用、外積のZ座標が正なら視線の右側、負なら左側にいる
                    CrossProductZ = cross(ReturnDict["V$vertex_num"], center_PDVec)[3]
                    is_RightSide::Bool = false
                    if CrossProductZ >= 0
                        is_RightSide = true
                    end
                    # 内積を取って走査範囲内か確認
                    DotProduct = dot(ReturnDict["V$vertex_num"], center_PDVec)
                    if DotProduct == NaN
                        print("WARNING::center and calculated polygon is in the samepoint!!")
                        is_calc_valid = false
                    elseif vertex_num == 1
                        # 右側にある場合は右側のcosと内積を比較し大きければ範囲内、左側は左側と
                        if (is_RightSide && DotProduct > cos(DetectionAngleRightLim)) || (!is_RightSide && DotProduct > cos(DetectionAngleLeftLim))
                            is_vert1_inner = true
                        else
                            if is_RightSide
                                is_vert1_outer_right = true
                            else
                                is_vert1_outer_left = true
                            end
                        end
                        println(ReturnDict["V$vertex_num"])
                        println("$(center_PDVec)")
                        println("Dot = $DotProduct")
                        println(cos(DetectionAngleRightLim))
                        println(cos(DetectionAngleLeftLim))
                        println("is_vert1_inner: $is_vert1_inner")
                        if is_RightSide
                            θ1 = -acos(DotProduct)
                        else
                            θ1 = acos(DotProduct)
                        end
                    elseif vertex_num == 2
                        if (is_RightSide && DotProduct > cos(DetectionAngleRightLim)) || (!is_RightSide && DotProduct > cos(DetectionAngleLeftLim))
                            is_vert2_inner = true
                        else
                            if is_RightSide
                                is_vert2_outer_right = true
                            else
                                is_vert2_outer_left = true
                            end
                        end
                        println(ReturnDict["V$vertex_num"])
                        println("$(center_PDVec)")
                        println("Dot = $DotProduct")
                        println(cos(DetectionAngleRightLim))
                        println(cos(DetectionAngleLeftLim))
                        println("is_vert2_inner: $is_vert2_inner")
                        if is_RightSide
                            θ2 = -acos(DotProduct)
                        else
                            θ2 = acos(DotProduct)
                        end
                    end
                end
                # 両方の頂点が走査範囲外ならここで走査終了
                # そうでないもののみ探知対象ピクセルを決定
                if (is_vert1_inner || is_vert2_inner || is_vert1_outer_left && is_vert2_outer_right) && is_calc_valid
                    # 各ピクセルに走査対象のポリゴンの属するIGO、法線ベクトルを与える
                    # Grid should be filled as CounterClockwise
                    # ここでDepthBufferingを行う
                    # 本来のグラフィクスパイプラインは、一度Bufferに描画対象ポリゴンを全て蓄える
                    # これは、メッシュの中には透明なものもあり、その影響を考えると
                    # この段階では全て格納して置いて後段のピクセルパイプラインに渡し、
                    # そこで各種計算とともに要らないものを捨てて最終的な絵(ポストエフェクト前)を出す
                    # という方法が好ましいから
                    # ViDARsのパイプラインはあくまで簡易的なメッシュ検出を目的としたもので、
                    # あくまでポリゴンは一番手前のもののみが「映る」ものと考えるが、
                    # 将来的な拡張と機能分離の観点からReturnDictにポリゴン毎にBufferを製造して格納するという実装を行った
                    grid_num1 = convert(Int64, (θ1 + DetectionAngleRightLim) ÷ dθ)
                    grid_num2 = convert(Int64, (θ2 + DetectionAngleRightLim) ÷ dθ)
                    println(grid_num1)
                    println(grid_num2)
                    if grid_num1 > grid_num2
                        grid_Llim = grid_num1
                        grid_Rlim = grid_num2
                        depth_Llim = ReturnDict["DisV1"]
                        depth_Rlim = ReturnDict["DisV2"]
                    else
                        grid_Llim = grid_num2
                        grid_Rlim = grid_num1
                        depth_Llim = ReturnDict["DisV2"]
                        depth_Rlim = ReturnDict["DisV1"]
                    end
                    # LinRangeでDepthを変更
                    println("----------------------------------------------")
                    println("grid_Llim: $grid_Llim")
                    println("grid_Rlim: $grid_Rlim")
                    if grid_Llim - grid_Rlim > 0
                        println(LinRange(depth_Rlim, depth_Llim, grid_Llim - grid_Rlim + 1))
                        DepthData = LinRange(depth_Rlim, depth_Llim, grid_Llim - grid_Rlim + 1)
                        ReturnDict["ScanningGrid"][max(grid_Rlim, 0):min(grid_Llim, θres)] = DepthData[max(0, grid_Rlim), min(grid_Llim, θres)]
                    else
                        println("No change to ScanningGrid")
                        println("----------------------------------------------")
                    end
                end
            end
            push!(NewAttributes, ReturnDict)
        end
        println(NewAttributes)
        IGO.Attributes = NewAttributes
    end
end

#= --- Comment Zone ---=#
#=
2Dに直す!!!!!
=#


#= --- TestCodes --- =#
tmpvertices = [
    0.0 0.0 0.0;
    1.0 0.0 3.0;
    1.0 1.0 -1.0;
    0.0 1.0 2.0;
    3.0 2.0 1.0;
    0.5 2.5 1.0;
    2.0 4.0 -2.0
]
tmpvertices_Z0 = [
    0.0 0.0 0.0;
    1.0 0.0 0.0;
    1.0 1.0 0.0;
    0.0 1.0 0.0;
    3.0 2.0 0.0;
    0.5 2.5 0.0;
    2.0 4.0 0.0
]
tmpfaces = [
    1 2 3;
    3 4 1;
    1 4 2;
    6 7 1;
    3 2 5;
    5 4 2
]
tmppos = [0.0, 0.0, 0.0]
tmpquarternion = MyNormalizedQuaternion(1, 0, 0, 0)
Mymesh = CollisionMesh(tmpvertices_Z0, tmpfaces)
tmpdict = Dict([("A", 1)])

Mygameobj = InGameObj(
    Pos=tmppos,
    PostureDirectionVector=[1.0, 0.0, 0.0],
    RotationQuaternion=tmpquarternion,
    CollisionMesh=Mymesh,
    AppearanceMesh=Nothing,
    Attributes=[tmpdict])

mesh(Mygameobj.CollisionMesh.Vertices, Mygameobj.CollisionMesh.Faces; color=:lightblue, shading=true)

igos = StructuredInGameObj([Mygameobj], Nothing)

polygons, num_faces = DecomposeMesh(Mygameobj.CollisionMesh)

NewAttributes = []

fig = Figure()
ax3d = Axis3(fig[1, 1], aspect=(1, 1, 1))
cam = ax3d.scene.camera

xlims!(ax3d, -3, 3)
ylims!(ax3d, -3, 3)
zlims!(ax3d, -3, 3)

center = [0.0, 0.0, 0.0]
center_PDVec = [1.0, 0.0, 0.0]

ViDARsLoop(center, center_PDVec, igos)

#igos.InGameObjArr[1].Attributes

mesh!(ax3d, igos.InGameObjArr[1].CollisionMesh.Vertices, igos.InGameObjArr[1].CollisionMesh.Faces; color=:lightblue)

xs = [igos.InGameObjArr[1].Attributes[1]["NormVec"][1], 0.0]
ys = [igos.InGameObjArr[1].Attributes[1]["NormVec"][2], 0.0]
zs = [igos.InGameObjArr[1].Attributes[1]["NormVec"][3], 0.0]

lines!(ax3d, xs, ys, zs; color=:orange)

display(fig)

#=
for i in 1:100
    @time a = MakeScanningGrid(i * 100, π / 2, π / 2)
end
# 0.001725 seconds (θres = 10000, alloc ∼2.5kB)

for i in 1:100
    @time a = MakeScanningGrid_static(i * 100, π / 2, π / 2)
end
# 0.000457 seconds (θres = 10000, alloc ∼1kB)
=#
#=
IGO = Mygameobj
polygons, num_faces = DecomposeMesh(IGO.CollisionMesh)
# Global coordinate Translation
# Parallel Translation
println(IGO.Pos)
println(IGO.RotationQuaternion)
IGO.Pos = [1.0, 2.0, 3.0]

polygons = polygons .+ IGO.Pos
# Rotation
for i in 1:num_faces
    for j in 1:3
        polygons[:, j, i] = Quaternion2Vector(RotateVectorbyQuaternion(
            IGO.RotationQuaternion, polygons[:, j, i]))
    end
end
=#

#Mygameobj.PostureDirectionVector = [1.0, 0.0, 0.0]
# いやこれ勝手にノルム1のベクトルを外から書き換えられるのまずいな