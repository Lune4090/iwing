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
        q.θ * vec.x - q.z * vec.y + q.y * vec.z + q.x * vec.θ,
        q.z * vec.x + q.θ * vec.y - q.x * vec.z + q.y * vec.θ,
        -q.y * vec.x + q.x * vec.y + q.θ * vec.z + q.z * vec.θ,
        -q.x * vec.x - q.y * vec.y - q.z * vec.z + q.θ * vec.θ
    )
end

@kwdef mutable struct CollisionMesh
    Vertices::Matrix{Float64} # ここでローカル座標が決まる
    Faces::Matrix{Int}
end

@kwdef mutable struct InGameObj
    Pos::Vector{Float64} # これがグローバル座標
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
end

# IGO間の関係性がメインのコンテナ
@kwdef mutable struct StructuredInGameObj
    InGameObjArr::Vector{InGameObj}
    InGameObjsRelation
end

# 2Dに直す!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

function CalcNormVec(polygon::Matrix{Float64})
    vec_a = polygon[1, :] - polygon[2, :]
    vec_b = polygon[3, :] - polygon[2, :]
    NormalVector = cross(vec_a, vec_b)
    return NormalVector / norm(NormalVector)
end

function CalcDistPoint2Polygon(eyepos::Vec{3,Float32}, polygon::Matrix{Float64})
    return norm(transpose(eyepos) - reduce(+, polygon; dims=1))
end

function Eyepos2Polygon(eyepos::Vector, polygon::Matrix{Float64})
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

function MakeScanningGrid_static(θres::Int, θRlim::Float64, θLlim::Float64; GridDim=1)
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
        ScanningGrid = Vector{Vector{Float64}}(undef, θres)
        step_θ = 1
        while step_θ <= θres
            #println("scanningstep: $step_θ")
            # Grid is filled as Clockwise
            # θ,cos(nθ),sin(nθ),cos(n+1θ),cos(n+1θ)
            ScanningGrid[step_θ] = [
                step_θ * dθ,
                cos(step_θ * dθ - θLlim), sin(step_θ * dθ - θLlim),
                cos((step_θ + 1) * dθ - θLlim), cos((step_θ + 1) * dθ - θLlim)
            ]
            step_θ += 1
        end
    end
    if GridDim == 2
        error("sorry, 2D grid is not ready...")
    end
    return ScanningGrid
end

function ProjectePolygons()

end


# 以下のループをR-bufferに書き直す(まあほぼ変わらない)
function ViDARsLoop(center::Vector{Float64}, RayCastedObject::StructuredInGameObj; CalcReflection=true)
    # 一旦愚直にforで実装
    for IGO in RayCastedObject.InGameObjArr
        polygons, num_faces = DecomposeMesh(IGO.CollisionMesh)
        NewAttributes = []
        # ブロードキャストできるはず
        for idx in 1:num_faces
            ReturnDict = Dict()
            ReturnDict["Dis"] = Eyepos2Polygon(center, polygons[:, :, idx])
            ReturnDict["NormVec"] = CalcNormVec(polygons[:, :, idx])
            #=
            if CalcReflection
              ReturnDict["Reflection"] =  CalcReflectionPoint2Plain(center::Point3f, NormVec)
            end
            if Calc
            =#
            push!(NewAttributes, ReturnDict)
        end
        println(NewAttributes)
        IGO.Attributes = NewAttributes
    end
end
function ViDARsLoop_tmp(cam::Camera, RayCastedObject::StructuredInGameObj; CalcReflection=true)
    # 一旦愚直にforで実装
    for IGO in RayCastedObject.InGameObjArr
        polygons, num_faces = DecomposeMesh(IGO.CollisionMesh)
        NewAttributes = []
        # ブロードキャストできるはず
        for idx in num_faces
            ReturnDict = Dict()
            cam_pos = to_value(cam.eyeposition)
            println("camera position:$cam_pos")

            ReturnDict["Dis"] = CalcDistPoint2Polygon(cam_pos, polygons[:, :, idx])
            ReturnDict["NormVec"] = CalcNormVec(polygons[:, :, idx])
            #=
            if CalcReflection
              ReturnDict["Reflection"] =  CalcReflectionPoint2Plain(center::Point3f, NormVec)
            end
            if Calc
            =#
            push!(NewAttributes, ReturnDict)
        end
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
    0.0 1.0 2.0
]
tmpfaces = [
    1 2 3;
    3 4 1
]
tmppos = [0.0, 0.0, 0.0]
tmpquarternion = MyNormalizedQuaternion(1, 0, 0, 0)
Mymesh = CollisionMesh(tmpvertices, tmpfaces)
tmpdict = Dict([("A", 1)])

Mygameobj = InGameObj(
    Pos=tmppos,
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

center = [8.0, 3.0, 1.0]

ViDARsLoop(center, igos)
println(igos.InGameObjArr[1].Attributes)
mesh!(ax3d, igos.InGameObjArr[1].CollisionMesh.Vertices, igos.InGameObjArr[1].CollisionMesh.Faces; color=:lightblue)
println(igos.InGameObjArr[1].Attributes[1]["NormVec"])
println(igos.InGameObjArr[1].CollisionMesh.Vertices[1, :])

xs = [igos.InGameObjArr[1].Attributes[1]["NormVec"][1], 0.0]
ys = [igos.InGameObjArr[1].Attributes[1]["NormVec"][2], 0.0]
zs = [igos.InGameObjArr[1].Attributes[1]["NormVec"][3], 0.0]

lines!(ax3d, xs, ys, zs; color=:orange)


display(fig)

for i in 1:100
    @time a = MakeScanningGrid(i * 100, π / 2, π / 2)
end
# 0.001725 seconds (θres = 10000, alloc ∼2.5kB)

for i in 1:100
    @time a = MakeScanningGrid_static(i * 100, π / 2, π / 2)
end
# 0.000457 seconds (θres = 10000, alloc ∼1kB)

println(igos.InGameObjArr[1].Attributes[1]["NormVec"])
dot(igos.InGameObjArr[1].Attributes[1]["NormVec"], [1, 0, 3])
dot(igos.InGameObjArr[1].Attributes[1]["NormVec"], [1, 1, -1])


