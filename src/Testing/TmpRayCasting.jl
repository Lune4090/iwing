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

MyQuaternion(x::Real,y::Real,z::Real,θ::Real)=MyQuaternion(convert(Float64, x),convert(Float64, y),convert(Float64, z),convert(Float64, θ))

MyQuaternion(vec::Vector) = MyQuaternion(vec[1], vec[2], vec[3], 0.)

mutable struct MyNormalizedQuaternion <: MyAbstractNormalizedQuaternion
  x::Float64
  y::Float64
  z::Float64
  θ::Float64
  MyNormalizedQuaternion(x,y,z,θ) = x^2+y^2+z^2 != 1 ? error("RotationQuaternion's Norm != 1!") : new(x,y,z,θ)
end

# Int->Float64 Conversion
MyNormalizedQuaternion(x::Int,y::Int,z::Int,θ::Int)=MyNormalizedQuaternion(convert(Float64, x),convert(Float64, y),convert(Float64, z),convert(Float64, θ))



#= --- Original Structures' functions --- =#
function RotateQuaternion(Rotator::MyAbstractQuaternion, Rotated::MyAbstractQuaternion)
  return MyQuaternion(
     Rotator.θ*Rotated.x - Rotator.z*Rotated.y + Rotator.y*Rotated.z + Rotator.x*Rotated.θ,
     Rotator.z*Rotated.x + Rotator.θ*Rotated.y - Rotator.x*Rotated.z + Rotator.y*Rotated.θ,
    -Rotator.y*Rotated.x + Rotator.x*Rotated.y + Rotator.θ*Rotated.z + Rotator.z*Rotated.θ,
    -Rotator.x*Rotated.x - Rotator.y*Rotated.y - Rotator.z*Rotated.z + Rotator.θ*Rotated.θ
    )
end

function RotateVectorbyQuaternion(q::MyAbstractQuaternion, vec::Vector)
  return MyQuaternion(
     q.θ*vec.x - q.z*vec.y + q.y*vec.z + q.x*vec.θ,
     q.z*vec.x + q.θ*vec.y - q.x*vec.z + q.y*vec.θ,
    -q.y*vec.x + q.x*vec.y + q.θ*vec.z + q.z*vec.θ,
    -q.x*vec.x - q.y*vec.y - q.z*vec.z + q.θ*vec.θ
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
@kwdef mutable struct  StructuredInGameObj
  InGameObjArr::Vector{InGameObj}
  InGameObjsRelation
end




function CalcNormVec(polygon::Matrix{Float64})
  vec_a = polygon[2,:] - polygon[1,:]
  vec_b = polygon[3,:] - polygon[1,:]
  println(vec_a)
  println(typeof(vec_a))
  VecNorm = cross(vec_a, vec_b)
  return VecNorm/norm(VecNorm)
end

# secne.camera.eyepositionがVEc{3, Float}なので
function CalcDistPoint2Polygon(eyepos::Vec{3, Float32}, polygon::Matrix{Float64})
  return norm(transpose(eyepos) - reduce(+, polygon; dims = 1))
end

# VerticesからFacesに面毎に格納された頂点番号を引数にpolygonを生成
# Facesがn*3のArrayであることを仮定してるから注意
function DecomposeMesh(mesh::CollisionMesh)
  num_faces = size(mesh.Faces, 1)
  polygons = Array{Float64}(undef, 3, 3, num_faces)
    for EachFace in 1:num_faces
      for EachVertex in 1:3
        polygons[EachVertex,:, EachFace] .= mesh.Vertices[mesh.Faces[EachFace, EachVertex], :]
        println(mesh.Faces[EachFace, EachVertex])
      end
    end
  return polygons, num_faces
end

function RayCastingIGOs(cam::Camera, RayCastedObject::StructuredInGameObj; CalcReflection = true)
  # 一旦愚直にforで実装
  for IGO in RayCastedObject.InGameObjArr
    polygons, num_faces = DecomposeMesh(IGO.CollisionMesh)
    NewAttributes = []
    # ブロードキャストできるはず
    for idx in num_faces
      ReturnDict = Dict()
      cam_pos = to_value(cam.eyeposition)
      println("camera position:$cam_pos")
    
      ReturnDict["Dis"] = CalcDistPoint2Polygon(cam_pos, polygons[:,:,idx])
      ReturnDict["NormVec"] = CalcNormVec(polygons[:,:,idx])
      #=
      if CalcReflection
        ReturnDict["Reflection"] =  CalcReflectionPoint2Plain(center::Point3f, NormVec)
      end
      =#
      push!(NewAttributes, ReturnDict)
    end
    IGO.Attributes = NewAttributes
  end
end



#= --- TestCodes --- =#
tmpvertices = [
    0.0 0.0 1.0;
    1.0 0.0 3.0;
    1.0 1.0 -1.0;
    0.0 1.0 2.0;
]
tmpfaces = [
    1 2 3;
    3 4 1;
]
tmppos = [0.0,0.0,0.0]
tmpquarternion = MyNormalizedQuaternion(1,0,0,0)
Mymesh = CollisionMesh(tmpvertices,tmpfaces)
tmpdict = Dict([("A", 1)])

Mygameobj = InGameObj(
  Pos= tmppos,
  RotationQuaternion= tmpquarternion,
  CollisionMesh= Mymesh,
  AppearanceMesh= Nothing,
  Attributes= [tmpdict])

mesh(Mygameobj.CollisionMesh.Vertices, Mygameobj.CollisionMesh.Faces; color = :lightblue, shading = true)

igos =StructuredInGameObj([Mygameobj], Nothing)

polygons , num_faces= DecomposeMesh(Mygameobj.CollisionMesh)

NewAttributes = []

fig = Figure()
ax3d = Axis3(fig[1,1])
cam = ax3d.scene.camera

RayCastingIGOs(cam, igos)
println(igos.InGameObjArr[1].Attributes)
mesh!(ax3d,igos.InGameObjArr[1].CollisionMesh.Vertices, igos.InGameObjArr[1].CollisionMesh.Faces ;color = :lightblue)
