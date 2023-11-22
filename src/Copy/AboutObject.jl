include("AboutMesh.jl")

#= --- Interface --- =#

@kwdef mutable struct MyGameObject
    Pos::Vector{Float64} # これがグローバル座標
    PostureDirectionVector::Vector{Float64}
    CollisionMesh::CollisionMesh
    AppearanceMesh::Any
    # AppearanceMeshは外部ツールで製作した.stl形式のオブジェクトのloadを想定
    # AppearanceMeshは触らない(Makie.meshscatter!()にmarkerとして渡すだけ)
    Attributes::Vector{Dict}
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
    function MyGameObject(
        Pos::Vector{Float64}, PostureDirectionVector::Vector{Float64},
        CollisionMesh::CollisionMesh, AppearanceMesh::Any,
        Attributes::Vector{Dict{Any,Any}})
        norm_PDVec = PostureDirectionVector[1]^2 + PostureDirectionVector[2]^2 + PostureDirectionVector[3]^2
        if norm_PDVec != 1
            error("PostureDirectionVector should be equal to 1 !!!")
        else
            new(Pos, PostureDirectionVector, CollisionMesh, AppearanceMesh, Attributes)
        end
    end
end

# 基本的に2D、3Dを問わずポリゴンは3次元座標系で保持される為流用可能
function TranslatePolygonLocal2Global!(mesh, GameObj::MyGameObject)
    decompose_mesh!(mesh)
    # Local2Global Coordinate Translation
    # Rotation (by Quarternion[λx*sin(θ/2),λy*sin(θ/2),λz*sin(θ/2),cos(θ/2)])
    RotateAxis = cross(mesh.PostureDirectionVector, GameObj.PostureDirectionVector)
    # 通常は回転軸を中心に回転
    if RotateAxis != [0.0, 0.0, 0.0]
        RotateAxis = RotateAxis / norm(RotateAxis)
        # CounterClockwiseを仮定
        RotateAngle = acos(mesh.PostureDirectionVector ⋅ GameObj.PostureDirectionVector)
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
    elseif mesh.PostureDirectionVector != GameObj.PostureDirectionVector
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
            mesh.polygons[:, vert, face] .= mesh.polygons[:, vert, face] .+ GameObj.Pos
        end
    end
end
