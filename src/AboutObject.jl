#= --- Comments--- =#
# In general, original structs or functions are minimal for easy maintanance
# Original structs are designed to be compatible with Makie & GeometryBasics

#= --- Structs --- =#

mutable struct MyGameObject
    position::Vector{Point3f}
    direction::QuatRotation
    # n, θ is against QuatRotation(1,0,0,0)
    # Cuz this game is look-down shooter, so always n = (0,0,1)
    collisionmesh <: AbstractMesh
    # Below attributes are kept in collisionmesh
    appearancemesh <: AbstractMesh
    attributes::Dict
    # some attributes(eg. reflection) will be calculated for each faces
    # KinecticProperty
    # 	Velocity
    # 	ArgVelocity
    # PassiveMaterialProperty
    # 	Roughness
    # 	ReflectionSpectrum
    # ActiveMaterialProperty
    # 	RadiationSpectrum
end

#= --- Functions --- =#

function generate_globalizedmesh(obj::T) where {T<:MyGameObject}
    # ParallelTranslation
    collisionverts, appearverts = obj.collisionmesh, obj.appearancemesh |> coordinates

    for i in eachindex(collisionverts)
        collisionverts[i] = collisionverts[i] |> (x -> obj.direction * x) |> (x + obj.position)
        appearverts[i] = appearverts[i] |> (x -> obj.direction * x) |> (x + obj.position)
    end
    return Mesh(collisionverts, faces(obj.collisionmesh)), Mesh(appearverts, faces(obj.appearancemesh))
end
