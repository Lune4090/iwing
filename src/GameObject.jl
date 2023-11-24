#= --- Comments--- =#
# In general, original structs or functions are minimal for easy maintanance
# Original structs are designed to be compatible with Makie & GeometryBasics

#= --- Structs --- =#

mutable struct MyGameObject
    position::Point3f
    direction::QuatRotation
    # n, θ is against QuatRotation(1,0,0,0)
    # Cuz this game is look-down shooter, so always n = (0,0,1)
    collisionmesh::T where {T<:AbstractMesh}
    # Below attributes are kept in collisionmesh
    appearancemesh::T where {T<:AbstractMesh}
    attributes::Dict
    # some attributes(eg. reflection) will be calculated for each faces
    # KinecticProperty
    # 	Velocity::Point3f(automatic calc)
    #	Size::Float32(0<1)
    # PassiveMaterialProperty
    # 	ReflectionSpectrum::Vector{Float32}
    #	Tranparency::Float32(0<1)
    collisionmesh_world::T where {T<:AbstractMesh}
    appearancemesh_world::T where {T<:AbstractMesh}
end

#= --- Functions --- =#

function generate_globalizedmesh2D!(obj::T) where {T<:MyGameObject}
    # ParallelTranslation
    collisionvert = obj.collisionmesh |> coordinates
    appearvert = obj.appearancemesh |> coordinates
    collisionvert = collisionvert .|> x -> obj.direction * x |> x -> x .+ obj.position |> Point3f
    appearvert = appearvert .|> x -> obj.direction * x |> x -> x .+ obj.position |> Point3f
    obj.collisionmesh_world = Mesh(collisionvert, faces(obj.collisionmesh))
    obj.appearancemesh_world = Mesh(appearvert, faces(obj.appearancemesh))
end

function object_mover!(obj::T, Movement::Point3f, FPS::Float32) where {T<:MyGameObject}
    obj.position += Movement
    obj.attributes["Velocity"] = Movement .* FPS
end
