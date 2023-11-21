#= --- Define Original structs --- =#

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

#= --- Define Original functions --- =#

# VerticesからFacesに面毎に格納された頂点番号を引数にpolygonを生成
# linearindexingによるアンローリングよりforの二重ループの方が速かったのでこのまま
function decompose_mesh!(mesh::MyAbstractMesh)
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