# 基本は「Pointを作ってVectorに入れる」ってことだと思う

p1 = Point(3, 1)
p2 = Point(1, 4)
p3 = Point(2, 5)

rect = Rect([p1, p2]) #頭良くて草
# Rect(0, 0, 1, 2) こういうのもできるけどあまり使わないか？

# Pointに限らず，GeometryPrimitiveにはメタデータつけれる
poi = meta(p1, city="Abuja", rainfall=1221.2)

meta(poi)
typeof(meta(poi))
poi.city
metafree(poi)
poi

# MultiPoint型というPrimitiveもある
multipoi = MultiPointMeta([p1], city="Abuja", rainfall=1221.2)

l1 = Line(p1, p2)
l2 = Line(p2, p3)

LineString([p1, p2, p3])

poly2 = Polygon([p1, p2, p3]) # なんかめっちゃ大量のtupleやらが出てきて怖い

p31 = Point3f(3.0, 1.0, 9.0)
p32 = Point3f(1.0, 4.0, 3.0)
p33 = Point3f(2.0, 5.0, 2.0)

poly3 = Polygon([p31, p32, p33])

#= Meshableは
const Meshable{Dim,T} = Union{
	Tesselation{Dim,T},
	Mesh{Dim,T},
	AbstractPolygon{Dim,T},
	GeometryPrimitive{Dim,T},
	AbstractVector{<:AbstractPoint{Dim,T}}
}
=#
# らしい

# FaceView()とかいう関数は"全く"使えない(ドキュメントに書いてある通りだと動かないし，そもそも実装が変)，issueも有ったし完全無視

circle = Circle(p1, 1.0)
mulpo = MultiPoint([p1, p2])
datatype = typeof(circle)
supertype(datatype)  # CircleやRectはHyperSphereやHyperRectのAlias，そしてそれらのSuperTypeはGeometryPrimitive
datatype2 = typeof(mulpo)
supertype(datatype2) # 一方でMultiPointのSuperTypeはVector


faceface = TriangleFace([1, 2, 2]) # TriangleFace<:AbstractNgonFace<:AbstractFace<:StaticVectorだった
faceface2 = TriangleFace([2, 1, 2])

# これで作れる
meshmesh = Mesh([p31, p32], [faceface])
# これで通る
@show result = RotX(π / 2) * coordinates(meshmesh)


#= --- Tesselation --- =#

res = Tesselation(rect, (8, 8))
sphere = Sphere(Point3f(0), 1)
m1 = GeometryBasics.mesh(sphere) # uses a default value for tesselation
m2 = GeometryBasics.mesh(Tesselation(sphere, 64)) # uses 64 for te

coordinates(meshmesh)
faces(meshmesh)

area(coordinates(meshmesh), faces(meshmesh)) # Meshのままぶちこめないのがちょいめんどい
area(coordinates(m2), faces(m2)) # Meshのままぶちこめないのがちょいめんどい


# 面積比較してみる
r = 1
sphere = Sphere(Point3f(0), r)
m1 = GeometryBasics.mesh(sphere) # uses a default value for tesselation
m2 = GeometryBasics.mesh(Tesselation(sphere, 64)) # uses 64 for te
# Tesselationはprimitiveのみに使用可能
# 但しMeshableなので上記の様にTesselationした後にMeshにすることもできる

m1area = area(coordinates(m1), faces(m1))
m2area = area(coordinates(m2), faces(m2))
analyticarea = r * r * 4π
Δarea_m1 = analyticarea - m1area
Δarea_m2 = analyticarea - m2area
# 1桁精度が上がる

using GLMakie
fig = Figure()
ax3 = Axis3(fig[1, 1])
meshscatter!(ax3, p31, marker=m1)

meta(m1)
#= 
全体として小文字関数の内exportされているものが少ない
恐らく使用されることを想定していない？
例えばmesh()という関数はあるが，MeshやUVMeshの方が使い勝手がいい
でも上記の様にmesh(Primitive)したいときはこれしかない

Polygonはあほみたいに内部が複雑でよーわからん
=#

p1 = Point3f(1.0, 0.0, 0.0)
p2 = Point3f(2.0, 2.0, 0.0)
tmpvertices2D = [p1, p2]
tmpfaces2D = [TriangleFace([1, 2, 2])]
meshmesh = Mesh(tmpvertices2D, tmpfaces2D)
meshmesh[1][1]


p31 = Point3f(3.0, 1.0, 9.0)
p32 = Point3f(1.0, 4.0, 3.0)

faceface = TriangleFace([1, 2, 1])
faceface2 = TriangleFace([2, 1, 2]) # TriangleFace<:AbstractNgonFace<:AbstractFace<:StaticVectorだった

meshmesh = Mesh([p31, p32], [faceface, faceface2])
meshmesh[1]
meshmesh[2]
using Rotations

QuatRotation(RotZ(π / 120) * QuatRotation(1, 0, 0, 0))

answer = RotZ(π / 120) * QuatRotation(1, 0, 0, 0) |> QuatRotation

a = QuatRotation(1, 0, 0, 0)

a = a |> x -> RotZ(π / 120) * x |> QuatRotation

# あるPoint3から別のPoint3fに

include("UtlityFuncs.jl")

pt = Point3f(1, 1, 0)
pt2 = Point3f(1, 0, 1)

pt * pt2

norm(pt)

rotation_angle(QuatRotation(1, 0, 0, 0))

rotation_angle(RotZ(π / 120) * QuatRotation(1, 0, 0, 0))

using GLMakie

fig = Figure()

polax = PolarAxis(fig[1, 1])

xs = 1:10

ys = xs .|> x -> x .* [1, 2]

xs = xs .|> x -> x .* [1, 1]

plt = map(zip, xs, ys)

scatter!(polax, [π / 3, π / 6, π / 6], [5, 1, 2])

using GLMakie, Rotations

QuatRotation(RotZ(0))