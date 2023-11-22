#= --- Define Original structs --- =#

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

#= --- Define Original functions --- =#

function rotate_quaternion(Rotator::MyAbstractQuaternion, Rotated::MyAbstractQuaternion)
    return MyQuaternion(
        Rotator.θ * Rotated.x - Rotator.z * Rotated.y + Rotator.y * Rotated.z + Rotator.x * Rotated.θ,
        Rotator.z * Rotated.x + Rotator.θ * Rotated.y - Rotator.x * Rotated.z + Rotator.y * Rotated.θ,
        -Rotator.y * Rotated.x + Rotator.x * Rotated.y + Rotator.θ * Rotated.z + Rotator.z * Rotated.θ,
        -Rotator.x * Rotated.x - Rotator.y * Rotated.y - Rotator.z * Rotated.z + Rotator.θ * Rotated.θ
    )
end

# https://qiita.com/kenjihiranabe/items/945232fbde58fab45681
function rotate_vector_by_quaternion(q::MyAbstractQuaternion, vec::Vector)
    tmp = MyQuaternion(
        q.θ * vec[1] - q.z * vec[2] + q.y * vec[3] + q.x * 1,
        q.z * vec[1] + q.θ * vec[2] - q.x * vec[3] + q.y * 1,
        -q.y * vec[1] + q.x * vec[2] + q.θ * vec[3] + q.z * 1,
        -q.x * vec[1] - q.y * vec[2] - q.z * vec[3] + q.θ * 1
    )
    return rotate_quaternion(tmp, MyQuaternion(-q.x, -q.y, -q.z, q.θ))

end

function quaternion2vector(q::MyQuaternion)
    v = Vector{Float64}(undef, 3)
    v[1] = q.x
    v[2] = q.y
    v[3] = q.z
    return v
end

function vector2quaternion(v::Vector{Real})
    return MyQuaternion(v[1], v[2], v[3], 1)
end
