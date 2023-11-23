#= --- Comments--- =#
# In general, original structs or functions are minimal for easy maintanance
# Original structs are designed to be compatible with Makie & GeometryBasics

#= --- Structs --- =#

#= --- Functions --- =#

"""
vector version is original
"""
function norm(vec::T) where {T<:AbstractPoint}
    return vec .|> (x -> x^2) |> sum |> sqrt
end

function decide_arg2D(Vec::Point3f)
    tmpVec = Vec / norm(Vec) # Normを1に戻す
    # グローバルのx軸を基準(θ = 0)とし、反時計回りを正とする
    # -π<θ<=πで走査範囲を定義
    θ = acos(tmpVec[1])
    # y座標が負なら負になる
    if tmpVec[2] < 0
        θ = -θ
    end
    return θ
end