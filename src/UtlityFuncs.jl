#= --- Comments--- =#
# In general, original structs or functions are minimal for easy maintanance
# Original structs are designed to be compatible with Makie & GeometryBasics

#= --- Structs --- =#

#= --- Functions --- =#

"""
vector version is original
"""
function decide_arg2D(Vec::Vector{Float64})
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

"""
rotation_axis(direction)>0 means rotation is counter-clockwise
"""
function decide_arg2D(direction::QuatRotation)
    return rotation_axis(direction) >= 0 ? rotation_angle(direction) : -rotation_angle(direction)
end