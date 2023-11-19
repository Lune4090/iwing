
# This module does not depend on other modules (except of Base)

#= --- Define Original structs --- =#

#= --- Define Original functions --- =#

function DecideArg2D(Vec::Vector{Float64})
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