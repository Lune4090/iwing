#= --- Comments--- =#
# In general, original structs or functions are minimal for easy maintanance
# Original structs are designed to be compatible with Makie & GeometryBasics

#= --- Structs --- =#

#= --- Functions --- =#

function ViDARsLoop(
    AllObjDict::Dict,
    center::Point3f, center_direction::QuatRotation,
    θRlim, θLlim, θres)

    #println("----------------------------------------------")
    #println("Start ViDARsLoop")
    # 描画情報を格納するScanningGridをθres個のDictを持ったVectorとして生成
    ScanningGrid = Vector{Vector}(undef, θres)
    dθ = (θLlim - θRlim) / θres
    # X軸正方向を基準に基準角を決定
    θstart = rotation_angle(center_direction)
    # 反時計回りに走査を行う
    for step_θ in 1:θres
        ScanningGrid[step_θ] = Vector{NamedTuple}(undef, length(AllObjDict))
        # https://risalc.info/src/line-plane-intersection-point.html
        ScanDir = Point3f(cos(step_θ * dθ + θstart + θRlim), sin(step_θ * dθ + θstart + θRlim), 0)
        for key in keys(AllObjDict)
            obj = AllObjDict[key]
            tmp = (dist=Inf, normv=nothing, num=nothing, objkey=nothing, attr=nothing)
            # ここからポリゴン毎の処理
            for face_num in eachindex(faces(obj.collisionmesh))
                polygon = obj.collisionmesh_world[face_num]
                V_v1_v2 = polygon[2] - polygon[1]
                # PlainEq: n ⋅ x = h (n: Normal vector)
                θpoly = decide_arg2D(V_v1_v2)
                n = Point3f(cos(θpoly + π / 2), sin(θpoly + π / 2), 0) # sign(π/2) doesn't matter cuz h also depends on it
                h = transpose(n) * polygon[1]
                t = Xpt_dist_calc(center, polygon, ScanDir, n, h)
                tmp.dist > t ? tmp = (dist=t, normv=n, num=face_num, objkey=key, attr=obj.attributes) : nothing
            end
            @show tmp
            ScanningGrid[step_θ][key] = tmp
        end
    end
    return ScanningGrid
end

function Xpt_dist_calc(center, polygon, ScanDir, n, h)
    dp = transpose(n) * ScanDir
    if dp != 0 # check scandir is parallel to plane or not
        t = (h - transpose(n) * center) / dp # Xpt = center + t*ScanDir
        if t >= 0 # check plane is directed on scandir
            Xpt = center + ScanDir * t
            V_v1_Xpt = (Xpt - polygon[1])
            V_v2_Xpt = (Xpt - polygon[2])
            if transpose(V_v1_Xpt) * V_v2_Xpt <= 0 # check Xpt is inside polygon
                return Float32(t)
            end
        end
    end
    return Float32(Inf)
end