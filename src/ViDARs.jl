#= --- GameLoop & Inner caller functions--- =#

function ViDARsLoop2D(
    AllObjDict::Dict,
    center::Vector{Float64}, center_PDVec::Vector{Float64};
    θRlim=π / 4,
    θLlim=π / 4,
    θres=1000,
    CalcReflection=false)

    #println("----------------------------------------------")
    #println("Start ViDARsLoop")
    # 描画情報を格納するScanningGridをθres個のDictを持ったVectorとして生成
    ScanningGrid = Vector{Dict{String,Any}}(undef, θres)
    # 走査偏角の決定
    dθ = (θRlim + θLlim) / θres
    # 2D走査の準備
    center_PDVec[3] = 0 # Z座標を強制的にゼロにする
    # X軸正方向を基準に基準角を決定
    θstart = DecideArg2D(center_PDVec)
    #println("----------------------------------------------")
    #println("Start Scanning")
    # 反時計回りに走査を行う
    for step_θ in 1:θres
        ScanningGrid[step_θ] = Dict("Dist" => Inf)
        # 計算はベクトルを利用
        # https://risalc.info/src/line-plane-intersection-point.html
        # 走査線の方向ベクトルScanDirを導出，走査線方程式はcenter + t*ScanDir (0<t)
        ScanDir = [cos(step_θ * dθ + θstart - θRlim), sin(step_θ * dθ + θstart - θRlim), 0]
        #println("----------------------------------------------")
        #println("Start Scanning Step : $step_θ/$θres")
        for key in keys(AllObjDict)
            GameObj = AllObjDict[key]
            #println("----------------------------------------------")
            #println("Start processing object : $key")
            # ここからポリゴン毎の処理
            for face_num in 1:GameObj.CollisionMesh.face_num
                #println("----------------------------------------------")
                #println("Start Processing polygon No.$face_num")
                ReturnDict = Dict{String,Any}()
                polygon = GameObj.CollisionMesh.polygons[:, :, face_num]
                # ポリゴンの頂点1から2へのベクトルをV_v1_v2，視点からポリゴン頂点1へのベクトルをV_c_v1とする
                V_v1_v2 = polygon[:, 2] - polygon[:, 1]
                V_c_v1 = polygon[:, 1] - center
                # 内積とノルムから角度ϕと距離hを導出，平面方程式n ⋅ x = hを導く
                θpoly = DecideArg2D(V_v1_v2)
                ϕ = DecideArg2D(-V_c_v1) - θpoly
                h = norm(V_c_v1) * sin(ϕ)
                n = [cos(θpoly - π / 2), sin(θpoly - π / 2), 0.0]
                # 交点座標は以下の通り
                tmp_dp = n ⋅ ScanDir
                is_parallel = tmp_dp == 0 ? true : false
                is_inner = false
                is_sameside = false
                if !is_parallel
                    t = (h - n ⋅ center) / tmp_dp
                    is_sameside = t >= 0 ? true : false
                    Xpt = center + ScanDir * t
                    # 交点の線内判定
                    # 両頂点からのベクトルとの内積の符号が異なる時のみ交点が線内となりその他の処理を行う
                    V_v1_Xpt = (Xpt - polygon[:, 1])
                    V_v2_Xpt = (Xpt - polygon[:, 2])
                    is_inner = dot(V_v1_Xpt, V_v2_Xpt) <= 0 ? true : false
                end
                if is_inner && is_sameside
                    ReturnDict["Dist"] = norm(Xpt - center)
                    if ScanningGrid[step_θ]["Dist"] > ReturnDict["Dist"]
                        ReturnDict["NormalVec"] = n
                        ReturnDict["FaceNum"] = face_num
                        ReturnDict["ObjectKey"] = key
                        ScanningGrid[step_θ] = ReturnDict
                    end
                end
            end
        end
    end
    return ScanningGrid
end
