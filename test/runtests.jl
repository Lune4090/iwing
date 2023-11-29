#= --- Comments--- =#
# Unit Tests have to be done for all of the function

#= --- Structs --- =#

#= --- Functions --- =#
using Test
include("../src/iwing.jl")

function TestGame()
    # 全てのオブジェクトは作られたのちにこのDictに追記されることで初めて名前と存在をゲームから認められる

    p1 = Point3f(1.0, 0.0, 0.0)
    p2 = Point3f(2.0, 2.0, 0.0)
    tmpvertices2D = [p1, p2]
    tmpfaces2D = [TriangleFace([1, 2, 1]), TriangleFace([2, 1, 2])]

    # 辞書に登録
    AllObjDict = Dict{Int,MyGameObject}()
    tmpdict = Dict("Normal" => 1.0f0, "IR" => 1.0f0)
    for i in 1:10
        AllObjDict[i] = MyGameObject(
            Point3f(2.0, 0.0, 0.0),
            QuatRotation(RotZ(0)),
            Mesh(tmpvertices2D, tmpfaces2D),
            Mesh(tmpvertices2D, tmpfaces2D),
            Dict("size" => 1.0f0, "reflectionspectrum" => tmpdict, "transparency" => 1.0f0),
            Mesh(tmpvertices2D, tmpfaces2D),
            Mesh(tmpvertices2D, tmpfaces2D))
    end

    # PlayerCharacter position (最終的にはObjとして引き渡す)
    player_pos = Point3f(0.0)
    player_direction = QuatRotation(RotZ(0))
    search_direction = QuatRotation(RotZ(0))
    chase_direction = QuatRotation(RotZ(0))

    savedata = (
        AllObjDict=AllObjDict,
        player_pos=player_pos,
        player_direction=player_direction,
        search_direction=search_direction,
        chase_direction=chase_direction
    )
    GameMain(savedata)
end

TestGame()

