#= --- Comments--- =#
# In general, original structs or functions are minimal for easy maintanance
# Original structs are designed to be compatible with Makie & GeometryBasics

# This file mainly aims to easy control for scripts with many expression variation
#=

このファイルは，「テキストとその見え方を渡したら，勝手に適切なオブジェクトを選んで返してくれる」
ことの実現が目的である．
例えば，"player name"というStringをどう表示するかは，実はかなり種類が考えられる．
1. Text!()を使い，テキストとして表示する
2. LaTeX表記でテキストとして表示する
3. フォントなりなんなりのpng画像を持ってきて，これを表示する
4. Mesh(MyGameObjectではない)として表示する
5. MyGameObjectとして表示する
...etc
これらを同一の関数で，引数を変えるだけで異なる型のオブジェクトを生成することがこのファイルの目的
生成された型は，Drawerの各関数に引き渡れた後多重ディスパッチで適当な関数を叩かれる

=#

using GLMakie
using FileIO

#= --- Structs --- =#

abstract type AbstractMetaText end

struct MetaTextasText <: AbstractMetaText
    Text::String
    Color::Union{Symbol,Float32}
end

#= --- Functions --- =#

function MetaTextgenerator(Text::String, generatedstruct::Symbol, kwargs...)
    return eval(generatedstruct, Text, kwargs)
end


brain = load(assetpath("brain.stl"))

gc = Figure()
ax3d = Axis3(gc[1, 1], title="Brain activation")
m = mesh!(
    ax3d,
    brain,
    color=[tri[1][2] for tri in brain for i in 1:3],
    colormap=Reverse(:magma),
)
Colorbar(gc[1, 2], m, label="BOLD level")

@show color = [tri[1][2] for tri in brain for i in 1:3]
