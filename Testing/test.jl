ab = 1

for i in 1:1
    nanka = Symbol("a" * "b")
    println(@eval Symbol("a" * "b"))
    println(@eval $(Symbol("a" * "b")))
    println(eval(Symbol("a" * "b")))
end

for i in 1:3
    nanka = Symbol("a" * "b")
    println(@eval nanka)
end

dump(Symbol("a" * "b"))
dump($(Symbol("a" * "b")))

using GLMakie

function gridlayout_base()
    fig = Figure()
    ax = Axis(fig[1, 1])
    obs = Observable(1)
    display(fig)
    while to_value(events(fig).window_open)
        @show fig[1, 2]
        ax2 = Axis(fig[1, 2])
        scatter!(ax2, 1, 1)
        obs[] = to_value(obs) + 1
        scatter!(ax, obs, 0)
        sleep(1)
        delete!(ax2)
        @show fig.layout
        sleep(1)
        trim!(fig.layout)
        @show fig.layout
        sleep(1)
    end
end

function boundingbox_base()
    fig = Figure()
    ax = Axis(fig[1, 1])
    obs = Observable(1)
    display(fig)
    while to_value(events(fig).window_open)
        @show fig[1, 2]
        ax2 = Axis(fig, bbox=BBox(100, 300, 100, 500))
        scatter!(ax2, 1, 1)
        obs[] = to_value(obs) + 1
        scatter!(ax, obs, 0)
        sleep(1)
        delete!(ax2)
        @show fig.layout
        sleep(1)
    end
end


gridlayout_base()

function f(Dic)
    Dic["1"] += Dic["1"]
    @show Dic["1"]
end

Dic = Dict("1" => 1)

for i in 1:5
    @show Dic
    sleep(2)
    f(Dic)
end


figfig = Figure()

axs = [Axis(figfig[1, i]) for i in 1:3]

scatters = map(axs) do ax
    [scatter!(ax, 0:0.1:10, x -> sin(x) + i) for i in 1:3]
end

typeof(scatters[1][1])

for i in 1:3
    delete!(axs[2], scatters[2][i])
end
empty!(axs[3])

figfig