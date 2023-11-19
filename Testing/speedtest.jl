using FileIO
using GLMakie

function tmpmain()
    brain = load(assetpath("brain.stl"))
    fig = Figure()
    ax3d = Axis3(fig[1, 1])
    display(fig)
    xs = cos.(1:0.5:20)
    ys = sin.(1:0.5:20)
    zs = LinRange(0, 3, length(xs))
    while to_value(events(fig).window_open)
        @time meshscatter!(ax3d, xs, ys, zs, markersize=0.001, color=zs, marker=brain)
    end
end

tmpmain()