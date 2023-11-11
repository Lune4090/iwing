# https://docs.makie.org/stable/tutorials/layout-tutorial/index.html#layout_tutorial

using GLMakie
using FileIO

# GLMakie initial setting
# https://docs.makie.org/stable/explanations/backends/glmakie/

GLMakie.activate!(; 
    focus_on_show = true,
    fullscreen = false,
    framerate = 60.0,
    decorated = true,
    title = "visualization"
    )

#= --------Make base Figure--------=#
f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),
resolution = (1000, 700))

#= --------Make nested grids--------=#
grid_a  = f[1,1] = GridLayout()
grid_b  = f[2,1] = GridLayout()
grid_cd = f[1:2,2] = GridLayout()
# A grid can be divided into 2 or more grids
grid_c  = grid_cd[1,1] = GridLayout()
grid_d  = grid_cd[2,1] = GridLayout()


#= -------Panel A-------- =#
axtop = Axis(grid_a[1, 1])
axmain = Axis(grid_a[2, 1], xlabel = "before", ylabel = "after")
axright = Axis(grid_a[2, 2])

linkyaxes!(axmain, axright)
linkxaxes!(axmain, axtop)

labels = ["treatement", "placebo", "control"]
data = randn(3, 100, 2) .+ [1, 3, 5]

#print(length(eachslice(data, dims=3)))
for (label, col) in zip(labels, eachslice(data, dims = 1))
    scatter!(axmain, col, label = label)
    density!(axtop, col[:, 1])
    density!(axright, col[:, 2], direction = :y)
end

# remove gap between plots and Axis
ylims!(axtop, low = 0)
xlims!(axright, low = 0)
# x/y ticks also can be changed
#axtop.xticks = 0:3:9

# legend
leg = Legend(grid_a[1, 2], axmain)
# hide axtop/axright decorations, modify leg size
hidedecorations!(axtop, grid = true)
hidedecorations!(axright, grid = false)
# leg can adjust its height to side axises
leg.tellheight = true
# remove gap between axtop/right and axmain
#colgap!(grid_a, 10)
#rowgap!(grid_a, 10)
# now add a label across top 2 elements
Label(
    grid_a[1, 1:2, Top()], "Stimulus ratings", valign = :bottom,
    font = :bold,
    padding = (0, 0, 5, 0)
    )
f


#= --------Panel B-------- =#
# create Axes by plotting into the right GridLayout slots
xs =LinRange(0.5, 6, 50)
ys =LinRange(0.5, 6, 50)
data1 = [sin(x^1.5)*cos(y^0.5) for x in xs, y in ys] .+ 0.1 .* randn.()
data2 = [sin(x^0.8)*cos(y^1.5) for x in xs, y in ys] .+ 0.1 .* randn.()

ax1, hm = contourf(
    grid_b[1, 1], xs, ys, data1,
    levels = 6
    )
ax1.title = "histological analysis"
contour!(ax1, xs, ys, data1, levels = 5, color = :black)
hidexdecorations!(ax1)

ax2, hm2 = contourf(
    grid_b[2, 1], xs, ys, data2,
    levels = 6
    )
contour!(ax2, xs, ys, data2, levels = 5, color = :black)

# add a colorbar
cb = Colorbar(grid_b[1:2, 2], hm, label = "cell group")
low, high = extrema(data1)
edges = range(low, high, length = 7)
centers = (edges[1:6] .+ edges[2:7]) .* 0.5
cb.ticks = (centers, string.(1:6))
# align colorbar to a plot in the same grid
# right=0 means setting cb's right side padding to 0
cb.alignmode = Mixed(right=0)
colgap!(grid_b, 10)
rowgap!(grid_b, 10)
f


#= --------Panel C-------- =#
brain = load(assetpath("brain.stl"))

ax3d = Axis3(grid_c[1, 1], title = "Brain activation")
m = mesh!(
    ax3d,
    brain,
    color = [tri[1][2] for tri in brain for i in 1:3],
    colormap = Reverse(:magma),
)
Colorbar(grid_c[1, 2], m, label = "BOLD level")

f


#= --------Panel D-------- =#
# 内包表記
axs = [Axis(grid_d[row, col]) for row in 1:3, col in 1:2]
hidedecorations!.(axs, grid = true, label = false)

for row in 1:3, col in 1:2
    xrange = col == 1 ? (0:0.1:6π) : (0:0.1:10π)

    eeg = [sum(sin(π*randn() + k*x) / k for k in 1:10)
    for x in xrange] .+ 0.1 .* randn.()

    lines!(axs[row, col], eeg, color = (:black, 0.5))
end

axs[3, 1].xlabel = "Day1"
axs[3, 2].xlabel = "Day2"

Label(
    grid_d[1, :, Top()], "EEG traces", 
    valign = :bottom,
    font = :bold,
    padding = (0, 0, 5, 0)
    )
rowgap!(grid_d, 10)
colgap!(grid_d, 10)

# set EEG labels
for (i, label) in enumerate(["sleep", "awake", "test"])
    Box(grid_d[i, 3], color = :gray90)
    Label(grid_d[i, 3], label, rotation = π/2, tellheight = false)
end

# second arg mean 2nd col gap (day2-labels)
colgap!(grid_d, 2, 0)

f


#= --------Subplot labels-------- =#
for (label, layout) in zip(["A", "B", "C", "D"], [grid_a,grid_b,grid_c,grid_d])
    Label(layout[1,1,TopLeft()], label,
        fontsize = 26,
        font = :bold,
        padding = (0, 5, 5, 0), # l,r,b,t
        halign = :right)
end

# final tweaking
# Auto(n<1) : grid_level auto align
colsize!(f.layout, 1, Auto(0.5))
rowsize!(f.layout, 1, Auto(1.5))

f