# https://docs.makie.org/stable/tutorials/basic-tutorial/

using GLMakie

bg = Figure(backgroundcolor = :orange, resolution = (800, 600))

#=
f   = Figure()
ax  = Axis(
    f[1, 1],
    title = "A Makie Axis",
    xlabel = "The x label",
    ylabel = "The y label"
)
x = range(0, 10, length = 100)
y= sin.(x)

# line plot
lines!(ax, x, y)
f

# scatter plot
scatter!(ax, x, y)
f
=#

#=
# Simple graph drawing
x = range(0, 10, length = 100)
y = sin.(x)

fig, axis, lineplot = lines(x, y)
println(typeof(fig))
println(typeof(axis))
println(typeof(lineplot))

lines(x, y)
scatter(x, y)
# note they don't need to call f::Figure

# more simple way
lines(0:0.5:10, cos)
=#

# layering plots
# must use lines!(), cuz only modifying existing fig
x = range(0, 10, length = 100)
f, ax, l1 = lines(x, sin)
# Axis keeps track of what plotted / color
l2 = lines!(ax, x, cos)
f

# Attributes
#=
colors https://juliagraphics.github.io/Colors.jl
    name    : :red
    hex str : "ffccbk"
    func    : RGBf(r,g,b), RGBAf(r,g,b,α)
    hex + α : (:red, 0.5)
=#
x = range(0, 10, length = 100)
f, ax, sc1 = scatter(
    x, sin, 
    color = RGBAf(0.8, 0.1, 0.2, 0.8), 
    markersize = 5
    )

sc2 = scatter!(
    ax, x, cos, 
    color = RGBAf(0.2, 0.1, 0.8, 0.4), 
    markersize = 15
    )

# Manupulate attribute by (plot.atribute = new value) 
sc1.marker = :utriangle
sc1.markersize = 20
#sc2.color = :transparent
sc2.markersize = 5
#sc2.strokewidth = 1
#sc2.strokecolor = :green

f

# Array Attributes
x = range(0, 10, length = 100)

# if color has numerical array, values maps to colormap
# colorrange can clip color range
scatter(
    x, sin,
    markersize = range(5,15,length=100),
    color = range(0,1,length=100),
    colormap = :thermal,
    colorrange = (0.33, 0.66)
    )
# using array of color is also possible
colors = repeat([:crimson, :dodgerblue, :slateblue1, :sienna1, :orchid1], 20)

scatter(x, sin, color = colors, markersize = 20)

# legend
x = range(0, 10, length=100)

lines(x, sin, color = :red, label = "sin")
lines!(x, cos, color = :blue, label = "cos")
axislegend()
current_figure()

# subplots
x = LinRange(0,10,100)
y = sin.(x)

fig = Figure()
lines(fig[1, 1], x, y, color = :red)
lines(fig[1, 2], x, y, color = :blue)
lines(fig[2, 1:2], x, y, color = :green)

fig

function invisibleAxis(pos)
    ax = Axis(
        pos,
        leftspinevisible = false,
        rightspinevisible = false,
        topspinevisible = false,
        bottomspinevisible = false,
        xgridvisible = false,
        ygridvisible = false,
        xticksvisible = false,
        yticksvisible = false,
        xticklabelsvisible = false,
        yticklabelsvisible = false
    )
    return ax
end

# make empty subplots
fig = Figure()
ax1 = Axis(fig[1, 1])
ax2 = Axis(fig[1, 2:3])
ax3 = invisibleAxis(fig[2, 1:3])
ax4 = Axis(fig[3, 1])
ax5 = Axis(fig[3, 2])
ax6 = Axis(fig[3, 3])

lines!(ax1, 0..10, exp2)
lines!(
    ax2, x, sqrt,
    markersize = 15,
    color = range(0,1,length=length(x)),
    colormap = :thermal,
    colorrange = (0.33, 0.66)
    )
lines!(ax3, x, y)
lines!(ax4, 0..10, sin)
hm = heatmap!(ax6, randn(20, 20))
Colorbar(fig[3, 2], hm)
fig


