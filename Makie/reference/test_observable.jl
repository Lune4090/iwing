using GLMakie
using Random

GLMakie.activate!(
    focus_on_show = true,
    fullscreen = false,
    framerate = 60.0,
    decorated = true,
    title = "Makie_test_observable"
)

x1 = Observable(1:10)
y1 = sin.(getindex(x1))

on(x1) do val
    println("x1 has a new value $val")
end

x1[]= 1:10

val1 = x1[]
val2 = to_value(x1)

val1 == val2 ? println("true") : println("false")


# Chaining observables by "lift"
f(x1) = x1.^2
y1 = lift(f, x1)

z1 = lift(y1) do val
    -val
end

@show x1[]
@show y1[]
@show z1[]

# @lift
# map(f, c...)
x2m = Observable(rand(100))
y2m = Observable(rand(100))
x2 = x2m
y2 = y2m
z2 = lift((x2,y2) -> x2 .+ y2, x2, y2)
z2m = @lift($x2m .+ $y2m)

println(z2 == z2m) # >false
println(to_value(z2) == to_value(z2m))# >true
# Objects are distincted even they have same values


#
xo = 0
yo = 0

ball = Observable([Point2f(xo, yo)])

fig = Figure();display(fig)
ax = Axis(fig[1,1])

scatter!(ax, ball)

xlims!(ax, -5,5)
ylims!(ax, -5,5)

# step-by-step calculation
function animationstep!(ball, xo, yo)
    xo += 0.3*(rand() -0.5)
    yo += 0.3*(rand() -0.5)
    ball[] = [Point2f(xo, yo)]
end

println(to_value(ball)[1][1])


# anim update
display(fig);sleep(3)
for i in 1:1000
    animationstep!(ball, to_value(ball)[1][1], to_value(ball)[1][2])
    xlims!(ax, -5-i*0.01,5+i*0.01)
    ylims!(ax, -5-i*0.01,5+i*0.01)
    sleep(0.001)
    if i%100 == 0
        println(ball)
    end
end
