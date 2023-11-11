# Tutorial movie
# https://www.youtube.com/watch?v=L-gyDvhjzGQ
# Source code
# https://gist.github.com/Datseris/4b9d25a3ddb3936d3b83d3037f8188dd


using DynamicalSystems
using OrdinaryDiffEq
using GLMakie
using DataStructures: CircularBuffer

# Observable is a mutable container{T::anytype}
👻 = Observable(1)

👀 = on(👻) do val
    println("Observable now has $val")
end

#👻[] = 5;

#=
Makie is a Observer and 
if any Observable(object) is changed,
then Makie can reflect this.
=#

# ex.
ox = 1:4
oy = Observable(rand(4))
lw = Observable(2)

fig, ax = lines(ox, oy; linewidth = lw)
ylims!(ax, 0, 1)

lw[] = 50
oy = rand(4)

# ex. double pendulum
# 1. Initialize simulation in a stepping manner

const L1 = 1.0
const L2 = 0.9
M = 2
u0 = [π/3, 0, 3π/4, -2]
dp = Systems.double_pendulum(u0; L1, L2)

diffeq = (alg = Tsit5(), adaptive = false, dt = 0.005)

integ = dp.integ

function xycoords(state)
    θ1 = state[1]
    θ2 = state[3]
    x1 = L1*sin(θ1)
    y1 = -L1*cos(θ1)
    x2 = x1 + L2*sin(θ2)
    y2 = y1 - L2*cos(θ2)
    return x1, x2, y1, y2
end

# step-by-step calculation
function progress_for_one_step!(integ)
    step!(integ)
    u = integ.u
    return xycoords(u)
end


# 2. Initialize "Observable" in the animation
# decide which var is static and dynamic

x1, x2, y1, y2 = xycoords(u0)
typeof(x1)
typeof(to_value(x1))
rod   = Observable([Point2f(0, 0), Point2f(x1, y1), Point2f(x2, y2)])
balls = Observable([Point2f(x1, y1), Point2f(x2, y2)])
# (Remember: the most optimal way to plot 2D things in Makie.jl is to
# give it a vector of `Point2f`, the coordinates for the plot)

# Here we have initialized two _different_ observables, because
# rods and balls will be plotted in a different manner (lines/scatter)

# Next is the observable for the tail
tail = 300 # length of plotted trajectory, in units of `dt`
# The circular buffer datastructure makes making stepping-based
# animations very intuitive
traj = CircularBuffer{Point2f}(tail)
fill!(traj, Point2f(x2, y2)) # add correct values to the circular buffer
traj = Observable(traj) # make it an observable


# %% 3. Plot the `Observable`s and any other static elements
# Before plotting we need to initialie a figure
fig = Figure(); display(fig)
# in my experience it leads to cleaner code if we first initialize 
# an axis and populate it accordingly.
ax = Axis(fig[1,1])

# Now we plot the observables _directly_! First the pendulum
lines!(ax, rod; linewidth = 4, color = :purple)
scatter!(ax, balls; marker = :circle, strokewidth = 2, 
    strokecolor = :purple,
    color = :black, markersize = [8, 12]
)

# then its trajectory, with a nice fadeout color
c = to_color(:purple)
tailcol = [RGBAf(c.r, c.g, c.b, (i/tail)^2) for i in 1:tail]
lines!(ax, traj; linewidth = 3, color = tailcol)

# We can also plot now any other static elements
ax.title = "double pendulum"
ax.aspect = DataAspect()
l = 1.05(L1+L2)
xlims!(ax, -l, l)
ylims!(ax, -l, 0.5l)

# %% 4. Create the "animation stepping function"
# Using the functions of step 1, we now define a function
# that updates the observables. Makie.jl understands observable
# updates and directly reflects this on the plotted elements.
function animstep!(integ, rod, balls, traj)
    x1,x2,y1,y2 = progress_for_one_step!(integ)
    rod[] = [Point2f(0, 0), Point2f(x1, y1), Point2f(x2, y2)]
    balls[] = [Point2f(x1, y1), Point2f(x2, y2)]
    push!(traj[], Point2f(x2, y2))
    traj[] = traj[] # <- important! Updating in-place the value of an
                    # `Observable` does not trigger an update!
end

# %% 5. Test it
for i in 1:1000
    animstep!(integ, rod, balls, traj)
    sleep(0.001)
end
