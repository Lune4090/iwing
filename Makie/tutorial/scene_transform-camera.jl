# https://docs.makie.org/stable/tutorials/scenes/

using GLMakie

#=
Scene do 4 things
1. holds a local theme, gets applied to all plot obj
2. manage camera, projection and transformation matrices
3. defince window size (Child scene area < Parent's one)
4. holds a reference to all window "events"
=#

GLMakie.activate!(
    focus_on_show = true,
    fullscreen = false,
    framerate = 60.0,
    decorated = true,
    title = "Makie_tutorial_scene"
)

scene = Scene(backgroundcolor = :gray)

subwindow = Scene(
    scene, # Scene constructor has parent arg
    px_area = Rect(100, 100, 200, 200), 
    clear = true,
    #backgroundcolor = :white
    )
println(typeof(subwindow)) # >Scene
println(subwindow)
#=
>Scene (200px, 200px):
  0 Plots
  0 Child Scenes
=#

cam3d!(subwindow)
meshscatter!(subwindow, rand(Point3f, 10), color=:orange)
center!(subwindow)
scene

subwindow.clear = false
relative_space = Makie.camrelative(subwindow)
# this draws a line at the scene window
lines!(relative_space, Rect(0, 0, 1, 1))
scene

#global i = 0

#=
while i<5
    sleep(2)
    lines!(relative_space, Rect(0, 0, 1, 1))
    scene
    global i += 1
end
=#

campixel!(scene)
w, h = size(scene)
#println("window size : $w*$h")
image!(scene, [sin(i/w) + cos(j/h) for i in 1:w, j in 1:h])
translate!(scene.plots[1], 0, 0, -10000)

scene

on(scene.events.mouseposition) do mousepos
    if ispressed(subwindow, Mouse.left & Keyboard.left_control)
        subwindow.px_area[] = Rect(Int.(mousepos)..., 200, 200)
    end
end


#= --------Interaction with Axis & Layouts-------- =#
#=
The axis contains a scene,
=#
figure, axis, plot_obj = scatter(1:4)
relative_projection = Makie.camrelative(axis.scene);
println(typeof(relative_projection)) # >Scene
println(relative_projection)
#=
Scene (740px, 541px):
0 Plots
0 Child Scenes
=#
scatter!(relative_projection, [Point2f(0.5)], color=:red)
# some 

# offset & text are in pixelspace
text!(relative_projection, "Hi", position=Point2f(0.5), offset=Vec2f(5))
lines!(relative_projection, Rect(0,0,1,1), color=:blue, linewidth=3)
figure


