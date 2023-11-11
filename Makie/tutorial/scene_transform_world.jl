# https://docs.makie.org/stable/tutorials/scenes/

using GLMakie

#=
like camera transformation, 
scene(world) transformation is also possible in Makie
through "Transformation" struct.
Scenes and plots has this struct and they are regarded
as "Makie.Transformable".
=#

#=
Makie.Transformable has 3 function
tranlate!
    relative : translate!(Accum, scene::Transformable, xyz...)
    relative : translate!(scene, xyz::VecTypes)
    absolute : translate!(scene, xyz::VecTypes)
rotate!
    relative : rotate!(Accum, scene, axis_rot)
    absolute : rotate!(t::Transformable, axis_rot::Quaternion)
    absolute : rotate!(t::Transformable, axis_rot::AbstracutFloat)
    absolute : rotate!(t::Transformable, axis_rot...)
scale! (all absolute)
    scale!(t::Transformable, x, y)
    scale!(t::Transformable, x, y, z)
    scale!(t::Transformable, xyz)
    scale!(t::Transformable, xyz...)
=#

scene = Scene()
cam3d!(scene)
sphere_plot = mesh!(scene, Sphere(Point3f(0), 0.5), color =:red)
scale!(scene, 0.5, 0.5, 0.5)
rotate!(scene, Vec3f(1, 0, 0), 0.5) # 0.5 rad around y-ax
scene

translate!(sphere_plot, Vec3f(0, 0, 1))
scene

parent = Scene()
cam3d!(parent)

# One can set the camera lookat and eyeposition, by getting the camera controls and using `update_cam!`
camc = cameracontrols(parent)
update_cam!(parent, camc, Vec3f(0, 8, 0), Vec3f(4.0, 0, 0))

s1 = Scene(parent, camera=parent.camera)
mesh!(s1, Rect3f(Vec3f(0, -0.1, -0.1), Vec3f(5, 0.2, 0.2)))
s2 = Scene(s1, camera=parent.camera)
mesh!(s2, Rect3f(Vec3f(0, -0.1, -0.1), Vec3f(5, 0.2, 0.2)), color=:red)
translate!(s2, 5, 0, 0)
s3 = Scene(s2, camera=parent.camera)
mesh!(s3, Rect3f(Vec3f(-0.2), Vec3f(0.4)), color=:blue)
translate!(s3, 5, 0, 0)
parent

# Now, rotate the "joints"
rotate!(s2, Vec3f(0, 1, 0), 0.5)
rotate!(s3, Vec3f(1, 0, 0), 0.5)
parent