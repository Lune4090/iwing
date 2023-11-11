using GLMakie
GLMakie.activate!(ssao = true)
ssao = Makie.SSAO(radius = 5, blur = 3)

fig = Figure()
pl = PointLight(Point3f(0), RGBf(20, 20, 20))
al = AmbientLight(RGBf(0.2, 0.2, 0.2))
lscene = LScene(fig[1, 1], show_axis = false, scenekw = (lights = [pl, al], backgroundcolor =:black, clear = true))
meshscatter!(lscene, randn(300, 3), color=:gray)

lscene.scene.ssao.bias[] = 0.025

display(fig)