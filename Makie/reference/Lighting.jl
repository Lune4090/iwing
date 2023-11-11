using GLMakie
GLMakie.activate!(ssao = true)
GLMakie.closeall() # close screen to refrect ssao setting

fig = Figure()

# https://ja.wikipedia.org/wiki/SSAO
ssao = Makie.SSAO(radius = 5.0, blur = 3)
ax = LScene(fig[1, 1], scenekw = (ssao = ssao,))
# SSAO attributes
ax.scene.ssao.bias[] = 0.025

box = Rect3(Point3f(-0.5), Vec3f(1))
positions = [Point3f(x, y, 2*rand()-1) for x in -5:5 for y in -5:5]
meshscatter!(ax, positions, marker=box, markersize=1, color=:lightblue, ssao=true)
fig