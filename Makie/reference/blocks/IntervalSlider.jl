using GLMakie

fig = Figure()

ax = Axis(fig[1, 1], limits = (0, 1, 0, 1))

hslider = IntervalSlider(fig[2, 1], range = LinRange(0, 1, 1000),
startvalues = (0.2, 0.8))

vslider = IntervalSlider(fig[1, 2], range = LinRange(0, 1, 1000),
startvalues = (0.4, 0.9), horizontal = false)

labeltext1 = lift(hslider.interval) do int
    string(round.(int, digits = 2))
end
Label(fig[3, 1], labeltext1, tellwidth = false)

labeltext2 = lift(vslider.interval) do int
    string(round.(int, digits = 2))
end
Label(fig[1, 3], labeltext2, tellheight= false, rotation = π/2)

points = rand(Point2f, 300)

# color points differently if they are within the 2 intervals
colors = lift(hslider.interval, vslider.interval) do h_int, v_int
    map(points) do p
        (h_int[1] < p[1] < h_int[2])&&(v_int[1] < p[2] < v_int[2])
    end
end

scatter!(points, color = colors, colormap = [:black, :orange], strokewidth = 0)

display(fig)