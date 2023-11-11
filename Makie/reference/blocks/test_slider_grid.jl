using GLMakie

fig = Figure()

ax = Axis(fig[1, 1])

sg = SliderGrid(
    fig[1, 2],
    (label = "Voltage", range=0:0.1:10, format = "{:.1f}V", startvalue = 5.3),
    (label = "Current", range = 0:0.1:20, format = "{:.1f}A", startvalue = 10.2),
    (label = "Resistance", range = 0:0.1:30, format = "{:.1f}Ω", startvalue = 15.9),
    width = 350,
    tellheight = false
)

slider_o = [s.value for s in sg.sliders]
bars = lift(slider_o...) do slvalues...
    [slvalues...]
end

barplot!(ax, bars, color = [:yellow, :orange, :red])
ylims!(ax, 0, 30)

display(fig)