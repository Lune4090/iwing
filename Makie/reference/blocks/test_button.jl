using GLMakie
fig = Figure()

ax = Axis(fig[1, 1])
fig[2, 1] = buttongrid = GridLayout(tellwidth = false)

counts = Observable([1, 4, 3, 7, 2])

buttonlabels = [lift(x -> "Count: $(x[i])", counts) for i in 1:5]

buttons = buttongrid[1, 1:5] = [Button(fig, label = l) for l in buttonlabels]

println(typeof(buttons[1].clicks)) # >Observable

for i in 1:5
    on(buttons[i].clicks) do n
        counts[][i] += 1
        notify(counts)
    end
end

barplot!(ax, counts, color = cgrad(:Spectral)[LinRange(0,1,5)])
ylims!(ax, 0, 20)

display(fig)