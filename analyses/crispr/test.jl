using CairoMakie

fig = Figure()
ax1 = Axis(fig[1,1])
ax2 = Axis(fig[1,1], yscale=log2, yaxisposition = :right,) # scale is where the problem occurs
ax3 = Axis(fig[1,1]) # scale is where the problem occurs

scatter!(ax1, rand(10), rand(10))
scatter!(ax2, rand(10), rand(10))

text!(ax1, ["this", "displays"]; position=[Point2(0.5, 0.5), Point2(0.25, 0.25)])
text!(ax3, ["I", "don't"]; position=[Point2(0.75, 0.75), Point2(0.35, 0.35)])

fig
