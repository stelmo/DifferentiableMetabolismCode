# function align_ylabels!(axs...)
#     max_extent = maximum(axs) do ax
#         tight_yticklabel_spacing!(ax)
#         ax.yticklabelspace[]
#     end
#     for ax in axs
#         ax.yticklabelspace = max_extent
#     end
#     return
# end

# f = Figure()
# ax1 = Axis(f[1, 1], ylabel = "Y Label")
# ax2 = Axis(f[2, 1], limits = (nothing, (1000, 9000)), ylabel = "Y Label")
# ax3 = Axis(f[3, 1], limits = (nothing, (0.0001, 0.0002)), ylabel = "Y Label")
# align_ylabels!(ax1, ax2, ax3)
# f

# using CairoMakie

# fig = Figure(resolution=(400, 300))

# ga = fig[1, 1] = GridLayout()
# gb12 = fig[1, 2] = GridLayout()
# gb1 = gb12[1, 1] = GridLayout()
# gb2 = gb12[2, 1] = GridLayout()

# ax1 = Axis(ga[1, 1])
# scatter!(ax1, rand(10), rand(10))

# ax2 = Axis(
#     gb1[1, 1],

# )

# for i in 1:5
#     lines!(
#         ax2,
#         1:10,
#         rand(10),
#         color = ColorSchemes.Paired_9[i],
#         label = string("label $i"),
#     )
# end
# gb1[1, 2] = Legend(
#     fig,
#     ax2,
#     "Legend",
#     framevisible = false,
# )

# # heatmap 
# ax3 = Axis(
#     gb2[1, 1],
#     # xscale=log10,

# )
# data = rand(649, 9)
# heatmap!(
#     ax3,
#     1:size(data, 1),
#     1:size(data, 2),
#     data,
# )

# fig


using CairoMakie

#: Plot
fig = Figure(
    resolution = (1200, 600), 
);

ga = fig[1, 1] = GridLayout()
gb12 = fig[1, 2] = GridLayout()
gb1 = gb12[1, 1] = GridLayout()
gb2 = gb12[2, 1] = GridLayout()

#: Losses
ax1 = Axis(
    ga[1, 1],
    xscale = log10,
);
for i in 1:5
    lines!(
        ax1,
        1:100_000,
        rand(100_000)
    )
end
hidexdecorations!(ax1, label = false, ticklabels = false, ticks = false)
hideydecorations!(ax1, label = false, ticklabels = false, ticks = false)

elms = [
    LineElement(color = :blue, linewidth = 3),
    LineElement(color = :black, linewidth = 2),
]

Legend(
    ga[1, 1],
    elms,
    ["Legend item", "Other conditions"],
    halign = :right,
    valign = :top,
    tellheight = false,
    tellwidth = false,
    margin = (10, 10, 10, 10),
)

#: kcats
ax2 = Axis(
    title = "Title",
    gb1[1, 1],
    xscale = log10,
    yscale = log10,
);

for i in 1:5
    lines!(
        ax2,
        1:100_000,
        rand(100_000),
        label = string(i),
    )
end
gb1[1, 2] = Legend(
    fig,
    ax2,
    "Enzyme",
    framevisible = false,
)
hidexdecorations!(ax, label = false, ticklabels = false, ticks = false)
hideydecorations!(ax, label = false, ticklabels = false, ticks = false)

#: derivatives
hmdata = rand(692, 9)
iters = float.([1:10; 20:10:100; 200:100:1000; 1250:250:154000])
ax3 = Axis(
    title = "title",
    gb2[1, 1],
    xscale = log10,
);

hm = heatmap!(
    ax2,
    xs,
    1:size(hmdata, 2),
    hmdata,
)

fig

using CairoMakie

fig = Figure(resolution = (300, 200),)
ax = Axis(fig[1,1], xscale=log10,)

xs = unique(float.([1:10; 20:10:100; 200:100:1000; 250:250:154000]))
ys = 1:9
hmdata = rand(length(xs), length(ys))
heatmap!(ax, xs, ys, hmdata)

fig