using JSON,
    ColorSchemes, CairoMakie, DataFrames, DataFramesMeta, Chain, CSV, COBREXA, Statistics

rdir = "linesearch"
losses_dir = filter(endswith("losses.csv"), readdir(joinpath("results", rdir)))
params_dir = filter(endswith("params.csv"), readdir(joinpath("results", rdir)))

dfl = DataFrame(Condition = String[], Loss = Float64[], Iteration = Int64[])
for dir in losses_dir
    append!(dfl, DataFrame(CSV.File(joinpath("results", rdir, dir))))
end

master_id = "WT1#B2"
pdirfile = master_id * "#params.csv"
dfp = DataFrame(CSV.File(joinpath("results", rdir, pdirfile)))
sort!(dfp, :Iteration)

#: load model
model =
    load_model(StandardModel, joinpath("model_construction", "model_files", "iML1515.json"))
subsys_df = DataFrame(
    KcatID = "k#" .* reactions(model),
    Subsystem = [model.reactions[rid].subsystem for rid in reactions(model)],
)


#: Plot
fig = Figure(resolution = (1200, 600), backgroundcolor = :transparent);

ga = fig[1, 1] = GridLayout()
gb12 = fig[1, 2] = GridLayout()
gb1 = gb12[1, 1] = GridLayout()
gb2 = gb12[2, 1] = GridLayout()

#: Losses
conditions = unique(dfl[!, :Condition])
losses_ax = Axis(
    ga[1, 1],
    # yscale=log10,
    xscale = log10,
    xlabel = "Iteration",
    ylabel = "Mean squared relative error (L)",
);
for condition in conditions
    dflt = @subset(dfl, :Condition .== condition)
    sort!(dflt, :Iteration)
    if condition != master_id
        color = ColorSchemes.Greys_4[3]
        lwidth = 2.0
        lines!(
            losses_ax,
            dflt[!, :Iteration], 
            dflt[!, :Loss],
            color = color,
            linewidth = lwidth,
        )
    end
end
for condition in conditions
    dflt = @subset(dfl, :Condition .== condition)
    sort!(dflt, :Iteration)
    if condition == master_id
        color = ColorSchemes.Dark2_3[1]
        lwidth = 3.0
        lines!(
            losses_ax,
            dflt[!, :Iteration],
            dflt[!, :Loss],
            color = color,
            linewidth = lwidth,
        )
    end
end

hidexdecorations!(losses_ax, label = false, ticklabels = false, ticks = false)
hideydecorations!(losses_ax, label = false, ticklabels = false, ticks = false)

elms = [
    LineElement(color = ColorSchemes.Dark2_3[1], linewidth = 3),
    LineElement(color = ColorSchemes.ColorSchemes.Greys_4[3], linewidth = 2),
]

Legend(
    ga[1, 1],
    elms,
    [master_id, "Other conditions"],
    halign = :right,
    valign = :top,
    tellheight = false,
    tellwidth = false,
    margin = (10, 10, 10, 10),
)

fig

#: kcats
subsystem = "Glycolysis/Gluconeogenesis"

subsystem_kcatids = @subset(subsys_df, :Subsystem .== subsystem)[!, :KcatID]
df = @subset dfp @byrow begin
    in(:KcatID, subsystem_kcatids)
end

ax = Axis(
    title = "Glycolysis for $master_id",
    gb1[1, 1],
    xscale = log10,
    yscale = log10,
    xlabel = "Iteration",
    ylabel = "Turnover number [1/s]",
);

# kcatids = unique(df[!, :KcatID])
kcatids = [ # downselected list
    "k#GAPD",
    "k#PFK",
    "k#PGM",
    "k#PPS",
    "k#TPI",
    "k#PGK",
    "k#FBA",
    "k#PGI",
    "k#ENO",
]
for (i, kcatid) in enumerate(kcatids)
    tdf = @subset(dfp, :KcatID .== kcatid)
    kcats = 1 / (3600 * 1e-6) .* tdf[!, :Kcat]
    lines!(
        ax,
        tdf[!, :Iteration],
        kcats,
        color = ColorSchemes.Paired_9[i],
        label = last(split(kcatid, "#")),
    )
end
gb1[1, 2] = Legend(
    fig,
    ax,
    "Enzyme",
    framevisible = false,
    #    orientation = :horizontal,
    #    tellheight = false,
    #    tellwidth = false,
)
hidexdecorations!(ax, label = false, ticklabels = false, ticks = false)
hideydecorations!(ax, label = false, ticklabels = false, ticks = false)
fig

#: derivatives
ax2 = Axis(
    title = "Glycolysis for $master_id",
    gb2[1, 1],
    xscale = log10,
    xlabel = "Iteration",
    ylabel = "Turnover number\nscaled derivative",
);

iters = unique(df[!, :Iteration])
hmdata = zeros(length(iters), length(kcatids))
psdlog(x) = sign(x) * log10(abs(x) + 1)
for (i, kcatid) in enumerate(kcatids)
    tdf = @subset(dfp, :KcatID .== kcatid)
    hmdata[:, i] .= psdlog.(tdf[!, :Derivative])
end

hm = heatmap!(
    ax2,
    iters,
    1:size(hmdata, 2),
    hmdata,
    colormap = Reverse(ColorSchemes.RdYlBu_9[2:8]),
    highclip = ColorSchemes.RdYlBu_9[1],
    lowclip = ColorSchemes.RdYlBu_9[9],
    colorrange = (-1.0, 1.0),
)

ax2.yticks = (1:size(hmdata, 2), last.(split.(kcatids, "#")))
Colorbar(gb2[1, 2], hm, label = "[AU]")

for (label, layout) in zip(["A", "B", "C"], [ga, gb1, gb2])
    Label(
        layout[1, 1, TopLeft()],
        label,
        textsize = 26,
        # font = noto_sans_bold,
        padding = (0, 5, 5, 0),
        halign = :right,
    )
end

fig

CairoMakie.FileIO.save(joinpath("..", "DifferentiableMetabolismPaper", "docs", "imgs", "gd_iteration_charact.pdf"), fig)
