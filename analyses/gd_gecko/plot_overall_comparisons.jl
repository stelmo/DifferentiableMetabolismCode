using ColorSchemes, CairoMakie, DataFrames, DataFramesMeta, Chain, CSV
using Measurements

kmax_brenda_df = DataFrame(CSV.File(joinpath("results", "gd_gecko", "kmax_brenda_df.csv")))
polish_df = DataFrame(CSV.File(joinpath("results", "gd_gecko", "polish_df.csv")))

#: Plot figure
fig = Figure(resolution = (1200, 600), backgroundcolor = :transparent);

ga = fig[1, 1] = GridLayout()
gb = fig[1, 2] = GridLayout()

#: Plot hold outs
xlabs = polish_df[!, :Condition]
rloss = polish_df[!, :Rloss] .± polish_df[!, :RlossSTD]
ploss = polish_df[!, :Ploss] .± polish_df[!, :PlossSTD]

frac_improv = (rloss .- ploss) ./ rloss

hold_ax = Axis(
    ga[1, 1],
    xlabel = "Condition",
    ylabel = "Fraction improved mean squared relative error",
)

barplot!(
    hold_ax,
    1:length(xlabs),
    Measurements.value.(frac_improv),
    color = ColorSchemes.Set2_3[1],
)

errorbars!(
    hold_ax,
    1:length(xlabs),
    Measurements.value.(frac_improv),
    Measurements.uncertainty.(frac_improv),
    Measurements.uncertainty.(frac_improv),
    whiskerwidth = 10,
)

hold_ax.xticks = (1:length(xlabs), xlabs)
hidexdecorations!(hold_ax, label = false, ticklabels = false, ticks = false)
hideydecorations!(hold_ax, label = false, ticklabels = false, ticks = false)

#: Plot BRENDA vs Polish

kb = kmax_brenda_df[!, :Kcat]
kgd = kmax_brenda_df[!, :Kmax]
usubs = unique(kmax_brenda_df[!, :SubID]) # unique subsystems
kss = [
    ColorSchemes.Set2_6[findfirst(x -> x == id, usubs)] for id in kmax_brenda_df[!, :SubID]
]

brenda_ax = Axis(
    gb[1, 1],
    yscale = log10,
    xscale = log10,
    xlabel = "BRENDA turnover number [1/s]",
    ylabel = "Improved turnover number [1/s]",
)
scatter!(brenda_ax, kb, kgd, color = kss)
lb = 10^-2
ub = 10^4
lines!(brenda_ax, [lb, ub], [lb, ub], color = ColorSchemes.Greys_9[3], linestyle = :dash)
xlims!(brenda_ax, lb, ub)
ylims!(brenda_ax, lb, ub)

hidexdecorations!(brenda_ax, label = false, ticklabels = false, ticks = false)
hideydecorations!(brenda_ax, label = false, ticklabels = false, ticks = false)

elms = [
    MarkerElement(
        color = ColorSchemes.Set2_6[findfirst(x -> x == id, usubs)],
        marker = '●',
        markersize = 15,
    ) for id in usubs
]

Legend(
    gb[1, 1],
    elms,
    usubs,
    "Metabolic module",
    halign = :right,
    valign = :bottom,
    tellheight = false,
    tellwidth = false,
    margin = (10, 10, 10, 10),
)

for (label, layout) in zip(["A", "B"], [ga, gb])
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
CairoMakie.FileIO.save(joinpath("..", "DifferentiableMetabolismPaper", "docs", "imgs", "kmax_comparison.pdf"), fig)


using GLM
df = DataFrame(Y = log10.(kb), X = log10.(kgd))
f = lm(@formula(Y ~ X), df)
r2(f)

using Statistics
mean(frac_improv)