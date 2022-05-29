using ColorSchemes, CairoMakie, DataFrames, DataFramesMeta, Chain, CSV, Statistics, COBREXA
const notnan(x) = !isnan(x)

model = load_model(StandardModel, joinpath("model_construction", "model_files", "iML1515.json"))
rid_subsystem = Dict{String, String}()
for rid in reactions(model)
    rid_subsystem[rid] = reaction_subsystem(model, rid)
end

df = DataFrame(CSV.File(joinpath("results", "gecko", "variability_sensitivities2.csv")))

@subset! df begin 
    abs.(:Sensitivity) .> 1e-6 # only analyze non-zero sensitivities
end

sdf = @combine groupby(df, :Parameter) begin 
    :Mean = mean(:Sensitivity)
    :Std = std(:Sensitivity)
    :StdMean = abs(std(:Sensitivity)/mean(:Sensitivity))
end

@subset! sdf begin 
    notnan.(:StdMean)
end

describe(sdf)

@rtransform! sdf begin
    :Subsystem = rid_subsystem[:Parameter]
end

ssdf = @combine groupby(sdf, :Subsystem) begin 
    :StdMean = mean(:StdMean)
    :StdStd = std(:StdMean)
end

@subset! ssdf begin 
    notnan.(:StdStd)
end

describe(ssdf)

fig = Figure(
    resolution = (600, 1000), 
    backgroundcolor = :transparent,
);


#: Plot hold outs
xlabs = ssdf[!, :Subsystem]
ys = ssdf[!, :StdMean]
yerr = ssdf[!, :StdStd]

ax = Axis(
    fig[1, 1],
    xlabel = "Metabolic module",
    ylabel = "Ratio of standard deviation to mean\n
    of all direct sensitivitie across\n
    all reactions in metabolic module",
    xticklabelrotation=-pi/2,
)

barplot!(
    ax,
    1:length(xlabs),
    ys,
    color = ColorSchemes.Set2_3[1],
)

# errorbars!(
#     ax,
#     1:length(xlabs),
#     ys,
#     yerr,
#     yerr,
#     whiskerwidth = 10,
# )

ax.xticks = (1:length(xlabs), xlabs)
hidexdecorations!(ax, label = false, ticklabels = false, ticks = false)
hideydecorations!(ax, label = false, ticklabels = false, ticks = false)

fig
CairoMakie.FileIO.save(joinpath("..", "DifferentiableMetabolismPaper", "docs", "imgs", "sensitivity_std.pdf"), fig)
