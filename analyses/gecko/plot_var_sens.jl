using ColorSchemes, CairoMakie, DataFrames, DataFramesMeta, Chain, CSV, Statistics, COBREXA
const notnan(x) = !isnan(x)

model = load_model(StandardModel, joinpath("model_construction", "model_files", "iML1515.json"))
rid_subsystem = Dict{String, String}()
for rid in reactions(model)
    rid_subsystem[rid] = reaction_subsystem(model, rid)
end

df = DataFrame(CSV.File(joinpath("results", "gecko", "variability_sensitivities.csv")))

@subset! df begin 
    abs.(:Sensitivity) .> 1e-3
end

sdf = @combine groupby(df, :Parameter) begin 
    :Mean = mean(:Sensitivity)
    :Std = std(:Sensitivity)
    :StdMean = abs(std(:Sensitivity)/mean(:Sensitivity))
end

@subset! sdf begin 
    notnan.(:StdMean)
end

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

fig = Figure();
ax = Axis(
    fig[1, 1], 
    # yscale=log10,
);

density!(sdf[!, :Std])

fig
