using JSON, DataFrames, DataFramesMeta, Chain, CSV, Statistics, COBREXA, CairoMakie, ColorSchemes, Colors

model = load_model(StandardModel, joinpath("model_construction", "model_files", "iML1515.json"))
subsys_df = DataFrame(
    KcatID = "k#" .* reactions(model),
    Subsystem = [model.reactions[rid].subsystem for rid in reactions(model)],
)

subsystem = "Glycolysis/Gluconeogenesis"
subsystem_kcatids = @subset(subsys_df, :Subsystem .== subsystem)[!, :KcatID]


#: Load data from random init 
rdir = "gd_random_init"
params_dir = filter(endswith("params.csv"), readdir(joinpath("results", rdir)))
maxiter = 0

rinit_dict = Dict{String, Dict{String, Dict{Int, Vector{Float64}}}}()

for dir in params_dir
    println(dir)
    random_init_df = DataFrame(CSV.File(joinpath("results", rdir, dir)))
    @select!(random_init_df, :Condition, :KcatID, :Kcat, :Iteration)
    @subset! random_init_df @byrow begin
        in(:KcatID, subsystem_kcatids)
    end
    dname = join(split(dir, "#")[1:2], "#")
    if !haskey(rinit_dict, dname)
        rinit_dict[dname] = Dict{String, Dict{Int, Vector{Float64}}}()
    end
    for gdf in groupby(random_init_df, :KcatID)
        sort!(gdf, [:Iteration])
        kname = String(first(gdf[!, :KcatID]))
        if haskey(rinit_dict[dname], kname)
            for (x, y) in zip(gdf[!, :Iteration], gdf[!, :Kcat])
                if haskey(rinit_dict[dname][kname], x)
                    push!(rinit_dict[dname][kname][x], y)
                else
                    rinit_dict[dname][kname][x] = [y]
                end
            end
        else
            rinit_dict[dname][kname] = Dict(x => [y] for (x, y) in zip(gdf[!, :Iteration], gdf[!, :Kcat]))
        end
    end
    maxiter = max(maxiter, maximum(random_init_df[!, :Iteration]))    
end
rinit_dict

#: load kcat from ML GD data
rdir = "linesearch"
params_dir = filter(startswith("WT"), filter(endswith("params.csv"), readdir(joinpath("results", rdir))))

ml_init_df = DataFrame(Condition = String[], KcatID = String[], Kmax = Float64[])

ml_dict = Dict{String, Dict{String, Dict{Int, Float64}}}()
for dir in params_dir
    df = DataFrame(CSV.File(joinpath("results", rdir, dir)))
    @subset! df @byrow begin
        in(:KcatID, subsystem_kcatids)
        :Iteration <= maxiter
    end
    dname = join(split(dir, "#")[1:2], "#")
    ml_dict[dname] = Dict{String, Dict{Int, Float64}}()
    for gdf in groupby(df, :KcatID)
        sort!(gdf, [:Iteration])
        kname = String(first(gdf[!, :KcatID]))
        ml_dict[dname][kname] = Dict(gdf[!, :Iteration] .=> gdf[!, :Kcat])
    end
end

ml_dict

#: Normalize
for (k, v) in rinit_dict
    for (kk, vv) in v 
        for (kkk, vvv) in vv 
            vvv ./= ml_dict[k][kk][kkk]
        end
    end
end
rinit_dict
#: Plot

p = subsystem_kcatids[2]
fig = Figure()
ax = Axis(fig[1,1], xscale=log10, yscale=log10)
css = ColorSchemes.Dark2_4

for (i, k) in enumerate(keys(rinit_dict))
    v = rinit_dict[k][p]
    for (kk, vv) in v 
        scatter!(ax, fill(kk, length(vv)), vv, color=css[i])
    end
end
fig



