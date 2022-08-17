using JSON, DataFrames, DataFramesMeta, Chain, CSV, Statistics, Colors, CairoMakie

#: Load random starting point GD kcats
rdir = "gd_random_init"
losses_dir = filter(endswith("losses.csv"), readdir(joinpath("results", rdir)))
params_dir = filter(endswith("params.csv"), readdir(joinpath("results", rdir)))

dfl = DataFrame(Condition = String[], Loss = Float64[], Iteration = Int64[])
for dir in losses_dir
    append!(dfl, DataFrame(CSV.File(joinpath("results", rdir, dir))))
end
master_ids = unique(dfl[!, :Condition])

cond_minlossiter = Dict()
for gdf in groupby(dfl, :Condition)
    losses = gdf[!, :Loss]
    idx = argmin(losses)
    cond_minlossiter[gdf[idx, :Condition]] = gdf[idx, :Iteration]
end

random_df = DataFrame(Condition = String[], KcatID = String[], Kcat = Float64[], Derivative = Float64[], Iteration = Float64[])
for master_id in master_ids
    try
        pdirfile = master_id * "#params.csv"
        dfp = DataFrame(CSV.File(joinpath("results", rdir, pdirfile)))

        df = @subset dfp @byrow begin 
            :Iteration == cond_minlossiter[master_id]
        end

        append!(
            random_df,
            df,
        )
    catch err
        println("failed on ", master_id)
    end
end

@transform!(random_df, :Kcat = 1 / (3600 * 1e-6) .* :Kcat)

@transform! random_df @byrow begin 
    :Condition = join(split(:Condition, "#")[1:2], "#")
end

@subset! random_df @byrow begin 
    abs(:Derivative) .< 1e-1 # only select values that are not changing much
end

rd_df = combine(groupby(random_df, :KcatID), :Kcat => mean => :Kcat_mean, :Kcat => std => :Kcat_std)

#: Load ML starting point GD kcats
rdir = "linesearch"
losses_dir = filter(endswith("losses.csv"), readdir(joinpath("results", rdir)))
params_dir = filter(endswith("params.csv"), readdir(joinpath("results", rdir)))

dfl = DataFrame(Condition = String[], Loss = Float64[], Iteration = Int64[])
for dir in losses_dir
    append!(dfl, DataFrame(CSV.File(joinpath("results", rdir, dir))))
end
master_ids = unique(dfl[!, :Condition])

cond_minlossiter = Dict()
for gdf in groupby(dfl, :Condition)
    losses = gdf[!, :Loss]
    idx = argmin(losses)
    cond_minlossiter[gdf[idx, :Condition]] = gdf[idx, :Iteration]
end

#: load kcat data
mlgd_df = DataFrame(Condition = String[], KcatID = String[], Kcat = Float64[], Derivative = Float64[], Iteration = Float64[])
for master_id in master_ids
    !startswith(master_id, "WT") && continue
    try
        pdirfile = master_id * "#params.csv"
        dfp = DataFrame(CSV.File(joinpath("results", rdir, pdirfile)))

        df = @subset dfp @byrow begin 
            :Iteration == cond_minlossiter[master_id]
        end

        append!(
            mlgd_df,
            df,
        )
    catch err
        println("failed on ", master_id)
    end
end

@transform!(mlgd_df, :Kcat = 1 / (3600 * 1e-6) .* :Kcat)

@subset! mlgd_df @byrow begin 
    abs(:Derivative) .< 1e-1 # only select values that are not changing much
end

_gd_df = combine(groupby(mlgd_df, :KcatID), :Kcat => mean => :Kcat_mean_gd, :Kcat => std => :Kcat_std_gd)
gd_df = @subset _gd_df @byrow begin
    !isnan(:Kcat_std_gd)
end 

#: Plot

df = innerjoin(rd_df, gd_df, on = :KcatID)

xs = log10.(df[!, :Kcat_mean]./df[!, :Kcat_mean_gd])

count(-0.5 .< xs .< 0.5)/length(xs)

fig = Figure(
    # resolution = (1200, 1200),
    backgroundcolor = :transparent,
)
ax = Axis(
    fig[1,1],
    xlabel = "Log ratio of random start kcat to ML start kcat",
    ylabel = "Density"
)

density!(ax, xs; npoints=50, color= :darkturquoise)
lines!(ax, [0,0], [0, 1.5], color = :tomato, linewidth=8)
lines!(ax, [-0.5,-0.5], [0, 1.], color = :mediumorchid2, linewidth=6)
lines!(ax, [0.5, 0.5], [0, 1.], color = :mediumorchid2, linewidth=6)

arrows!(ax, [0,0],[1,1],[0.45,-0.45],[0.0,0], color=:mediumorchid2, linewidth=6,arrowsize=16)

hidexdecorations!(ax, label = false, ticklabels = false, ticks = false)
hideydecorations!(ax, label = false, ticklabels = false, ticks = false)

Legend(
    fig[1,1],
    [
        [PolyElement(color = :mediumorchid2, strokecolor = :mediumorchid2, strokewidth = 1),],
        [PolyElement(color = :tomato, strokecolor = :tomato, strokewidth = 1),],
    ],
    [
        "Optimal random start iterate within 3\nfold of ML start optimal iterate",
        "Optimal random start iterate not changed\nwrt ML start optimal start iterate",
    ],
    tellheight = false,
    tellwidth = false,
    halign = :right, 
    valign = :top,
    margin = (10, 10, 10, 10),
)

fig

CairoMakie.FileIO.save(joinpath("..", "DifferentiableMetabolismPaper", "docs", "imgs", "random_start_vs_ml_start_ratio.pdf"), fig)