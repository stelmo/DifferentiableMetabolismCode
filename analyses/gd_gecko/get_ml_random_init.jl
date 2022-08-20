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

#: Plot first 

df = innerjoin(rd_df, gd_df, on = :KcatID)

xs = log10.(df[!, :Kcat_mean]./df[!, :Kcat_mean_gd])

count(-0.5 .< xs .< 0.5)/length(xs)

fig = Figure(resolution = (1200, 600), backgroundcolor = :transparent);

ga = fig[1, 1] = GridLayout()
gb = fig[1, 2] = GridLayout()

ax = Axis(
    ga[1,1],
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
    ga[1,1],
    [
        [PolyElement(color = :mediumorchid2, strokecolor = :mediumorchid2, strokewidth = 1),],
        [PolyElement(color = :tomato, strokecolor = :tomato, strokewidth = 1),],
    ],
    [
        "Optimal random start iterate\nwithin 3 fold of ML\nstart optimal iterate",
        "Optimal random start iterate\nnot changed wrt ML start\noptimal start iterate",
    ],
    tellheight = false,
    tellwidth = false,
    halign = :right, 
    valign = :top,
    margin = (10, 10, 10, 10),
)




#: Plot the next dataset
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
kbest_df = DataFrame(Condition = String[], KcatID = String[], Kcat = Float64[], Derivative = Float64[], Iteration = Float64[])
for master_id in master_ids
    try
        pdirfile = master_id * "#params.csv"
        dfp = DataFrame(CSV.File(joinpath("results", rdir, pdirfile)))

        df = @subset dfp @byrow begin 
            :Iteration == cond_minlossiter[master_id]
        end

        append!(
            kbest_df,
            df,
        )
    catch err
        println("failed on ", master_id)
    end
end

until_digit(x) = begin
    if startswith(x, "WT")
        return "WT"
    elseif startswith(x, "pgi")
        return "pgi"
    elseif startswith(x, "tpi")
        return "tpi"
    elseif startswith(x, "pts")
        return "pts"
    elseif startswith(x, "sdh")
        return "sdh"
    else
        throw(error("no match"))
    end
end

gd_df = transform(
    kbest_df,
    :Condition => x -> until_digit.(x),
    :KcatID => x -> last.(split.(x, "#")),
)
@select!(gd_df, :Kcat, :KcatID_function, :Condition_function)
rename!(gd_df, Dict("KcatID_function" => :KcatID, "Condition_function" => :Condition))
@transform!(gd_df, :Kcat = 1 / (3600 * 1e-6) .* :Kcat)

df = combine(
    groupby(gd_df, :KcatID),
    :Kcat => x -> std(x)/mean(x),
)

@subset! df @byrow begin 
    !isnan(:Kcat_function)
end

xs = df[!, :Kcat_function]

ax2 = Axis(
    gb[1,1],
    xlabel = "Ratio of standard deviation to mean of the\nfinal kcat estimates for each enzyme across all conditions",
    ylabel = "Density"
)

density!(ax2, xs; color= :darkturquoise, boundary=(0.0, 10.0))

hidexdecorations!(ax2, label = false, ticklabels = false, ticks = false)
hideydecorations!(ax2, label = false, ticklabels = false, ticks = false)

xlims!(ax2, -0.5, 6) # sigh density wraps around :/ 

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

CairoMakie.FileIO.save(joinpath("..", "DifferentiableMetabolismPaper", "docs", "imgs", "random_start_vs_ml_start_ratio.pdf"), fig)