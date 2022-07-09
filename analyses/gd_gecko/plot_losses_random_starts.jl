using JSON, ColorSchemes, CairoMakie, DataFrames, DataFramesMeta, Chain, CSV, COBREXA, Statistics, Colors

rdir = "linesearch"
losses_dir = filter(startswith("WT"), filter(endswith("losses.csv"), readdir(joinpath("results", rdir))))
dfl = DataFrame(Condition = String[], Loss = Float64[], Iteration = Int64[])
for dir in losses_dir
    append!(dfl, DataFrame(CSV.File(joinpath("results", rdir, dir))))
end

# no regularization, but random starts 
rdir2 = "gd_random_init"
losses_dir2 = filter(endswith("losses.csv"), readdir(joinpath("results", rdir2)))
dfl2 = DataFrame(Condition = String[], Loss = Float64[], Iteration = Int64[])
for dir in losses_dir2
    append!(dfl2, DataFrame(CSV.File(joinpath("results", rdir2, dir))))
end 

# # λ = 0.01
# rdir3 = "gd_random_init_reg"
# losses_dir3 = filter(endswith("losses.csv"), readdir(joinpath("results", rdir3)))
# dfl3 = DataFrame(Condition = String[], Loss = Float64[], Iteration = Int64[])
# for dir in losses_dir3
#     append!(dfl3, DataFrame(CSV.File(joinpath("results", rdir3, dir))))
# end 

# # λ = 0.1
# rdir4 = "gd_random_init_reg2"
# losses_dir4 = filter(endswith("losses.csv"), readdir(joinpath("results", rdir4)))
# dfl4 = DataFrame(Condition = String[], Loss = Float64[], Iteration = Int64[])
# for dir in losses_dir4
#     append!(dfl4, DataFrame(CSV.File(joinpath("results", rdir4, dir))))
# end 

#: Plot
fig = Figure(
    # resolution = (1200, 600), 
    # backgroundcolor = :transparent,
);
 
#: Losses
losses_ax = Axis(
    fig[1, 1],
    yscale=log10,
    xscale = log10,
    xlabel = "Iteration",
    ylabel = "Mean squared relative error (L)",
);
for condition in unique(dfl[!, :Condition])
    dflt = @subset(dfl, :Condition .== condition)
    sort!(dflt, :Iteration)
    color = :aqua
    lwidth = 2.0
    lines!(
        losses_ax,
        dflt[!, :Iteration], 
        dflt[!, :Loss],
        color = color,
        linewidth = lwidth,
    )
end

for condition in unique(dfl2[!, :Condition])
    dflt = @subset(dfl2, :Condition .== condition)
    sort!(dflt, :Iteration)
    color = :tomato
    lwidth = 2.0
    lines!(
        losses_ax,
        dflt[!, :Iteration], 
        dflt[!, :Loss],
        color = color,
        linewidth = lwidth,
    )
end

# for condition in unique(dfl3[!, :Condition])
#     dflt = @subset(dfl3, :Condition .== condition)
#     sort!(dflt, :Iteration)
#     color = :lightgreen
#     lwidth = 2.0
#     lines!(
#         losses_ax,
#         dflt[!, :Iteration], 
#         dflt[!, :Loss],
#         color = color,
#         linewidth = lwidth,
#     )
# end

# for condition in unique(dfl4[!, :Condition])
#     dflt = @subset(dfl4, :Condition .== condition)
#     sort!(dflt, :Iteration)
#     color = :gold
#     lwidth = 2.0
#     lines!(
#         losses_ax,
#         dflt[!, :Iteration], 
#         dflt[!, :Loss],
#         color = color,
#         linewidth = lwidth,
#     )
# end


hidexdecorations!(losses_ax, label = false, ticklabels = false, ticks = false)
hideydecorations!(losses_ax, label = false, ticklabels = false, ticks = false)

elms = [
    LineElement(color = :aqua, linewidth = 2),
    LineElement(color = :tomato, linewidth = 2),
]

Legend(
    fig[1,1],
    elms,
    ["ML starting points", "Random starting points"],
    halign = :right,
    valign = :top,
    tellheight = false,
    tellwidth = false,
    margin = (10, 10, 10, 10),
)

fig 
CairoMakie.FileIO.save(joinpath("..", "DifferentiableMetabolismPaper", "docs", "imgs", "random_starts.pdf"), fig)
