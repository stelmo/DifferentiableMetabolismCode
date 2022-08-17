using JSON, DataFrames, DataFramesMeta, Chain, CSV, Statistics

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

gd_df = transform(
    combine(groupby(kbest_df, :KcatID), :Kcat => maximum => :Kmax),
    :KcatID => x -> last.(split.(x, "#")),
)

@select!(gd_df, :Kmax, :KcatID_function)
rename!(gd_df, Dict("KcatID_function" => :KcatID))
@transform!(gd_df, :Kmax = 1 / (3600 * 1e-6) .* :Kmax)

#: load ML kcats 
ml_kcats = JSON.parsefile(joinpath("data", "kcats", "heckmann2020", "e_coli_kcats.json"))
ml_df = DataFrame(Kcat=collect(values(ml_kcats)), KcatID=collect(keys(ml_kcats)))
df = leftjoin(ml_df, gd_df, on=:KcatID) # units = 1/s

# number changed
@subset df @byrow begin 
    !ismissing(:Kcat) && ! (abs(:Kcat - :Kmax) < 0.1) # substantially changed
end

rename!(df, Dict("Kcat" => "Kcat_(1/s)_Heckmann2020", "Kmax" => "Kcat_(1/s)_Wilken2022"))


CSV.write(
    joinpath("..", "DifferentiableMetabolismPaper", "data", "supplementary_dataset_1.csv"),
    df,
)
