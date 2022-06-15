using JSON, DataFrames, DataFramesMeta, Chain, CSV, Statistics

#: load kcat data
rdir = "linesearch"
params_dir = filter(endswith("params.csv"), readdir(joinpath("results", rdir)))

kmax_df = DataFrame(Condition = String[], KcatID = String[], Kmax = Float64[])
for dir in params_dir
    try
        gdf = groupby(DataFrame(CSV.File(joinpath("results", rdir, dir))), :KcatID)
        df = combine(gdf, :Kcat => maximum => :Kmax)
        cond = join(split(dir, "#")[1:2], "#")
        n = size(df, 1)
        append!(
            kmax_df,
            DataFrame(
                Condition = fill(cond, n),
                KcatID = df[!, :KcatID],
                Kmax = df[!, :Kmax],
            ),
        )
    catch err
        println("failed on ", dir)
    end
end

gd_df = transform(
    combine(groupby(kmax_df, :KcatID), :Kmax => maximum => :Kmax),
    :KcatID => x -> last.(split.(x, "#")),
)

@select!(gd_df, :Kmax, :KcatID_function)
rename!(gd_df, Dict("KcatID_function" => :KcatID))
@transform!(gd_df, :Kmax = 1 / (3600 * 1e-6) .* :Kmax)

#: load ML kcats 
ml_kcats = JSON.parsefile(joinpath("data", "kcats", "heckmann2020", "e_coli_kcats.json"))

ml_df = DataFrame(Kcat=collect(values(ml_kcats)), KcatID=collect(keys(ml_kcats)))

df = leftjoin(ml_df, gd_df, on=:KcatID) # units = 1/s
rename!(df, Dict("Kcat" => "Kcat_(1/s)_Heckmann2020", "Kmax" => "Kcat_(1/s)_Wilken2022"))
CSV.write(
    joinpath("..", "DifferentiableMetabolismPaper", "data", "supplementary_dataset_1.csv"),
    df,
)
