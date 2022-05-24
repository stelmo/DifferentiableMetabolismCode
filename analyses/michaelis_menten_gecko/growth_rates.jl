using DataFrames, DataFramesMeta, Chain, CSV, CairoMakie, ColorSchemes, LsqFit

growth_df = DataFrame(CSV.File(joinpath("data", "crispr", "growth_crispr.csv")))

function get_mu(df, col)
    ts = df[!, "Time (h)"]
    od = df[!, col]
    m(t, p) = p[1] * exp.(p[2] * t)
    p0 = [0.5, 0.5]
    f0 = findfirst(od .≈ -1.0)
    idx = isnothing(f0) ? length(od) : f0
    curve_fit(m, ts[1:idx], od[1:idx], p0)
end

df = DataFrame(Condition = String[], Growth = Float64[], GrowthSTD = Float64[])
for col in filter(!startswith("Time"), names(growth_df))
    fit = get_mu(growth_df, col)
    cov = estimate_covar(fit)
    p = fit.param
    push!(df, (col, p[2], cov[2, 2]))
end
df


CSV.write(joinpath("results", "michaelis_menten_gecko", "crispr_growth.csv"), df)
