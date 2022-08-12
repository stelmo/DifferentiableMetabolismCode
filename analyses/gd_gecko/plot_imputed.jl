using JSON, ColorSchemes, CairoMakie, DataFrames, DataFramesMeta, Chain, CSV, COBREXA, Statistics, GLM, Colors

#: load kcat data from GD
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

#: load kcat data from heckmann
heckmann_df = DataFrame(CSV.File(joinpath("results", "gd_gecko", "heckmann_df.csv")))

#: Load E. coli data from BRENDA
brenda_df = DataFrame(CSV.File(joinpath("data", "gecko2", "source_data_1.csv")))
ecoli_df = @subset brenda_df @byrow begin
    occursin("escherichia coli", :Organism) 
end
mean_ecoli_df = transform(
    combine(groupby(ecoli_df, :ECnumber), :Kcat => mean => :Kcat),
    :ECnumber => x -> last.(split.(x, "EC")),
)
select!(mean_ecoli_df, :Kcat, :ECnumber_function)
rename!(mean_ecoli_df, Dict("ECnumber_function" => :EC))

#: Load model 
base_load_path = joinpath("model_construction", "processed_models_files", "ecoli")
model = load_model(StandardModel, joinpath(base_load_path, "iml1515_fixed.json"))

#: load measurements
data_dir = joinpath("data", "kcats", "heckmann2020")
proteindata = JSON.parsefile(joinpath(data_dir, "proteindata.json"))
measured_genes = [collect(keys(x)) for x in values(proteindata)]
umgs = unique(reduce(vcat, measured_genes)) # unique measured genes 
mrids = String[]
for rid in reactions(model)
    isnothing(reaction_gene_association(model, rid)) && continue
    for grr in reaction_gene_association(model, rid)
        all(gid in umgs for gid in grr) && push!(mrids, rid)
        break
    end
end 
umrids = unique(mrids)

imputed_rids = setdiff(reactions(model), umrids)
imputed_rids_with_gd = DataFrame(KcatID = intersect(imputed_rids, gd_df[!, "KcatID"]))

heckmann_df_unmeasured = innerjoin(heckmann_df, imputed_rids_with_gd; on = :KcatID)
rename!(heckmann_df_unmeasured, Dict("Kmax" => :KmaxML))
gd_df_unmeasured = innerjoin(gd_df, imputed_rids_with_gd; on = :KcatID)
rename!(gd_df_unmeasured, Dict("Kmax" => :KmaxGD))

gd_heckmann = innerjoin(gd_df_unmeasured, heckmann_df_unmeasured; on =:KcatID)

#: connect brenda data to reactions 
rid_ecs = Dict()
for rid in gd_heckmann[!, :KcatID]
    annos = reaction_annotations(model, String(rid))
    if haskey(annos, "ec-code")
        rid_ecs[rid] = annos["ec-code"]
    end
end

ec_kcat = Dict(String(ec) => kcat for (ec, kcat) in zip(mean_ecoli_df[!, :EC], mean_ecoli_df[!, :Kcat]))

rid_kcat = Dict()
for (rid, ecs) in rid_ecs
    x = Float64[]
    for ec in ecs 
        haskey(ec_kcat, ec) && push!(x, ec_kcat[ec])
    end
    !isempty(x) && (rid_kcat[rid] = mean(x))
end

cleaned_brenda_df = DataFrame(KcatID = string.(collect(keys(rid_kcat))), Kcat = float.(collect(values(rid_kcat))))

all_data = innerjoin(gd_heckmann, cleaned_brenda_df; on = :KcatID)

#ml
df = DataFrame(Y = log10.(all_data[!, :KmaxML]), X = log10.(all_data[!, :Kcat]))
f = lm(@formula(Y ~ X), df)
r2(f)

#gd
df = DataFrame(Y = log10.(all_data[!, :KmaxGD]), X = log10.(all_data[!, :Kcat]))
g = lm(@formula(Y ~ X), df)
r2(g)

#: Plot figure
fig = Figure(
    # resolution = (1200, 1200),
    backgroundcolor = :transparent,
);
ax = Axis(
    fig[1, 1],
    yscale = log10,
    xscale = log10,
    xlabel = "BRENDA turnover number [1/s]",
    ylabel = "Estimated turnover number [1/s]",
)
scatter!(ax, all_data[!, :Kcat], all_data[!, :KmaxML], color=:tomato, label= "ML")
scatter!(ax, all_data[!, :Kcat], all_data[!, :KmaxGD], color=:lightskyblue, label = "GD")

bmin, bmax = extrema(all_data[!, :Kcat])
b_ex = DataFrame(X=log10.([bmin, bmax]))
f_ys = predict(f, b_ex)
lines!(ax, [bmin, bmax], 10 .^(f_ys), color=:tomato)

g_ys = predict(g, b_ex)
lines!(ax, [bmin, bmax], 10 .^(g_ys), color=:lightskyblue)


hidexdecorations!(ax, ticks = false, ticklabels = false, label = false)
hideydecorations!(ax, ticks = false, ticklabels = false, label = false)
axislegend("Estimation method", position = :rb)
fig

CairoMakie.FileIO.save(joinpath("..", "DifferentiableMetabolismPaper", "docs", "imgs", "unmeasured_proteome_kcats.pdf"), fig)
