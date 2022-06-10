using ColorSchemes, CairoMakie, DataFrames, DataFramesMeta, Chain, CSV
using Measurements, COBREXA, GLM, Statistics

dirs = [ # order of paper
    "purB.csv"
    "pyrF.csv"
    "aroA.csv" 
    "purC.csv"
    "ptsH.csv"
    "pfkA.csv"
    "gnd.csv"
    "metE.csv"
    "sdh.csv"
    "adk.csv"
    "eno.csv"
    "pgi.csv"
    "cysH.csv"
    "fbaA.csv"
    "gdhA.csv"
    "pykF.csv"
    "tpiA.csv"
    "gltA.csv"
    # "glmS.csv"
    # "zwf.csv"
    # "ppc.csv"
    # "carA.csv"
    # "dxs.csv"
    # "gapA.csv"
    # "prs.csv"
    # "ilvC.csv"
    # "pykA.csv"
    # "pfkB.csv"
    # "pck.csv"
]

filter!(x -> x in readdir(joinpath("results", "crispr")), dirs)

metabolite_df = DataFrame(CSV.File(joinpath("data", "crispr", "metabolite_data.csv")))
rename!(metabolite_df, Dict("KEGG number" => :Kegg))

keggids = String.(metabolite_df[!, :Kegg])
model = load_model(StandardModel, joinpath("model_construction", "model_files", "iML1515.json"))
kegg_mid_lu = Dict()
for mid in metabolites(model)
    keggs = get(model.metabolites[mid].annotations, "kegg.compound", [])
    for kegg in keggs
        kegg_mid_lu[kegg] = mid[1:end-2]
    end
end
kegg_mid_lu["C01094"] = "g6p"

_lu(x) = get(kegg_mid_lu, x, x) 
@transform!(metabolite_df, :Metabolite = _lu.(:Kegg))

#: Plot figure
fig = Figure(
    # resolution = (1600, 600),
    # backgroundcolor=:transparent,
);

ax = Axis(
    fig[1, 1],
    ylabel = "Scaled sensitivity of biomass growth rate\nto metabolite concentration",
    xlabel = "Measured metabolite concentration fold change",
    xscale = log10,
)

fdir(x, y) = x == 1 ? -1 : y == 1 ? 1 : 0
all_df = DataFrame(Metabolite=String[], Sensitivity=Float64[], FoldChange=Float64[], KD=String[])
for res in dirs

    df = DataFrame(CSV.File(joinpath("results", "crispr", res)))
    @transform!(
        df,
        :Metabolite = last.(split.(:Metabolite, "#")),
        :SubOrProdorNot = fdir.(:Substrate, :Product),
    )
    @eachrow! df begin 
        :Metabolite = :Metabolite[1:end-2]
    end
    @subset! df begin
        :Metabolite .!= "h"
        :Metabolite .!= "pi"
        :Metabolite .!= "h2o"
        :Metabolite .!= "co2"
        # :Metabolite .!= "atp"
        # :Metabolite .!= "adp"
        :SubOrProdorNot .!= 0 
    end

    #: plot measured data
    coln = Symbol(first(split(res, ".")))
    tdf = @select(metabolite_df, :Metabolite, $coln)
    dropmissing!(tdf)

    jdf = innerjoin(df, tdf, on=:Metabolite)
    n = size(jdf, 1)
    append!(all_df, DataFrame(Metabolite=jdf[!, :Metabolite], Sensitivity=jdf[!, :Sensitivity], FoldChange = jdf[!, coln], KD=fill(string(coln), n)))
    if n > 0
        scatter!(ax, jdf[!, coln], jdf[!, :Sensitivity], label=string(coln))
    end
end
fig[1, 2] = Legend(fig, ax, "Knockdowns", framevisible = false)
fig

CairoMakie.FileIO.save(joinpath("..", "DifferentiableMetabolismPaper", "docs", "imgs", "crispr2.pdf"), fig)


all_df
using GLM

f = lm(@formula(Sensitivity ~ FoldChange), all_df)
r2(f)