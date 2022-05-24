using ColorSchemes, CairoMakie, DataFrames, DataFramesMeta, Chain, CSV
using Measurements, COBREXA

dirs = [
    "purB.csv"
    "pyrF.csv"
    "aroA.csv"
    "purC.csv"
    "gnd.csv"
    "metE.csv"
    "eno.csv"
    "cysH.csv"
    "glmS.csv"
    "zwf.csv"
    "gapA.csv"
    "ilvC.csv"
]

knockdown_df = transform(
    DataFrame(CSV.File(joinpath("data", "crispr", "target_gene.csv"))),
    :Reaction => ByRow(x -> rstrip.(split(x, "; "))) => :Reaction,
)
knockdown_df = flatten(knockdown_df, :Reaction)

metabolite_df = DataFrame(CSV.File(joinpath("data", "crispr", "metabolite_data.csv")))
rename!(metabolite_df, Dict("KEGG number" => :Kegg))

keggids = String.(metabolite_df[!, :Kegg])
model =
    load_model(StandardModel, joinpath("model_construction", "model_files", "iML1515.json"))
kegg_mid_lu = Dict()
for mid in metabolites(model)
    keggs = get(model.metabolites[mid].annotations, "kegg.compound", [])
    for kegg in keggs
        kegg_mid_lu[kegg] = mid[1:end-2]
    end
end

#: Plot figure
fig = Figure(
    resolution = (1200, 600),
    backgroundcolor=:transparent,
);

ax = Axis(
    fig[1, 1],
    ylabel = "Scaled sensitivity of biomass growth rate\nto metabolite concentration",
    xlabel = "Gene knocked down",
)

ax2 = Axis(
    fig[1, 1],
    ylabel = "Measured metabolite concentration fold change",
    yaxisposition = :right,
    yscale = log2,
)
hlines!(ax, 0.0, color = ColorSchemes.Greys_9[3])
vlines!(ax, 1:3:(length(dirs)*3), color = ColorSchemes.Greys_9[3], linestyle = :dash)

fdir(x, y) = x == 1 ? -1 : y == 1 ? 1 : 0
for (resi, res) in zip(1:3:length(dirs)*3, dirs)
    # (resi, res) = first(zip(1:3:length(dirs)*3, dirs))
    df = DataFrame(CSV.File(joinpath("results", "crispr", res)))
    @transform!(
        df,
        :Metabolite = last.(split.(:Metabolite, "#")),
        :SubOrProdorNot = fdir.(:Substrate, :Product),
    )
    @subset! df begin
        :Metabolite .!= "h_c"
        :Metabolite .!= "pi_c"
        :Metabolite .!= "h2o_c"
        :Metabolite .!= "h_p"
        :Metabolite .!= "pi_p"
        :Metabolite .!= "h2o_p"
        :Metabolite .!= "h_e"
        :Metabolite .!= "pi_e"
        :Metabolite .!= "h2o_e"
    end

    sub_df = unique(@subset(df, :SubOrProdorNot .== -1), :Metabolite)
    substrate_ids = @. first(split(sub_df[!, :Metabolite], "_"))
    substrate_sens = sub_df[!, :Sensitivity]

    prod_df = unique(@subset(df, :SubOrProdorNot .== 1), :Metabolite)
    product_ids = @. first(split(prod_df[!, :Metabolite], "_"))
    product_sens = prod_df[!, :Sensitivity]

    no_sub_df = unique(@subset(df, :SubOrProdorNot .== 0), :Metabolite)
    no_sub_sens = no_sub_df[!, :Sensitivity]

    xs = (resi - 0.25) .+ 0.01 .* randn(length(no_sub_sens))
    scatter!(
        ax,
        xs,
        no_sub_sens,
        color = ColorSchemes.Set2_6[2],
        marker = :rect,
        markersize = 5,
    )


    sub_xs = (resi - 0.25) .+ 0.01 .* randn(length(substrate_sens))
    scatter!(
        ax,
        sub_xs,
        substrate_sens,
        color = ColorSchemes.Set2_6[1],
        marker = :rect,
        markersize = 10,
    )

    prod_xs = (resi - 0.25) .+ 0.01 .* randn(length(product_sens))
    scatter!(
        ax,
        prod_xs,
        product_sens,
        color = ColorSchemes.Set2_6[3],
        marker = :rect,
        markersize = 10,
    )

    mids = [substrate_ids; product_ids]
    midys = [substrate_sens; product_sens]
    midxs = [sub_xs; prod_xs]

    idxs = sortperm(midys)
    mids = mids[idxs]
    midys = midys[idxs]
    midxs = midxs[idxs]
    has_changed = true
    while has_changed
        has_changed = false
        for i = 1:length(mids)-1
            if abs(midys[i+1] - midys[i]) < 0.03
                midys[i] -= 0.03
                has_changed = true
            end
        end
    end
    text!(
        ax,
        mids;
        position = [Point2f(x, y) for (x, y) in zip(midxs, midys)],
        textsize = 10.0,
        align = (:right, :baseline),
    )

    #: plot measured data
    coln = Symbol(first(split(res, ".")))

    tdf = @select(metabolite_df, :Kegg, $coln)
    dropmissing!(tdf)

    mids = [get(kegg_mid_lu, x, x) for x in tdf[!, :Kegg]]
    l2fc = tdf[!, coln]

    sub_idxs = findall(x -> x in substrate_ids, mids)
    prod_idxs = findall(x -> x in product_ids, mids)
    nonsub_idxs = findall(x -> x ∉ substrate_ids && x ∉ product_ids, mids)

    xs = (resi + 0.25) .+ 0.01 .* randn(length(l2fc))
    scatter!(
        ax2,
        xs[nonsub_idxs],
        l2fc[nonsub_idxs],
        color = ColorSchemes.Set2_6[5],
        markersize = 5,
        marker = :circle,
    )
    scatter!(
        ax2,
        xs[sub_idxs],
        l2fc[sub_idxs],
        color = ColorSchemes.Set2_6[4],
        markersize = 10,
        marker = :circle,
    )
    scatter!(
        ax2,
        xs[prod_idxs],
        l2fc[prod_idxs],
        color = ColorSchemes.Set2_6[6],
        markersize = 10,
        marker = :circle,
    )
    text!(
        ax2,
        mids[[prod_idxs; sub_idxs]];
        position = [
            Point2f(x, y) for
            (x, y) in zip(xs[[prod_idxs; sub_idxs]], l2fc[[prod_idxs; sub_idxs]])
        ],
        textsize = 10.0,
        align = (:left, :baseline),
    )
end

hidexdecorations!(ax, label = false, ticklabels = false, ticks = false)
hideydecorations!(ax, label = false, ticklabels = false, ticks = false)
ax.xticks = (1:3:(length(dirs)*3), first.(split.(dirs, ".")))

hidespines!(ax2)
hidexdecorations!(ax2)
hideydecorations!(ax2, label = false, ticklabels = false, ticks = false)
ylims!(ax2, 2^-10, 2^10)
ylims!(ax, -1.1, 1.1)
linkxaxes!(ax, ax2)

elem1 = [MarkerElement(color = ColorSchemes.Set2_6[1], marker = :rect, markersize = 10)]
elem2 = [MarkerElement(color = ColorSchemes.Set2_6[2], marker = :rect, markersize = 10)]
elem3 = [MarkerElement(color = ColorSchemes.Set2_6[3], marker = :rect, markersize = 10)]
elem4 = [MarkerElement(color = ColorSchemes.Set2_6[4], marker = :circle, markersize = 10)]
elem5 = [MarkerElement(color = ColorSchemes.Set2_6[5], marker = :circle, markersize = 10)]
elem6 = [MarkerElement(color = ColorSchemes.Set2_6[6], marker = :circle, markersize = 10)]

Legend(
    fig[1, 1],
    [[elem1, elem2, elem3], [elem4, elem5, elem6]],
    [
        ["Substrate", "Other metabolite", "Product"],
        ["Substrate", "Other metabolite", "Product"],
    ],
    ["Model", "Experiment"],
    patchsize = (10, 10),
    # rowgap = 10,
    margin = (10, 10, 10, 10),
    halign = :left,
    valign = :bottom,
    tellwidth = false,
    tellheight = false,
    orientation = :horizontal,
    nbanks = 3,
)


fig

CairoMakie.FileIO.save(joinpath("..", "DifferentiableMetabolismPaper", "docs", "imgs", "crispr.pdf"), fig)
