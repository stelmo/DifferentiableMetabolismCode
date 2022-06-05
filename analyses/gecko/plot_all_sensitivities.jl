using CairoMakie, DataFrames, CSV, DataFramesMeta, Statistics, COBREXA, ColorSchemes

subsys_map = [
    "Nucleotide Salvage Pathway" "Nucleotide"
    "Tyrosine, Tryptophan, and Phenylalanine Metabolism" "Amino Acid"
    "Purine and Pyrimidine Biosynthesis" "Nucleotide"
    "Arginine and Proline Metabolism" "Amino Acid"
    "Citric Acid Cycle" "Carbon"
    "Cofactor and Prosthetic Group Biosynthesis" "Cofactor"
    "Pyruvate Metabolism" "Carbon"
    "Cysteine Metabolism" "Amino Acid"
    "Alternate Carbon Metabolism" "Carbon"
    "Anaplerotic Reactions" "Carbon"
    "Alanine and Aspartate Metabolism" "Amino Acid"
    "Glycolysis/Gluconeogenesis" "Carbon"
    "Extracellular exchange" "Membrane"
    "Valine, Leucine, and Isoleucine Metabolism" "Amino Acid"
    "Pentose Phosphate Pathway" "Carbon"
    "Oxidative Phosphorylation" "Carbon"
    "Lipopolysaccharide Biosynthesis / Recycling" "Membrane"
    "Membrane Lipid Metabolism" "Membrane"
    "Unassigned" "Other"
    "Methylglyoxal Metabolism" "Carbon"
    "tRNA Charging" "Other"
    "Threonine and Lysine Metabolism" "Amino Acid"
    "Folate Metabolism" "Cofactor"
    "Glycine and Serine Metabolism" "Amino Acid"
    "Glyoxylate Metabolism" "Carbon"
    "Cell Envelope Biosynthesis" "Membrane"
    "Methionine Metabolism" "Amino Acid"
    "Nitrogen Metabolism" "Amino Acid"
    "Histidine Metabolism" "Amino Acid"
    "Inorganic Ion Transport and Metabolism" "Cofactor"
    "Glutamate Metabolism" "Amino Acid"
    "Glycerophospholipid Metabolism" "Membrane"
    "Intracellular demand" "Other"
    "Transport, Outer Membrane Porin" "Membrane"
    "Transport, Inner Membrane" "Membrane"
    "Transport, Outer Membrane" "Membrane"
    "Murein Recycling" "Membrane"
    "Biomass and maintenance functions" "Other"
    "Murein Biosynthesis" "Membrane"
    "Metabolite Repair" "Other"
]
subs_lu = Dict(subsys_map[:,1] .=> subsys_map[:, 2])

model = load_model(StandardModel, joinpath("model_construction", "model_files", "iML1515.json"))
ss = Dict(rid => subs_lu[model.reactions[rid].subsystem] for rid in reactions(model))

df = DataFrame(CSV.File(joinpath("results", "gecko", "aerobic_sensitivity.csv")))

gd = groupby(df, :Parameter)
param_df = @combine(gd, :Ave = mean(:Sensitivity), :Max = maximum(:Sensitivity), :Min = minimum(:Sensitivity))

usubs = unique(values(subs_lu)) # unique subsystems
ssids = [ss[rid] for rid in param_df[!, :Parameter]]
kss = [
    ColorSchemes.Set2_6[findfirst(x -> x == id, usubs)] for id in ssids
]

ps = param_df[!, :Parameter]
maxfluxsens = @. abs(param_df[!, :Max])

fig = Figure(backgroundcolor = :transparent)
ax = Axis(
    fig[1,1],
    xlabel="Parameters",
    ylabel="Maximum absolute flux sensitivity",
    yscale=log10,
)

scatter!(ax, 1:length(ps), maxfluxsens, color=kss)

hidexdecorations!(ax, label = false)
hideydecorations!(ax, label = false, ticklabels = false, ticks = false)

elms = [
    MarkerElement(
        color = ColorSchemes.Set2_6[findfirst(x -> x == id, usubs)],
        marker = '●',
        markersize = 15,
    ) for id in usubs
]

Legend(
    fig[1, 1],
    elms,
    usubs,
    "Metabolic module",
    halign = :right,
    valign = :bottom,
    tellheight = false,
    tellwidth = false,
    margin = (10, 10, 10, 10),
)

fig
CairoMakie.FileIO.save(joinpath("..", "DifferentiableMetabolismPaper", "docs", "imgs", "average_sens.pdf"), fig)

