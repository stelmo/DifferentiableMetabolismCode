using ColorSchemes, CairoMakie, DataFrames, DataFramesMeta, Chain, CSV
using Measurements, COBREXA

davidi_df = DataFrame(CSV.File(joinpath("results", "gd_gecko", "davidi_df.csv")))
heckmann_df = DataFrame(CSV.File(joinpath("results", "gd_gecko", "heckmann_df.csv")))
polished_df = DataFrame(CSV.File(joinpath("results", "gd_gecko", "polish_kmaxs_df.csv")))
@rtransform!(polished_df, :KcatID = last(split(:KcatID,"#")))

vs_brenda = @select!(innerjoin(polished_df, davidi_df, on=:KcatID, makeunique=true), :KcatID, :Kcat, :Kmax) 
vs_davidi = @select!(innerjoin(polished_df, davidi_df, on=:KcatID, makeunique=true), :KcatID, :Kmax, :Kmax_1)
vs_heckmann = innerjoin(polished_df, heckmann_df, on=:KcatID, makeunique=true)
brenda_heckman = innerjoin(@select(davidi_df, :KcatID, :Kcat), heckmann_df, on=:KcatID, makeunique=true)

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
rid_subsystem = Dict{String, String}()
for rid in reactions(model)
    rid_subsystem[rid] = get(subs_lu, reaction_subsystem(model, rid), "Other")
end

_rid_subsystem(x) = get(rid_subsystem, x, "Other")
@transform!(vs_brenda, :Log10Diff = log10.(:Kmax) .- log10.(:Kcat), :FracDiff = (:Kmax .- :Kcat)./:Kcat, :Subsystem = _rid_subsystem.(:KcatID))
@transform!(vs_davidi, :Log10Diff = log10.(:Kmax) .- log10.(:Kmax_1), :FracDiff = (:Kmax .- :Kmax_1)./:Kmax_1, :Subsystem = _rid_subsystem.(:KcatID))
@transform!(vs_heckmann, :Log10Diff = log10.(:Kmax) .- log10.(:Kmax_1), :FracDiff = (:Kmax .- :Kmax_1)./:Kmax_1, :Subsystem = _rid_subsystem.(:KcatID))
@transform!(brenda_heckman, :Subsystem = _rid_subsystem.(:KcatID))

usubs = unique(vs_brenda[!, :Subsystem]) # unique subsystems

#: Plot figure
fig = Figure(
    resolution = (1200, 1200)
    # backgroundcolor = :transparent,
);

ga = fig[1, 1] = GridLayout()
gb = fig[1, 2] = GridLayout()
gc = fig[2, 1] = GridLayout()
gd = fig[2, 2] = GridLayout()

lb = 10^-2
ub = 10^4

#: Brenda 
kb = vs_brenda[!, :Kcat]
kb_imp = vs_brenda[!, :Kmax]

kss = [
    ColorSchemes.Set2_6[findfirst(x -> x == id, usubs)] for id in vs_brenda[!, :Subsystem]
]

brenda_ax = Axis(
    ga[1, 1],
    yscale = log10,
    xscale = log10,
    xlabel = "BRENDA turnover number [1/s]",
    ylabel = "Improved turnover number [1/s]",
)
scatter!(brenda_ax, kb, kb_imp, color = kss)

lines!(brenda_ax, [lb, ub], [lb, ub], color = ColorSchemes.Greys_9[3], linestyle = :dash)
xlims!(brenda_ax, lb, ub)
ylims!(brenda_ax, lb, ub)

hidexdecorations!(brenda_ax, label = false, ticklabels = false, ticks = false)
hideydecorations!(brenda_ax, label = false, ticklabels = false, ticks = false)

elms = [
    MarkerElement(
        color = ColorSchemes.Set2_6[findfirst(x -> x == id, usubs)],
        marker = '●',
        markersize = 15,
    ) for id in usubs
]

Legend(
    ga[1, 1],
    elms,
    usubs,
    "Metabolic module",
    halign = :right,
    valign = :bottom,
    tellheight = false,
    tellwidth = false,
    margin = (10, 10, 10, 10),
)

#: Davidi 
kb = vs_davidi[!, :Kmax_1]
kb_imp = vs_davidi[!, :Kmax]
# usubs = unique(vs_davidi[!, :Subsystem]) # unique subsystems
kss = [
    ColorSchemes.Set2_6[findfirst(x -> x == id, usubs)] for id in vs_davidi[!, :Subsystem]
]

davidi_ax = Axis(
    gb[1, 1],
    yscale = log10,
    xscale = log10,
    xlabel = "Davidi et. al. turnover number [1/s]",
    ylabel = "Improved turnover number [1/s]",
)
scatter!(davidi_ax, kb, kb_imp, color = kss)
lines!(davidi_ax, [lb, ub], [lb, ub], color = ColorSchemes.Greys_9[3], linestyle = :dash)
xlims!(davidi_ax, lb, ub)
ylims!(davidi_ax, lb, ub)

hidexdecorations!(davidi_ax, label = false, ticklabels = false, ticks = false)
hideydecorations!(davidi_ax, label = false, ticklabels = false, ticks = false)

elms = [
    MarkerElement(
        color = ColorSchemes.Set2_6[findfirst(x -> x == id, usubs)],
        marker = '●',
        markersize = 15,
    ) for id in usubs
]

Legend(
    gb[1, 1],
    elms,
    usubs,
    "Metabolic module",
    halign = :right,
    valign = :bottom,
    tellheight = false,
    tellwidth = false,
    margin = (10, 10, 10, 10),
)

#: Heckmann
kb = vs_heckmann[!, :Kmax_1]
kb_imp = vs_heckmann[!, :Kmax]
# usubs = unique(vs_heckmann[!, :Subsystem]) # unique subsystems
kss = [
    ColorSchemes.Set2_6[findfirst(x -> x == id, usubs)] for id in vs_heckmann[!, :Subsystem]
]

heckmann_ax = Axis(
    gc[1, 1],
    yscale = log10,
    xscale = log10,
    xlabel = "Heckmann et. al. turnover number [1/s]",
    ylabel = "Improved turnover number [1/s]",
)
scatter!(heckmann_ax, kb, kb_imp, color = kss)
lines!(heckmann_ax, [lb, ub], [lb, ub], color = ColorSchemes.Greys_9[3], linestyle = :dash)
xlims!(heckmann_ax, lb, ub)
ylims!(heckmann_ax, lb, ub)

hidexdecorations!(heckmann_ax, label = false, ticklabels = false, ticks = false)
hideydecorations!(heckmann_ax, label = false, ticklabels = false, ticks = false)

elms = [
    MarkerElement(
        color = ColorSchemes.Set2_6[findfirst(x -> x == id, usubs)],
        marker = '●',
        markersize = 15,
    ) for id in usubs
]

Legend(
    gc[1, 1],
    elms,
    usubs,
    "Metabolic module",
    halign = :right,
    valign = :bottom,
    tellheight = false,
    tellwidth = false,
    margin = (10, 10, 10, 10),
)

#: Heckman vs Brenda 
kb = brenda_heckman[!, :Kcat]
kb_heck = brenda_heckman[!, :Kmax]
# usubs = unique(brenda_heckman[!, :Subsystem]) # unique subsystems
kss = [
    ColorSchemes.Set2_6[findfirst(x -> x == id, usubs)] for id in brenda_heckman[!, :Subsystem]
]

brenda_heckman_ax = Axis(
    gd[1, 1],
    yscale = log10,
    xscale = log10,
    xlabel = "BRENDA turnover number [1/s]",
    ylabel = "Heckmann et. al. turnover number [1/s]",
)
scatter!(brenda_heckman_ax, kb, kb_heck, color = kss)

lines!(brenda_heckman_ax, [lb, ub], [lb, ub], color = ColorSchemes.Greys_9[3], linestyle = :dash)
xlims!(brenda_heckman_ax, lb, ub)
ylims!(brenda_heckman_ax, lb, ub)

hidexdecorations!(brenda_heckman_ax, label = false, ticklabels = false, ticks = false)
hideydecorations!(brenda_heckman_ax, label = false, ticklabels = false, ticks = false)

elms = [
    MarkerElement(
        color = ColorSchemes.Set2_6[findfirst(x -> x == id, usubs)],
        marker = '●',
        markersize = 15,
    ) for id in usubs
]

Legend(
    gd[1, 1],
    elms,
    usubs,
    "Metabolic module",
    halign = :right,
    valign = :bottom,
    tellheight = false,
    tellwidth = false,
    margin = (10, 10, 10, 10),
)

#: get R2 of brenda vs heckmann
using GLM, DataFrames
df = DataFrame(Y = log10.(kb_heck), X = log10.(kb))
f = lm(@formula(Y ~ X), df)
r2(f)

#: Add layout

for (label, layout) in zip(["A", "B", "C", "D"], [ga, gb, gc, gd])
    Label(
        layout[1, 1, TopLeft()],
        label,
        textsize = 26,
        padding = (0, 5, 5, 0),
        halign = :right,
    )
end

fig

CairoMakie.FileIO.save(joinpath("..", "DifferentiableMetabolismPaper", "docs", "imgs", "kmax_brenda_davidi_heckmann.pdf"), fig)