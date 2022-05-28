using JSON, ColorSchemes, CairoMakie, DataFrames, DataFramesMeta, Chain, CSV, COBREXA

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

#: load brenda data
davidi_df = DataFrame(CSV.File(joinpath("data", "kcats", "davidi2016", "brenda_kcats.csv")))
rename!(
    davidi_df,
    Dict(
        "ReactionID" => :KcatID,
    ),
)
brenda_df = @chain davidi_df begin
    @transform(:KcatID = "k#" .* first.(split.(:KcatID, "_")))
    @select(:KcatID, :Kcat)
end

#: load model
model = load_model(StandardModel, joinpath("model_construction", "model_files", "iML1515.json"))
subsys_df = DataFrame(
    KcatID = "k#" .* reactions(model),
    Subsystem = [model.reactions[rid].subsystem for rid in reactions(model)],
)

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

df_sys_map = DataFrame(Subsystem = subsys_map[:, 1], SubID = subsys_map[:, 2])
subsys_df = innerjoin(subsys_df, df_sys_map, on = :Subsystem)
@select!(subsys_df, :KcatID, :SubID)

#: find kmax overall
kmax_overall_df = combine(groupby(kmax_df, :KcatID), :Kmax => maximum => :Kmax)
@transform!(kmax_overall_df, :Kmax = 1 / (3600 * 1e-6) .* :Kmax)
CSV.write(joinpath("results", "gd_gecko", "polish_kmaxs_df.csv"), kmax_overall_df)

kmax_brenda_df = innerjoin(kmax_overall_df, brenda_df, on = :KcatID)
kmax_brenda_df = innerjoin(kmax_brenda_df, subsys_df, on = :KcatID)

CSV.write(joinpath("results", "gd_gecko", "kmax_brenda_df.csv"), kmax_brenda_df)
CSV.write(joinpath("results", "gd_gecko", "davidi_df.csv"), davidi_df)

base_load_path = joinpath("model_construction", "processed_models_files", "ecoli")
reaction_kcats = JSON.parsefile(joinpath(base_load_path, "reaction_kcats.json"))
rids = collect(keys(reaction_kcats))
heck = [first(vs)[1] for vs in values(reaction_kcats)] .* 1/3600 * 1e6
CSV.write(joinpath("results", "gd_gecko", "heckmann_df.csv"), 
    DataFrame(KcatID=rids,Kmax=heck),
)


