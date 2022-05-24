using CSV, DataFrames, DataFramesMeta, COBREXA, JSON

df = DataFrame(CSV.File(joinpath("data", "proteome_ecoli", "protein_mass_separated.csv")))

base_load_path = joinpath("model_construction", "processed_models_files", "ecoli")
model = load_model(StandardModel, joinpath(base_load_path, "iml1515_fixed.json"))
reaction_protein_stoichiometry =
    JSON.parsefile(joinpath(base_load_path, "protein_stoichiometry.json"))

gname_gid = Dict()
for gid in genes(model)
    gn = model.genes[gid].name
    if !isnothing(gn)
        gname_gid[gn] = gid
    end
end

gn_gm = Dict(gname_gid[k] => v for (k, v) in zip(df[!, :Gene], df[!, :Glucose_fg_cell]) if haskey(gname_gid, k))
# gn_gm = Dict(gname_gid[k] => v for (k, v) in zip(df[!, :Gene], df[!, :Chemostat_mu_0_5]) if haskey(gname_gid, k))


rid_mass = Dict()
for rid in reactions(model)
    grrs = reaction_gene_association(model, rid)
    isnothing(grrs) && continue
    total_mass = 0
    for (grr, sts) in zip(grrs, reaction_protein_stoichiometry[rid])
        for (gid, st) in zip(grr, sts)
            total_mass += st * get(gn_gm, gid, 0.0)
        end
    end
    rid_mass[rid] = total_mass
end

ks = collect(keys(rid_mass))
vs = collect(values(rid_mass))

open(joinpath("analyses", "gecko", "rid_mass.json"), "w") do io 
    JSON.print(io, Dict(ks .=> vs))
end

open(joinpath("analyses", "gecko", "gid_mass.json"), "w") do io 
    JSON.print(io, gn_gm)
end
