using COBREXA, JSON
include(joinpath("model_construction", "preprocess.jl"))
using .Preprocess

#: load data from papers and uniprot
proteome_data_loc = joinpath("data", "proteome", "e_coli.json")
proteome_data = JSON.parsefile(proteome_data_loc)
kcat_data_loc = joinpath("data", "kcats", "heckmann2020", "e_coli_kcats.json")
kcat_data = JSON.parsefile(kcat_data_loc)

#: load model
model =
    load_model(StandardModel, joinpath("model_construction", "model_files", "iML1515.json"))

#! increase protein concentrations by changing units
#! convert 1/s to M/h 
scale = 3600 * 1e-6

organism = "e_coli"

Preprocess.rm_spontaneous!(model; spont_rxn_genes = ["s0001"])
Preprocess.add_kcats!(model, kcat_data; transporter_kcat = 65.0, average_kcat = 25.0)
protein_stoichiometry =
    Preprocess.get_protein_stoichiometry!(model, proteome_data, organism)
protein_masses = Preprocess.get_protein_masses(model, proteome_data)
reaction_kcats = Preprocess.get_reaction_kcats(model, kcat_data, scale)

#! remove really low abundance metabolites from bof -  numerical precision issues
for (k, v) in model.reactions["BIOMASS_Ec_iML1515_core_75p37M"].metabolites
    if abs(v) < 1e-4 # cutoff
        println("Deleted ", k, " with coefficient ", v)
        delete!(model.reactions["BIOMASS_Ec_iML1515_core_75p37M"].metabolites, k)
    end
end

#! delete other biomass reaction
delete!(model.reactions, "BIOMASS_Ec_iML1515_WT_75p37M")

Preprocess.save_all(
    "ecoli",
    "iml1515_fixed",
    model,
    protein_masses,
    reaction_kcats,
    protein_stoichiometry,
)
