using COBREXA, JSON
include(joinpath("model_construction", "preprocess.jl"))
using .Preprocess

#: load data from papers and uniprot
proteome_data_loc = joinpath("data", "proteome", "e_coli.json")
proteome_data = JSON.parsefile(proteome_data_loc)
kcat_data_loc = joinpath("data", "kcats", "heckmann2020", "e_coli_kcats.json")
kcat_data = JSON.parsefile(kcat_data_loc)

#: add kcats for electron transport chain
kcat_data["ATPS4r"] = kcat_data["ATPS4rpp"]
kcat_data["CYTBD"] = kcat_data["CYTBDpp"]
kcat_data["NADH16"] = kcat_data["NADH17pp"]
kcat_data["FRD7"] = kcat_data["FRD2"]
kcat_data["THD2"] = kcat_data["THD2pp"]

#: load model
model_loc = joinpath("model_construction", "model_files", "e_coli_core.json")
model = load_model(StandardModel, model_loc)

#! increase protein concentrations by changing units
#! convert 1/s to k/h 
scale = 3600 * 1e-3

organism = "e_coli"

Preprocess.rm_spontaneous!(model; spont_rxn_genes = ["s0001"])
Preprocess.add_kcats!(model, kcat_data; transporter_kcat = 65.0, average_kcat = 25.0)
protein_stoichiometry =
    Preprocess.get_protein_stoichiometry!(model, proteome_data, organism)
protein_masses = Preprocess.get_protein_masses(model, proteome_data)
reaction_kcats = Preprocess.get_reaction_kcats(model, kcat_data, scale)

Preprocess.save_all(
    "ecoli_core",
    "e_coli_core_fixed",
    model,
    protein_masses,
    reaction_kcats,
    protein_stoichiometry,
)


