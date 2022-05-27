using COBREXA, DifferentiableMetabolism
using CPLEX, JSON, SparseArrays, Statistics, Serialization, SparseArrays, LinearAlgebra
using DataFrames, DataFramesMeta, Chain, CSV, CairoMakie, ColorSchemes

include(joinpath("analyses", "thermo.jl"))
using .Thermo

#: load growth rates
growth_df =
    DataFrame(CSV.File(joinpath("results", "michaelis_menten_gecko", "crispr_growth.csv")))

#: import model and data files
base_load_path = joinpath("model_construction", "processed_models_files", "ecoli")
model = load_model(StandardModel, joinpath(base_load_path, "iml1515_fixed.json"))
reaction_kcats = JSON.parsefile(joinpath(base_load_path, "reaction_kcats.json"))
reaction_protein_stoichiometry =
    JSON.parsefile(joinpath(base_load_path, "protein_stoichiometry.json"))
gene_product_molar_mass = Dict(
    k => v for (k, v) in JSON.parsefile(joinpath(base_load_path, "protein_masses.json"))
)

model.reactions["EX_glc__D_e"].lb = -1000.0 # unconstrain because enzyme constraints take over

#: setup gecko data 
reaction_isozymes = Dict(
    rid => [
        Isozyme(
            Dict(
                k => v for (k, v) in zip(
                    reaction_gene_association(model, rid)[i],
                    reaction_protein_stoichiometry[rid][i],
                )
            ),
            reaction_kcats[rid][i][1],
            reaction_kcats[rid][i][2],
        ) for i = 1:length(reaction_kcats[rid])
    ] for rid in keys(reaction_kcats)
)
for rid in keys(reaction_isozymes)
    reaction_isozymes[rid] = [
        argmax(
            smoment_isozyme_speed(x -> gene_product_molar_mass[x]),
            reaction_isozymes[rid],
        ),
    ]
end

gene_product_bounds(gid) = (0, 100_000.0)
gene_product_mass_group_bound = Dict("uncategorized" => 600_000.0) # more than enough, minimized later

#: knockdown data
#=
Under parsimonious enzyme usage conditions, these reactions do not carry flux in the model.
=#
zero_flux = [
    "CBPS",
    "PAPSR2",
    "FBA3",
    "PPCK",
    "PFK_3",
    "PFK_2",
    "PRPPS",
    "DHAPT",
    "PYK",
    "PYK",
    "PYK4",
    "PYK2",
    "PYK3",
    "PYK6",
    "SUCDi",
]
knockdown_df = transform(
    DataFrame(CSV.File(joinpath("data", "crispr", "target_gene.csv"))),
    :Reaction => ByRow(x -> rstrip.(split(x, "; "))) => :Reaction,
)
knockdown_df = flatten(knockdown_df, :Reaction)
knockdown_df = @subset knockdown_df @byrow begin
    :Reaction ∉ zero_flux
end

#: Get reference and KD ids 
ref_cond = "pyrF"
kd_cond = ref_cond * " + aTC"
df = @subset knockdown_df begin
    :Target .== ref_cond
end
rid_kd = String(strip(df[1, :Reaction]))

μ_ref = first(growth_df[growth_df[!, :Condition].==ref_cond, :Growth])
μ_ko = first(growth_df[growth_df[!, :Condition].==kd_cond, :Growth])
target_gene = first(knockdown_df[knockdown_df[!, :Target].==ref_cond, :Gene])

gm = make_gecko_model(
    model;
    reaction_isozymes,
    gene_product_bounds,
    gene_product_molar_mass,
    gene_product_mass_group_bound,
)

opt_model = flux_balance_analysis(
    gm,
    CPLEX.Optimizer;
    modifications = [
        change_optimizer_attribute("CPXPARAM_Emphasis_Numerical", 1),
        change_optimizer_attribute("CPX_PARAM_SCAIND", 1),
        COBREXA.silence,
        change_objective(genes(gm); sense = COBREXA.MIN_SENSE),
        change_constraint("BIOMASS_Ec_iML1515_core_75p37M"; lb = μ_ref, ub = 1000),
    ],
)
rfluxes_ref = flux_dict(gm, opt_model)
gpconcs_ref = gene_product_dict(gm, opt_model)
flux_summary(rfluxes_ref)
enzyme_mass_ref = gene_product_mass_group_dict(gm, opt_model)["uncategorized"]

#: Get knockdown condition
kd_factor = 5
gene_product_bounds_kd(x) =
    x == target_gene ? (0.0, gpconcs_ref[x] / kd_factor) : (0.0, 100_000.0)
gene_product_mass_group_bound_kd = Dict("uncategorized" => 600_000.0) # more than enough, protein limitation limits everything

gm_kd = make_gecko_model(
    model;
    reaction_isozymes,
    gene_product_bounds = gene_product_bounds_kd,
    gene_product_molar_mass,
    gene_product_mass_group_bound = gene_product_mass_group_bound_kd,
)

opt_model = flux_balance_analysis(
    gm_kd,
    CPLEX.Optimizer;
    modifications = [
        change_optimizer_attribute("CPXPARAM_Emphasis_Numerical", 1),
        change_optimizer_attribute("CPX_PARAM_SCAIND", 1),
        COBREXA.silence,
    ],
)
rfluxes_kd = flux_dict(gm_kd, opt_model)
gpconcs_kd = gene_product_dict(gm_kd, opt_model)
flux_summary(rfluxes_kd)
enzyme_mass_kd = gene_product_mass_group_dict(gm_kd, opt_model)["uncategorized"]

pmodel = prune_model(model, rfluxes_kd)

loopless_sol = flux_balance_analysis_dict(
    pmodel,
    CPLEX.Optimizer;
    modifications = [add_loopless_constraints()],
)

pmodel = prune_model(pmodel, loopless_sol)


#: load km data
_rid_km = JSON.parsefile(joinpath("data", "kms", "kroll2021", "kmdata_iml1515.json"))
rid_km = Dict{String,Dict{String,Float64}}()
for (rid, td) in _rid_km
    rid ∉ reactions(pmodel) && continue
    isnothing(reaction_gene_association(pmodel, rid)) && continue
    isempty(reaction_gene_association(pmodel, rid)) && continue

    rid_km[rid] = Dict{String,Float64}()
    for (mid, km) in td
        rid_km[rid][mid] = km * 1e-3 # convert to Molar from mM
    end
end

del_list = [
    "SUCDi"
    "IGPDH"
    "IG3PS"
    "DHDPRy"
    "DMPPS"
    "IPDPS"
    "DHDPS"
    "ECOAH7"
    "DHAD2"
    "DHAD1"
    "ACACT7r"
    "MTHFR2"
    "PERD"
    "E4PD"
    # "CS"
    # "GND"
    # "PPC"
]
delete!.(Ref(rid_km), del_list) # trial and error

fix_list = [
    "CS"
    "GND"
    "PPC"
    # "DHAD1"
]

brenda_ave = JSON.parsefile(joinpath("data", "brenda", "parsed_brenda.json"))
for rid in fix_list
    for (mid, km) in rid_km[rid]
        rid_km[rid][mid] = get(brenda_ave[rid], mid, km)
    end
end

_km_concs = Dict{String,Vector{Float64}}()
for (rid, d) in rid_km
    rs = reaction_stoichiometry(pmodel, rid)
    for (mid, km) in d
        if loopless_sol[rid] > 0 && rs[mid] < 0
            _km_concs[mid] = push!(get(_km_concs, mid, Float64[]), km)
        elseif loopless_sol[rid] < 0 && rs[mid] > 0
            _km_concs[mid] = push!(get(_km_concs, mid, Float64[]), km)
        end
    end
end
km_concs = Dict(k => mean(v) for (k, v) in _km_concs)

#: get reasonable concentrations
rid_dg0 = Dict(
    rid => float(v) for (rid, v) in
    JSON.parsefile(joinpath("data", "thermodynamics", "iML1515_thermo.json")) if
    rid in reactions(pmodel) &&
    !isnothing(reaction_gene_association(pmodel, rid)) &&
    !isempty(reaction_gene_association(pmodel, rid))
)

subsysts = [
    "Citric Acid Cycle"
    "Glycolysis/Gluconeogenesis"
    "Pentose Phosphate Pathway"
    "Pyruvate Metabolism"
    "Anaplerotic Reactions"
    "Alternate Carbon Metabolism"
    "Oxidative Phosphorylation"
    "Valine, Leucine, and Isoleucine Metabolism"
    "Cysteine Metabolism"
    "Glycine and Serine Metabolism"
    "Threonine and Lysine Metabolism"
    "Methionine Metabolism"
    "Histidine Metabolism"
    "Tyrosine, Tryptophan, and Phenylalanine Metabolism"
    "Arginine and Proline Metabolism"
    "Alanine and Aspartate Metabolism"
    "Glutamate Metabolism"
    "Nitrogen Metabolism"
    "Purine and Pyrimidine Biosynthesis"
    "Lipopolysaccharide Biosynthesis / Recycling"
    "Membrane Lipid Metabolism"
    "Glycerophospholipid Metabolism"
    "Murein Biosynthesis"
    "Folate Metabolism"
    "Inorganic Ion Transport and Metabolism"
]
for rid in keys(rid_dg0)
    if model.reactions[rid].subsystem ∉ subsysts || rid ∉ reactions(pmodel)
        delete!(rid_dg0, rid)
    end
end

mmdf = Thermo.max_min_driving_force3(
    pmodel,
    rid_dg0,
    km_concs,
    2.0,
    CPLEX.Optimizer;
    flux_solution = loopless_sol,
    proton_ids = ["h_c", "h_e", "h_p"],
    water_ids = ["h2o_c", "h2o_e", "h2o_p"],
    concentration_lb = mid -> get(km_concs, mid, 1e-6) / 1000,
    concentration_ub = mid -> get(km_concs, mid, 0.1 / 10) * 10,
    ignore_reaction_ids = ["Htex", "H2Otex", "H2Otpp"],
)

mmdf.concentrations

#: Differentiate but take into account kinetics
pgm_kd = make_gecko_model(
    pmodel;
    reaction_isozymes,
    gene_product_bounds = gene_product_bounds_kd,
    gene_product_molar_mass,
    gene_product_mass_group_bound = gene_product_mass_group_bound_kd,
)

opt_model = flux_balance_analysis(
    pgm_kd,
    CPLEX.Optimizer;
    modifications = [
        change_optimizer_attribute("CPXPARAM_Emphasis_Numerical", 1),
        change_optimizer_attribute("CPX_PARAM_SCAIND", 1),
        COBREXA.silence,
    ],
)

rfluxes_kd = flux_dict(pgm_kd, opt_model)
gpconcs_kd = gene_product_dict(pgm_kd, opt_model)
flux_summary(rfluxes_kd)

rid_enzyme = Dict(
    rid => isozyme_to_enzyme(first(isozyme_vec), gene_product_molar_mass) for
    (rid, isozyme_vec) in reaction_isozymes
)

diffmodel = with_parameters(
    pgm_kd,
    rid_enzyme,
    rid_dg0,
    rid_km,
    mmdf.concentrations;
    scale_equality = true,
    scale_inequality = true,
    ignore_reaction_ids = ["Htex", "H2Otex", "H2Otpp"],
    ignore_metabolite_ids = ["h2o_c", "h2o_e", "h2o_p", "h_c", "h_e", "h_p"],
)

#!
ignore_reaction_ids = ["Htex", "H2Otex", "H2Otpp"]
ignore_metabolite_ids = ["h2o_c", "h2o_e", "h2o_p", "h_c", "h_e", "h_p"]

sats = Float64[]
rids = String[]
dgs = Float64[]
for mangled_rid in reactions(pgm_kd)
    rid = String(first(split(mangled_rid, "#")))
    sat = DifferentiableMetabolism._saturation(
        pgm_kd,
        rid_enzyme,
        rid_km,
        rid,
        mangled_rid,
        diffmodel.θ;
        ignore_reaction_ids,
        ignore_metabolite_ids,
    )
    dg = DifferentiableMetabolism._dg(
        pgm_kd,
        rid_enzyme,
        rid_dg0,
        rid,
        mangled_rid,
        diffmodel.θ;
        ignore_reaction_ids,
        ignore_metabolite_ids,
    )
    push!(sats, sat)
    push!(dgs, dg)
    push!(rids, rid)
end

idxs = sortperm(dgs .* sats)
x = [rids[idxs] sats[idxs] dgs[idxs] sats[idxs] .* dgs[idxs]]

#: check if solveable
x = zeros(length(diffmodel.var_ids));
ν = zeros(length(diffmodel.d(diffmodel.θ)));
λ = zeros(length(diffmodel.h(diffmodel.θ)));
modifications = [
    change_optimizer_attribute("CPXPARAM_Emphasis_Numerical", 1),
    change_optimizer_attribute("CPX_PARAM_SCAIND", 1),
    # COBREXA.silence,
];

DifferentiableMetabolism._solve_model!(x, ν, λ, diffmodel, CPLEX.Optimizer; modifications)

