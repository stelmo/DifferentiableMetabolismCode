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

gene_product_bounds(gid) = (0, 9_000.0)
gene_product_mass_group_bound = Dict("uncategorized" => 500_000.0) # more than enough, minimized later

#: knockdown data
knockdown_df = DataFrame(CSV.File(joinpath("data", "crispr", "target_gene.csv")))

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
    x == target_gene ? (0.0, gpconcs_ref[x] / kd_factor) : (0.0, 9_000.0)
gene_product_mass_group_bound_kd = Dict("uncategorized" => 500_000.0) # more than enough, protein limitation limits everything

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
    "ACACT7r"
    "MTHFR2"
]
delete!.(Ref(rid_km), del_list) # trial and error

fix_list = [
    "CS"
    "GND"
    "PPC"
    "DHAD1"
]

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
    # "Cell Envelope Biosynthesis" #!
    # "Cofactor and Prosthetic Group Biosynthesis" #!
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
    # Dict{String,Dict{String,Float64}}(),
    rid_km,
    mmdf.concentrations;
    scale_equality = true,
    scale_inequality = true,
    ignore_reaction_ids = ["Htex", "H2Otex", "H2Otpp"],
    ignore_metabolite_ids = ["h2o_c", "h2o_e", "h2o_p", "h_c", "h_e", "h_p"],
)

mat = diffmodel.E(diffmodel.θ)
extrema(mat)
check_scaling(diffmodel; verbose = true)

nvars = length(diffmodel.c(diffmodel.θ))
update_Q!(diffmodel, x -> spdiagm(1e-9 .* ones(nvars))) # make sure can differentiate

x, dx = differentiate(
    diffmodel,
    CPLEX.Optimizer;
    modifications = [
        change_optimizer_attribute("CPXPARAM_Emphasis_Numerical", 1),
        change_optimizer_attribute("CPX_PARAM_SCAIND", 1),
        COBREXA.silence,
    ],
    use_analytic = false,
    scale_output = true,
)

rfs = Dict(fluxes(pgm_kd) .=> reaction_flux(pgm_kd)' * x)
rfs["BIOMASS_Ec_iML1515_core_75p37M"]
rfs[rid_kd]

biomass_idx = findfirst(startswith("BIOMASS_"), diffmodel.var_ids)
conc_idxs = findall(startswith("c#"), diffmodel.param_ids)
conc_ids = diffmodel.param_ids[conc_idxs]

conc_sens = dx[biomass_idx, conc_idxs]
conc_sens .= abs.(conc_sens) .+ minimum(abs, conc_sens)
# conc_sens .= conc_sens ./ norm(conc_sens)

rs = reaction_stoichiometry(pgm_kd, rid_kd)
mids = "c#" .* keys(rs)
conc_xs = randn(length(conc_sens))

met_idxs = [findfirst(x -> x == mid, conc_ids) for mid in mids]
filter!(x -> conc_sens[x] > 0.0, met_idxs)

big_idxs = filter(x -> x ∉ met_idxs, sortperm(abs.(conc_sens), rev = true))[1:300]

[conc_ids[met_idxs] conc_sens[met_idxs]]
[conc_ids[big_idxs] conc_sens[big_idxs]]

fig = Figure();
ax = Axis(
    fig[1, 1],
    # yscale=Makie.pseudolog10,
    yscale = log10,
);

# scatter!(ax, conc_xs, conc_sens)

scatter!(ax, conc_xs[met_idxs], conc_sens[met_idxs], color = :red)
text!(
    ax,
    last.(split.(conc_ids[met_idxs], "#"));
    position = [Point2f(x, y) for (x, y) in zip(conc_xs[met_idxs], conc_sens[met_idxs])],
)

scatter!(ax, conc_xs[big_idxs], conc_sens[big_idxs], color = :green)
text!(
    ax,
    last.(split.(conc_ids[big_idxs], "#"));
    position = [Point2f(x, y) for (x, y) in zip(conc_xs[big_idxs], conc_sens[big_idxs])],
)

fig


#!
idx = findfirst(x -> x == "c#nadh_c", diffmodel.param_ids)
diffmodel.θ[idx] *= 2.0
x = zeros(length(diffmodel.var_ids))
ν = zeros(length(diffmodel.d(diffmodel.θ)))
λ = zeros(length(diffmodel.h(diffmodel.θ)))
modifications = [
    change_optimizer_attribute("CPXPARAM_Emphasis_Numerical", 1),
    change_optimizer_attribute("CPX_PARAM_SCAIND", 1),
    # COBREXA.silence,
]

DifferentiableMetabolism._solve_model!(x, ν, λ, diffmodel, CPLEX.Optimizer; modifications)
x[biomass_idx]
(x[biomass_idx] - rfs["BIOMASS_Ec_iML1515_core_75p37M"]) /
rfs["BIOMASS_Ec_iML1515_core_75p37M"]

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

#!
x = zeros(length(diffmodel.var_ids))
ν = zeros(length(diffmodel.d(diffmodel.θ)))
λ = zeros(length(diffmodel.h(diffmodel.θ)))
nA = length(x) + length(ν) + length(λ)
A = spzeros(nA, nA)
B = zeros(nA, length(diffmodel.param_ids))
dx = zeros(size(A, 1), size(B, 2))
modifications = [
    change_optimizer_attribute("CPXPARAM_Emphasis_Numerical", 1),
    change_optimizer_attribute("CPX_PARAM_SCAIND", 1),
    # COBREXA.silence,
]

DifferentiableMetabolism._solve_model!(x, ν, λ, diffmodel, CPLEX.Optimizer; modifications)

A .= diffmodel.analytic_var_derivs(x, ν, λ, diffmodel.θ)
B .= diffmodel.analytic_par_derivs(x, ν, λ, diffmodel.θ)


DifferentiableMetabolism._differentiate_kkt!(x, ν, λ, A, B, dx, diffmodel; use_analytic)

differentiate!(
    x,
    ν,
    λ,
    A,
    B,
    dx,
    diffmodel,
    CPLEX.Optimizer;
    use_analytic = true,
    modifications,
)

spB = sparse(B[1:length(diffmodel.var_ids), :])
spB[:, met_idxs[1]]

m = last(split(mids[1], "#"))
mmdf.concentrations[m]
idx = findfirst(x -> x == mids[1], diffmodel.param_ids)
diffmodel.θ[idx]
