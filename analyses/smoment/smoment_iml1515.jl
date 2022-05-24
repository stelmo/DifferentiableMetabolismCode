using COBREXA, DifferentiableMetabolism, CPLEX, JSON, CairoMakie

#: import model and data files
base_load_path = joinpath("model_construction", "processed_models_files", "ecoli")
model = load_model(StandardModel, joinpath(base_load_path, "iml1515_fixed.json"))
reaction_kcats = JSON.parsefile(joinpath(base_load_path, "reaction_kcats.json"))
reaction_protein_stoichiometry =
    JSON.parsefile(joinpath(base_load_path, "protein_stoichiometry.json"))
gene_product_molar_mass = Dict(
    k => v for (k, v) in JSON.parsefile(joinpath(base_load_path, "protein_masses.json"))
)
gene_product_molar_mass_func(gid) = gene_product_molar_mass[gid]

model.reactions["EX_glc__D_e"].lb = -1000.0 # unconstrain because enzyme constraints take over

#: setup smoment data
reaction_isozyme = Dict{String,Isozyme}()
for rid in keys(reaction_kcats)
    isovec = [
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
    ]
    reaction_isozyme[rid] =
        argmax(smoment_isozyme_speed(gene_product_molar_mass_func), isovec)
end

total_enzyme_capacity = 500.0

#: run smoment to prune the model
smm = make_smoment_model(
    model;
    reaction_isozyme,
    gene_product_molar_mass,
    total_enzyme_capacity,
)

rfluxes = flux_balance_analysis_dict(smm, CPLEX.Optimizer)

pmodel = prune_model(model, rfluxes)

#: construct a pruned GeckoModel
pruned_reaction_isozyme = Dict{String,Isozyme}()
for rid in reactions(pmodel)
    !haskey(reaction_isozyme, rid) && continue
    isozyme = reaction_isozyme[rid]
    if abs(rfluxes[rid]) >= 1e-6
        pruned_reaction_isozyme[rid] = isozyme
    end
end

psmm = make_smoment_model(
    pmodel;
    reaction_isozyme = pruned_reaction_isozyme,
    gene_product_molar_mass,
    total_enzyme_capacity,
)

rfluxes = flux_balance_analysis_dict(smm, CPLEX.Optimizer)

#: differentiate model
rid_enzyme = Dict(
    rid => isozyme_to_enzyme(isozyme, gene_product_molar_mass) for
    (rid, isozyme) in pruned_reaction_isozyme
)

diffmodel = with_parameters(psmm, rid_enzyme)

x, dx = differentiate(diffmodel, CPLEX.Optimizer; use_analytic = false)

#: plot results
include(joinpath("analyses", "visualize.jl"))
using .Visualize

rids, gids, mids = Visualize.get_core_metabolism(pmodel)

#: Plot fluxes
fig = Figure(
    resolution = (1000, 1000),
    # backgroundcolor=:transparent,
);

heatmap_ax = Axis(fig[1, 1], xticklabelrotation = -pi / 2);

Visualize.plot_heatmap(
    diffmodel,
    dx,
    heatmap_ax,
    fig[1, 2];
    ignore_vids = filter(x -> first(split(x, "#")) ∉ rids, diffmodel.var_ids),
    ignore_pids = filter(x -> last(split(x, "#")) ∉ rids, diffmodel.param_ids),
    n_σ = 1.0,
    atol = 1e-6,
)
fig
