using COBREXA, DifferentiableMetabolism, CPLEX, JSON, CairoMakie

#: import model and data files
base_load_path = joinpath("model_construction", "processed_models_files", "ecoli_core")
model = load_model(StandardModel, joinpath(base_load_path, "e_coli_core_fixed.json"))
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
gene_product_bounds(gid) = (0, 100.0)
gene_product_mass_group_bound = Dict("uncategorized" => 500.0)

#: run gecko to prune the model
gm = make_gecko_model(
    model;
    reaction_isozymes,
    gene_product_bounds,
    gene_product_molar_mass,
    gene_product_mass_group_bound,
)

opt_model = flux_balance_analysis(gm, CPLEX.Optimizer)
fluxes = flux_dict(gm, opt_model)
gpconcs = gene_product_dict(gm, opt_model)


fig = Figure(resolution = (1000, 1000), backgroundcolor = :transparent);

reaction_mass = Dict{String,Float64}()
for (k, v) in fluxes
    !startswith(k, "b") && continue
    for (rid, enzyme) in reaction_isozymes
        if k in keys(enzyme.gene_product_count)
            reaction_mass[rid] = get(reaction_mass, rid, 0.0) + v
        end
    end
end

Visualize.plot_escher_viz(
    Axis(fig[1, 1]),
    fig[1, 2],
    rids,
    fluxes,
    reaction_mass,
    joinpath("data", "maps", "simplified_core_iml1515_map.json"),
)
fig
CairoMakie.FileIO.save(joinpath("..", "DifferentiableMetabolismPaper", "docs", "imgs", "gecko_map_iml1515.pdf"), fig)

# rfluxes = flux_dict(gm, opt_model)
# gpconcs = gene_product_dict(gm, opt_model)

# pmodel = prune_model(model, rfluxes)

# #: construct a pruned GeckoModel
# reaction_isozyme = Dict{String,Vector{Isozyme}}()
# for rid in reactions(pmodel)
#     for isozyme in get(reaction_isozymes, rid, [])
#         if all([get(gpconcs, gid, 0.0) for gid in keys(isozyme.gene_product_count)] .> 0)
#             reaction_isozyme[rid] = [isozyme]
#             break
#         end
#     end
# end

# pgm = make_gecko_model(
#     pmodel;
#     reaction_isozymes = reaction_isozyme,
#     gene_product_bounds,
#     gene_product_molar_mass,
#     gene_product_mass_group_bound,
# )

# opt_model = flux_balance_analysis(pgm, CPLEX.Optimizer)
# rfluxes = flux_dict(pgm, opt_model)
# gpconcs = gene_product_dict(pgm, opt_model)

# #: differentiate model
# rid_enzyme = Dict(
#     rid => isozyme_to_enzyme(first(isozyme_vec), gene_product_molar_mass) for
#     (rid, isozyme_vec) in reaction_isozyme
# )

# diffmodel = with_parameters(pgm, rid_enzyme)

# x, _dx = differentiate(diffmodel, CPLEX.Optimizer)
# dx = _dx[1:length(diffmodel.var_ids), :]

# #: plot results
# include(joinpath("analyses", "visualize.jl"))
# using .Visualize

# #: Plot fluxes
# fig = Figure(
#     resolution = (1000, 1000),
#     # backgroundcolor=:transparent,
# );

# heatmap_ax = Axis(fig[1, 1], xticklabelrotation = -pi / 2);

# Visualize.plot_heatmap(
#     diffmodel,
#     dx,
#     heatmap_ax,
#     fig[1, 2];
#     ignore_vids = filter(startswith("b"), diffmodel.var_ids),
#     ignore_pids = [],
#     n_σ = 1.0,
#     atol = 1e-6,
# )
# fig

# #: Plot gene products
# fig = Figure(
#     resolution = (1000, 1000),
#     # backgroundcolor=:transparent,
# );

# heatmap_ax = Axis(fig[1, 1], xticklabelrotation = -pi / 2);

# Visualize.plot_heatmap(
#     diffmodel,
#     dx,
#     heatmap_ax,
#     fig[1, 2];
#     ignore_vids = filter(!startswith("b"), diffmodel.var_ids),
#     ignore_pids = [],
#     n_σ = 1.0,
#     atol = 1e-6,
# )
# fig

# #! compare to flux control coefficients
# fcc = Dict()
# for rid in keys(reaction_isozyme)

#     original_isozyme = first(reaction_isozyme[rid])

#     new_isozyme = Isozyme(
#         original_isozyme.gene_product_count,
#         1.001 * original_isozyme.kcat_forward,
#         1.001 * original_isozyme.kcat_reverse,
#     )

#     reaction_isozyme[rid] = [new_isozyme]

#     pgm = make_gecko_model(
#         pmodel;
#         reaction_isozymes = reaction_isozyme,
#         gene_product_bounds,
#         gene_product_molar_mass,
#         gene_product_mass_group_bound,
#     )

#     opt_model = flux_balance_analysis(pgm, CPLEX.Optimizer)
#     rfluxes = flux_dict(pgm, opt_model)
#     gpconcs = gene_product_dict(pgm, opt_model)

#     fcc[rid] = (rfluxes, gpconcs)

#     reaction_isozyme[rid] = [original_isozyme]
# end

# pgm = make_gecko_model(
#     pmodel;
#     reaction_isozymes = reaction_isozyme,
#     gene_product_bounds,
#     gene_product_molar_mass,
#     gene_product_mass_group_bound,
# )

# opt_model = flux_balance_analysis(pgm, CPLEX.Optimizer)
# rfluxes = flux_dict(pgm, opt_model)
# gpconcs = gene_product_dict(pgm, opt_model)

# fcc_dx = similar(dx)
# for (i, pid) in enumerate(diffmodel.param_ids)
#     base_param_id = last(split(pid, "#"))
#     f, g = fcc[base_param_id]
#     vids = [first(split(x, "#")) for x in diffmodel.var_ids]
#     fvec = [-1000.0 * (1.0 - f[k] / rfluxes[k]) for k in filter(!startswith("b"), vids)]
#     gvec = [-1000.0 * (1.0 - g[k] / gpconcs[k]) for k in filter(startswith("b"), vids)]

#     fcc_dx[:, i] .= [fvec; gvec]
# end

# fig = Figure(
#     resolution = (1000, 1000),
#     # backgroundcolor=:transparent,
# );

# ax = Axis(fig[1, 1], xlabel = "FCC", ylabel = "MCA");

# scatter!(ax, vec(fcc_dx), vec(dx))

# hidexdecorations!(ax, ticks = false, ticklabels = false, label = false)
# hideydecorations!(ax, ticks = false, ticklabels = false, label = false)
# fig
