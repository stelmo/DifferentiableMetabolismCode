using COBREXA, DifferentiableMetabolism, CPLEX, JSON

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

modifications = [
    change_optimizer_attribute("CPXPARAM_Emphasis_Numerical", 1),
    change_optimizer_attribute("CPX_PARAM_SCAIND", 1),
    change_optimizer_attribute("CPX_PARAM_THREADS", 1),
    COBREXA.silence,
]

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
gene_product_bounds(gid) = (0, 100_000.0)
gene_product_mass_group_bound = Dict("uncategorized" => 320_000.0)

#: run gecko to prune the model
gm = make_gecko_model(
    model;
    reaction_isozymes,
    gene_product_bounds,
    gene_product_molar_mass,
    gene_product_mass_group_bound,
)

opt_model = flux_balance_analysis(gm, CPLEX.Optimizer; modifications)
rfluxes = flux_dict(gm, opt_model)
gpconcs = gene_product_dict(gm, opt_model)
flux_summary(rfluxes)

pmodel = prune_model(model, rfluxes)

#: construct a pruned GeckoModel
reaction_isozyme = Dict{String,Vector{Isozyme}}()
for rid in filter(x -> haskey(reaction_isozymes, x), reactions(pmodel))
    reaction_isozyme[rid] = [
        argmax(
            smoment_isozyme_speed(x -> gene_product_molar_mass[x]),
            reaction_isozymes[rid],
        ),
    ]
end

pgm = make_gecko_model(
    pmodel;
    reaction_isozymes = reaction_isozyme,
    gene_product_bounds,
    gene_product_molar_mass,
    gene_product_mass_group_bound,
)

opt_model = flux_balance_analysis(pgm, CPLEX.Optimizer; modifications)
rfluxes = flux_dict(pgm, opt_model)
gpconcs = gene_product_dict(pgm, opt_model)
flux_summary(rfluxes)
mass_groups = gene_product_mass_group_dict(pgm, opt_model)

#: differentiate model
rid_enzyme = Dict(
    rid => isozyme_to_enzyme(first(isozyme_vec), gene_product_molar_mass) for
    (rid, isozyme_vec) in reaction_isozyme
)

diffmodel = with_parameters(pgm, rid_enzyme)

x, _dx = differentiate(diffmodel, CPLEX.Optimizer; modifications)
dx = _dx[1:length(diffmodel.var_ids), :]

sol = Dict(string(first(split(k, "#"))) => v for (k, v) in zip(diffmodel.var_ids, x))
flux_summary(sol)

using CSV, DataFrames

vids = diffmodel.var_ids 
pids = diffmodel.param_ids
rs = String[]
ps = String[]
d = Float64[]
for (i, vid) in enumerate(vids)
    v = first(split(vid, "#"))
    for (j, pid) in enumerate(pids)
        p = last(split(pid, "#")) 
        push!(rs, v)
        push!(ps, p)
        push!(d, dx[i, j])
    end
end
df = DataFrame(Reaction=rs, Parameter=ps, Sensitivity=d)

CSV.write(joinpath("results", "gecko", "aerobic_sensitivity.csv"), df)

#: Plot
using CairoMakie, ColorSchemes, Escher 

include(joinpath("analyses", "visualize.jl"))
using .Visualize

rids, gids, mids = Visualize.get_core_metabolism(pmodel, "data/maps/simplified_core_iml1515_map.json")
# filter!(!startswith("EX_"), rids)

#: plot map
fig = Figure(resolution = (1000, 1000), backgroundcolor = :transparent);

fluxes = Dict(k => v for (k, v) in sol if !startswith(k, "b"))
reaction_mass = Dict{String,Float64}()
for (k, v) in sol
    !startswith(k, "b") && continue
    for (rid, enzyme) in rid_enzyme
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
    joinpath("data", "maps", "simplified_core_iml1515_map2.json");
    rts = 8,
    mts = 8,
)
fig
CairoMakie.FileIO.save(joinpath("..", "DifferentiableMetabolismPaper", "docs", "imgs", "gecko_map_iml1515.pdf"), fig)

#: plot sensitivities

fig = Figure(resolution = (2000, 1000), backgroundcolor = :transparent);

ga = fig[1, 1] = GridLayout()
gb = fig[1, 2] = GridLayout()

heatmap_ax_fluxes = Axis(
    ga[1, 1],
    xticklabelrotation = -pi / 2,
    xlabel = "Reaction\n(Response variable)",
    ylabel = "Turnover number\n(Perturbed parameter)",
);

ylabs = Visualize.plot_heatmap(
    diffmodel,
    dx,
    heatmap_ax_fluxes,
    ga[1, 2];
    ignore_vids = filter(
        x -> first(split(x, "#")) ∉ rids || startswith(x, "b"),
        diffmodel.var_ids,
    ),
    ignore_pids = filter(x -> last(split(x, "#")) ∉ rids, diffmodel.param_ids),
    n_σ = 1.0,
    atol = 1e-6,
)

fig

heatmap_ax_genes = Axis(
    gb[1, 1],
    xticklabelrotation = -pi / 2,
    xlabel = "Gene product\n(Response variable)",
    xlabelpadding = 70,
);

Visualize.plot_heatmap(
    diffmodel,
    dx,
    heatmap_ax_genes,
    gb[1, 2];
    ignore_vids = filter(x -> x ∉ gids, diffmodel.var_ids),
    ignore_pids = filter(x -> last(split(x, "#")) ∉ rids, diffmodel.param_ids),
    n_σ = 1.0,
    atol = 1e-6,
)
fig

for (label, layout) in zip(["A", "B"], [ga, gb])
    Label(
        layout[1, 1, TopLeft()],
        label,
        textsize = 26,
        # font = noto_sans_bold,
        padding = (0, 5, 5, 0),
        halign = :right,
    )
end

fig

CairoMakie.FileIO.save(joinpath("..", "DifferentiableMetabolismPaper", "docs", "imgs", "gecko_sens_iml1515.pdf"), fig)

#: Plot proteome masses vs sensitivities
gid_mass = JSON.parsefile(joinpath("results", "proteomic_mass", "gid_mass.json"))
rid_mass = JSON.parsefile(joinpath("results", "proteomic_mass", "rid_mass.json"))

fig = Figure(backgroundcolor = :transparent);
ax = Axis(
    fig[1,1],
    yscale=log10,
    xscale=log10,
    ylabel="Mass fraction of enzyme in reaction (fg/cell)",
    xlabel="Scaled flux control coefficient of reaction",
)

rids = [(last(split(x, "#")), i) for (i, x) in enumerate(diffmodel.param_ids)]
for (k, v) in rid_mass
    v == 0 && delete!(rid_mass, k)
    # remove transporters, poor proteomic resolution
    if contains(model.reactions[k].name, "transport") || contains(model.reactions[k].name, "Transport") 
        delete!(rid_mass, k)
    end 
end
filter!(x -> haskey(rid_mass, x[1]), rids)
idxs = [x[2] for x in rids]
rids = [x[1] for x in rids]
rs = filter(!startswith("b"), first.(split.(diffmodel.var_ids, "#")))
nrxns = length(rs)
rlu = Dict(r => first(indexin([r], rs)) for r in rids)
cf = [dx[rlu[r], i] for (r, i) in zip(rids, idxs)]
ms = [rid_mass[r] for r in rids]

scatter!(ax, cf, ms)
hidexdecorations!(ax, ticks = false, ticklabels = false, label = false)
hideydecorations!(ax, ticks = false, ticklabels = false, label = false)
fig
CairoMakie.FileIO.save(joinpath("..", "DifferentiableMetabolismPaper", "docs", "imgs", "gecko_flux_control_vs_proteome.pdf"), fig)

using GLM, DataFrames
df = DataFrame(Y = log10.(ms), X = log10.(cf))
f = lm(@formula(Y ~ X), df)
r2(f)

#! compare to flux control coefficients
fcc = Dict()
for rid in keys(reaction_isozyme)

    original_isozyme = first(reaction_isozyme[rid])

    new_isozyme = Isozyme(
        original_isozyme.gene_product_count,
        1.001 * original_isozyme.kcat_forward,
        1.001 * original_isozyme.kcat_reverse,
    )

    reaction_isozyme[rid] = [new_isozyme]

    pgm = make_gecko_model(
        pmodel;
        reaction_isozymes = reaction_isozyme,
        gene_product_bounds,
        gene_product_molar_mass,
        gene_product_mass_group_bound,
    )

    opt_model = flux_balance_analysis(pgm, CPLEX.Optimizer)
    rfluxes = flux_dict(pgm, opt_model)
    gpconcs = gene_product_dict(pgm, opt_model)

    fcc[rid] = (rfluxes, gpconcs)

    reaction_isozyme[rid] = [original_isozyme]
end

pgm = make_gecko_model(
    pmodel;
    reaction_isozymes = reaction_isozyme,
    gene_product_bounds,
    gene_product_molar_mass,
    gene_product_mass_group_bound,
)

opt_model = flux_balance_analysis(pgm, CPLEX.Optimizer)
rfluxes = flux_dict(pgm, opt_model)
gpconcs = gene_product_dict(pgm, opt_model)

fcc_dx = similar(dx)
for (i, pid) in enumerate(diffmodel.param_ids)
    base_param_id = last(split(pid, "#"))
    f, g = fcc[base_param_id]
    vids = [first(split(x, "#")) for x in diffmodel.var_ids]
    fvec = [-1000.0 * (1.0 - f[k] / rfluxes[k]) for k in filter(!startswith("b"), vids)]
    gvec = [-1000.0 * (1.0 - g[k] / gpconcs[k]) for k in filter(startswith("b"), vids)]

    fcc_dx[:, i] .= [fvec; gvec]
end

fig = Figure(
    resolution = (1000, 1000),
    # backgroundcolor=:transparent,
);

ax = Axis(
    fig[1, 1],
    xlabel = "Finite difference flux control coefficient",
    ylabel = "Derivative flux control coefficient",
);

scatter!(ax, vec(fcc_dx), vec(dx))

hidexdecorations!(ax, ticks = false, ticklabels = false, label = false)
hideydecorations!(ax, ticks = false, ticklabels = false, label = false)
fig

CairoMakie.FileIO.save(
    joinpath("..", "DifferentiableMetabolismPaper", "docs", "imgs", "gecko_sens_vs_flux_control_coefficients.pdf"),
    fig,
)
