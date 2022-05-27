using COBREXA, DifferentiableMetabolism, CPLEX, JSON, CairoMakie, ColorSchemes

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
gene_product_mass_group_bound = Dict("uncategorized" => 500_000.0)

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
bmid = "BIOMASS_Ec_iML1515_core_75p37M"
μ = rfluxes[bmid]

fva_mins, fva_maxs = flux_variability_analysis_dict(
    gm,
    CPLEX.Optimizer;
    bounds = objective_bounds(0.9),
    modifications,
)

model.reactions[bmid].lb = μ * 0.9

all_dxs = []
all_xs = []
for rid in keys(fva_maxs)
    for bnd in [fva_maxs[rid][rid], fva_mins[rid][rid]]
        try 
            dcmodel = deepcopy(model)
            if bnd < 0 
                dcmodel.reactions[rid].lb = bnd * 1.01
                dcmodel.reactions[rid].ub = bnd * 0.99    
            else
                dcmodel.reactions[rid].ub = bnd * 1.01
                dcmodel.reactions[rid].lb = bnd * 0.99
            end
            
            gm = make_gecko_model(
                dcmodel;
                reaction_isozymes,
                gene_product_bounds,
                gene_product_molar_mass,
                gene_product_mass_group_bound,
            )

            opt_model = flux_balance_analysis(gm, CPLEX.Optimizer; modifications)
            rfluxes = flux_dict(gm, opt_model)

            pmodel = prune_model(dcmodel, rfluxes)

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

            #: differentiate model
            rid_enzyme = Dict(
                rid => isozyme_to_enzyme(first(isozyme_vec), gene_product_molar_mass) for
                (rid, isozyme_vec) in reaction_isozyme
            )

            diffmodel = with_parameters(pgm, rid_enzyme)

            x, _dx = differentiate(diffmodel, CPLEX.Optimizer; modifications)
            dx = _dx[1:length(diffmodel.var_ids), :]
            push!(all_xs, [diffmodel.var_ids, x]) 
            push!(all_dxs, [diffmodel.var_ids, diffmodel.param_ids, dx])
        catch
            println("Error on reaction: ", rid)
        end
    end
end

using DataFrames, CSV, Serialization

serialize(joinpath("results", "gecko", "var_sens2.jls"), all_dxs)
serialize(joinpath("results", "gecko", "vars2.jls"), all_xs)


df = DataFrame(Reaction=String[], Parameter=String[], Sensitivity=Float64[])

for (vids, pids, dx) in all_dxs
    # (vids, pids, dx) = first(all_dxs)
    rs = String[]
    ps = String[]
    d = Float64[]
    for (i, vid) in enumerate(vids)
        v = first(split(vid, "#"))
        for (j, pid) in enumerate(pids)
            p = last(split(pid, "#"))
            if p == v  
                push!(rs, v)
                push!(ps, p)
                push!(d, dx[i, j])
            end
        end
    end
    append!(
        df,
        DataFrame(Reaction=rs, Parameter=ps, Sensitivity=d)
    )
end

CSV.write(joinpath("results", "gecko", "variability_sensitivities2.csv"), df)


# all_dxs = deserialize(joinpath("results", "gecko", "var_sens2.jls"))
# df = DataFrame(Reaction=String[], Parameter=String[], Sensitivity=Float64[])

# for (vids, pids, dx) in all_dxs
#     # (vids, pids, dx) = first(all_dxs)
#     rs = String[]
#     ps = String[]
#     d = Float64[]
#     for (i, vid) in enumerate(vids)
#         v = first(split(vid, "#"))
#         for (j, pid) in enumerate(pids)
#             p = last(split(pid, "#"))
#             if p == v  
#                 push!(rs, v)
#                 push!(ps, p)
#                 push!(d, dx[i, j])
#             end
#         end
#     end
#     append!(
#         df,
#         DataFrame(Reaction=rs, Parameter=ps, Sensitivity=d)
#     )
# end