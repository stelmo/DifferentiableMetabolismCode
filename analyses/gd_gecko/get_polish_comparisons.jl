using JSON, ColorSchemes, CairoMakie, DataFrames, DataFramesMeta, Chain, CSV, Statistics
using COBREXA, JuMP, CPLEX, Gurobi
include(joinpath("analyses", "subgradient_descent.jl"))
include(joinpath("analyses", "data_constants.jl"))
using .SubgradientDescent
using .DataConstants

scale_factor = 1e-6

#: load kcat data
rdir = "linesearch_new"
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

#: load model and data files
base_load_path = joinpath("model_construction", "processed_models_files", "ecoli")
model =
    load_model(StandardModel, joinpath("model_construction", "model_files", "iML1515.json"))
reaction_kcats = JSON.parsefile(joinpath(base_load_path, "reaction_kcats.json"))
reaction_protein_stoichiometry =
    JSON.parsefile(joinpath(base_load_path, "protein_stoichiometry.json"))
gene_product_molar_mass = Dict(
    k => v for (k, v) in JSON.parsefile(joinpath(base_load_path, "protein_masses.json"))
)

model.reactions["EX_glc__D_e"].lb = -1000.0 # unconstrain because enzyme constraints take over

#: load measured data 
data_dir = joinpath("data", "kcats", "heckmann2020")
proteindata = JSON.parsefile(joinpath(data_dir, "proteindata.json"))
fluxdata = JSON.parsefile(joinpath(data_dir, "fluxdata.json"))

#: load condition specific data
master_kos = DataConstants.master_kos
protein_upper_bound = DataConstants.protein_upper_bound

#: setup gecko data 
ref_reaction_isozymes = Dict(
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
gene_product_bounds(gid) = (0, 99_000.0)

#: make opt func 
function run_model(
    optimizer,
    model,
    measured_fluxes,
    measured_proteins;
    scale_factor = 1e-6,
    modifications = [],
)

    m, n = size(stoichiometry(model))
    xl, xu = bounds(model)

    n_rxns = n_reactions(model) - n_genes(model)
    R = reaction_flux(model)[1:n_rxns, :]
    F = fluxes(model)
    fvals = [first(get(measured_fluxes, x, [missing])) for x in F]
    fvals = [ismissing(x) || abs(x) < 1e-3 ? missing : x for x in fvals]

    G = genes(model)
    pvals = [scale_factor * get(measured_proteins, x, missing) for x in G]
    pvals = [ismissing(x) || abs(x) < 1e-3 ? missing : x for x in pvals]
    mfact =
        length(pvals) + length(fvals) - count(ismissing.(pvals)) - count(ismissing.(fvals))
    opt_model = JuMP.Model(optimizer)
    x = @variable(opt_model, x[i = 1:n])

    Rx = R' * x[1:n_rxns]
    #: means squared relative deviation
    @objective(
        opt_model,
        Min,
        1 / mfact * (
            sum(
                (1 - Rx[i] / xhat)^2 for (i, xhat) in enumerate(fvals) if !ismissing(xhat)
            ) + sum(
                (1 - x[i+n_rxns] / xhat)^2 for
                (i, xhat) in enumerate(pvals) if !ismissing(xhat)
            )
        )
    )

    @constraint(opt_model, mb, stoichiometry(model) * x .== balance(model)) # mass balance
    @constraint(opt_model, lbs, xl .<= x) # lower bounds
    @constraint(opt_model, ubs, x .<= xu) # upper bounds

    C = coupling(model) # empty if no coupling
    cl, cu = coupling_bounds(model)
    @constraint(opt_model, c_lbs, cl .<= C * x) # coupling lower bounds
    @constraint(opt_model, c_ubs, C * x .<= cu) # coupling upper bounds

    for mod in modifications
        mod(nothing, opt_model)
    end

    optimize!(opt_model)

    if termination_status(opt_model) ∉ [JuMP.OPTIMAL, JuMP.LOCALLY_SOLVED]
        throw(DomainError(termination_status(opt_model), " model not solved optimally!"))
    end

    return objective_value(opt_model)
end


#! condition specific analysis 
master_ids = unique(kmax_df[!, :Condition])
ref_losses = Float64[]
polish_losses = Float64[]

for master_id in master_ids
    condid = first(split(master_id, "#"))

    holdout_df = @subset kmax_df @byrow begin
        !startswith(:Condition, condid)
    end

    holdout_df = transform(
        combine(groupby(holdout_df, :KcatID), :Kmax => maximum => :Kmax),
        :KcatID => x -> last.(split.(x, "#")),
    )
    @select!(holdout_df, :Kmax, :KcatID_function)
    rename!(holdout_df, Dict("KcatID_function" => :KcatID))
    polished_rid_kcat = Dict(holdout_df[!, :KcatID] .=> holdout_df[!, :Kmax])

    condition_reaction_isozymes = Dict(
        rid => [
            Isozyme(
                Dict(
                    k => v for (k, v) in zip(
                        reaction_gene_association(model, rid)[i],
                        reaction_protein_stoichiometry[rid][i],
                    )
                ),
                get(polished_rid_kcat, rid, reaction_kcats[rid][i][1]),
                get(polished_rid_kcat, rid, reaction_kcats[rid][i][2]),
            ) for i = 1:length(reaction_kcats[rid])
        ] for rid in keys(reaction_kcats)
    )

    cond_model = deepcopy(model)

    ko_id = first(split(master_id, "#"))
    if !isempty(master_kos[ko_id])
        delete!.(Ref(cond_model.reactions), master_kos[ko_id])
        delete!.(Ref(ref_reaction_isozymes), master_kos[ko_id])
        delete!.(Ref(condition_reaction_isozymes), master_kos[ko_id])
    end

    gene_product_mass_group_bound = Dict("uncategorized" => protein_upper_bound[master_id])
    # gene_product_mass_group_bound = Dict("uncategorized" => 320_000.0) # use the constraint they use in paper! 

    #: modifications for solver
    # modifications = [
    #     change_optimizer_attribute("CPXPARAM_Emphasis_Numerical", 1),
    #     change_optimizer_attribute("CPX_PARAM_SCAIND", 1),
    #     change_optimizer_attribute("CPX_PARAM_THREADS", 1),
    #     # COBREXA.silence,
    # ]

    modifications = [change_optimizer_attribute("NumericFocus", 3), COBREXA.silence]

    #: run reference gecko
    ref_gm = make_gecko_model(
        cond_model;
        reaction_isozymes = ref_reaction_isozymes,
        gene_product_bounds,
        gene_product_molar_mass,
        gene_product_mass_group_bound,
    )
    ref_loss = run_model(
        Gurobi.Optimizer,
        ref_gm,
        fluxdata[master_id],
        proteindata[master_id];
        scale_factor,
        modifications,
    )
    push!(ref_losses, ref_loss)

    #: run polished model
    polish_gm = make_gecko_model(
        cond_model;
        reaction_isozymes = condition_reaction_isozymes,
        gene_product_bounds,
        gene_product_molar_mass,
        gene_product_mass_group_bound,
    )
    polish_loss = run_model(
        Gurobi.Optimizer,
        polish_gm,
        fluxdata[master_id],
        proteindata[master_id];
        scale_factor,
        modifications,
    )
    push!(polish_losses, polish_loss)
end

df = DataFrame(Condition = master_ids, Rloss = ref_losses, Ploss = polish_losses)
df = @transform df @byrow begin
    :Condition = :Condition[1:findfirst(isnumeric, :Condition)-1]
end

CSV.write(
    joinpath("results", "gd_gecko", "polish_df.csv"),
    df,
)
