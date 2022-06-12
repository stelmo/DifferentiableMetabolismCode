using COBREXA, DifferentiableMetabolism
using CPLEX, JSON, SparseArrays, Statistics, Serialization, SparseArrays, LinearAlgebra
using DataFrames, DataFramesMeta, Chain, CSV, CairoMakie, ColorSchemes

#: now find better estimates
include(joinpath("analyses", "metabolite_estimates.jl"))
using .MetaboliteEstimates

function get_metabolite_sensitivities(ref_cond, kd_factor)
    try
        #: import model and data files
        base_load_path = joinpath("model_construction", "processed_models_files", "ecoli")
        model = load_model(StandardModel, joinpath(base_load_path, "iml1515_fixed.json"))
        reaction_kcats = JSON.parsefile(joinpath(base_load_path, "reaction_kcats.json"))
        reaction_protein_stoichiometry =
            JSON.parsefile(joinpath(base_load_path, "protein_stoichiometry.json"))
        gene_product_molar_mass = Dict(
            k => v for
            (k, v) in JSON.parsefile(joinpath(base_load_path, "protein_masses.json"))
        )

        model.reactions["EX_glc__D_e"].lb = -1000.0 # unconstrain because enzyme constraints take over
        model.reactions["ATPM"].lb = 0.0 # remove because only interested in derivatives
        model.reactions["GLUDy"].grr = [["b1761"]] # only this isozyme is expressed
        reaction_kcats["GLUDy"] = [first(reaction_kcats["GLUDy"])] 

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

        #: set protein bounds
        gene_product_bounds(gid) = (0, 200_000.0)
        gene_product_mass_group_bound = Dict("uncategorized" => 320_000.0)

        #: knockdown data
        knockdown_df = transform(
            DataFrame(CSV.File(joinpath("data", "crispr", "target_gene_rid.csv"))),
            :Reaction => ByRow(x -> rstrip.(split(x, "; "))) => :Reaction,
        )
        knockdown_df = flatten(knockdown_df, :Reaction)

        #: Get KD gene and reaction ids 
        df = @subset knockdown_df begin
            :Target .== ref_cond
        end
        target_reaction = first(String.(strip.(df[!, :Reaction])))
        target_gene = String(first(knockdown_df[knockdown_df[!, :Target].==ref_cond, :Gene]))

        #: Get a reference solution that uses the gene and associated reaction
        gm = make_gecko_model(
            model;
            reaction_isozymes,
            gene_product_bounds,
            gene_product_molar_mass,
            gene_product_mass_group_bound,
        )

        biomass_max = flux_balance_analysis(
            gm,
            CPLEX.Optimizer;
            modifications = [
                change_optimizer_attribute("CPXPARAM_Emphasis_Numerical", 1),
                change_optimizer_attribute("CPX_PARAM_SCAIND", 1),
                COBREXA.silence,
            ],
        )
        μ_ref = flux_dict(gm, biomass_max)["BIOMASS_Ec_iML1515_core_75p37M"]

        #: ensure gene is used, except if isozyme exists then just min total enzyme 
        modifications = [
                change_optimizer_attribute("CPXPARAM_Emphasis_Numerical", 1),
                change_optimizer_attribute("CPX_PARAM_SCAIND", 1),
                change_constraint("BIOMASS_Ec_iML1515_core_75p37M", lb = 0.9*μ_ref),
                COBREXA.silence,
            ]
        if target_gene ∉ ["b1854", "b1723", "b3403"] 
            push!(modifications, change_objective(target_gene; sense = COBREXA.MAX_SENSE))
        else
            push!(modifications, change_objective(target_gene; sense = COBREXA.MIN_SENSE))    
        end
        
        opt_model = flux_balance_analysis(
            gm,
            CPLEX.Optimizer;
            modifications=modifications,
        )
        tg_ref = gene_product_dict(gm, opt_model)[target_gene]
        # println(ref_cond, "\t", tg_ref)
        rfluxes_ref = flux_dict(gm, opt_model)

        pmodel = prune_model(model, rfluxes_ref)
      
        loopless_sol = flux_balance_analysis_dict(
            pmodel,
            CPLEX.Optimizer;
            modifications = [
                add_loopless_constraints(), 
                COBREXA.silence,
            ],
        )

        pmodel = prune_model(pmodel, loopless_sol)
        
        #: get dg0 and km constants
        rid_dg0, rid_km, km_concs = MetaboliteEstimates.get_parameters(pmodel, loopless_sol)
        
        #: get reasonable intracellular metabolite concentrations 
        minlims = MetaboliteEstimates.estimate_metabolite_concentrations(pmodel, rid_dg0, rid_km, km_concs, loopless_sol)
        # # ks = collect(keys(minlims.saturation))
        # # vs = collect(values(minlims.saturation))
        # # idxs = sortperm(vs)
        # # [ks[idxs] vs[idxs]]
        # # minimum(vs) < 0.01 && error("Too limited!")

        # #: get rid_enzyme
        for rid in keys(reaction_isozymes)
            reaction_isozymes[rid] = [
                    argmax(
                        smoment_isozyme_speed(x -> gene_product_molar_mass[x]),
                        reaction_isozymes[rid],
                    ),
                ]
        end
        rid_enzyme = Dict(
            rid => isozyme_to_enzyme(first(isozyme_vec), gene_product_molar_mass) for
            (rid, isozyme_vec) in reaction_isozymes
        )

        #: Simulate WT to get reference growth rate
        pgm_nokd = make_gecko_model(
            pmodel;
            reaction_isozymes,
            gene_product_bounds = gene_product_bounds,
            gene_product_molar_mass,
            gene_product_mass_group_bound = gene_product_mass_group_bound,
        )

        opt_model = flux_balance_analysis(
            pgm_nokd,
            CPLEX.Optimizer;
            modifications = [
                change_optimizer_attribute("CPXPARAM_Emphasis_Numerical", 1),
                change_optimizer_attribute("CPX_PARAM_SCAIND", 1),
                COBREXA.silence,
            ],
        )

        rfluxes_nokd = flux_dict(pgm_nokd, opt_model)
        gpconcs_nokd = gene_product_dict(pgm_nokd, opt_model)
        gpconcs_nokd[target_gene]
        flux_summary(rfluxes_nokd)

        diffmodel = with_parameters(
            pgm_nokd,
            rid_enzyme;
            rid_dg0,
            rid_km,
            mid_concentration = minlims.concentrations,
            scale_equality = true,
            scale_inequality = true,
            ignore_reaction_ids = ["Htex", "H2Otex", "H2Otpp"],
            ignore_metabolite_ids = ["h2o_c", "h2o_e", "h2o_p", "h_c", "h_e", "h_p"],
        )

        nvars = length(diffmodel.c(diffmodel.θ))
        update_Q!(diffmodel, x -> spdiagm(1e-9 .* ones(nvars))) # make sure can differentiate

        x = zeros(length(diffmodel.var_ids))
        ν = zeros(length(diffmodel.d(diffmodel.θ)))
        λ = zeros(length(diffmodel.h(diffmodel.θ)))
        DifferentiableMetabolism._solve_model!(
            x,
            ν,
            λ,
            diffmodel,
            CPLEX.Optimizer;
            modifications = [
                change_optimizer_attribute("CPXPARAM_Emphasis_Numerical", 1),
                change_optimizer_attribute("CPX_PARAM_SCAIND", 1),
                COBREXA.silence,
            ]
        )
        
        nokdrfs = Dict(fluxes(pgm_nokd) .=> reaction_flux(pgm_nokd)' * x)
        nokdgcs = Dict(diffmodel.var_ids[idx] => x[idx] for idx in findall(startswith("b"), diffmodel.var_ids))
        flux_summary(nokdrfs)
        nokdgcs[target_gene]
        # nokdrfs[target_reaction]

        #: limit enzyme 
        target_gene_prot_req = nokdgcs[target_gene] / kd_factor
        gene_product_bounds_kd(x) = x == target_gene ? (0.0, target_gene_prot_req) : (0.0, 200_000.0)
        gene_product_mass_group_bound_kd = Dict("uncategorized" => 320_000.0)

        #: Differentiate KD
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
        gpconcs_kd[target_gene]
        flux_summary(rfluxes_kd)

        diffmodel = with_parameters(
            pgm_kd,
            rid_enzyme;
            rid_dg0,
            rid_km,
            mid_concentration = minlims.concentrations,
            scale_equality = true,
            scale_inequality = true,
            ignore_reaction_ids = ["Htex", "H2Otex", "H2Otpp"],
            ignore_metabolite_ids = ["h2o_c", "h2o_e", "h2o_p", "h_c", "h_e", "h_p"],
        )

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
        kdrfs = Dict(fluxes(pgm_kd) .=> reaction_flux(pgm_kd)' * x)
        kdgcs = Dict(diffmodel.var_ids[idx] => x[idx] for idx in findall(startswith("b"), diffmodel.var_ids))
        flux_summary(kdrfs)
        # kdrfs[target_reaction]
        kdgcs[target_gene]

        #: get all affected reactions 
        affected_rids = String[]
        for rid in reactions(pmodel)
            grrs = reaction_gene_association(pmodel, rid)
            isnothing(grrs) && continue 
            for grr in grrs 
                if target_gene in grr 
                    push!(affected_rids, rid)
                    break
                end
            end
        end
        unique!(affected_rids)

        #: get data
        biomass_idx = findfirst(startswith("BIOMASS_"), diffmodel.var_ids)
        conc_idxs = findall(startswith("c#"), diffmodel.param_ids)
        conc_ids = diffmodel.param_ids[conc_idxs]

        conc_sens = dx[biomass_idx, conc_idxs]

        rmets = String[]
        rsubs = Float64[]
        rprods = Float64[]
        for affected_rid in affected_rids
            if abs(get(kdrfs, affected_rid, 0)) > 0
                for (k, v) in reaction_stoichiometry(pgm_kd, affected_rid)
                    k in rmets && continue
                    if v * kdrfs[affected_rid] > 0 # product 
                        push!(rmets, k)
                        push!(rsubs, 0.0)
                        push!(rprods, 1.0)
                    else # substrate
                        push!(rmets, k)
                        push!(rsubs, 1.0)
                        push!(rprods, 0.0)
                    end
                end
            end
        end
        met_idxs = [findfirst(x -> x == mid, conc_ids) for mid in "c#" .* rmets]

        big_idxs = filter(x -> x ∉ met_idxs, sortperm(abs.(conc_sens), rev = true))[1:25]

        mids = [conc_ids[met_idxs]; conc_ids[big_idxs]]
        sens = [conc_sens[met_idxs]; conc_sens[big_idxs]]
        subs = [rsubs; fill(0.0, length(big_idxs))]
        prods = [rprods; fill(0.0, length(big_idxs))]
        kd_growth = fill(kdrfs["BIOMASS_Ec_iML1515_core_75p37M"], length(subs))
        ref_growth = fill(nokdrfs["BIOMASS_Ec_iML1515_core_75p37M"], length(subs))
        df = DataFrame(
            Metabolite = mids,
            Sensitivity = sens,
            Substrate = subs,
            Product = prods,
            RefGrowth = ref_growth,
            KDGrowth = kd_growth,
        )
        CSV.write(
            joinpath("results", "crispr", ref_cond * ".csv"),
            df,
        )
    catch err
        println("Error on: ", ref_cond)
        println(err)
    end
end

targets = [
    "adk"
    "aroA"
    "carA"
    "cysH"
    "dxs"
    "eno"
    "fbaA"
    "gapA"
    "gdhA"
    "glmS"
    "gltA"
    "gnd"
    "ilvC"
    "metE"
    "pck"
    "pfkA"
    "pfkB"
    "pgi"
    "ppc"
    "prs"
    "ptsH"
    "purB"
    "purC"
    "pykA"
    "pykF"
    "pyrF"
    "sdhC"
    "tpiA"
    "zwf"
]

for target in targets
    get_metabolite_sensitivities(target, 5)
end
