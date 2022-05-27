using COBREXA, DifferentiableMetabolism
using CPLEX, JSON, SparseArrays, Statistics, Serialization, SparseArrays, LinearAlgebra
using DataFrames, DataFramesMeta, Chain, CSV, CairoMakie, ColorSchemes

#: now find better estimates
include(joinpath("analyses", "metabolite_estimates.jl"))
using .MetaboliteEstimates
using JuMP, KNITRO

function get_metabolite_sensitivities(ref_cond, biomass_relax=1.0)
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
            if rid == "GLUDy" # the slower homohexamer is expressed! 
                reaction_isozymes[rid] = [first(reaction_isozymes[rid])]
            else
                reaction_isozymes[rid] = [
                    argmax(
                        smoment_isozyme_speed(x -> gene_product_molar_mass[x]),
                        reaction_isozymes[rid],
                    ),
                ]
            end
        end

        gene_product_bounds(gid) = (0, 200_000.0)
        gene_product_mass_group_bound = Dict("uncategorized" => 650_000.0) # more than enough, minimized later

        #: knockdown data
        knockdown_df = transform(
            DataFrame(CSV.File(joinpath("data", "crispr", "target_gene.csv"))),
            :Reaction => ByRow(x -> rstrip.(split(x, "; "))) => :Reaction,
        )
        knockdown_df = flatten(knockdown_df, :Reaction)

        #: Get reference and KD ids 
        df = @subset knockdown_df begin
            :Target .== ref_cond
        end
        rid_kds = String.(strip.(df[!, :Reaction]))
        target_gene = String(first(knockdown_df[knockdown_df[!, :Target].==ref_cond, :Gene]))

        gm = make_gecko_model(
            model;
            reaction_isozymes,
            gene_product_bounds,
            gene_product_molar_mass,
            gene_product_mass_group_bound,
        )

        rfluxes_ref = flux_balance_analysis_dict(
            gm,
            CPLEX.Optimizer;
            modifications = [
                change_optimizer_attribute("CPXPARAM_Emphasis_Numerical", 1),
                change_optimizer_attribute("CPX_PARAM_SCAIND", 1),
                COBREXA.silence,
            ],
        )
        μ_ref = rfluxes_ref["BIOMASS_Ec_iML1515_core_75p37M"]

        opt_model = flux_balance_analysis(
            gm,
            CPLEX.Optimizer;
            modifications = [
                change_optimizer_attribute("CPXPARAM_Emphasis_Numerical", 1),
                change_optimizer_attribute("CPX_PARAM_SCAIND", 1),
                COBREXA.silence,
                change_objective(target_gene; sense = COBREXA.MAX_SENSE),
                change_constraint("BIOMASS_Ec_iML1515_core_75p37M"; lb = biomass_relax*μ_ref, ub = 1000),
            ],
        )
        rfluxes_ref = flux_dict(gm, opt_model)
        gpconcs_ref = gene_product_dict(gm, opt_model)
        tg_ref = gpconcs_ref[target_gene]
        flux_summary(rfluxes_ref)

        opt_model = flux_balance_analysis(
            gm,
            CPLEX.Optimizer;
            modifications = [
                change_optimizer_attribute("CPXPARAM_Emphasis_Numerical", 1),
                change_optimizer_attribute("CPX_PARAM_SCAIND", 1),
                COBREXA.silence,
                change_objective(genes(gm); sense = COBREXA.MIN_SENSE),
                change_constraint("BIOMASS_Ec_iML1515_core_75p37M"; lb = biomass_relax*μ_ref, ub = 1000),
                change_constraint(target_gene; lb = tg_ref, ub = 200_000),
            ],
        )
        rfluxes_ref = flux_dict(gm, opt_model)
        gpconcs_ref = gene_product_dict(gm, opt_model)
        flux_summary(rfluxes_ref)
        gpconcs_ref[target_gene]

        #: Get knockdown condition
        kd_factor = 5
        target_gene_prot_req = gpconcs_ref[target_gene] / kd_factor
        gene_product_bounds_kd(x) =
            x == target_gene ? (target_gene_prot_req*0.75, target_gene_prot_req*1.25) : (0.0, 200_000.0)
        gene_product_mass_group_bound_kd = Dict("uncategorized" => 650_000.0) # more than enough, protein limitation limits everything

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

        pmodel = prune_model(model, rfluxes_kd)

        loopless_mods = [
            add_loopless_constraints(), 
            COBREXA.silence,
        ]

        # ensure reaction of interest is not lost
        for rid_kd in rid_kds
            fff = rfluxes_kd[rid_kd]
            if fff > 1e-3
                push!(
                    loopless_mods,
                    change_constraint(rid_kd; lb=0.9*fff, ub=1000)
                )
            elseif fff < -1e-3
                push!(
                    loopless_mods,
                    change_constraint(rid_kd; lb=-1000, ub=0.9*fff)
                )
            end
        end
      
        loopless_sol = flux_balance_analysis_dict(
            pmodel,
            CPLEX.Optimizer;
            modifications = loopless_mods,
        )
        pmodel = prune_model(pmodel, loopless_sol)

        #: load km data
        _rid_km =
            JSON.parsefile(joinpath("data", "kms", "kroll2021", "kmdata_iml1515.json"))
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
            "ACOAD6f"
            "ANPRT"
            "CPPPGO"
            "ACOAD5f"
            "ACOAD3f"
            "DHORD2"
            "ACOAD7f"
            "AHCYSNS"
            "HACD2"
            "HACD4"
            "HACD7"
            "HACD6"
            "HACD1"
            "HACD5"
            "HACD3"
        ]
        delete!.(Ref(rid_km), del_list) # trial and error

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

        mmdf = MetaboliteEstimates.get_thermo_and_saturation_concentrations(
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
            modifications = [COBREXA.silence],
        )

        minlims = MetaboliteEstimates.min_limitation(
            pmodel,
            rid_dg0,
            rid_km,
            mmdf.concentrations,
            2.0,
            optimizer_with_attributes(KNITRO.Optimizer,"honorbnds" => 1);
            flux_solution = loopless_sol,
            proton_ids = ["h_c", "h_e", "h_p"],
            water_ids = ["h2o_c", "h2o_e", "h2o_p"],
            concentration_lb = mid -> get(km_concs, mid, 1e-6) / 1000,
            concentration_ub = mid -> get(km_concs, mid, 0.1 / 10) * 10,
            ignore_reaction_ids = ["Htex", "H2Otex", "H2Otpp"],
            modifications = [COBREXA.silence],
            ignore_metabolite_ids = ["h_c", "h_e", "h_p", "h2o_c", "h2o_e", "h2o_p"],
        )

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

        #: Plot 
        rfs = Dict(fluxes(pgm_kd) .=> reaction_flux(pgm_kd)' * x)
        flux_summary(rfs)

        biomass_idx = findfirst(startswith("BIOMASS_"), diffmodel.var_ids)
        conc_idxs = findall(startswith("c#"), diffmodel.param_ids)
        conc_ids = diffmodel.param_ids[conc_idxs]

        conc_sens = dx[biomass_idx, conc_idxs]

        rmets = String[]
        rsubs = Float64[]
        rprods = Float64[]
        for rid_kd in rid_kds
            if abs(get(rfs, rid_kd, 0)) > 0
                for (k, v) in reaction_stoichiometry(pgm_kd, rid_kd)
                    k in rmets && continue
                    if v * rfs[rid_kd] > 0 # product 
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
        df = DataFrame(
            Metabolite = mids,
            Sensitivity = sens,
            Substrate = subs,
            Product = prods,
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
    get_metabolite_sensitivities(target)
end


