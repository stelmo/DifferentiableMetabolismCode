module MetaboliteEstimates

using COBREXA
using JuMP #max driving force of metabolic reactions
using CPLEX, KNITRO, DifferentiableMetabolism
using JSON, Statistics, SparseArrays

function get_thermo_and_saturation_concentrations(
    model::MetabolicModel,
    reaction_standard_gibbs_free_energies::Dict{String,Float64},
    mid_concs::Dict{String,Float64},
    mmdf_lb::Float64,
    optimizer;
    flux_solution::Dict{String,Float64} = Dict{String,Float64}(),
    proton_ids::Vector{String} = ["h_c", "h_e"],
    water_ids::Vector{String} = ["h2o_c", "h2o_e"],
    constant_concentrations::Dict{String,Float64} = Dict{String,Float64}(),
    concentration_ratios::Dict{Tuple{String,String},Float64} = Dict{
        Tuple{String,String},
        Float64,
    }(),
    concentration_lb = x -> 1e-9,
    concentration_ub = x -> 100e-3,
    T = 298.15, # Kelvin
    R = 8.31446261815324e-3, # kJ/K/mol
    small_flux_tol = 1e-6,
    ignore_reaction_ids = [],
    modifications = [],
)
    opt_model = Model(optimizer)

    @variables opt_model begin
        mmdf
        logcs[1:n_metabolites(model)]
        dgrs[1:n_reactions(model)]
    end

    # set proton log concentration to zero so that it won't impact any calculations (biothermodynamics assumption)
    proton_idxs = Int.(indexin(proton_ids, metabolites(model)))
    for idx in proton_idxs
        JuMP.fix(logcs[idx], 0.0)
    end

    # set water concentrations to zero (aqueous condition assumptions)
    water_idxs = Int.(indexin(water_ids, metabolites(model)))
    for idx in water_idxs
        JuMP.fix(logcs[idx], 0.0)
    end

    # only consider reactions with supplied thermodynamic data AND have a flux bigger than
    # small_flux_tol => finds a thermodynamic profile that explains flux_solution
    active_rids = filter(
        rid ->
            haskey(reaction_standard_gibbs_free_energies, rid) &&
                abs(get(flux_solution, rid, small_flux_tol / 2)) > small_flux_tol &&
                !(rid in ignore_reaction_ids),
        reactions(model),
    )
    active_ridxs = Int.(indexin(active_rids, reactions(model)))

    # give dummy dG0 for reactions that don't have data
    dg0s =
        [get(reaction_standard_gibbs_free_energies, rid, 0.0) for rid in reactions(model)]

    S = stoichiometry(model)

    @constraint(opt_model, dgrs .== dg0s .+ (R * T) * S' * logcs)

    # thermodynamics should correspond to the fluxes
    flux_signs = [sign(get(flux_solution, rid, 1.0)) for rid in reactions(model)]

    # only constrain reactions that have thermo data
    @constraint(opt_model, dgrs[active_ridxs] .* flux_signs[active_ridxs] .<= 0)

    # add the absolute bounds
    missing_mets =
        [mid for mid in keys(constant_concentrations) if !(mid in metabolites(model))]
    !isempty(missing_mets) &&
        throw(DomainError(missing_mets, "metabolite(s) not found in model."))
    for (midx, mid) in enumerate(metabolites(model)) # idx in opt_model (missing ignore_metabolites)
        midx in water_idxs && continue
        midx in proton_idxs && continue
        if haskey(constant_concentrations, mid)
            JuMP.fix(logcs[midx], log(constant_concentrations[mid]))
        else
            # this metabolite needs bounds
            @constraint(
                opt_model,
                log(concentration_lb(mid)) <= logcs[midx] <= log(concentration_ub(mid))
            )
        end
    end

    # add the relative bounds
    for ((mid1, mid2), val) in concentration_ratios
        idxs = indexin([mid1, mid2], metabolites(model)) # TODO: this is not performant
        any(isnothing.(idxs)) &&
            throw(DomainError((mid1, mid2), "metabolite pair not found in model."))
        @constraint(opt_model, logcs[idxs[1]] == log(val) + logcs[idxs[2]])
    end

    @constraint(opt_model, mmdf_lb .<= -dgrs[active_ridxs] .* flux_signs[active_ridxs])

    mid_midx_lu = Dict(mid => midx for (midx, mid) in enumerate(metabolites(model)))

    @objective(
        opt_model,
        Min,
        sum((1 - logcs[mid_midx_lu[mid]] / log(c))^2 for (mid, c) in mid_concs)
    )

    # apply the modifications, if any
    for mod in modifications
        mod(model, opt_model)
    end

    optimize!(opt_model)
    is_solved(opt_model) || return nothing

    return (
        dg_reactions = Dict(
            rid => value(opt_model[:dgrs][i]) for (i, rid) in enumerate(reactions(model))
        ),
        concentrations = Dict(
            mid => exp(value(opt_model[:logcs][i])) for
            (i, mid) in enumerate(metabolites(model))
        ),
    )
end

function min_limitation(
    model::MetabolicModel,
    reaction_standard_gibbs_free_energies::Dict{String,Float64},
    reaction_michaelis::Dict{String,Dict{String, Float64}},
    mid_concs::Dict{String,Float64},
    mmdf_lb::Float64,
    optimizer;
    flux_solution::Dict{String,Float64} = Dict{String,Float64}(),
    proton_ids::Vector{String} = ["h_c", "h_e"],
    water_ids::Vector{String} = ["h2o_c", "h2o_e"],
    constant_concentrations::Dict{String,Float64} = Dict{String,Float64}(),
    concentration_ratios::Dict{Tuple{String,String},Float64} = Dict{
        Tuple{String,String},
        Float64,
    }(),
    concentration_lb = x -> 1e-9,
    concentration_ub = x -> 100e-3,
    T = 298.15, # Kelvin
    R = 8.31446261815324e-3, # kJ/K/mol
    ignore_reaction_ids = [],
    ignore_metabolite_ids = [],
    modifications = [],
)

    opt_model = Model(optimizer)

    @variables opt_model begin
        mmdf
        logcs[1:n_metabolites(model)]
        dgrs[1:n_reactions(model)]
        sats[1:n_reactions(model)]
        minlim
    end

    # set proton log concentration to zero so that it won't impact any calculations (biothermodynamics assumption)
    proton_idxs = Int.(indexin(proton_ids, metabolites(model)))
    for idx in proton_idxs
        JuMP.fix(logcs[idx], 0.0)
    end

    # set water concentrations to zero (aqueous condition assumptions)
    water_idxs = Int.(indexin(water_ids, metabolites(model)))
    for idx in water_idxs
        JuMP.fix(logcs[idx], 0.0)
    end

    # only consider reactions with supplied thermodynamic data AND have a flux bigger than
    # small_flux_tol => finds a thermodynamic profile that explains flux_solution
    thermo_active_rids = filter(
        rid -> haskey(reaction_standard_gibbs_free_energies, rid) && !(rid in ignore_reaction_ids), reactions(model),
    )
    thermo_active_ridxs = Int.(indexin(thermo_active_rids, reactions(model)))

    # give dummy dG0 for reactions that don't have data
    dg0s = [get(reaction_standard_gibbs_free_energies, rid, 0.0) for rid in reactions(model)]

    S = stoichiometry(model)

    @constraint(opt_model, dgrs .== dg0s .+ (R * T) * S' * logcs)

    # thermodynamics should correspond to the fluxes
    flux_signs = [sign(get(flux_solution, rid, 1.0)) for rid in reactions(model)]

    # only constrain reactions that have thermo data
    @constraint(opt_model, dgrs[thermo_active_ridxs] .* flux_signs[thermo_active_ridxs] .<= 0)

    # add the absolute bounds
    missing_mets =
        [mid for mid in keys(constant_concentrations) if !(mid in metabolites(model))]
    !isempty(missing_mets) &&
        throw(DomainError(missing_mets, "metabolite(s) not found in model."))
    for (midx, mid) in enumerate(metabolites(model)) # idx in opt_model (missing ignore_metabolites)
        midx in water_idxs && continue
        midx in proton_idxs && continue
        if haskey(constant_concentrations, mid)
            JuMP.fix(logcs[midx], log(constant_concentrations[mid]))
        else
            # this metabolite needs bounds
            @constraint(
                opt_model,
                log(concentration_lb(mid)) <= logcs[midx] <= log(concentration_ub(mid))
            )
        end
    end

    # add the relative bounds
    for ((mid1, mid2), val) in concentration_ratios
        idxs = indexin([mid1, mid2], metabolites(model)) # TODO: this is not performant
        any(isnothing.(idxs)) &&
            throw(DomainError((mid1, mid2), "metabolite pair not found in model."))
        @constraint(opt_model, logcs[idxs[1]] == log(val) + logcs[idxs[2]])
    end

    @constraint(opt_model, mmdf_lb .<= -dgrs[thermo_active_ridxs] .* flux_signs[thermo_active_ridxs])

    sat_active_rids = filter(
        rid -> haskey(reaction_michaelis, rid) && !(rid in ignore_reaction_ids), reactions(model),
    )

    mid_midx_lu = Dict(mid => midx for (midx, mid) in enumerate(metabolites(model)))

    for (ridx, rid) in enumerate(reactions(model))
        if rid ∉ sat_active_rids
            JuMP.fix(sats[ridx], 1.0)
            continue
        end
        (!haskey(reaction_michaelis, rid) || rid in ignore_reaction_ids) && continue
        _rs = reaction_stoichiometry(model, rid)
        rs = Dict(k => v for (k, v) in _rs if k ∉ ignore_metabolite_ids)
        stoich = values(rs)
        mids = collect(keys(rs))
        midxs = [mid_midx_lu[mid] for mid in mids]
        is_forward = flux_solution[rid] > 0 ? true : false
        
        s_term = @NLexpression(
            opt_model,
            prod(
                (exp(logcs[midx]) / reaction_michaelis[rid][mid])^(-1*nu) for
                (nu, midx, mid) in zip(stoich, midxs, mids) if nu < 0
            )
        )
        p_term = @NLexpression(
            opt_model,
            prod(
                (exp(logcs[midx]) / reaction_michaelis[rid][mid])^nu for
                (nu, midx, mid) in zip(stoich, midxs, mids) if nu > 0
            )
        )

        if is_forward
            @NLconstraint(
                opt_model, 
                sats[ridx] == s_term / (1.0 + s_term + p_term),
            )
        else
            @NLconstraint(
                opt_model, 
                sats[ridx] == p_term / (1.0 + s_term + p_term),
            )
        end
    end

    # overall constraints
    for i in eachindex(reactions(model))
        @NLconstraint(opt_model, minlim <= sats[i])
    end

    # start values 
    for (midx, mid) in enumerate(metabolites(model))
        set_start_value(logcs[midx], mid_concs[mid])
    end

    @objective(
        opt_model,
        Max,
        minlim
    )

    # apply the modifications, if any
    for mod in modifications
        mod(model, opt_model)
    end

    optimize!(opt_model)
    println("Min limitation = ", objective_value(opt_model))
    is_solved(opt_model) || return nothing

    return (
        dg_reactions = Dict(
            rid => value(opt_model[:dgrs][i]) for (i, rid) in enumerate(reactions(model))
        ),
        concentrations = Dict(
            mid => exp(value(opt_model[:logcs][i])) for
            (i, mid) in enumerate(metabolites(model))
        ),
        saturation = Dict(
            rid => value(opt_model[:sats][i]) for (i, rid) in enumerate(reactions(model))
        )
    )
end

function estimate_metabolite_concentrations(pmodel, rid_dg0, rid_km, km_concs, loopless_sol)

    mmdf = max_min_driving_force(
        pmodel,
        rid_dg0,
        CPLEX.Optimizer;
        flux_solution = loopless_sol,
        proton_ids = ["h_c", "h_e", "h_p"],
        water_ids = ["h2o_c", "h2o_e", "h2o_p"],
        concentration_lb = 1e-9,
        concentration_ub = 0.1,
        ignore_reaction_ids = ["Htex", "H2Otex", "H2Otpp"],
        modifications = [COBREXA.silence],
    )

    thermo_relax = 3
    kmconcs = MetaboliteEstimates.get_thermo_and_saturation_concentrations(
        pmodel,
        rid_dg0,
        km_concs,
        mmdf.mmdf/thermo_relax,
        CPLEX.Optimizer;
        flux_solution = loopless_sol,
        proton_ids = ["h_c", "h_e", "h_p"],
        water_ids = ["h2o_c", "h2o_e", "h2o_p"],
        concentration_lb = x -> 1e-9,
        concentration_ub = x -> 0.1,
        ignore_reaction_ids = ["Htex", "H2Otex", "H2Otpp"],
        modifications = [COBREXA.silence],
    )

    minlims = MetaboliteEstimates.min_limitation(
        pmodel,
        rid_dg0,
        rid_km,
        kmconcs.concentrations,
        mmdf.mmdf/thermo_relax,
        optimizer_with_attributes(KNITRO.Optimizer,"honorbnds" => 1);
        flux_solution = loopless_sol,
        proton_ids = ["h_c", "h_e", "h_p"],
        water_ids = ["h2o_c", "h2o_e", "h2o_p"],
        concentration_lb = mid -> get(kmconcs.concentrations, mid, 1e-9) / 10,
        concentration_ub = mid -> get(kmconcs.concentrations, mid, 0.1) * 10,
        ignore_reaction_ids = ["Htex", "H2Otex", "H2Otpp"],
        modifications = [COBREXA.silence],
        ignore_metabolite_ids = ["h_c", "h_e", "h_p", "h2o_c", "h2o_e", "h2o_p"],
    )

    return minlims 
end

function get_parameters(pmodel, loopless_sol)
    #: some kms make everything infeasible, remove these
    skip_kms = [
        "UPPDC1"
        "ACOAD6f"
        "ACOAD1fr"
        "ACOAD5f"
        "MTHFR2"
        "FADRx"
        "ACOAD3f"
    ]
    
    #: load km data
    _rid_km = JSON.parsefile(joinpath("data", "kms", "kroll2021", "kmdata_iml1515.json"))
    rid_km = Dict{String,Dict{String,Float64}}()
    for (rid, td) in _rid_km
        rid ∉ reactions(pmodel) && continue
        rid in skip_kms && continue
        isnothing(reaction_gene_association(pmodel, rid)) && continue
        isempty(reaction_gene_association(pmodel, rid)) && continue

        rid_km[rid] = Dict{String,Float64}()
        for (mid, km) in td
            rid_km[rid][mid] = km * 1e-3 # convert to Molar from mM
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
        if pmodel.reactions[rid].subsystem ∉ subsysts || rid ∉ reactions(pmodel)
            delete!(rid_dg0, rid)
        end
    end

    return rid_dg0, rid_km, km_concs
end

function simulate_wt(
    pmodel, 
    reaction_isozymes, 
    gene_product_bounds, 
    gene_product_molar_mass, 
    gene_product_mass_group_bound, 
    target_gene,
    target_reaction,
    rid_enzyme,
    rid_dg0,
    rid_km,
    mid_concentration,
)
    #: Simulate WT to get reference growth rate
    pgm_nokd = make_gecko_model(
        pmodel;
        reaction_isozymes,
        gene_product_bounds,
        gene_product_molar_mass,
        gene_product_mass_group_bound,
    )

    # opt_model = flux_balance_analysis(
    #     pgm_nokd,
    #     CPLEX.Optimizer;
    #     modifications = [
    #         change_optimizer_attribute("CPXPARAM_Emphasis_Numerical", 1),
    #         change_optimizer_attribute("CPX_PARAM_SCAIND", 1),
    #         COBREXA.silence,
    #     ],
    # )

    # rfluxes_nokd = flux_dict(pgm_nokd, opt_model)
    # gpconcs_nokd = gene_product_dict(pgm_nokd, opt_model)
    # gpconcs_nokd[target_gene]
    # flux_summary(rfluxes_nokd)

    diffmodel = with_parameters(
        pgm_nokd,
        rid_enzyme;
        rid_dg0,
        rid_km,
        mid_concentration,
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
    # flux_summary(nokdrfs)
    # nokdgcs[target_gene]
    # nokdrfs[target_reaction]
    return nokdgcs[target_gene], nokdrfs["BIOMASS_Ec_iML1515_core_75p37M"]
end

end # module
