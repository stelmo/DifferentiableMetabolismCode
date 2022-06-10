module MetaboliteEstimates

using COBREXA
using JuMP #max driving force of metabolic reactions

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

end # module
