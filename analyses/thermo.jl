module Thermo

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
    println(objective_value(opt_model))
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

end # module
