using COBREXA, JSON, Statistics
using eQuilibrator, Measurements, Unitful

# Conditions
pH_cytosol = 7.6
pH_noncytosol = 7.0
Δψ = -0.14u"V"
ionic_strength_noncytosol = 0.25u"M"
pMg_noncytosol = 3.0

# Cytosol environment
cytosol = eQuilibrator.Equilibrator(pH = pH_cytosol, ionic_strength = 250u"mM")

# Model - get all reactions and metabolites
model = load_model(StandardModel, joinpath("data", "models", "iML1515.json"))

# Find ΔG'⁰ all internal reactions
rid_dg0 = Dict{String,Float64}()
for rid in reactions(model)
    try
        println(rid)
        looks_like_biomass_reaction(rid) && continue
        looks_like_exchange_reaction(rid) && continue
        rid == "ATPM" && continue
        startswith(rid, "DM_") && continue

        r_stoich = reaction_stoichiometry(model, rid) # NB: this overwrites the model
        if all(endswith(k, "_c") for k in keys(r_stoich)) # internal reactions
            rid_dg0[rid] =
                ustrip(
                    u"kJ/mol",
                    standard_dg_prime(
                        cytosol,
                        bigg(stoichiometry_string(r_stoich; format_id = x -> x[1:end-2])),
                    ),
                ).val
        else # intercompartmental reactions

            non_cyto_compartment = first(
                filter(x -> x != "c", [last(split(mid, "_")) for mid in keys(r_stoich)]),
            )
            if haskey(r_stoich, "h_c") && haskey(r_stoich, "h_$non_cyto_compartment")
                r_stoich["h_c"] =
                    sign(r_stoich["h_c"]) * abs(r_stoich["h_$non_cyto_compartment"]) # number of protons that cross the membrane (assume no acid base reaction happens outside & true number)
            end
            if haskey(r_stoich, "h_p") && haskey(r_stoich, "h_e")
                r_stoich["h_p"] = sign(r_stoich["h_p"]) * abs(r_stoich["h_e"]) # number of protons that cross the membrane (assume no acid base reaction happens outside & true number)
            end

            cyto_str = Dict(k => v for (k, v) in r_stoich if endswith(k, "_c"))
            extra_str =
                Dict(k => v for (k, v) in r_stoich if endswith(k, "_$non_cyto_compartment"))

            cyto_rxn = bigg(stoichiometry_string(cyto_str; format_id = x -> x[1:end-2]))
            extra_rxn = bigg(stoichiometry_string(extra_str; format_id = x -> x[1:end-2]))

            dg = multi_compartmental_standard_dg_prime(
                cytosol,
                cyto_rxn,
                extra_rxn;
                potential_difference = Δψ,
                pH_outer = pH_noncytosol,
                pMg_outer = pMg_noncytosol,
                ionic_strength_outer = ionic_strength_noncytosol,
            )
            rid_dg0[rid] = ustrip(u"kJ/mol", dg).val
        end
    catch
        println("Error with ", rid)
    end
end

rid_dg0_pruned = Dict(k => v for (k, v) in rid_dg0 if abs(v) < 100)

open(joinpath("data", "thermodynamics", "iML1515_thermo.json"), "w") do io
    JSON.print(io, rid_dg0_pruned)
end
