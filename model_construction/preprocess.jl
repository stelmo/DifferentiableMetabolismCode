module Preprocess

using COBREXA, JSON, Statistics

"""
Check if reaction has a gene reaction rule assigned to it.
"""
function has_reaction_grr(model, rid)
    grr = reaction_gene_association(model, rid)
    if isnothing(grr) || isempty(grr) || isnothing(first(grr)) || isempty(first(grr))
        return false
    end
    return true
end

"""
If reaction spontaneous, remove all grrs for simplicity  
"""
function rm_spontaneous!(model; spont_rxn_genes = ["s0001"])
    for rid in reactions(model)
        if has_reaction_grr(model, rid)
            grrs = reaction_gene_association(model, rid)
            if any(any(in.(spont_rxn_genes, Ref(grr))) for grr in grrs)
                model.reactions[rid].grr = nothing
            end
        end
    end
    delete!.(Ref(model.genes), spont_rxn_genes)
    return nothing
end

"""
Add kcats to 
1) transporters if missing and has grr
2) physiological reactions without a kcat get an average one if grr present
"""
function add_kcats!(
    model,
    kcat_data;
    transporter_kcat = 65.0,
    average_kcat = 25.0,
    block_brackets = false,
)
    for rid in reactions(model)
        !has_reaction_grr(model, rid) && continue
        rxn = reaction_stoichiometry(model, rid)
        if block_brackets
            suffixes = [split(met, "[")[end] for met in keys(rxn)]
        else
            suffixes = [split(met, "_")[end] for met in keys(rxn)]
        end
        if length(unique(suffixes)) > 1
            if !haskey(kcat_data, rid)
                println("Transporter ", rid, " assigned kcat.")
                kcat_data[rid] = transporter_kcat
            end
        end
    end

    for rid in reactions(model)
        !has_reaction_grr(model, rid) && continue
        if !haskey(kcat_data, rid)
            kcat_data[rid] = average_kcat
            println("Assigned average kcat to: ", rid)
        elseif haskey(kcat_data, rid)
            continue
        else
            println("Rid not assigned a kcat: ", rid)
        end
    end

    return nothing
end

"""
Save all data for future use
"""
function save_all(
    foldername,
    modelname,
    model,
    protein_masses,
    reaction_kcats,
    protein_stoichiometry,
)
    save_loc = joinpath("model_construction", "processed_models_files", foldername)
    for f in readdir(save_loc)
        rm(joinpath(save_loc, f))
    end

    save_model(model, joinpath(save_loc, "$modelname.json"))

    open(joinpath(save_loc, "protein_masses.json"), "w") do io
        JSON.print(io, protein_masses)
    end

    open(joinpath(save_loc, "reaction_kcats.json"), "w") do io
        JSON.print(io, reaction_kcats)
    end

    open(joinpath(save_loc, "protein_stoichiometry.json"), "w") do io
        JSON.print(io, protein_stoichiometry)
    end

    return nothing
end

"""
Make reaction => kcat map    
"""
function get_reaction_kcats(model, kcat_data, scale)
    reaction_kcats = Dict()
    for rid in reactions(model)
        if has_reaction_grr(model, rid) # all rxns with grr should be assigned kcat by now
            reaction_kcats[rid] = [
                fill(scale * kcat_data[rid], 2) for
                _ in reaction_gene_association(model, rid)
            ]
        end
    end
    return reaction_kcats
end

"""
Make gid => molar mass map
"""
function get_protein_masses(model, proteome_data; replace_dash = false)
    protein_masses = Dict()
    avg_mass = mean([x[1] for x in values(proteome_data)])
    for gid in genes(model)
        if replace_dash
            gid_id = replace(gid, "-" => "")
        else
            gid_id = gid
        end
        if haskey(proteome_data, gid)
            protein_masses[gid_id] = proteome_data[gid][1]
        else
            protein_masses[gid_id] = avg_mass
        end
    end
    return protein_masses
end

"""
This must be done first because it changes the GRR structure of the model as well.
In a set of complexes, it gets rid of low confidence complexes if any complex in the 
set can be found in the ComplexPortal.
"""
function get_protein_stoichiometry!(model, proteome_data, organism)

    #: protein stoich map, infer from uniprot
    mer_map = [
        "Homotetramer" 4
        "Homodimer" 2
        "Homotrimer" 3
        "Homohexamer" 6
        "Homopentamer" 5
        "Homodecamer" 10
        "Homooctamer" 8
        "Homoheptamer" 7
        "Homododecamer" 12
        "Homomonomer" 1
        "Monomer" 1
    ]

    #: infer protein stoichiometry from uniprot annotations
    protein_stoichiometry = Dict()
    for rid in reactions(model)
        if has_reaction_grr(model, rid)
            grrs = reaction_gene_association(model, rid)
            smer = []
            for grr in grrs
                if length(grr) == 1 # only assign homomers
                    gid = grr[1]
                    if haskey(proteome_data, gid)
                        mer = proteome_data[gid][2]
                        if any(startswith.(Ref(mer), mer_map[:, 1]))
                            idx = findfirst(x -> startswith(mer, x), mer_map[:, 1])
                            push!(smer, [mer_map[idx, 2]])
                        else
                            push!(smer, [1.0])
                        end
                    else # no data
                        push!(smer, [1.0])
                    end
                else # assume complexes have uni-stoichiometry, manually fix later
                    push!(smer, fill(1.0, length(grr)))
                end
            end
            isempty(smer[1]) && continue
            protein_stoichiometry[rid] = smer
        end
    end

    #:fix complex stoichiometry, use ComplexPortal database
    complex = JSON.parsefile(joinpath("data", "proteome", organism * "_complex.json"))
    for rid in reactions(model)
        !has_reaction_grr(model, rid) && continue
        grrs = reaction_gene_association(model, rid)
        length(first(grrs)) == 1 && continue # skip monomers

        accurate_complex_idxs = Int[]
        for (i, grr) in enumerate(grrs)
            stoichs = []
            for v in values(complex)
                if length(intersect(collect(keys(v)), grr)) == length(grr) &&
                   length(intersect(collect(keys(v)), grr)) == length(v)
                    push!(stoichs, v)
                end
            end
            if length(stoichs) > 1
                @warn("Uncertain which complex to choose for reaction", rid)
            elseif length(stoichs) == 1
                stoich = first(stoichs)
                protein_stoichiometry[rid][i] = [
                    get(stoich, gid, 1.0) == 0 ? 1.0 : get(stoich, gid, 1.0) for gid in grr
                ]
                push!(accurate_complex_idxs, i)
            end
        end

        #: remove grrs for complexes not found in the database (conservative)
        if !isempty(accurate_complex_idxs)
            rem_idxs = filter(x -> x ∉ accurate_complex_idxs, 1:length(grrs))
            deleteat!(model.reactions[rid].grr, rem_idxs)
            deleteat!(protein_stoichiometry[rid], rem_idxs)
        end
    end

    return protein_stoichiometry
end

end
