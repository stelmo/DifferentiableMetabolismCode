using COBREXA, JSON, FASTX
using HTTP

model = load_model(StandardModel, joinpath("model_construction", "model_files", "iML1515.json"))

function get_reaction(model, bn)
    rxns = String[]
    for rid in reactions(model)
        gss = reaction_gene_association(model, rid)
        isnothing(gss) && continue
        for gs in gss
            for bnum in gs
                if String(bnum) == String(bn)
                    push!(rxns, rid)
                end
            end
        end
    end
    return rxns
end

rid_id_bid = Dict()
open(joinpath("data", "proteome_sequences", "ecoli_seqs.tab")) do io
    firstline = true
    for ln in eachline(io)
        firstline && (firstline = false; continue)
        prts = split(ln, "\t")
        id = prts[1]
        bnum = split(prts[3], " ")
        if !isempty(bnum) && startswith(first(bnum), "b") 
            bnum = String(first(bnum))
            rxns = get_reaction(model, bnum)
            isempty(rxns) && continue
            for rxn in rxns
                if haskey(rid_id_bid, rxn)
                    push!(rid_id_bid[rxn], (bnum, id))
                else
                    rid_id_bid[rxn] = [(bnum, id)]
                end
            end
        end
    end
end

for (k, v) in rid_id_bid
    rid_id_bid[k] = unique(rid_id_bid[k])
end


seqs = Dict{String,String}()
reader = FASTA.Reader(open(joinpath("data", "proteome_sequences", "ecoli_seqs.fasta"), "r"))
for record in reader
    id = split(identifier(record), "|")[2]
    seqs[id] = string(sequence(record))
end
close(reader)

# "succoa_c": -1.0, VNOYUJKHFWYWIR-ITIYDSSPSA-I
# "suchms_c": 1.0 GNISQJGXJIDKDJ-YFKPBYRVSA-M

mid_inchikey = Dict()
for (i, mid) in enumerate(metabolites(model))
    println(i, "out of ", n_metabolites(model))
    inchikeys = get(model.metabolites[mid].annotations, "inchi_key", nothing)
    isnothing(inchikeys) && continue
    inchikey = first(inchikeys)
    # println(inchikey)
    r = HTTP.request(
        "GET",
        "http://www.chemspider.com/InChI.asmx/InChIKeyToInChI?inchi_key=$inchikey",
    )
    b = String(r.body)
    _sidx = findfirst("InChI", b)
    isnothing(_sidx) && continue
    sidx = last(_sidx) + 2
    eidx = first(findlast("</string>", b)) - 1
    mid_inchikey[mid] = String(b[sidx:eidx])
end

list = []
num_rs_correct = 0
for (rid, v) in rid_id_bid
    temp_list = []
    for (bnum, id) in v
        enzyme = get(seqs, id, nothing)
        isnothing(enzyme) && continue
        rs = reaction_stoichiometry(model, rid)
        for mid in keys(rs)
            !haskey(mid_inchikey, mid) && break
            inchikey = mid_inchikey[mid]
            push!(temp_list, (rid, mid, inchikey, enzyme))
        end
        if length(rs) == length(temp_list)
            num_rs_correct += 1
            append!(list, temp_list)
        end
    end
end
println(num_rs_correct/n_reactions(model))

modeldata = Dict("rxnlist" => list)

open(joinpath("data", "proteome_sequences", "rid_mid_inchi_enzyme.json"), "w") do io
    JSON.print(io, modeldata)
end

#: smiles 
using CSV, DataFrames

df = DataFrame(CSV.File(joinpath("data", "kms", "kroll2021", "chebi.tsv")))
chebi_smile = Dict(zip(df[!, :id], df[!, :smiles]))

mid_smile = Dict()
for (i, mid) in enumerate(metabolites(model))
    ma = metabolite_annotations(model, mid)
    if haskey(ma, "chebi")
        for chebi in ma["chebi"]
            if haskey(chebi_smile, chebi)
                mid_smile[mid] = chebi_smile[chebi]
                break
            end
        end
    end
end

list = []
num_rs_correct = 0
for (rid, v) in rid_id_bid
    temp_list = []
    for (bnum, id) in v
        enzyme = get(seqs, id, nothing)
        isnothing(enzyme) && continue
        rs = reaction_stoichiometry(model, rid)
        for mid in keys(rs)
            !haskey(mid_smile, mid) && break
            smile = mid_smile[mid]
            push!(temp_list, (rid, mid, smile, enzyme))
        end
        if length(rs) == length(temp_list)
            num_rs_correct += 1
            append!(list, temp_list)
        end
    end
end
println(num_rs_correct/n_reactions(model))

modeldata = Dict("rxnlist" => list)

open(joinpath("data", "proteome_sequences", "rid_mid_smile_enzyme.json"), "w") do io
    JSON.print(io, modeldata)
end
