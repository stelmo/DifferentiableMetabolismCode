using COBREXA, JSON, FASTX
using HTTP

model =
    load_model(StandardModel, joinpath("model_construction", "model_files", "iML1515.json"))

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
        if !isempty(bnum) && startswith(first(bnum), "b") && prts[4] != ""
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

list = []
for (rid, v) in rid_id_bid
    for (bnum, id) in v
        enzyme = get(seqs, id, nothing)
        isnothing(enzyme) && continue
        rs = reaction_stoichiometry(model, rid)
        for mid in keys(rs)
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

            push!(list, (rid, mid, String(b[sidx:eidx]), enzyme))
        end
    end
end

modeldata = Dict("rxnlist" => list)

open(joinpath("data", "proteome_sequences", "rid_mid_inchi_enzyme.json"), "w") do io
    JSON.print(io, modeldata)
end
