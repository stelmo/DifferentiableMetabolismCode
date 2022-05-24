using COBREXA, Dates, ReadableRegex

function download_models()
    base_model_loc = joinpath("model_construction", "model_files")

    for f in readdir(base_model_loc)
        rm(joinpath(base_model_loc, f))
    end

    yeast_model_url = "https://raw.githubusercontent.com/SysBioChalmers/yeast-GEM/main/model/yeast-GEM.mat"
    yeast_model_path = joinpath(base_model_loc, "yeast-GEM.mat")
    download(yeast_model_url, yeast_model_path)

    model = load_model(yeast_model_path)
    #: need to fix this model's identifiers
    grrs = String[]
    regex = exactly(1, "x(") * maybe(one_or_more(DIGIT)) * exactly(1, ")")
    for grr in model.mat["rules"][:]
        for m in eachmatch(regex, grr)
            gidx = parse(Int64, foldl(replace, ["x(" => "", ")" => ""]; init = m.match))
            grr = foldl(replace, [m.match => genes(model)[gidx], "-" => ""]; init = grr)
        end
        push!(grrs, grr)
    end
    model.mat["grRules"] = grrs

    gs = [replace(g, "-" => "") for g in model.mat["genes"][:]]

    model.mat["genes"] = gs

    rid_names = Dict(
        model.mat["rxns"] .=> [
            string(a, " (", b, ")") for
            (a, b) in zip(model.mat["rxnNames"], model.mat["rxnBiGGID"])
        ],
    )
    mid_names = Dict(
        model.mat["mets"] .=> [
            string(a, " (", b, ")") for
            (a, b) in zip(model.mat["metNames"], model.mat["metBiGGID"])
        ],
    )

    jmodel = convert(JSONModel, model)
    #: fix various things
    for rxn in jmodel.json["reactions"]
        rxn["name"] = rid_names[rxn["id"]]
    end
    for met in jmodel.json["metabolites"]
        met["name"] = mid_names[met["id"]]
    end

    save_model(jmodel, joinpath(base_model_loc, "yeast-GEM.json"))

    core_ecoli_url = "http://bigg.ucsd.edu/static/models/e_coli_core.json"
    core_ecoli_path = joinpath(base_model_loc, "e_coli_core.json")
    download(core_ecoli_url, core_ecoli_path)

    iML1515_url = "http://bigg.ucsd.edu/static/models/iML1515.json"
    iML1515_path = joinpath(base_model_loc, "iML1515.json")
    download(iML1515_url, iML1515_path)

    log_file_loc = joinpath(base_model_loc, "log.txt")
    open(log_file_loc, "w") do io
        write(io, "Model files downloaded from urls on ", string(now()), "\n")
        write(io, "Consensus yeast model: ", yeast_model_url, "\n")
        write(io, "E. coli iML1514 model: ", iML1515_url, "\n")
        write(io, "E. coli core model: ", core_ecoli_url, "\n")
    end

    return nothing
end

download_models()
