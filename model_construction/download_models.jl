using COBREXA, Dates, ReadableRegex

function download_models()
    base_model_loc = joinpath("model_construction", "model_files")

    for f in readdir(base_model_loc)
        rm(joinpath(base_model_loc, f))
    end

    core_ecoli_url = "http://bigg.ucsd.edu/static/models/e_coli_core.json"
    core_ecoli_path = joinpath(base_model_loc, "e_coli_core.json")
    download(core_ecoli_url, core_ecoli_path)

    iML1515_url = "http://bigg.ucsd.edu/static/models/iML1515.json"
    iML1515_path = joinpath(base_model_loc, "iML1515.json")
    download(iML1515_url, iML1515_path)

    log_file_loc = joinpath(base_model_loc, "log.txt")
    open(log_file_loc, "w") do io
        write(io, "Model files downloaded from urls on ", string(now()), "\n")
        write(io, "E. coli iML1514 model: ", iML1515_url, "\n")
        write(io, "E. coli core model: ", core_ecoli_url, "\n")
    end

    return nothing
end

download_models()
