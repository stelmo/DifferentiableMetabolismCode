using JSON

function get_kcats_ecoli()
    kcat_data = Dict()
    open(joinpath("data", "kcats", "heckmann2020", "Dataset_S1C_turnonver_n.csv")) do io
        firstline = true
        for ln in eachline(io)
            firstline && (firstline = false; continue)
            prts = split(ln, ",")
            kcat_data[prts[1]] = parse(Float64, prts[2:end][3])
        end
    end

    open(joinpath("data", "kcats", "heckmann2020", "e_coli_kcats.json"), "w") do io
        JSON.print(io, kcat_data)
    end

    return nothing
end

get_kcats_ecoli()

function get_kcats_yeast()
    kcat_data = Dict()
    open(joinpath("data", "kcats", "chen2021", "turnover_number.csv")) do io
        firstline = true
        for ln in eachline(io)
            firstline && (firstline = false; continue)
            prts = split(ln, "\t")
            id = join(split(prts[1], "_")[1:2], "_")
            kcat_data[id] = parse(Float64, prts[2])
        end
    end

    open(joinpath("data", "kcats", "chen2021", "s_cere_kcats.json"), "w") do io
        JSON.print(io, kcat_data)
    end

    return nothing
end

get_kcats_yeast()
