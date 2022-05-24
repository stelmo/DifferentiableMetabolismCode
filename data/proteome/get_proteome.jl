using JSON

function get_proteomic_data(organism)
    proteome_data = Dict()
    open(joinpath("data", "proteome", "$organism.tab")) do io
        firstline = true
        for ln in eachline(io)
            firstline && (firstline = false; continue)
            prts = split(ln, "\t")
            prts[3] == "" && continue
            bnum = first(split(prts[3], " "))
            mass = parse(Float64, replace(prts[2], "," => "")) / 1000.0 # in Da originally
            if prts[4] == ""
                subunit = ""
            else
                subunit =
                    split(replace(replace(prts[4], "SUBUNIT: " => ""), "." => ""), " ")
                subunit = subunit[1]
                if !contains(subunit, "mer")
                    subunit = ""
                end
            end
            proteome_data[bnum] = (mass, subunit)
        end
    end

    open(joinpath("data", "proteome", "$organism.json"), "w") do io
        JSON.print(io, proteome_data)
    end

    return nothing
end

get_proteomic_data("e_coli")
get_proteomic_data("s_cere")
