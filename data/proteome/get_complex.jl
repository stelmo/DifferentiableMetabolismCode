using JSON

function get_complexes(organism, complexportalnum)

    pid_bid = Dict()
    open(joinpath("data", "proteome", "$organism.tab")) do io
        firstline = true
        for ln in eachline(io)
            firstline && (firstline = false; continue)
            prts = split(ln, "\t")
            prts[3] == "" && continue
            pid = prts[1]
            bnum = first(split(prts[3], " "))
            pid_bid[pid] = bnum
        end
    end

    data = Dict()
    open(joinpath("data", "proteome", complexportalnum * ".tsv")) do io
        firstline = true
        for ln in eachline(io)
            firstline && (firstline = false; continue)
            prts = split(ln, "\t")
            id = prts[1]
            raw_stoich = prts[5]
            stoich = Dict{String,Float64}()
            for elem in split(raw_stoich, "|")
                if startswith(elem, "P") || startswith(elem, "Q")
                    pid_stoich = split(elem, "(")
                    pid = first(split(first(pid_stoich), "-"))
                    stoich_num = first(split(last(pid_stoich), ")"))
                    !haskey(pid_bid, pid) && continue
                    bid = pid_bid[pid]
                    stoich[bid] = parse(Float64, stoich_num)
                end
            end
            isempty(stoich) && continue
            data[id] = stoich
        end
    end

    open(joinpath("data", "proteome", "$(organism)_complex.json"), "w") do io
        JSON.print(io, data)
    end

    return nothing
end

get_complexes("e_coli", "83333")
get_complexes("s_cere", "559292")
