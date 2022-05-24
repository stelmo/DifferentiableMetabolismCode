using JSON

"""
format of data
proteindata[sample_id§bio_rep][gid] = fmol_per_gDW_calib_avg
fluxdata[sample_id][rid] = mfa_sampling_ave
"""
function get_proteomics_fluxomics()

    proteindata = Dict{String,Dict{String,Float64}}()
    fluxdata = Dict{String,Dict{String,Tuple{Float64,Float64,Float64}}}()
    open(joinpath("data", "kcats", "heckmann2020", "Dataset_S1A_protein_ab.csv")) do io
        firstline = true
        for ln in eachline(io)
            firstline && (firstline = false; continue)
            prts = split(ln, ",")
            rid = prts[1]
            sample_id = prts[3]
            bnum = prts[5]
            bio_rep = prts[6]
            fmol = parse(Float64, prts[7])
            mfa = prts[9]
            mfa_lb = prts[11]
            mfa_ub = prts[12]

            master_id = sample_id * "#" * bio_rep
            t = get(proteindata, master_id, Dict{String,Float64}())
            t[bnum] = fmol
            proteindata[master_id] = t

            if mfa != "NA"
                if endswith(rid, "_b")
                    rid = rid[1:end-2]
                    dir = -1.0
                else
                    dir = 1.0
                end
                t = get(fluxdata, master_id, Dict{String,Tuple{Float64,Float64,Float64}}())
                t[rid] = (
                    dir * parse(Float64, mfa),
                    parse(Float64, mfa_lb),
                    parse(Float64, mfa_ub),
                )
                fluxdata[master_id] = t
            end
        end
    end

    open(joinpath("data", "kcats", "heckmann2020", "fluxdata.json"), "w") do io
        JSON.print(io, fluxdata)
    end
    open(joinpath("data", "kcats", "heckmann2020", "proteindata.json"), "w") do io
        JSON.print(io, proteindata)
    end

    return nothing
end

get_proteomics_fluxomics()
