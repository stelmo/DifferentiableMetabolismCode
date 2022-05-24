using JSON, COBREXA, DataFrames, CSV, Statistics, DataFramesMeta

dfr = DataFrame(CSV.File(joinpath("data", "kms", "kroll2021", "kmdata.csv")))
rename!(dfr, Dict("KM [mM]" => :Km))
@select!(dfr, :Km, :RIDS, :MIDS)

model =
    load_model(StandardModel, joinpath("model_construction", "model_files", "iML1515.json"))

d = Dict()
for df in groupby(dfr, :RIDS)
    rid = String(df[1, :RIDS])
    rs = keys(reaction_stoichiometry(model, rid))
    td = Dict()
    for i = 1:size(df, 1)
        td[df[i, :MIDS]] = df[i, :Km]
    end

    ks = keys(td)
    if all(in.(ks, Ref(rs))) && all(in.(rs, Ref(ks)))
        d[rid] = td
    else
        @warn("Missing metabolites for $rid")
    end
end

open(joinpath("data", "kms", "kroll2021", "kmdata_iml1515.json"), "w") do io
    JSON.print(io, d)
end

