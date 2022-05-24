using JSON, DataFrames, DataFramesMeta, CSV

df = DataFrame(CSV.File(joinpath("data", "metabolites", "bennet2016", "metabolites.csv")))

d = Dict()
for i = 1:size(df, 1)
    mid = df[i, :MetaboliteID]
    conc = df[i, :ConcentrationMolar]
    d[mid] = conc
end

open(joinpath("data", "metabolites", "bennet2016", "metabolites.json"), "w") do io
    JSON.print(io, d)
end
