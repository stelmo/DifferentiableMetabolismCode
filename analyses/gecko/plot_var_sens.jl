using ColorSchemes, CairoMakie, DataFrames, DataFramesMeta, Chain, CSV, Statistics

df = DataFrame(CSV.File(joinpath("results", "gecko", "variability_sensitivities2.csv")))


enzymes = [
    "ATPS4rpp"
    "GLUDy"
    "CYTBO3_4pp"
    "GAPD"
    "PFK"
    "PGM"
    "PPS"
    "TPI"
    "PGK"
    "FBA"
    "PGI"
    "ENO"
]

rdf = @subset df @byrow begin 
    :Reaction in enzymes
end

fig = Figure();
ax = Axis(
    fig[1, 1], 
    # yticks = ((1:12) ./ 4,  reverse(months)),
    # yscale=log10,
);

for (i, enzyme) in enumerate(enzymes)
    rrdf = @subset rdf @byrow begin 
        :Reaction == enzyme
    end 
    x = rrdf[!, :Sensitivity]
    d = density!(
        ax, 
        x.*1_000, 
        offset = i,
        color = :x, 
        # colormap = :thermal, 
        # colorrange = (-5, 5),
        # strokewidth = 1, 
        # strokecolor = :black,
    )
end
fig

mean(x)
std(x)

rrdf = @subset rdf @byrow begin 
    :Reaction == "FBA"
end 
x = rrdf[!, :Sensitivity]
d = density(
    x.*10000, 
    color = :x, 
    # colormap = :thermal, 
    # colorrange = (-5, 5),
    # strokewidth = 1, 
    # strokecolor = :black,
)