using COBREXA, DifferentiableMetabolism
using CPLEX, JSON, SparseArrays, Statistics, Polynomials, CairoMakie

include(joinpath("analyses", "subgradient_descent.jl"))
using .SubgradientDescent

include(joinpath("analyses", "data_constants.jl"))
using .DataConstants

conds = [
    "WT1",
    "WT2",
    "pgi1",
    "pgi2",
    "pgi3",
    "pgi4",
    "pgi5",
    "pgi6",
    "pgi7",
    "pgi8",
    "pts1",
    "pts2",
    "pts3",
    "pts4",
    "sdh1",
    "sdh2",
    "sdh3",
    "tpi1",
    "tpi2",
    "tpi3",
    "tpi4",
]

rescale_factor = 1e-6

protein_upper_bound = Dict{String,Float64}()
for cond in conds
    for experi in ["B1", "B2"]
        master_id = "$cond#$experi"
        try
            pmodel,
            reaction_isozymes,
            gene_product_bounds,
            gene_product_molar_mass,
            pruned_protein_data,
            pruned_flux_data = SubgradientDescent.make_ecoli_pruned_model(
                master_id,
                rescale_factor;
                master_kos = DataConstants.master_kos,
            )

            masses = [0.0]
            gid_added = String[]
            for (rid, isozymes) in reaction_isozymes # these are metabolic enzymes
                isozyme = first(isozymes)
                gids = collect(keys(isozyme.gene_product_count))
                if !any(in.(gids, Ref(gid_added))) &&
                   all(haskey.(Ref(pruned_protein_data), gids))
                    mass = sum([
                        pruned_protein_data[gid] * gene_product_molar_mass[gid] for
                        gid in gids
                    ])
                    push!(masses, last(masses) + mass)
                    append!(gid_added, gids)
                end
            end
            masses = masses[2:end]

            linfit = fit(1:length(masses), masses, 1)

            metabolic_mass = linfit(length(reaction_isozymes))

            protein_upper_bound[master_id] = metabolic_mass
        catch err
            println("Failed on: ", master_id)
        end

    end

end
protein_upper_bound
master_id = "WT1#B1"
