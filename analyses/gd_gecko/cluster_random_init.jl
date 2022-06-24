using COBREXA, DifferentiableMetabolism
using CPLEX, JSON, SparseArrays, Statistics
using ArgParse

if Sys.iswindows()
    include(joinpath("analyses", "subgradient_descent.jl"))
    include(joinpath("analyses", "data_constants.jl"))
else
    include(joinpath("..", "subgradient_descent.jl"))
    include(joinpath("..", "data_constants.jl"))
end
using .SubgradientDescent
using .DataConstants

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "master_id"
        required = true
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    master_id = parsed_args["master_id"]
    
    master_id_plus = master_id * "#" * join(rand(["A","B", "C", "D", "E", "F"], 8),"")
    savedir = "/nfs/qtbscratch/stelmo/differentiable-metabolism-2022-random-init/$master_id_plus"
    # savedir = "results"
    !isdir(savedir) && mkdir(savedir)

    rescale_factor = 1e-6

    #: import pruned model 
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

    theoretical_protein_bound = DataConstants.protein_upper_bound[master_id]
    gene_product_mass_group_bound = Dict("uncategorized" => theoretical_protein_bound) # overestimate, pruned later

    #: pruned gecko model
    pgm = make_gecko_model(
        pmodel;
        reaction_isozymes,
        gene_product_bounds,
        gene_product_molar_mass,
        gene_product_mass_group_bound,
    )

    #: make differentiable model
    rid_enzyme = Dict(
        rid => isozyme_to_enzyme(first(isozyme_vec), gene_product_molar_mass) for
        (rid, isozyme_vec) in reaction_isozymes
    )

    #: randomize the kcatskk
    for v in values(rid_enzyme)
        v.kcat *= 10^(2*rand() - 1)
    end

    #: make differentiable model
    diffmodel =
        with_parameters(pgm, rid_enzyme; scale_equality = true, scale_inequality = true)
    
    open(joinpath(savedir, "var_params.json"), "w") do io 
        d = Dict("var_ids" => diffmodel.var_ids, "param_ids" => diffmodel.param_ids)
        JSON.print(io, d)
    end

    _Q, _c, loss_offset = SubgradientDescent.qp_objective_measured(
        reactions(pgm),
        genes(pgm),
        pruned_flux_data,
        pruned_protein_data;
        vtol = 1e-3,
        etol = 1e-3,
        reg = 1e-3,
    )

    update_c!(diffmodel, _ -> _c)
    update_Q!(diffmodel, _ -> _Q)

    cap(θ, dθ; max_frac_change = 0.5, lb = 0.0001, ub = 100.0) = begin
        t = sign(dθ) * min(max_frac_change, abs(dθ / θ)) * θ
        (θ - t <= lb) || (θ - t >= ub) ? 0.0 : t
    end

    # check_scaling(diffmodel; verbose=true)

    modifications = [
        change_optimizer_attribute("CPXPARAM_Emphasis_Numerical", 1),
        change_optimizer_attribute("CPX_PARAM_SCAIND", 1),
        change_optimizer_attribute("CPX_PARAM_THREADS", 1),
        COBREXA.silence,
    ]

    SubgradientDescent.linear_seach_descent(
        diffmodel;
        niters = 10_000_000,
        η_start = 0.1,
        ls_maxiter = 25,
        ls_decr = 0.5,
        modifications,
        verbose = false,
        savedir,
        cap = cap,
        save_on_iter = 250,
        loss_offset,
    )

    nothing
end

main()
