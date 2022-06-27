module SubgradientDescent

using COBREXA, DifferentiableMetabolism
using JuMP, CPLEX
using JSON
using SparseArrays, Statistics, LinearAlgebra

function make_ecoli_pruned_model(
    master_id,
    rescale_factor;
    gene_bound_factor = 5,
    base_load_path = joinpath("model_construction", "processed_models_files", "ecoli"),
    data_dir = joinpath("data", "kcats", "heckmann2020"),
    master_kos = Dict(),
)

    #: import model and data files
    model = load_model(StandardModel, joinpath(base_load_path, "iml1515_fixed.json"))
    reaction_kcats = JSON.parsefile(joinpath(base_load_path, "reaction_kcats.json"))
    reaction_protein_stoichiometry =
        JSON.parsefile(joinpath(base_load_path, "protein_stoichiometry.json"))
    gene_product_molar_mass = Dict(
        k => v for
        (k, v) in JSON.parsefile(joinpath(base_load_path, "protein_masses.json"))
    )

    #: unconstrain because enzyme constraints take over
    model.reactions["EX_glc__D_e"].lb = -1000.0

    #: load measured protein data
    proteindata = JSON.parsefile(joinpath(data_dir, "proteindata.json"))
    fluxdata = JSON.parsefile(joinpath(data_dir, "fluxdata.json"))

    #: load reaction isozyme data
    _reaction_isozymes = Dict(
        rid => [
            Isozyme(
                Dict(
                    k => v for (k, v) in zip(
                        reaction_gene_association(model, rid)[i],
                        reaction_protein_stoichiometry[rid][i],
                    )
                ),
                reaction_kcats[rid][i][1],
                reaction_kcats[rid][i][2],
            ) for i = 1:length(reaction_kcats[rid])
        ] for rid in keys(reaction_kcats)
    )

    #: get only the expressed grr 
    reaction_isozymes = Dict{String,Vector{Isozyme}}()
    for rid in filter(x -> haskey(_reaction_isozymes, x), reactions(model))
        total_prot_isozymes = Float64[]
        for isozyme in _reaction_isozymes[rid]
            push!(
                total_prot_isozymes,
                sum([
                    get(proteindata[master_id], gid, 0.0) for
                    gid in keys(isozyme.gene_product_count)
                ]),
            )
        end
        if all(total_prot_isozymes .≈ 0.0)
            total_prot_isozymes = Float64[]
            for isozyme in _reaction_isozymes[rid]
                push!(
                    total_prot_isozymes,
                    sum([
                        get(proteindata[master_id], gid, 1.0) for
                        gid in keys(isozyme.gene_product_count)
                    ]),
                )
            end
            reaction_isozymes[rid] = [_reaction_isozymes[rid][argmax(total_prot_isozymes)]]
        else
            reaction_isozymes[rid] = [_reaction_isozymes[rid][argmax(total_prot_isozymes)]]
        end
    end

    ko_id = first(split(master_id, "#"))
    if !isempty(master_kos[ko_id])
        delete!.(Ref(model.reactions), master_kos[ko_id])
        delete!.(Ref(reaction_isozymes), master_kos[ko_id])
    end

    modifications = [
        change_optimizer_attribute("CPXPARAM_Emphasis_Numerical", 1),
        change_optimizer_attribute("CPX_PARAM_SCAIND", 1),
        change_optimizer_attribute("CPX_PARAM_THREADS", 1),
        COBREXA.silence,
    ]

    #: setup gecko data, use observed genes to find all reactions that should be active
    gene_product_mass_group_bound = Dict("uncategorized" => 1000_000.0) # overestimate, pruned later

    gene_product_bounds = SubgradientDescent.get_ecoli_gene_bounds(
        model,
        reaction_isozymes,
        gene_product_molar_mass,
        gene_product_mass_group_bound,
        proteindata,
        master_id,
        rescale_factor;
        gene_bound_factor,
        modifications,
    )

    gm = make_gecko_model(
        model;
        reaction_isozymes,
        gene_product_bounds,
        gene_product_molar_mass,
        gene_product_mass_group_bound,
    )

    opt_model = flux_balance_analysis(gm, CPLEX.Optimizer; modifications)
    rfluxes = flux_dict(gm, opt_model)
    gpconcs = gene_product_dict(gm, opt_model)

    #: prune model
    pmodel = prune_model(model, rfluxes)

    rs_to_del = filter(x -> x ∉ reactions(pmodel), keys(reaction_isozymes))
    delete!.(Ref(reaction_isozymes), rs_to_del)
    gs_to_del = filter(x -> x ∉ genes(pmodel), keys(gene_product_molar_mass))
    delete!.(Ref(gene_product_molar_mass), gs_to_del)

    gene_product_bounds = Dict(gid => (0.0, 1000.0) for gid in genes(pmodel))

    pproteindata = Dict(
        k => v * rescale_factor for (k, v) in proteindata[master_id] if k in genes(pmodel)
    )
    pfluxdata =
        Dict(k => abs(first(v)) for (k, v) in fluxdata[master_id] if k in reactions(pmodel))

    return pmodel,
    reaction_isozymes,
    gene_product_bounds,
    gene_product_molar_mass,
    pproteindata,
    pfluxdata
end

"""
Get bounds for genes that match experimental measurements
"""
function get_ecoli_gene_bounds(
    model,
    reaction_isozymes,
    gene_product_molar_mass,
    gene_product_mass_group_bound,
    proteindata,
    master_id,
    rescale_factor;
    gene_bound_factor = 5.0,
    modifications,
)
    gene_product_bounds = Dict(gid => (0.0, 9000.0) for gid in genes(model))

    for gid in genes(model)
        v = get(proteindata[master_id], gid, 0.0)
        if rescale_factor * v > 1e-1

            gene_product_bounds[gid] = (
                rescale_factor * v / gene_bound_factor,
                rescale_factor * v * gene_bound_factor,
            )

            #: run gecko to prune the model
            gm = make_gecko_model(
                model;
                reaction_isozymes,
                gene_product_bounds,
                gene_product_molar_mass,
                gene_product_mass_group_bound,
            )

            opt_model = flux_balance_analysis(gm, CPLEX.Optimizer; modifications)
            rfluxes = flux_dict(gm, opt_model)

            if isnothing(rfluxes)
                gene_product_bounds[gid] = (0.0, 9000.0)
            end
        end
    end
    return gene_product_bounds
end

function qp_objective_measured(
    rids,
    gids,
    obs_v_dict,
    obs_e_dict;
    vtol = 1e-3,
    etol = 1e-3,
    reg = 1e-3,
)
    n_vars = length(rids)
    reaction_ids = filter(x -> x ∉ gids, rids)
    c = zeros(n_vars) #  linear component
    q = zeros(n_vars) #  quadratic component
    n = 0 # count number of loss variables

    # flux
    for (i, rid) in enumerate(reaction_ids)
        unmangled_rid = first(split(rid, "#"))
        if !haskey(obs_v_dict, unmangled_rid) || abs(obs_v_dict[unmangled_rid]) < vtol
            q[i] = reg
        else
            c[i] = -1.0 / abs(obs_v_dict[unmangled_rid]) #! fluxes are positive in model
            q[i] = 1.0 / obs_v_dict[unmangled_rid]^2
            n += 1
        end
    end

    # protein
    k = length(reaction_ids)
    for (i, gid) in enumerate(gids)
        if !haskey(obs_e_dict, gid) || abs(obs_e_dict[gid]) < etol
            q[k+i] = reg
        else
            c[k+i] = -1.0 / obs_e_dict[gid]
            q[k+i] = 1.0 / obs_e_dict[gid]^2
            n += 1
        end
    end

    return spdiagm(q), sparse(c), n
end

function save_data(savedir, iter, θ, dLdθ, loss, nvars)
    open(joinpath(savedir, "iter_" * string(iter) * ".json"), "w") do io
        JSON.print(io, Dict("θ" => θ, "dx" => dLdθ, "loss" => loss, "nvars" => nvars))
    end
end

function linear_seach_descent(
    diffmodel;
    niters = 10,
    η_start = 0.01,
    ls_maxiter = 20,
    ls_decr = 0.5,
    modifications = [],
    verbose = false,
    savedir = "results",
    cap = (θ, dθ) -> dθ,
    save_on_iter = 100,
    loss_offset = 1,
    reg_kcats = 0,
)

    #: make reference problem
    rf = ReferenceFunctions(diffmodel)

    #: rescale equality constraints
    rescale_factors = scaling_factor(rf.E(diffmodel.θ), rf.d(diffmodel.θ))
    update_E!(diffmodel, θ -> rescale_factors .* rf.E(θ))
    update_d!(diffmodel, θ -> rescale_factors .* rf.d(θ))

    x = zeros(rf.nx)
    ν = zeros(rf.neq)
    λ = zeros(rf.nineq)
    A = spzeros(rf.nx + rf.neq + rf.nineq, rf.nx + rf.neq + rf.nineq)
    B = zeros(rf.nx + rf.neq + rf.nineq, rf.nθ)
    dxdθ = zeros(rf.nx + rf.neq + rf.nineq, rf.nθ)

    var_deriv, par_deriv = make_derivatives_with_equality_scaling(rf) # very time consuming
    diffmodel.analytic_var_derivs = (x, ν, λ, θ) -> var_deriv(x, ν, λ, θ, rescale_factors)
    diffmodel.analytic_par_derivs = (x, ν, λ, θ) -> par_deriv(x, ν, λ, θ, rescale_factors)

    #: setup gradient descent
    dLdθ = zeros(length(diffmodel.θ))
    θ = zeros(length(diffmodel.θ))
    dLdx = zeros(size(diffmodel.Q(diffmodel.θ), 1))

    loss_offset = 0.5 * loss_offset

    save_iters_extra = [
        collect(1:10)
        collect(20:10:100)
        collect(200:100:1000)
    ]
    for iter = 1:niters

        differentiate!(
            x,
            ν,
            λ,
            A,
            B,
            dxdθ,
            diffmodel,
            CPLEX.Optimizer;
            modifications,
            use_analytic = true,
            scale_output = false,
        )

        loss =
            0.5 * x' * diffmodel.Q(diffmodel.θ) * x +
            diffmodel.c(diffmodel.θ)' * x +
            loss_offset

        dLdx .= diffmodel.Q(diffmodel.θ) * x + diffmodel.c(diffmodel.θ)

        dLdθ .= dxdθ[1:rf.nx, :]' * dLdx .+ reg_kcats .* diffmodel.θ 

        #: line search gradient descent
        θ .= diffmodel.θ
        η = η_start # for memory to display later
        for ls_iter = 1:ls_maxiter
            η = η_start * ls_decr^(ls_iter - 1)
            diffmodel.θ .= θ - map(cap, θ, η .* dLdθ)
            current_loss = solve_model(loss_offset, diffmodel; modifications)
            if !isnothing(current_loss)
                if current_loss < loss
                    break
                end
            end
        end

        @assert(all(diffmodel.θ .> 0))

        #: rescale if necessary
        ((best_scaling, current_scaling), _, _) = check_scaling(diffmodel; verbose = false)
        if current_scaling - best_scaling > 1.5
            rescale_factors = scaling_factor(rf.E(diffmodel.θ), rf.d(diffmodel.θ))
            update_E!(diffmodel, θ -> rescale_factors .* rf.E(θ))
            update_d!(diffmodel, θ -> rescale_factors .* rf.d(θ))

            diffmodel.analytic_var_derivs =
                (x, ν, λ, θ) -> var_deriv(x, ν, λ, θ, rescale_factors)
            diffmodel.analytic_par_derivs =
                (x, ν, λ, θ) -> par_deriv(x, ν, λ, θ, rescale_factors)
        end

        #: save data
        (iter % save_on_iter == 0 || iter in save_iters_extra) &&
            save_data(savedir, iter, diffmodel.θ, dLdθ, loss, loss_offset)

        #: print loss
        verbose && println(
            "Iteration: ",
            iter,
            " with loss: ",
            round(loss, digits = 2),
            " with scaling = ",
            round(current_scaling, digits = 2),
            " with η = ",
            η,
            " with loss variables = ",
            loss_offset,
        )
    end
    return nothing
end

function solve_model(loss_offset, diffmodel::DifferentiableModel; modifications = [])
    #: forward pass, solve the optimization problem
    opt_model = JuMP.Model(CPLEX.Optimizer)
    set_silent(opt_model)
    @variable(opt_model, z[1:length(diffmodel.var_ids)])

    if all(diffmodel.Q(diffmodel.θ) .== 0) # is LP
        @objective(opt_model, Min, diffmodel.c(diffmodel.θ)' * z)
    else
        @objective(
            opt_model,
            Min,
            0.5 * z' * diffmodel.Q(diffmodel.θ) * z + diffmodel.c(diffmodel.θ)' * z
        )
    end

    @constraint(opt_model, eq, diffmodel.E(diffmodel.θ) * z .== diffmodel.d(diffmodel.θ))
    @constraint(opt_model, ineq, diffmodel.M(diffmodel.θ) * z .<= diffmodel.h(diffmodel.θ))

    # apply the modifications
    for mod in modifications
        mod(nothing, opt_model)
    end

    optimize!(opt_model)

    if termination_status(opt_model) ∉ [JuMP.OPTIMAL, JuMP.LOCALLY_SOLVED]
        nothing
    else
        x = value.(opt_model[:z])
        0.5 * x' * diffmodel.Q(diffmodel.θ) * x +
        diffmodel.c(diffmodel.θ)' * x +
        loss_offset
    end
end

end

