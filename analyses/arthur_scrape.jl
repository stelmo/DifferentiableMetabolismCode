using JSON, DataFrames, DifferentiableMetabolism, COBREXA, CSV

function get_data_out_as_df(results_dir, save_dir)
    results_name = last(split(results_dir, "/"))
    files = filter(startswith("iter_"), readdir(joinpath(results_dir)))
    diffmodel = JSON.parsefile(joinpath(results_dir, "var_params.json"))

    param_ids = diffmodel["param_ids"]
    n_params = length(param_ids)
    n_files = length(files)

    losses = zeros(n_files)
    iter_nums = zeros(Int64, n_files)

    df_params = DataFrame(
        Condition = String[],
        KcatID = String[],
        Kcat = Float64[],
        Derivative = Float64[],
        Iteration = Float64[],
    )

    for (i, file) in enumerate(files)
        j = JSON.parsefile(joinpath(results_dir, file))
        iter_num = Int(parse(Int, last(split(first(split(file, ".")), "_"))))
        iter_nums[i] = iter_num
        losses[i] = convert(Float64, j["loss"])/(convert(Float64, j["nvars"])*2)

        append!(
            df_params,
            DataFrame(
                Condition = fill(results_name, n_params),
                KcatID = param_ids,
                Kcat = convert.(Float64, j["θ"]),
                Derivative = convert.(Float64, j["dx"]),
                Iteration = fill(iter_num, n_params),
            ),
        )
    end

    df_losses = DataFrame(
        Condition = fill(results_name, n_files),
        Loss = losses,
        Iteration = iter_nums,
    )

    CSV.write(joinpath(save_dir, "$results_name#losses.csv"), df_losses)
    CSV.write(joinpath(save_dir, "$results_name#params.csv"), df_params)

    nothing
end

function get_all_dfs(results_dir, save_dir)
    for dir in readdir(results_dir)
        get_data_out_as_df(joinpath(results_dir, dir), save_dir)
    end
end


# function get_iter_loss(dir)
#     iter_counts = Int[]
#     losses = Float64[]
#     for file in readdir(dir)
#         !startswith(file, "iter_") && continue
#         push!(iter_counts, parse(Int64, last(split(first(split(file, ".")), "_"))))
#         push!(losses, JSON.parsefile(joinpath(dir, file))["loss"])
#     end
#     pidxs = sortperm(iter_counts)
#     return iter_counts[pidxs], losses[pidxs]
# end

# function plot_loss(dir; kwargs...)
#     iter_counts, losses = get_iter_loss(dir)
#     return lineplot(
#         iter_counts,
#         losses;
#         xscale = log10,
#         yscale = log10,
#         xlabel = "Iterations",
#         ylabel = "Loss",
#         kwargs...,
#     )
# end

# function plot_loss!(plt, dir, name)
#     iter_counts, losses = get_iter_loss(dir)
#     plt = lineplot!(plt, iter_counts, losses; name) # inherits some kwargs
#     return plt
# end

# function plot_all_losses(master_dir, prefix; kwargs...)
#     plt = nothing
#     for dir in readdir(master_dir)
#         !startswith(dir, prefix) && continue
#         if isnothing(plt)
#             plt = plot_loss(joinpath(master_dir, dir); name = dir, kwargs...)
#         else
#             plt = plot_loss!(plt, joinpath(master_dir, dir), dir)
#         end
#     end
#     return println(plt)
# end

# function get_kcats()

#     base_dir = "ecoli_kcats"
#     largest_kcats = Dict{String,Float64}()
#     for cond in readdir(base_dir)
#         println(cond)
#         jf = JSON.parsefile(joinpath(base_dir, cond, "problem_data_$cond.json"))
#         kcats = jf["kcats"]
#         largest_n = -1
#         for iter in filter(startswith("iter_"), readdir(joinpath(base_dir, cond)))
#             n = Int(parse(Int, replace(split(iter, ".")[1], "iter_" => "")))
#             if n > largest_n
#                 largest_n = n
#             end
#         end
#         jf = JSON.parsefile(joinpath(base_dir, cond, "iter_$largest_n.json"))
#         for (k, v) in Dict(kcats .=> jf["θ"][1:end-1])
#             if haskey(largest_kcats, k)
#                 if largest_kcats[k] < v
#                     largest_kcats[k] = v
#                 end
#             else
#                 largest_kcats[k] = v
#             end
#         end
#     end

#     open("grad_desc_kcats.json", "w") do io
#         JSON.print(io, largest_kcats)
#     end

# end
