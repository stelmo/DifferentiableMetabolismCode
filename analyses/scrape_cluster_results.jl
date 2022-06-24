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