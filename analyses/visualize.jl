module Visualize

using CairoMakie, COBREXA, ColorSchemes, JSON, Statistics
using Escher, COBREXA

pseudolog(x) = sign(x) * log10(abs(x) + 1)

function plot_heatmap(
    diffmodel,
    dx,
    heatmap_ax,
    color_bar_ax;
    ignore_vids = [],
    ignore_pids = [],
    transform_data = x -> x,
    n_σ = 1.0,
    atol = 1e-6,
)
    vids =
        isempty(ignore_vids) ? diffmodel.var_ids :
        filter(x -> x ∉ ignore_vids, diffmodel.var_ids)
    pids =
        isempty(ignore_pids) ? diffmodel.param_ids :
        filter(x -> x ∉ ignore_pids, diffmodel.param_ids)
    vidxs = Int.(indexin(vids, diffmodel.var_ids))
    pidxs = Int.(indexin(pids, diffmodel.param_ids))

    ndx = transform_data.(dx[:, pidxs][vidxs, :]) # select heatmap

    neg_and_pos = maximum(ndx) > atol && minimum(ndx) < -atol
    μ = mean(ndx)
    σ = std(ndx)
    _lclip = μ - n_σ * σ
    lclip = !neg_and_pos && _lclip < 0.0 ? 0.0 : _lclip
    hclip = μ + n_σ * σ

    if neg_and_pos
        colormap = Reverse(ColorSchemes.RdYlBu_9[2:8])
        highclip = ColorSchemes.RdYlBu_9[1]
        lowclip = ColorSchemes.RdYlBu_9[9]
    else
        colormap = ColorSchemes.YlOrRd_9[2:8]
        highclip = ColorSchemes.YlOrRd_9[9]
        lowclip = ColorSchemes.YlOrRd_9[1]
    end

    hm = heatmap!(
        heatmap_ax,
        1:size(ndx, 1),
        1:size(ndx, 2),
        ndx,
        colormap = colormap,
        colorrange = (lclip, hclip),
        highclip = highclip,
        lowclip = lowclip,
    )

    ylabs = [first(split(x, "#")) for x in vids]
    heatmap_ax.xticks = (1:length(vids), ylabs)
    heatmap_ax.yticks = (1:length(pids), [last(split(x, "#")) for x in pids])

    Colorbar(color_bar_ax, hm, label = "Scaled sensitivity", labelsize=20)

    return ylabs
end

function get_core_metabolism(model, core_met_loc)
    coremet = JSON.parsefile(core_met_loc)
    core_rxns = [rxn["bigg_id"] for rxn in values(coremet[2]["reactions"])]
    rxns = intersect(core_rxns, reactions(model))
    gs = String[]
    mets = String[]
    for rid in rxns
        !has_reaction_grr(model, rid) && continue
        for g in reaction_gene_association(model, rid)
            append!(gs, g)
        end
        append!(mets, collect(keys(reaction_stoichiometry(model, rid))))
    end
    return rxns, unique(gs), unique(mets)
end

"""
Check if reaction has a gene reaction rule assigned to it.
"""
function has_reaction_grr(model, rid)
    grr = reaction_gene_association(model, rid)
    if isnothing(grr) || isempty(grr) || isnothing(first(grr)) || isempty(first(grr))
        return false
    end
    return true
end

function plot_gradients(results_dir)

    dirs = readdir(results_dir)
    j = JSON.parsefile(joinpath(results_dir, first(dirs)))
    dx_mat = zeros(length(j["dx"]), length(dirs))
    for dir in dirs
        iter_num = parse(Int, last(split(first(split(dir, ".")), "_")))
        j = JSON.parsefile(joinpath(results_dir, dir))
        dx_mat[:, iter_num] .= j["dx"]
    end

    fig = Figure(resolution = (1000, 4000))
    ax = Axis(
        fig[1, 1],
        # yscale=log10,
        xscale = log10,
    )
    heatmap!(
        ax,
        1:size(dx_mat, 2),
        1:size(dx_mat, 1),
        dx_mat',
        colormap = Reverse(ColorSchemes.RdYlBu_9[2:8]),
        highclip = ColorSchemes.RdYlBu_9[1],
        lowclip = ColorSchemes.RdYlBu_9[9],
        colorrange = (-1.0, 1.0),
    )
    fig

    CairoMakie.FileIO.save("dx.pdf", fig)

end

function plot_escher_viz(ax1, ax2, rids, sol, reaction_mass, escher_loc; mts=4, rts=4)

    rd = Dict(k => v > 0 ? :f : :b for (k, v) in sol)

    maxmass = maximum(abs.(values(reaction_mass)))
    minmass = minimum(abs.(values(reaction_mass)))
    width_interp(x) = 2 + 5 * (abs(x) - minmass) / (maxmass - minmass) # widths between 2 and 5
    re = Dict(k => width_interp(v) for (k, v) in reaction_mass) # map reaction id to reaction edge width  

    # Find min and max absolute fluxes for normalization
    fluxes = [sol[x] for x in rids if haskey(sol, x)]
    maxflux = maximum(log.(abs.(values(fluxes))))
    minflux = minimum(log.(abs.(values(fluxes))))

    # Scale color of reaction edges to fluxes (manually binned)
    rng = range(0, stop = 1.0, length = 9)
    bins = reverse([a .. b for (a, b) in zip(rng[1:end-1], rng[2:end])])
    indx(y) = findfirst(x -> y in x, bins)

    color_interp(x) = begin
        normed_x = (log(abs(x)) - minflux) / (maxflux - minflux)
        ColorSchemes.RdYlBu_9[indx(normed_x)]
    end
    rc = Dict(k => color_interp(v) for (k, v) in sol if k in rids) # map reaction id to reaction edge color

    # Normal Makie plotting features all work (escherplot is a full recipe)
    escherplot!(
        ax1,
        escher_loc;
        reaction_edge_colors = rc,
        reaction_show_text = true,
        metabolite_show_text = true,
        annotation_show_text = true,
        # reaction_directions = rd,
        reaction_edge_widths = re,
        reaction_text_size = rts,
        metabolite_text_size = mts,
    )
    hidexdecorations!(ax1)
    hideydecorations!(ax1)

    Colorbar(
        ax2,
        limits = (0, 1),
        colormap = reverse(ColorSchemes.RdYlBu_9),
        label = "Relative flux",
    )
    nothing
end

end # module
