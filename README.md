# DifferentiableMetabolismCode
 
The code in this repo can be used to reproduce all of the results in the [associated paper](https://www.sciencedirect.com/science/article/pii/S1096717622001173)
Interrogating the effect of enzyme kinetics on metabolism using differentiable constraint-based models
(preprint [here](https://www.biorxiv.org/content/10.1101/2022.07.11.499575v1)).

There are three main results in the paper:
1. Sensitivity analysis of an enzyme constrained GECKO model
2. Gradient descent on measurements to fit enzyme turnover numbers
3. Sensitivity analysis on a model that incorporates full Michaelis-Menten kinetics

Note, the Julia environment used for the analyses can be recreated by activating
the `Project.toml` file, as described
[here](https://pkgdocs.julialang.org/v1/environments/). The directory `data`
contains all the data (turnover numbers, Michaelis constants, experimental data,
etc.), as well as the scripts used to parse them into the format used in this
work. The directory `model_construction` contains scripts used to download the
metabolic models, as well as add stoichiometric data for enzyme subunits, enzyme
turnover numbers, and enzyme masses to these models (called `fixed_models` in the 
context of this work).

If anything is unclear, please file an issue on this repo.
## 1. Sensitivity analysis of a GECKO model 
The directory `analyses/gecko` contains all the scripts used to generate the
associated results. In particular, the script `gecko_iml1515.jl`
can be used to run the sensitivity analysis of the flux and concentration
predictions to enzyme turnover numbers. Note, it will attempt to save the images
created in the script to another folder not contained in this repo (the repo
used to write the paper), so comment those image file saving lines out.  
## 2. Gradient descent to fit turnover numbers 
The directory `analyses/gd_gecko` contains all the scripts used to generate the
associated results. The script `cluster.jl` can be used to run the gradient
descent algorithm for a specific `master_id`, which is the ID of the respective
experimental data. Valid IDs are taken from the data source (see the paper), but
include `WT1#B1`, `WT2#B2`etc.
## 3. Sensitivity analysis of a Michaelis-Menten constrained model
The directory `analyses/crispr` contains all the scripts used to generate the
associated results. In particular, the file `michaelis_menten.jl` can be used to
generate the results, and `plot_mm_all.jl` will plot all the sensitivities. In
contrast, `plot_mm_few.jl` can be used to plot a representative selection of the
results (which is easier to  inspect).
