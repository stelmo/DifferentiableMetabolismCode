
# _concs = JSON.parsefile(joinpath("data", "metabolites", "bennet2016", "metabolites.json"))
# concs = Dict(k => v for (k, v) in _concs if k in metabolites(pmodel))
# concs["pi_c"] = 20e-3
#from comments
concs = Dict(
    "ppcoa_c" => 5.3e-6,
    "gdp_c" => 0.00068,
    "lys__L_c" => 0.00041,
    "nadh_c" => 8.3e-5,
    "phpyr_c" => 9.0e-5,
    "his__L_c" => 6.8e-5,
    "pep_c" => 0.00018,
    "gmp_c" => 2.4e-5,
    "5drib_c" => 0.0003,
    "amp_c" => 0.00028,
    "succoa_c" => 0.00023,
    "ctp_c" => 0.0027,
    "dttp_c" => 0.0046,
    "utp_c" => 0.0083,
    "pi_c" => 0.02,
    "udp_c" => 0.0018,
    "uacgam_c" => 0.0092,
    "gtp_c" => 0.0049,
    "quln_c" => 1.2e-5,
    "dtdp_c" => 0.00038,
    "amet_c" => 0.00018,
    "tyr__L_c" => 2.9e-5,
    "accoa_c" => 0.00061,
    "fum_c" => 0.00012,
    "dhor__S_c" => 1.2e-5,
    "thr__L_c" => 0.00018,
    "succ_c" => 0.00057,
    "ribflv_c" => 1.9e-5,
    "dctp_c" => 3.5e-5,
    "fad_c" => 0.00017,
    "atp_c" => 0.0096,
    "dhap_c" => 0.00037,
    "pro__L_c" => 0.00039,
    "akg_c" => 0.00044,
    "skm_c" => 1.4e-5,
    "phe__L_c" => 1.8e-5,
    "acon_C_c" => 1.6e-5,
    "val__L_c" => 0.004,
    "glyc3p_c" => 4.9e-5,
    "citr__L_c" => 0.0014,
    "4hbz_c" => 5.2e-5,
    "datp_c" => 1.6e-5,
    "cmp_c" => 0.00036,
    "uri_c" => 0.0021,
    "glu__L_c" => 0.096,
    "prpp_c" => 0.00026,
    "6pgl_c" => 0.001,
    "cit_c" => 0.002,
    "coa_c" => 0.0014,
    "nad_c" => 0.0026,
    # "mal__L_c"  => 0.0017, # infeasible
    "gam6p_c" => 0.0012,
    "6pgc_c" => 0.0038,
    "acgam1p_c" => 8.2e-5,
    "histd_c" => 1.3e-5,
    "ser__L_c" => 6.8e-5,
    "malcoa_c" => 3.5e-5,
    "arg__L_c" => 0.00057,
    "fmn_c" => 5.4e-5,
    "gln__L_c" => 0.0038,
    "trp__L_c" => 1.2e-5,
    "nadph_c" => 0.00012,
    "acorn_c" => 4.3e-5,
    "cbasp_c" => 0.00059,
    "ala__L_c" => 0.0026,
    "ade_c" => 1.5e-6,
    "asn__L_c" => 0.00051,
    "hcys__L_c" => 0.00037,
    "aps_c" => 6.6e-6,
    "imp_c" => 0.00027,
    "anth_c" => 3.5e-6,
    "nadp_c" => 2.1e-6,
    "adp_c" => 0.00056,
    "3pg_c" => 0.0015,
    "met__L_c" => 0.00015,
)

mmdf = max_min_driving_force(
    pmodel,
    rid_dg0,
    CPLEX.Optimizer;
    flux_solution = loopless_sol,
    proton_ids = ["h_c", "h_e", "h_p"],
    water_ids = ["h2o_c", "h2o_e", "h2o_p"],
    constant_concentrations = concs,
    concentration_lb = 1e-7,
    concentration_ub = 0.1,
    ignore_reaction_ids = ["Htex", "H2Otex", "H2Otpp"],
)
max_driving_force = mmdf.mmdf

idxs = Int.(indexin(keys(rid_dg0), reactions(pmodel)))

swap_obj(model, optmodel) = begin
    @constraint(optmodel, optmodel[:mmdf] >= max_driving_force)
    @objective(optmodel, Min, sum(optmodel[:dgrs][i]^2 for i in idxs))
end

mmdf = max_min_driving_force(
    pmodel,
    rid_dg0,
    CPLEX.Optimizer;
    flux_solution = loopless_sol,
    proton_ids = ["h_c", "h_e", "h_p"],
    water_ids = ["h2o_c", "h2o_e", "h2o_p"],
    constant_concentrations = concs,
    concentration_lb = 1e-7,
    concentration_ub = 0.1,
    ignore_reaction_ids = ["Htex", "H2Otex", "H2Otpp"],
    modifications = [swap_obj],
)
mmdf
mmdf.concentrations["mal__L_c"]
mmdf.concentrations
