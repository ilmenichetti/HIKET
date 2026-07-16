# run_config.R — shared MCMC run settings for all HIKET models.
# Edit here to change settings across all five models simultaneously.

N_PLOTS_TEST <- NA    # NA = full dataset; set to e.g. 20L for quick tests
N_CHAINS     <- 5L
N_ITER       <- 50000L
N_BURNIN     <- 5000L
N_LOG        <- 200L

# C5 (2026-07-16): down-weight the systematically-suspect 1985 (VMI8) SOC
# campaign by inflating its observation SD. The 1985->2006 mineral-soil jump is a
# measurement artefact (D3: paired 414 plots, organic stable but mineral ~doubles
# incl. deep 20-40cm), so rather than drop VMI8 we distrust it: the likelihood SD
# for every 1985 observation is multiplied by SIGMA_1985_INFL. Fixed factor, NOT
# calibrated (a 1985 offset is non-identifiable vs sigma_init). Keyed on the
# campaign YEAR (== 1985), not is_first (not every plot has VMI8). ~2x puts the
# ~30% (0.36 log) 1985-low discrepancy inside 1 SD of the base sigma_obs (~0.47).
# Sensitivity: try 1.5 / 2 / 3.
SIGMA_1985_INFL <- 2.0