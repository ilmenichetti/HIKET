# run_config.R — shared MCMC run settings for all HIKET models.
# Edit here to change settings across all five models simultaneously.

N_PLOTS_TEST <- NA    # NA = full dataset; set to e.g. 20L for quick tests
N_CHAINS     <- 5L
N_ITER       <- 50000L
N_BURNIN     <- 5000L
N_LOG        <- 200L