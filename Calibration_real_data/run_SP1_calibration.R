# =============================================================================
# run_SP1_calibration.R
#
# CALIBRATION STAGE for the SP1 (Single-Pool) model -- HIKET intercomparison.
#
# SP1 is the simplest model in the comparison suite: one carbon pool, one
# decay rate, no inter-pool transfers. It serves as the low-complexity anchor
# of the SP1 -> TP2 -> Yasso07 complexity spectrum. Its free parameter set
# (6 parameters) vs Yasso07's (20) quantifies how much structural complexity
# is needed to explain Finnish NFI SOC dynamics.
#
# MODEL:
#   dC/dt = J(t) - alpha * xi(t) * C(t)
#   xi(t) = exp(4*beta1*T_mean + beta2*(4*T_mean^2+2*T_amp^2)) * (1-exp(gamma*P))
#   All litter cohorts (NWL + FWL + CWL, all AWEN fractions) collapsed to
#   a single total flux J. Pure R implementation; no Fortran.
#
# FREE PARAMETERS (6):
#   alpha       : decay rate (yr^-1). Log transform.
#   beta1       : temperature response, linear. Log transform (enforced > 0).
#   beta2       : temperature response, quadratic. Unconstrained.
#   gamma       : precipitation response. Unconstrained (sign enforced by prior).
#   sigma_init  : initial-condition uncertainty. Log transform.
#   sigma_input : global litter scaling factor. Log transform.
#
# WHAT DIFFERS FROM run_Yasso07_calibration.R:
#   - No Fortran (.so); model functions are pure R (sp1_wrapper.R).
#   - param_spec has no stick-breaking groups (no compositional constraints).
#   - assemble_model_params is a pass-through (all parameters are free; none
#     are fixed from the literature -- SP1 has no published defaults to fix).
#   - inputs_by_plot carries a single J_total column, not the 12 AWEN columns.
#   - litter_means carries J_total_mean (scalar), not nwl/fwl/cwl means.
#   - HIGHLIGHT covers all 6 free parameters; no GROUP1/GROUP2 split needed.
#
# SHARED ARCHITECTURE (identical to all HIKET models):
#   - Data loading, filtering, obs_meta, obs_cv / sigma_obs_fixed computation.
#   - Prior predictive matching for sigma_ppm.
#   - make_likelihood() / run_mcmc_chains() / run_diagnostics() engine calls.
#   - Output file layout: runs/, diagnostics/SP1/, model_inputs/.
# =============================================================================

source("./Calibration_real_data/calibration_engine.R")
source("./Model_functions_real_data/Decomposition_functions/SimpleModels/sp1_wrapper.R")
source("./Model_functions_real_data/Decomposition_functions/Yasso/yasso07_wrapper.R") #loading the xi functions

library(dplyr)


# =============================================================================
# 0.  Configuration
# =============================================================================

MODEL_NAME <- "SP1"
RUN_ID     <- format(Sys.time(), "%Y%m%d_%H%M%S")

N_PLOTS_TEST <- NA    # NA = all calibration-ready plots
N_CHAINS     <- 4L
N_ITER       <- 2000L
N_BURNIN     <- 200L
N_LOG        <- 20L

# Production (Roihu): uncomment below
# N_ITER   <- 50000L
# N_BURNIN <- 10000L
# N_LOG    <- 1000L

CORES_PER_CHAIN <- if (grepl("puhti|mahti", Sys.info()[["nodename"]])) {
  parallelly::availableCores()
} else {
  parallel::detectCores() - 1L
}
# CORES_PER_CHAIN <- 90L   # Roihu production

# Steady-state climate window: same value as Yasso07 for consistency.
# All models in the intercomparison use STEADY_STATE_YEARS = 20 so that any
# difference in the steady-state init is attributable to model structure,
# not to different window choices.
STEADY_STATE_YEARS <- 20L

DIR_LOGS   <- "./Calibration_real_data/progress_logs/"
DIR_DIAG   <- file.path("./Calibration_real_data/diagnostics", MODEL_NAME)
DIR_RUNS   <- "./Calibration_real_data/runs/"
DIR_INPUTS <- "./Data/model_inputs/"
DIR_PLOTS  <- "./Calibration_real_data/plots/"
for (d in c(DIR_LOGS, DIR_DIAG, DIR_RUNS, DIR_INPUTS, DIR_PLOTS))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

run_config <- list(
  MODEL_NAME      = MODEL_NAME,
  RUN_ID          = RUN_ID,
  DIR_LOGS        = DIR_LOGS,
  DIR_DIAG        = DIR_DIAG,
  DIR_RUNS        = DIR_RUNS,
  DIR_INPUTS      = DIR_INPUTS,
  CORES_PER_CHAIN = CORES_PER_CHAIN
)

mcmc_settings <- list(
  N_CHAINS = N_CHAINS, N_ITER = N_ITER,
  N_BURNIN = N_BURNIN, N_LOG  = N_LOG
)

message("=============================================================")
message(sprintf("  %s Calibration  |  Run: %s", MODEL_NAME, RUN_ID))
message("=============================================================\n")


# =============================================================================
# 1.  Parameter specification
# =============================================================================
# SP1 has no fixed parameters (unlike Yasso07 which fixes 6 decomposition
# rates and p_H from the literature). All 4 model parameters and 2 nuisance
# parameters are free.
#
# No stick-breaking groups: SP1 has no compositional constraints (no fractions
# that must sum to a budget). Only log and unconstrained transforms are used.
#
# beta1 enforced > 0 (log): in Finnish conditions a positive linear temperature
# response is physically unambiguous. beta2 and gamma are unconstrained; their
# sign is controlled by prior placement, not by the transform.
# =============================================================================

param_spec <- list(
  list(names = "alpha",       type = "log"),           # decay rate; strictly positive
  list(names = "beta1",       type = "log"),           # T linear response; strictly positive
  list(names = "beta2",       type = "unconstrained"), # T quadratic; can be negative
  list(names = "gamma",       type = "unconstrained"), # precip response; sign via prior
  list(names = "sigma_init",  type = "log"),           # init uncertainty; strictly positive
  list(names = "sigma_input", type = "log")            # litter scaling; strictly positive
)

transforms       <- build_transforms(param_spec)
to_original      <- transforms$to_original
to_unconstrained <- transforms$to_unconstrained
log_jacobian     <- transforms$log_jacobian
FREE_NAMES       <- transforms$param_names
N_FREE           <- transforms$n_params

# Prior centres in PHYSICAL space.
# alpha: prior centre at MRT ~11 years at xi=1, consistent with a mixed-pool
#   average for Finnish mineral upland forest. Computed from rough target:
#   C* = J/(alpha*xi) ~ 70 tC/ha at J~2.4 tC/ha/yr, xi~0.35 -> alpha ~ 0.098.
# beta1, beta2, gamma: Yasso07 published defaults (Tuomi et al. 2009, 2011).
#   Using the same starting point as Yasso07 so any post-calibration difference
#   reflects data-driven adaptation, not prior placement differences.
# sigma_init, sigma_input: weakly informative, same as Yasso07.
free_defaults <- c(
  alpha       = 0.09,      # MRT ~11 yr at xi=1; see note above
  beta1       = 0.095,     # Yasso07 published default
  beta2       = -0.00014,  # Yasso07 published default (near zero)
  gamma       = -1.21,     # Yasso07 published default
  sigma_init  = 0.10,      # ~10% initial-condition uncertainty
  sigma_input = 1.00       # no a priori litter scaling
)

best_x <- to_unconstrained(free_defaults)

stopifnot(
  max(abs(to_original(best_x)[FREE_NAMES] - free_defaults[FREE_NAMES])) < 1e-10)
message("Transform round-trip: OK")


# =============================================================================
# 2.  Parameter assembly
# =============================================================================
# SP1 has no fixed parameters and no Fortran routine requiring a specific
# parameter ordering. assemble_model_params simply passes through the free
# parameter vector. It is kept as an explicit function (rather than being
# skipped) so the engine interface is identical across all models.

assemble_model_params <- function(p_free) {
  p_free    # all parameters are free; no fixed values to merge
}


# =============================================================================
# 3.  Data preparation
# =============================================================================
# Data loading and filtering are identical to run_Yasso07_calibration.R.
# The only differences are:
#   (a) inputs_by_plot carries J_total (scalar per year) rather than 12 AWEN
#       columns. The mapping collapses all litter cohorts after the Yasso07
#       mapping, so the same input_raw_monthly.csv is used.
#   (b) litter_means carries J_total_mean (scalar) rather than nwl/fwl/cwl means.
#   (c) climate_by_plot is identical to Yasso07 (same xi formula, same climate
#       columns: temp_mean, temp_amplitude, precip).
# =============================================================================

message("Loading data...")

input_raw <- read.csv("./Data/model_inputs/input_raw_monthly.csv")
site_raw  <- read.csv("./Data/model_inputs/site_raw.csv")

litter_cols <- grep("^C_nwl|^C_fwl|^C_cwl", names(input_raw), value = TRUE)

# Filter (a): calib_ready plots only
calib_plots <- site_raw$plot_id[site_raw$calib_ready]
message(sprintf("Calibration-ready plots: %d", length(calib_plots)))
input_calib <- input_raw[input_raw$plot_id %in% calib_plots, ]

# Filter (b): drop plots with NA climate
n_climate_na <- sum(is.na(input_calib$temp_air))
if (n_climate_na > 0) {
  message(sprintf("WARNING: %d rows with NA climate. Filtering...", n_climate_na))
  plots_with_na <- unique(input_calib$plot_id[is.na(input_calib$temp_air)])
  calib_plots   <- setdiff(calib_plots, plots_with_na)
  input_calib   <- input_raw[input_raw$plot_id %in% calib_plots, ]
  message(sprintf("After climate filter: %d plots", length(calib_plots)))
}

# Filter (c): drop plots with NA in any litter column
n_lit_na <- sum(is.na(input_calib[, litter_cols]))
if (n_lit_na > 0) {
  message(sprintf("WARNING: %d NA in litter columns. Filtering...", n_lit_na))
  bad_rows <- !complete.cases(input_calib[, litter_cols])
  plots_with_lit_na <- unique(input_calib$plot_id[bad_rows])
  calib_plots <- setdiff(calib_plots, plots_with_lit_na)
  input_calib <- input_calib[input_calib$plot_id %in% calib_plots, ]
  message(sprintf("After litter filter: %d plots", length(calib_plots)))
}

# Climate mapping: same as Yasso07 (temp_mean, temp_amplitude, precip per year).
# SP1's xi formula uses the same three climate columns.
source("./Model_functions_real_data/input_compatibility_layer.R")
Yasso07_climate <- map_climate_yasso07(input_calib)

# Input mapping: collapse all 12 AWEN litter columns to a single total J.
# This is done AFTER the Yasso07 mapping so we reuse the existing compatibility
# layer without modification.
Yasso07_inputs_raw <- map_inputs_yasso07(input_calib)
LITTER_COLS_AWEN   <- grep("^(nwl|fwl|cwl)_", names(Yasso07_inputs_raw), value = TRUE)
SP1_inputs <- Yasso07_inputs_raw
SP1_inputs$J_total <- rowSums(Yasso07_inputs_raw[, LITTER_COLS_AWEN])
SP1_inputs <- SP1_inputs[, c("plot_id", "year", "J_total")]

plots_real <- as.character(unique(Yasso07_climate$plot_id))
message(sprintf("Plots after mapping: %d", length(plots_real)))

if (!is.na(N_PLOTS_TEST) && N_PLOTS_TEST < length(plots_real)) {
  set.seed(42)
  plots <- sort(sample(plots_real, N_PLOTS_TEST))
} else {
  plots <- plots_real
}
plots <- as.character(plots)

SOC_obs_all <- input_calib %>%
  filter(!is.na(soc_obs_tCha)) %>%
  select(plot_id, year, month, soc_obs_tCha) %>%
  arrange(plot_id, year, month) %>%
  group_by(plot_id) %>%
  mutate(obs_rank = row_number()) %>%
  ungroup()

obs_cv          <- sd(SOC_obs_all$soc_obs_tCha) / mean(SOC_obs_all$soc_obs_tCha)
sigma_obs_fixed <- obs_cv
message(sprintf("SOC obs: %d  |  obs CV: %.3f  |  sigma_obs fixed: %.3f",
                nrow(SOC_obs_all), obs_cv, sigma_obs_fixed))

climate_by_plot <- split(Yasso07_climate, Yasso07_climate$plot_id)
inputs_by_plot  <- split(SP1_inputs, SP1_inputs$plot_id)

# litter_means: scalar J_total_mean per plot, averaged over the steady-state
# window. Same window width (STEADY_STATE_YEARS = 20) as all other models.
litter_means <- lapply(plots_real, function(pid) {
  inp    <- SP1_inputs[as.character(SP1_inputs$plot_id) == pid, ]
  n_ss   <- min(STEADY_STATE_YEARS, nrow(inp))
  inp_ss <- inp[seq_len(n_ss), ]
  list(J_total_mean = mean(inp_ss$J_total))
})
names(litter_means) <- plots_real

obs_meta <- lapply(plots_real, function(pid) {
  clim     <- Yasso07_climate[as.character(Yasso07_climate$plot_id) == pid, ]
  obs_plot <- SOC_obs_all[as.character(SOC_obs_all$plot_id) == pid, ]
  list(
    idx      = match(obs_plot$year, clim$year),
    soc_obs  = obs_plot$soc_obs_tCha,
    is_first = obs_plot$obs_rank == 1L
  )
})
names(obs_meta) <- plots_real

plot_info <- lapply(plots_real, function(pid) {
  row <- site_raw[as.character(site_raw$plot_id) == pid, ]
  as.list(row[1, , drop = FALSE])
})
names(plot_info) <- plots_real

run_config$n_plots <- length(plots)


# =============================================================================
# 4.  Model interface wrappers
# =============================================================================
# Glue functions adapting sp1_wrapper.R functions to the engine interface.
# The xi functions delegate to compute_xi_yasso07 / compute_xi_mean_yasso07
# (sourced from yasso07_wrapper.R) -- same formula and data flow as Yasso07,
# independently calibrated parameters.
# =============================================================================

compute_xi_sp1_engine <- function(clim, model_params) {
  compute_xi_yasso07(
    temp_mean = clim$temp_mean,
    temp_amp  = clim$temp_amplitude,
    precip    = clim$precip,
    beta1     = model_params["beta1"],
    beta2     = model_params["beta2"],
    gamma     = model_params["gamma"]
  )
}

compute_xi_mean_sp1_engine <- function(clim_ss, model_params) {
  compute_xi_mean_yasso07(
    clim_ss = clim_ss,
    beta1   = model_params["beta1"],
    beta2   = model_params["beta2"],
    gamma   = model_params["gamma"]
  )
}

steady_state_sp1_engine <- function(model_params, lm, xi_mean) {
  sp1_steady_state(model_params, lm, xi_mean)
}

sp1_run_engine <- function(inputs, model_params, C_init, xi_array) {
  sp1_run(inputs, model_params, C_init, xi_array)
}


# =============================================================================
# 5.  Prior specification
# =============================================================================
# alpha: sigma_ppm = 0.5 on log scale (exp(log(0.09) +/- 2*0.5) = [0.033, 0.24]).
#   Wide enough to span MRTs from ~4 to 30 years.
# beta1: sigma_ppm = 0.20 (log space), matching Yasso07 for comparability.
# beta2, gamma: sigma_ppm matching Yasso07 (tight to keep xi well-behaved).
# sigma_init, sigma_input: same as Yasso07.
# =============================================================================

sigma_ppm <- setNames(rep(1.0, N_FREE), FREE_NAMES)
sigma_ppm["alpha"]       <- 0.50   # log space; wide prior on MRT
sigma_ppm["beta1"]       <- 0.20   # log space; Tuomi et al. 2008/2009
sigma_ppm["beta2"]       <- 0.05   # physical space; keep quadratic bounded
sigma_ppm["gamma"]       <- 0.30   # physical space; Finnish precip range
sigma_ppm["sigma_init"]  <- 0.50
sigma_ppm["sigma_input"] <- 0.50

stopifnot(length(sigma_ppm) == N_FREE, all(sigma_ppm > 0),
          all(names(sigma_ppm) == FREE_NAMES))

message("\nPrior widths (unconstrained space):")
print(round(sigma_ppm, 3))


# =============================================================================
# 6.  Pre-MCMC sanity check
# =============================================================================

message("\n-------------------------------------------------------------")
message("Pre-MCMC sanity: forward run at defaults...")

set.seed(2025)
sanity_pids   <- sort(sample(plots, min(4, length(plots))))
sanity_params <- assemble_model_params(free_defaults)

sanity_results <- lapply(sanity_pids, function(pid) {
  clim   <- climate_by_plot[[pid]]
  inputs <- inputs_by_plot[[pid]]
  lm     <- litter_means[[pid]]
  meta   <- obs_meta[[pid]]

  xi_array  <- compute_xi_sp1_engine(clim, sanity_params)
  n_ss      <- min(STEADY_STATE_YEARS, nrow(clim))
  xi_for_ss <- compute_xi_mean_sp1_engine(
    clim[seq_len(n_ss), , drop = FALSE], sanity_params)
  C_init    <- steady_state_sp1_engine(sanity_params, lm, xi_for_ss)
  run_out   <- sp1_run_engine(inputs, sanity_params, C_init, xi_array)

  list(pid = pid, run_out = run_out, meta = meta, xi_for_ss = xi_for_ss,
       C_init_total = sum(C_init))
})

forward_run_pass <- TRUE
for (sr in sanity_results) {
  hat <- sr$run_out$total_soc[sr$meta$idx]
  if (any(!is.finite(hat)) || any(hat <= 0)) { forward_run_pass <- FALSE; break }
  if (length(sr$meta$soc_obs) > 0) {
    ratio <- hat / sr$meta$soc_obs
    if (any(ratio < 0.1 | ratio > 10, na.rm = TRUE)) {
      forward_run_pass <- FALSE
      message(sprintf("  WARNING: plot %s out of order-of-magnitude at defaults",
                      sr$pid))
    }
  }
}

PX_PER_IN <- 150L
sanity_png <- file.path(DIR_DIAG,
                        sprintf("%s_forward_at_defaults_%s.png", MODEL_NAME, RUN_ID))
png(sanity_png, width = 10L * PX_PER_IN, height = 8L * PX_PER_IN, res = PX_PER_IN)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for (sr in sanity_results) {
  yr  <- sr$run_out$year
  soc <- sr$run_out$total_soc
  obs_yr  <- yr[sr$meta$idx]
  obs_soc <- sr$meta$soc_obs
  ylim_r <- range(c(soc, obs_soc), na.rm = TRUE)
  plot(yr, soc, type = "l", col = "steelblue", lwd = 2, ylim = ylim_r,
       xlab = "Year", ylab = "SOC (tC/ha)",
       main = sprintf("Plot %s (defaults)", sr$pid))
  if (length(obs_soc) > 0)
    points(obs_yr, obs_soc, pch = 19, col = "firebrick", cex = 1.3)
  legend("topleft",
         legend = sprintf("xi(mean clim) = %.3f", sr$xi_for_ss),
         bty = "n", cex = 0.8)
}
mtext(sprintf("%s forward run at defaults  |  %s", MODEL_NAME, RUN_ID),
      side = 3, outer = TRUE, line = -1.5, cex = 0.95, font = 2)
dev.off()
message(sprintf("Sanity plot: %s", sanity_png))
message(sprintf("Forward-run sanity: %s",
                if (forward_run_pass) "PASS" else "WARN -- inspect plot"))



# =============================================================================
# 7.  Save inputs
# =============================================================================

save_inputs(
  inputs_list = list(
    climate_by_plot    = climate_by_plot,
    inputs_by_plot     = inputs_by_plot,
    litter_means       = litter_means,
    obs_meta           = obs_meta,
    plot_info          = plot_info,
    plots              = plots,
    plots_real         = plots_real,
    obs_cv             = obs_cv,
    sigma_obs_fixed    = sigma_obs_fixed,
    sigma_ppm          = sigma_ppm,
    free_defaults      = free_defaults,
    best_x             = best_x,
    STEADY_STATE_YEARS = STEADY_STATE_YEARS
  ),
  run_config = run_config
)


# =============================================================================
# 8.  Build likelihood and prior
# =============================================================================

ll_fn <- make_likelihood(
  n_cores         = CORES_PER_CHAIN,
  to_original     = to_original,
  log_jacobian    = log_jacobian,
  assemble_params = assemble_model_params,
  compute_xi      = compute_xi_sp1_engine,
  compute_xi_mean = compute_xi_mean_sp1_engine,
  steady_state    = steady_state_sp1_engine,
  run_model       = sp1_run_engine,
  sigma_obs_fixed = sigma_obs_fixed,
  plots           = plots,
  climate_by_plot = climate_by_plot,
  inputs_by_plot  = inputs_by_plot,
  litter_means    = litter_means,
  obs_meta        = obs_meta,
  steady_state_n  = STEADY_STATE_YEARS
)

ll_check <- ll_fn(best_x)
message(sprintf("\nLikelihood at defaults: %.3f", ll_check))
if (!is.finite(ll_check)) stop("FAIL: likelihood not finite at defaults.")
message("Likelihood check: PASS\n")

pre_mcmc_sanity <- list(
  ll_at_defaults     = ll_check,
  n_plots            = length(plots),
  n_obs              = nrow(SOC_obs_all),
  steady_state_years = STEADY_STATE_YEARS,
  forward_run_pass   = forward_run_pass
)

prior <- createPrior(
  density = function(x) sum(dnorm(x, mean = best_x, sd = sigma_ppm, log = TRUE)),
  sampler = function(n = 1) {
    m <- matrix(0, nrow = n, ncol = N_FREE)
    for (j in seq_len(N_FREE))
      m[, j] <- rnorm(n, mean = best_x[j], sd = sigma_ppm[j] * 0.1)
    m
  },
  lower = rep(-Inf, N_FREE),
  upper = rep( Inf, N_FREE)
)


# =============================================================================
# 9.  MCMC
# =============================================================================

t_run <- system.time({
  mcmc_out <- run_mcmc_chains(
    ll_fn         = ll_fn,
    prior         = prior,
    best_x        = best_x,
    param_names   = FREE_NAMES,
    mcmc_settings = mcmc_settings,
    run_config    = run_config)
})[["elapsed"]]

chain_results <- mcmc_out$chain_results
chain_health  <- mcmc_out$chain_health

message(sprintf("\nAll chains complete. Wallclock: %.1f min\n", t_run / 60))


# =============================================================================
# 11.  Diagnostics, save, report
# =============================================================================
# With only 6 free parameters, all fit comfortably in one marginal PNG
# without group splitting. No HIGHLIGHT subset needed: plot all.

HIGHLIGHT <- FREE_NAMES   # all 6 for SP1

diag_out <- run_diagnostics(
  chain_results    = chain_results,
  to_original      = to_original,
  param_names      = FREE_NAMES,
  best_x           = best_x,
  sigma_ppm        = sigma_ppm,
  run_config       = run_config,
  mcmc_settings    = mcmc_settings,
  highlight_params = HIGHLIGHT,
  group1_params    = NULL,    # no split needed; 6 params fit in one grid
  group2_params    = NULL
)

save_results(
  chain_results   = chain_results,
  all_phys        = diag_out$all_phys,
  run_config      = run_config,
  mcmc_settings   = mcmc_settings,
  sigma_ppm       = sigma_ppm,
  sigma_obs_fixed = sigma_obs_fixed,
  t_elapsed       = t_run,
  gr              = diag_out$gr,
  ess_vals        = diag_out$ess_vals
)

write_calibration_report(
  run_config      = run_config,
  mcmc_settings   = mcmc_settings,
  pre_mcmc_sanity = pre_mcmc_sanity,
  chain_health    = chain_health,
  diag_out        = diag_out,
  param_names     = FREE_NAMES,
  best_x          = best_x,
  sigma_ppm       = sigma_ppm,
  sigma_obs_fixed = sigma_obs_fixed,
  t_run           = t_run
)

message("\n=============================================================")
message(sprintf("  %s Calibration Complete", MODEL_NAME))
message(sprintf("  Run ID:    %s", RUN_ID))
message(sprintf("  Wallclock: %.1f min", t_run / 60))
message(sprintf("  Plots:     %d", length(plots)))
message(sprintf("  Free params: %d", N_FREE))
message("=============================================================\n")
message("Next: run_SP1_predictive.R with RUN_ID=", RUN_ID)
