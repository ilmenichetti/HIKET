# =============================================================================
# run_TP3_transient_calibration.R
#
# TP3 calibration with transient pre-initialization (1917 pre-run).
# Structurally identical to run_TP2_transient_calibration.R with these changes:
#
#   1. Sources tp3_wrapper_transient.R (three-pool ASH cascade).
#   2. Three decomposition rates (alpha_A, alpha_S, alpha_H) and two transfer
#      fractions (p_S, p_H) are free, giving 10 free parameters total:
#
#        alpha_A     : Active pool decomposition rate             [log-transform]
#        alpha_S     : Slow pool decomposition rate               [log-transform]
#        alpha_H     : Humus pool decomposition rate              [log-transform]
#        p_S         : Fraction of A decomp flux entering S       [logit-transform]
#        p_H         : Fraction of S decomp flux entering H       [logit-transform]
#        beta1       : Yasso07 climate param (temp effect)        [log-transform]
#        beta2       : Yasso07 climate param (temp sq. effect)    [unconstrained]
#        gamma       : Yasso07 climate param (precip effect)      [unconstrained]
#        sigma_init  : Pre-run litter multiplier (1917)           [log-transform]
#        sigma_input : Global litter multiplier                   [log-transform]
#
#   3. Prior centres taken from TP2 where parameters are shared;
#      alpha_S = 0.10 (intermediate between alpha_A and alpha_H);
#      p_S = 0.028 (matches TP2's p_H for the first cascade step);
#      p_H = 0.50  (neutral centre for the second cascade step).
#
# Everything else — engine interface, litter_means, MCMC loop, diagnostics,
# save logic — is structurally identical to run_TP2_transient_calibration.R.
# =============================================================================

source("./Calibration_real_data_transient/calibration_engine_transient.R")
source("./Model_functions_real_data_transient/Decomposition_functions/Yasso/yasso07_wrapper_transient.R")
# No dyn.load: only the pure-R xi functions are used (compute_xi_yasso07,
# compute_xi_mean_yasso07), not yasso07_run which requires the Fortran .so
source("./Model_functions_real_data_transient/Decomposition_functions/SimpleModels/tp3_wrapper_transient.R")
source("./Model_functions_real_data_transient/input_compatibility_layer.R")

library(dplyr)


# =============================================================================
# 0.  Configuration
# =============================================================================

MODEL_NAME <- "TP3"
RUN_ID     <- format(Sys.time(), "%Y%m%d_%H%M%S")

# Read the calibration configuration (N_PLOTS_TEST, N_CHAINS, N_ITER, N_BURNIN, N_LOG)
source("./Calibration_real_data_transient/calib_config.R")

CORES_PER_CHAIN <- if (grepl("puhti|mahti", Sys.info()[["nodename"]])) {
  parallelly::availableCores()
} else {
  parallel::detectCores() - 1L
}

STEADY_STATE_YEARS <- 20L    # climate window for xi_mean (unchanged from TP2)
PREINIT_YEAR       <- 1917L  # first year of pre-run (Finland independence)
N_PREINIT_SMOOTH   <- 5L     # years averaged for J_t0_mean endpoint

DIR_LOGS   <- "./Calibration_real_data_transient/progress_logs/"
DIR_DIAG   <- file.path("./Calibration_real_data_transient/diagnostics", MODEL_NAME)
DIR_RUNS   <- "./Calibration_real_data_transient/runs/"
DIR_INPUTS <- "./Data/model_inputs/"
for (d in c(DIR_LOGS, DIR_DIAG, DIR_RUNS, DIR_INPUTS))
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
message(sprintf("  %s Transient Calibration  |  Run: %s", MODEL_NAME, RUN_ID))
message(sprintf("  Pre-run: %d -> %d (%d years, smoothing = %d yr)",
                PREINIT_YEAR, 1985L, 1985L - PREINIT_YEAR, N_PREINIT_SMOOTH))
message("=============================================================\n")


# =============================================================================
# 1.  Parameter specification
# =============================================================================
# TP3 has 10 free parameters: the 8 from TP2 plus alpha_S (Slow pool rate) and
# p_H (S-to-H transfer fraction). sigma_init and sigma_input are identical in
# role to TP2.

param_spec <- list(
  list(names = "alpha_A",     type = "log"),
  list(names = "alpha_S",     type = "log"),
  list(names = "alpha_H",     type = "log"),
  list(names = "p_S",         type = "logit"),
  list(names = "p_H",         type = "logit"),
  list(names = "beta1",       type = "log"),
  list(names = "beta2",       type = "unconstrained"),
  list(names = "gamma",       type = "unconstrained"),
  list(names = "sigma_init",  type = "log"),
  list(names = "sigma_input", type = "log")
)

transforms       <- build_transforms(param_spec)
to_original      <- transforms$to_original
to_unconstrained <- transforms$to_unconstrained
log_jacobian     <- transforms$log_jacobian
FREE_NAMES       <- transforms$param_names
N_FREE           <- transforms$n_params

# Prior centres:
#   alpha_A, alpha_H, beta1/beta2/gamma, sigma_init, sigma_input: identical to TP2.
#   alpha_S = 0.10: intermediate between alpha_A (0.73) and alpha_H (0.0015).
#   p_S     = 0.028: matches TP2's p_H (A-to-H fraction), applied here to A-to-S.
#   p_H     = 0.50: neutral centre for the S-to-H transfer fraction.
free_defaults <- c(
  alpha_A     = 0.73,
  alpha_S     = 0.10,
  alpha_H     = 0.0015,
  p_S         = 0.028,
  p_H         = 0.50,
  beta1       = 0.095,
  beta2       = -0.00014,
  gamma       = -1.21,
  sigma_init  = 1.00,
  sigma_input = 1.00
)

best_x <- to_unconstrained(free_defaults)

stopifnot(
  max(abs(to_original(best_x)[FREE_NAMES] - free_defaults[FREE_NAMES])) < 1e-10)
message("Transform round-trip: OK")


# =============================================================================
# 2.  Parameter assembly (pass-through; all TP3 params are free)
# =============================================================================

assemble_model_params <- function(p_free) p_free


# =============================================================================
# 3.  Data preparation
# =============================================================================
# Identical to TP2: total litter J_total aggregated from AWEN fractions;
# climate via Yasso07 mapping; litter_means extended for transient initialisation.

message("Loading data...")

input_raw <- read.csv("./Data/model_inputs/input_raw_monthly.csv")
site_raw  <- read.csv("./Data/model_inputs/site_raw.csv")
litter_cols <- grep("^C_nwl|^C_fwl|^C_cwl", names(input_raw), value = TRUE)

calib_plots <- site_raw$plot_id[site_raw$calib_ready]
message(sprintf("Calibration-ready plots: %d", length(calib_plots)))
input_calib <- input_raw[input_raw$plot_id %in% calib_plots, ]

n_climate_na <- sum(is.na(input_calib$temp_air))
if (n_climate_na > 0) {
  plots_with_na <- unique(input_calib$plot_id[is.na(input_calib$temp_air)])
  calib_plots   <- setdiff(calib_plots, plots_with_na)
  input_calib   <- input_raw[input_raw$plot_id %in% calib_plots, ]
  message(sprintf("After climate filter: %d plots", length(calib_plots)))
}

n_lit_na <- sum(is.na(input_calib[, litter_cols]))
if (n_lit_na > 0) {
  bad_rows <- !complete.cases(input_calib[, litter_cols])
  plots_with_lit_na <- unique(input_calib$plot_id[bad_rows])
  calib_plots <- setdiff(calib_plots, plots_with_lit_na)
  input_calib <- input_calib[input_calib$plot_id %in% calib_plots, ]
  message(sprintf("After litter filter: %d plots", length(calib_plots)))
}

Yasso07_climate    <- map_climate_yasso07(input_calib)
Yasso07_inputs_raw <- map_inputs_yasso07(input_calib)
LITTER_COLS_AWEN   <- grep("^(nwl|fwl|cwl)_", names(Yasso07_inputs_raw), value = TRUE)
TP3_inputs         <- Yasso07_inputs_raw
TP3_inputs$J_total <- rowSums(Yasso07_inputs_raw[, LITTER_COLS_AWEN])
TP3_inputs         <- TP3_inputs[, c("plot_id", "year", "J_total")]

plots_real <- as.character(unique(Yasso07_climate$plot_id))
message(sprintf("Plots after mapping: %d", length(plots_real)))

if (!is.na(N_PLOTS_TEST) && N_PLOTS_TEST < length(plots_real)) {
  set.seed(42)
  plots <- sort(sample(plots_real, N_PLOTS_TEST))
} else {
  plots <- plots_real
}
plots <- as.character(plots)

# --- Holdout split: rename, flag, reassign (is_holdout from site_raw.csv) ---
plots_all     <- plots
is_holdout_fl <- as.logical(
  site_raw$is_holdout[match(plots_all, as.character(site_raw$plot_id))])
is_holdout_fl[is.na(is_holdout_fl)] <- FALSE
plots         <- plots_all[!is_holdout_fl]   # calibration set
holdout_plots <- plots_all[ is_holdout_fl]   # held out from likelihood
message(sprintf("Holdout split: %d calibration | %d holdout",
                length(plots), length(holdout_plots)))

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
inputs_by_plot  <- split(TP3_inputs,      TP3_inputs$plot_id)

# litter_means: identical structure to TP2 (same three fields).
litter_means <- lapply(plots_real, function(pid) {
  inp  <- TP3_inputs[as.character(TP3_inputs$plot_id) == pid, ]
  n_ss <- min(STEADY_STATE_YEARS, nrow(inp))
  n_t0 <- min(N_PREINIT_SMOOTH,   nrow(inp))
  list(
    J_total_mean = mean(inp$J_total[seq_len(n_ss)]),
    J_full_mean  = mean(inp$J_total),
    J_t0_mean    = mean(inp$J_total[seq_len(n_t0)])
  )
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
# Same xi computation as TP2 (Yasso07-style, single scalar applied to all pools
# except H). Transient init via tp3_transient_init.

compute_xi_tp3_engine <- function(clim, model_params) {
  compute_xi_yasso07(
    temp_mean = clim$temp_mean,
    temp_amp  = clim$temp_amplitude,
    precip    = clim$precip,
    beta1     = model_params["beta1"],
    beta2     = model_params["beta2"],
    gamma     = model_params["gamma"]
  )
}

compute_xi_mean_tp3_engine <- function(clim_ss, model_params) {
  compute_xi_mean_yasso07(
    clim_ss = clim_ss,
    beta1   = model_params["beta1"],
    beta2   = model_params["beta2"],
    gamma   = model_params["gamma"]
  )
}

steady_state_tp3_engine <- function(model_params, lm, xi_mean) {
  tp3_transient_init(model_params, lm, xi_mean)
}

tp3_run_engine <- function(inputs, model_params, C_init, xi_array) {
  tp3_run(inputs, model_params, C_init, xi_array)
}


# =============================================================================
# 5.  Prior specification
# =============================================================================
# sigma_ppm widths follow TP2 for shared parameters.
# alpha_S and p_H receive the same widths as their analogues (alpha_A/alpha_H
# and p_S respectively).
# 95% CI for sigma_init: exp(log(1.0) ± 2*0.5) = [0.37, 2.72] (same as TP2).

sigma_ppm <- setNames(rep(1.0, N_FREE), FREE_NAMES)
sigma_ppm["alpha_A"]     <- 0.50
sigma_ppm["alpha_S"]     <- 0.50
sigma_ppm["alpha_H"]     <- 0.50
sigma_ppm["p_S"]         <- 1.00
sigma_ppm["p_H"]         <- 1.00
sigma_ppm["beta1"]       <- 0.20
sigma_ppm["beta2"]       <- 0.05
sigma_ppm["gamma"]       <- 0.30
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

  xi_array  <- compute_xi_tp3_engine(clim, sanity_params)
  n_ss      <- min(STEADY_STATE_YEARS, nrow(clim))
  xi_for_ss <- compute_xi_mean_tp3_engine(
    clim[seq_len(n_ss), , drop = FALSE], sanity_params)
  C_init    <- steady_state_tp3_engine(sanity_params, lm, xi_for_ss)
  run_out   <- tp3_run_engine(inputs, sanity_params, C_init, xi_array)

  list(pid = pid, run_out = run_out, meta = obs_meta[[pid]],
       xi_for_ss = xi_for_ss, C_init_total = sum(C_init))
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
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
for (sr in sanity_results) {
  yr      <- sr$run_out$year
  soc     <- sr$run_out$total_soc
  obs_yr  <- yr[sr$meta$idx]
  obs_soc <- sr$meta$soc_obs
  ylim_r  <- range(c(soc, obs_soc), na.rm = TRUE)
  plot(yr, soc, type = "l", col = "steelblue", lwd = 2, ylim = ylim_r,
       xlab = "Year", ylab = "SOC (tC/ha)",
       main = sprintf("Plot %s (defaults)", sr$pid))
  if (length(obs_soc) > 0)
    points(obs_yr, obs_soc, pch = 19, col = "firebrick", cex = 1.3)
  legend("topleft",
         legend = sprintf("xi=%.3f  sigma_init=%.2f",
                          sr$xi_for_ss, free_defaults["sigma_init"]),
         bty = "n", cex = 0.8)
}
mtext(sprintf("%s transient forward run at defaults  |  %s", MODEL_NAME, RUN_ID),
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
    holdout_plots      = holdout_plots,
    obs_cv             = obs_cv,
    sigma_obs_fixed    = sigma_obs_fixed,
    sigma_ppm          = sigma_ppm,
    free_defaults      = free_defaults,
    best_x             = best_x,
    STEADY_STATE_YEARS = STEADY_STATE_YEARS,
    PREINIT_YEAR       = PREINIT_YEAR,
    N_PREINIT_SMOOTH   = N_PREINIT_SMOOTH
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
  compute_xi      = compute_xi_tp3_engine,
  compute_xi_mean = compute_xi_mean_tp3_engine,
  steady_state    = steady_state_tp3_engine,
  run_model       = tp3_run_engine,
  sigma_obs_fixed = sigma_obs_fixed,
  plots           = plots,
  climate_by_plot = climate_by_plot,
  inputs_by_plot  = inputs_by_plot,
  litter_means    = litter_means,
  obs_meta        = obs_meta,
  steady_state_n  = STEADY_STATE_YEARS,
  transient_init  = TRUE
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
# 10.  Diagnostics, save, report
# =============================================================================

HIGHLIGHT <- FREE_NAMES
GROUP1    <- c("alpha_A", "alpha_S", "alpha_H", "p_S", "p_H")
GROUP2    <- c("beta1", "beta2", "gamma", "sigma_init", "sigma_input")

diag_out <- run_diagnostics(
  chain_results    = chain_results,
  to_original      = to_original,
  param_names      = FREE_NAMES,
  best_x           = best_x,
  sigma_ppm        = sigma_ppm,
  run_config       = run_config,
  mcmc_settings    = mcmc_settings,
  highlight_params = HIGHLIGHT,
  group1_params    = GROUP1,
  group2_params    = GROUP2
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
message(sprintf("  %s Transient Calibration Complete", MODEL_NAME))
message(sprintf("  Run ID:    %s", RUN_ID))
message(sprintf("  Wallclock: %.1f min", t_run / 60))
message(sprintf("  Pre-run:   %d yr  |  sigma_init centre: %.2f",
                1985L - PREINIT_YEAR, free_defaults["sigma_init"]))
message("=============================================================\n")
message("Next: run_TP3_transient_predictive.R with RUN_ID=", RUN_ID)
