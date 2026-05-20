# =============================================================================
# run_TP2_calibration.R
#
# CALIBRATION STAGE for the TP2 (Two-Pool, ICBM-style) model -- HIKET.
#
# TP2 occupies the middle position in the HIKET complexity spectrum:
#   SP1 (1 pool, 6 params) -> TP2 (2 pools, 8 params) -> Yasso07 (5 pools, 20 params)
#
# MODEL (ICBM topology):
#   dA/dt = J(t) - alpha_A * xi(t) * A(t)
#   dH/dt = p_H * alpha_A * xi(t) * A(t) - alpha_H * xi(t) * H(t)
#
#   All litter inputs enter pool A (fast/active). H (humus/slow) receives
#   only transferred mass from A. Same topology as ICBM (Andren & Katterer
#   1997) but parameterised with Yasso07 notation (alpha_A, alpha_H, p_H)
#   so parameters are directly interpretable alongside Yasso07 posteriors.
#
#   xi formula: identical to SP1 and Yasso07 -- same functional form,
#   independently calibrated. Pure R implementation; no Fortran.
#
# FREE PARAMETERS (8):
#   alpha_A    : A pool decay rate (yr^-1). Log transform.
#   alpha_H    : H pool decay rate (yr^-1). Log transform.
#   p_H        : humification fraction A -> H. In (0,1), logit transform.
#   beta1      : T response, linear. Log transform.
#   beta2      : T response, quadratic. Unconstrained.
#   gamma      : precip response. Unconstrained.
#   sigma_init : init uncertainty. Log transform.
#   sigma_input: global litter scaling. Log transform.
#
# WHAT DIFFERS FROM run_Yasso07_calibration.R:
#   - No Fortran; pure R (tp2_wrapper.R).
#   - param_spec: log + logit transforms only (no stick-break).
#   - p_H calibrated as a free parameter (in Yasso07 it is fixed at 0.028).
#   - inputs_by_plot / litter_means use J_total (same as SP1).
#   - GROUP1 = structural params (alpha_A, alpha_H, p_H).
#     GROUP2 = climate + nuisance.
# =============================================================================

source("./Calibration_real_data/calibration_engine.R")
source("./Model_functions_real_data/Decomposition_functions/SimpleModels/tp2_wrapper.R")
source("./Model_functions_real_data/Decomposition_functions/Yasso/yasso07_wrapper.R")
source("./Model_functions_real_data/Decomposition_functions/Yasso/yasso07_wrapper.R") #loading the xi functions

library(dplyr)


# =============================================================================
# 0.  Configuration
# =============================================================================

MODEL_NAME <- "TP2"
RUN_ID     <- format(Sys.time(), "%Y%m%d_%H%M%S")

N_PLOTS_TEST <- NA
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
# alpha_A, alpha_H: log-transformed (strictly positive rates). Note that
#   alpha_A is a FREE parameter here, unlike Yasso07 where alpha_A is fixed.
#   This is intentional: SP1 and TP2 have no published calibrations to fix
#   from, so everything is free. Comparing TP2's posterior alpha_A to Yasso07's
#   fixed alpha_A is itself a structural question of interest.
#
# p_H: logit-transformed so it stays in (0,1). Unlike the stick-breaking
#   groups in Yasso07 (which enforce sum of fractions <= budget), p_H here is
#   a single proportion with no budget constraint -- it is the only transfer
#   coefficient in the model.
# =============================================================================

param_spec <- list(
  list(names = "alpha_A",     type = "log"),           # A pool decay rate; > 0
  list(names = "alpha_H",     type = "log"),           # H pool decay rate; > 0
  list(names = "p_H",         type = "logit"),         # humification fraction; (0,1)
  list(names = "beta1",       type = "log"),           # T linear; > 0
  list(names = "beta2",       type = "unconstrained"), # T quadratic
  list(names = "gamma",       type = "unconstrained"), # precip response
  list(names = "sigma_init",  type = "log"),
  list(names = "sigma_input", type = "log")
)

transforms       <- build_transforms(param_spec)
to_original      <- transforms$to_original
to_unconstrained <- transforms$to_unconstrained
log_jacobian     <- transforms$log_jacobian
FREE_NAMES       <- transforms$param_names
N_FREE           <- transforms$n_params

# Prior centres in PHYSICAL space.
# alpha_A, alpha_H, p_H: Yasso07 published values used as structural anchors.
#   In Yasso07, alpha_A = 0.73, alpha_H = 0.0015, p_H = 0.028. Using these
#   as TP2 starting points is physically motivated (these are well-validated
#   pool timescales for forest litter / humus). The sampler is free to move
#   away from them given the Finnish data.
# beta1, beta2, gamma: Yasso07 published defaults, same rationale as SP1.
free_defaults <- c(
  alpha_A     = 0.73,      # Yasso07 published (active pool rate, ~1.4 yr MRT)
  alpha_H     = 0.0015,    # Yasso07 published (humus pool rate, ~667 yr MRT)
  p_H         = 0.028,     # Yasso07 published (humification fraction)
  beta1       = 0.095,
  beta2       = -0.00014,
  gamma       = -1.21,
  sigma_init  = 0.10,
  sigma_input = 1.00
)

best_x <- to_unconstrained(free_defaults)

stopifnot(
  max(abs(to_original(best_x)[FREE_NAMES] - free_defaults[FREE_NAMES])) < 1e-10)
message("Transform round-trip: OK")


# =============================================================================
# 2.  Parameter assembly
# =============================================================================

assemble_model_params <- function(p_free) {
  p_free    # all parameters free; no fixed values to merge
}


# =============================================================================
# 3.  Data preparation
# =============================================================================
# Identical to SP1: climate uses map_climate_yasso07(); inputs collapse all
# AWEN litter to J_total. Only litter_means$J_total_mean is needed for
# tp2_steady_state().

message("Loading data...")

input_raw <- read.csv("./Data/model_inputs/input_raw_monthly.csv")
site_raw  <- read.csv("./Data/model_inputs/site_raw.csv")

litter_cols <- grep("^C_nwl|^C_fwl|^C_cwl", names(input_raw), value = TRUE)

calib_plots <- site_raw$plot_id[site_raw$calib_ready]
message(sprintf("Calibration-ready plots: %d", length(calib_plots)))
input_calib <- input_raw[input_raw$plot_id %in% calib_plots, ]

n_climate_na <- sum(is.na(input_calib$temp_air))
if (n_climate_na > 0) {
  message(sprintf("WARNING: %d rows with NA climate. Filtering...", n_climate_na))
  plots_with_na <- unique(input_calib$plot_id[is.na(input_calib$temp_air)])
  calib_plots   <- setdiff(calib_plots, plots_with_na)
  input_calib   <- input_raw[input_raw$plot_id %in% calib_plots, ]
  message(sprintf("After climate filter: %d plots", length(calib_plots)))
}

n_lit_na <- sum(is.na(input_calib[, litter_cols]))
if (n_lit_na > 0) {
  message(sprintf("WARNING: %d NA in litter columns. Filtering...", n_lit_na))
  bad_rows <- !complete.cases(input_calib[, litter_cols])
  plots_with_lit_na <- unique(input_calib$plot_id[bad_rows])
  calib_plots <- setdiff(calib_plots, plots_with_lit_na)
  input_calib <- input_calib[input_calib$plot_id %in% calib_plots, ]
  message(sprintf("After litter filter: %d plots", length(calib_plots)))
}

source("./Model_functions_real_data/input_compatibility_layer.R")
Yasso07_climate    <- map_climate_yasso07(input_calib)
Yasso07_inputs_raw <- map_inputs_yasso07(input_calib)
LITTER_COLS_AWEN   <- grep("^(nwl|fwl|cwl)_", names(Yasso07_inputs_raw), value = TRUE)
TP2_inputs         <- Yasso07_inputs_raw
TP2_inputs$J_total <- rowSums(Yasso07_inputs_raw[, LITTER_COLS_AWEN])
TP2_inputs         <- TP2_inputs[, c("plot_id", "year", "J_total")]

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
inputs_by_plot  <- split(TP2_inputs, TP2_inputs$plot_id)

litter_means <- lapply(plots_real, function(pid) {
  inp    <- TP2_inputs[as.character(TP2_inputs$plot_id) == pid, ]
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
# xi functions delegate to compute_xi_yasso07 / compute_xi_mean_yasso07
# (sourced from yasso07_wrapper.R) -- same formula and data flow as Yasso07,
# independently calibrated parameters.
# =============================================================================
compute_xi_tp2_engine <- function(clim, model_params) {
  compute_xi_yasso07(
    temp_mean = clim$temp_mean,
    temp_amp  = clim$temp_amplitude,
    precip    = clim$precip,
    beta1     = model_params["beta1"],
    beta2     = model_params["beta2"],
    gamma     = model_params["gamma"]
  )
}
compute_xi_mean_tp2_engine <- function(clim_ss, model_params) {
  compute_xi_mean_yasso07(
    clim_ss = clim_ss,
    beta1   = model_params["beta1"],
    beta2   = model_params["beta2"],
    gamma   = model_params["gamma"]
  )
}
steady_state_tp2_engine <- function(model_params, lm, xi_mean) {
  tp2_steady_state(model_params, lm, xi_mean)
}
tp2_run_engine <- function(inputs, model_params, C_init, xi_array) {
  tp2_run(inputs, model_params, C_init, xi_array)
}

# =============================================================================
# 5.  Prior specification
# =============================================================================
# alpha_A: log space. Default 0.73 (Yasso07 active pool). sigma_ppm = 0.5
#   gives 95% CI exp(log(0.73) +/- 1) = [0.27, 1.98] -- allows up to ~4 yr MRT.
# alpha_H: log space. Default 0.0015 (Yasso07 humus). sigma_ppm = 0.5
#   gives 95% CI [0.00055, 0.0041] -- 240 to 1800 yr MRT. Tight enough to
#   prevent the humus pool from behaving like a fast pool.
# p_H: logit space. Default 0.028. sigma_ppm = 1.0 gives physical CI
#   [inv_logit(logit(0.028)-2), inv_logit(logit(0.028)+2)] ~ [0.005, 0.14].
#   Covers a range of reasonable humification fractions.

sigma_ppm <- setNames(rep(1.0, N_FREE), FREE_NAMES)
sigma_ppm["alpha_A"]     <- 0.50   # log space; active pool MRT
sigma_ppm["alpha_H"]     <- 0.50   # log space; humus pool MRT
sigma_ppm["p_H"]         <- 1.00   # logit space; humification fraction
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
  meta   <- obs_meta[[pid]]

  xi_array  <- compute_xi_tp2_engine(clim, sanity_params)
  n_ss      <- min(STEADY_STATE_YEARS, nrow(clim))
  xi_for_ss <- compute_xi_mean_tp2_engine(
    clim[seq_len(n_ss), , drop = FALSE], sanity_params)
  C_init    <- steady_state_tp2_engine(sanity_params, lm, xi_for_ss)
  run_out   <- tp2_run_engine(inputs, sanity_params, C_init, xi_array)

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
  yr  <- sr$run_out$year;  soc <- sr$run_out$total_soc
  obs_yr  <- yr[sr$meta$idx];  obs_soc <- sr$meta$soc_obs
  ylim_r <- range(c(soc, obs_soc), na.rm = TRUE)
  plot(yr, soc, type = "l", col = "steelblue", lwd = 2, ylim = ylim_r,
       xlab = "Year", ylab = "SOC (tC/ha)",
       main = sprintf("Plot %s (defaults)", sr$pid))
  if (length(obs_soc) > 0)
    points(obs_yr, obs_soc, pch = 19, col = "firebrick", cex = 1.3)
  legend("topleft",
         legend = c(sprintf("xi = %.3f", sr$xi_for_ss),
                    sprintf("A* = %.1f  H* = %.1f", sr$C_init_total - sr$run_out$H[1],
                            sr$run_out$H[1])),
         bty = "n", cex = 0.75)
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
  compute_xi      = compute_xi_tp2_engine,
  compute_xi_mean = compute_xi_mean_tp2_engine,
  steady_state    = steady_state_tp2_engine,
  run_model       = tp2_run_engine,
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
# 9  MCMC
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
# GROUP1: structural pool parameters (the ones that change from SP1 to Yasso07).
# GROUP2: climate response + nuisance (shared parameterisation across models).
# Splitting into two PNGs keeps panels large enough to read.

HIGHLIGHT <- FREE_NAMES   # all 8 in trace plots
GROUP1    <- c("alpha_A", "alpha_H", "p_H")          # structural
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
message(sprintf("  %s Calibration Complete", MODEL_NAME))
message(sprintf("  Run ID:    %s", RUN_ID))
message(sprintf("  Wallclock: %.1f min", t_run / 60))
message(sprintf("  Plots:     %d", length(plots)))
message(sprintf("  Free params: %d", N_FREE))
message("=============================================================\n")
message("Next: run_TP2_predictive.R with RUN_ID=", RUN_ID)
