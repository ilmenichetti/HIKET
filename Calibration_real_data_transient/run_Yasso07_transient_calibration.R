# =============================================================================
# run_Yasso07_transient_calibration.R
#
# Yasso07 calibration with transient pre-initialization (1917 pre-run).
# Differs from run_Yasso07_calibration.R in the following sections only:
#
#   1. Sources calibration_engine_transient.R and yasso07_wrapper_transient.R.
#   2. PREINIT_YEAR, N_PREINIT_SMOOTH added to config.
#   3. assemble_model_params: includes sigma_init in model_params so that
#      yasso07_transient_init can read it (sigma_init was previously only
#      extracted in the likelihood for sd_vec; it is now a model parameter).
#   4. litter_means: gains nwl/fwl/cwl_full_mean and nwl/fwl/cwl_t0_mean fields.
#   5. free_defaults["sigma_init"] = 1.00 (not 0.10).
#   6. steady_state engine wrapper calls yasso07_transient_init.
#   7. make_likelihood(..., transient_init = TRUE).
# =============================================================================

source("./Calibration_real_data_transient/calibration_engine_transient.R")
source("./Model_functions_real_data_transient/Decomposition_functions/Yasso/yasso07_wrapper_transient.R")
source("./Model_functions_real_data_transient/input_compatibility_layer.R")

dyn.load("./Model_functions_real_data_transient/Decomposition_functions/Yasso/yasso07.so")

library(dplyr)


# =============================================================================
# 0.  Configuration
# =============================================================================

MODEL_NAME       <- "Yasso07"
RUN_ID           <- format(Sys.time(), "%Y%m%d_%H%M%S")

# read the calibration configuration (N_PLOTS_TEST, N_CHAINS, N_ITER, N_BURNIN, N_LOG)
source("./Calibration_real_data_transient/calib_config.R")


CORES_PER_CHAIN <- if (grepl("puhti|mahti", Sys.info()[["nodename"]])) {
  parallelly::availableCores()
} else {
  parallel::detectCores() - 1L
}

STEADY_STATE_YEARS <- 20L
PREINIT_YEAR       <- 1917L
N_PREINIT_SMOOTH   <- 5L

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
# 1.  Parameter specification  (identical to non-transient Yasso07)
# =============================================================================

p_default        <- YASSO07_DEFAULT_PARAMS
FIXED_RATE_NAMES <- c("alpha_A","alpha_W","alpha_E","alpha_N","p_H","alpha_H")
fixed_rates      <- p_default[FIXED_RATE_NAMES]
B                <- 1.0 - fixed_rates["p_H"]

FRAC_COL_A <- c("p_AW","p_AE","p_AN")
FRAC_COL_W <- c("p_WA","p_WE","p_WN")
FRAC_COL_E <- c("p_EA","p_EW","p_EN")
FRAC_COL_N <- c("p_NA","p_NW","p_NE")

param_spec <- list(
  list(names = FRAC_COL_A,    type = "stick_break",   budget = B),
  list(names = FRAC_COL_W,    type = "stick_break",   budget = B),
  list(names = FRAC_COL_E,    type = "stick_break",   budget = B),
  list(names = FRAC_COL_N,    type = "stick_break",   budget = B),
  list(names = "beta1",       type = "log"),
  list(names = "beta2",       type = "unconstrained"),
  list(names = "gamma",       type = "unconstrained"),
  list(names = "delta1",      type = "unconstrained"),
  list(names = "delta2",      type = "log"),
  list(names = "r",           type = "log"),
  list(names = "sigma_init",  type = "log"),
  list(names = "sigma_input", type = "log")
)

transforms       <- build_transforms(param_spec)
to_original      <- transforms$to_original
to_unconstrained <- transforms$to_unconstrained
log_jacobian     <- transforms$log_jacobian
FREE_NAMES       <- transforms$param_names
N_FREE           <- transforms$n_params

MODEL_FREE_NAMES <- FREE_NAMES[!FREE_NAMES %in% c("sigma_init", "sigma_input")]

free_defaults <- c(
  p_default[MODEL_FREE_NAMES],
  sigma_init  = 1.00,   # CHANGED: prior centre = 1.0 (same productivity as today)
  sigma_input = 1.00
)
free_defaults["r"] <- abs(free_defaults["r"])

best_x <- to_unconstrained(free_defaults)
stopifnot(
  max(abs(to_original(best_x)[FREE_NAMES] - free_defaults[FREE_NAMES])) < 1e-10)
message("Transform round-trip: OK")


# =============================================================================
# 2.  Parameter assembly
# =============================================================================
# CHANGED: sigma_init added to model_params so yasso07_transient_init can
# read it. In the original script sigma_init was only in p_free (extracted
# in the likelihood for sd_vec). In the transient version it must enter the
# model to scale the 1917 litter anchor.

assemble_model_params <- function(p_free) {
  yasso_params <- c(fixed_rates, p_free[MODEL_FREE_NAMES])
  yasso_params <- yasso_params[names(p_default)]
  c(yasso_params,
    sigma_input = unname(p_free["sigma_input"]),
    sigma_init  = unname(p_free["sigma_init"]))   # NEW
}


# =============================================================================
# 3.  Data preparation
# =============================================================================

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
  bad_rows          <- !complete.cases(input_calib[, litter_cols])
  plots_with_lit_na <- unique(input_calib$plot_id[bad_rows])
  calib_plots       <- setdiff(calib_plots, plots_with_lit_na)
  input_calib       <- input_calib[input_calib$plot_id %in% calib_plots, ]
  message(sprintf("After litter filter: %d plots", length(calib_plots)))
}

Yasso07_climate <- map_climate_yasso07(input_calib)
Yasso07_inputs  <- map_inputs_yasso07(input_calib)

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
message(sprintf("SOC observations: %d  |  obs CV: %.3f  |  sigma_obs fixed: %.3f",
                nrow(SOC_obs_all), obs_cv, sigma_obs_fixed))

climate_by_plot <- split(Yasso07_climate, Yasso07_climate$plot_id)
inputs_by_plot  <- split(Yasso07_inputs,  Yasso07_inputs$plot_id)

AWEN_NWL <- c("nwl_A","nwl_W","nwl_E","nwl_N")
AWEN_FWL <- c("fwl_A","fwl_W","fwl_E","fwl_N")
AWEN_CWL <- c("cwl_A","cwl_W","cwl_E","cwl_N")

# litter_means: EXTENDED for transient initialization.
#   Original fields: nwl_mean, fwl_mean, cwl_mean (steady-state window means)
#   New fields:
#     *_full_mean : mean over full time series (1917 anchor)
#     *_t0_mean   : mean over first N_PREINIT_SMOOTH years (1985 endpoint)
litter_means <- lapply(plots_real, function(pid) {
  inp    <- Yasso07_inputs[as.character(Yasso07_inputs$plot_id) == pid, ]
  n_ss   <- min(STEADY_STATE_YEARS,  nrow(inp))
  n_t0   <- min(N_PREINIT_SMOOTH,    nrow(inp))
  inp_ss <- inp[seq_len(n_ss), ]
  inp_t0 <- inp[seq_len(n_t0), ]
  list(
    # Original: steady-state window means
    nwl_mean      = colMeans(inp_ss[, AWEN_NWL]),
    fwl_mean      = colMeans(inp_ss[, AWEN_FWL]),
    cwl_mean      = colMeans(inp_ss[, AWEN_CWL]),
    # New: full-series means for 1917 anchor
    nwl_full_mean = colMeans(inp[, AWEN_NWL]),
    fwl_full_mean = colMeans(inp[, AWEN_FWL]),
    cwl_full_mean = colMeans(inp[, AWEN_CWL]),
    # New: first-N-year means for 1985 endpoint
    nwl_t0_mean   = colMeans(inp_t0[, AWEN_NWL]),
    fwl_t0_mean   = colMeans(inp_t0[, AWEN_FWL]),
    cwl_t0_mean   = colMeans(inp_t0[, AWEN_CWL])
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

compute_xi_yasso07_engine <- function(clim, model_params) {
  compute_xi_yasso07(
    temp_mean = clim$temp_mean,
    temp_amp  = clim$temp_amplitude,
    precip    = clim$precip,
    beta1     = model_params["beta1"],
    beta2     = model_params["beta2"],
    gamma     = model_params["gamma"]
  )
}

compute_xi_mean_yasso07_engine <- function(clim_ss, model_params) {
  compute_xi_mean_yasso07(
    clim_ss = clim_ss,
    beta1   = model_params["beta1"],
    beta2   = model_params["beta2"],
    gamma   = model_params["gamma"]
  )
}

# CHANGED: yasso07_transient_init replaces yasso07_steady_state
steady_state_yasso07_engine <- function(model_params, lm, xi_mean) {
  yasso07_transient_init(model_params, lm, xi_mean)
}

yasso07_run_engine <- function(inputs, model_params, C_init, xi_array) {
  si <- model_params["sigma_input"]
  input_scaled        <- inputs
  input_scaled[, c("nwl_A","nwl_W","nwl_E","nwl_N",
                   "fwl_A","fwl_W","fwl_E","fwl_N",
                   "cwl_A","cwl_W","cwl_E","cwl_N")] <-
    inputs[, c("nwl_A","nwl_W","nwl_E","nwl_N",
               "fwl_A","fwl_W","fwl_E","fwl_N",
               "cwl_A","cwl_W","cwl_E","cwl_N")] * si
  yasso07_run(
    input_df = input_scaled,
    params   = model_params[YASSO07_PARAM_NAMES],
    C_init   = C_init,
    xi_array = xi_array
  )
}


# =============================================================================
# 5.  Prior specification
# =============================================================================

sigma_ppm <- setNames(rep(1.0, N_FREE), FREE_NAMES)
sigma_ppm["beta1"]       <- 0.20
sigma_ppm["beta2"]       <- 0.05
sigma_ppm["gamma"]       <- 0.30
sigma_ppm["delta1"]      <- 0.15
sigma_ppm["delta2"]      <- 0.10
sigma_ppm["r"]           <- 0.015
sigma_ppm["sigma_init"]  <- 0.50   # physical 95% CI ~[0.37, 2.72] around 1.0
sigma_ppm["sigma_input"] <- 0.50

stopifnot(length(sigma_ppm) == N_FREE, all(sigma_ppm > 0),
          all(names(sigma_ppm) == FREE_NAMES))

message("\nPrior widths (unconstrained space):")
print(round(sigma_ppm[c("beta1","beta2","gamma","delta1","delta2","r",
                         "sigma_init","sigma_input")], 4))


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

  xi_array  <- compute_xi_yasso07_engine(clim, sanity_params)
  n_ss      <- min(STEADY_STATE_YEARS, nrow(clim))
  xi_for_ss <- compute_xi_mean_yasso07_engine(
    clim[seq_len(n_ss), , drop = FALSE], sanity_params)
  C_init    <- steady_state_yasso07_engine(sanity_params, lm, xi_for_ss)
  run_out   <- yasso07_run_engine(inputs, sanity_params, C_init, xi_array)

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
  compute_xi      = compute_xi_yasso07_engine,
  compute_xi_mean = compute_xi_mean_yasso07_engine,
  steady_state    = steady_state_yasso07_engine,
  run_model       = yasso07_run_engine,
  sigma_obs_fixed = sigma_obs_fixed,
  plots           = plots,
  climate_by_plot = climate_by_plot,
  inputs_by_plot  = inputs_by_plot,
  litter_means    = litter_means,
  obs_meta        = obs_meta,
  steady_state_n  = STEADY_STATE_YEARS,
  transient_init  = TRUE              # NEW
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

# Thermodynamic recycling constraint: carbon flow from a pool back toward
# faster pools must be a minority flux (< 0.5). Column A has no constraint
# as it is the fastest labile pool. Yasso07 uses the full three-term version
# since p_NW and p_EW are free (not structurally fixed as in Yasso15/20).
check_recycling_fractions <- function(p_phys) {
  (p_phys["p_WA"]                   < 0.5) &&
    (p_phys["p_EA"] + p_phys["p_EW"] < 0.5) &&
    (p_phys["p_NA"] + p_phys["p_NW"] < 0.5)
}

prior <- createPrior(
  density = function(x) {
    p_phys <- to_original(x)
    if (!check_recycling_fractions(p_phys)) return(-Inf)
    sum(dnorm(x, mean = best_x, sd = sigma_ppm, log = TRUE))
  },
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

HIGHLIGHT <- c("beta1","beta2","gamma","delta1","delta2","r",
               "sigma_init","sigma_input")
GROUP1    <- c(FRAC_COL_A, FRAC_COL_W, FRAC_COL_E, FRAC_COL_N)
GROUP2    <- HIGHLIGHT

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
message("Next: run_Yasso07_transient_predictive.R with RUN_ID=", RUN_ID)
