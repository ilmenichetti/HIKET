# =============================================================================
# run_Yasso20_calibration.R
#
# CALIBRATION STAGE — Yasso20 / HIKET pipeline.
#
# Structure mirrors run_Yasso15_calibration.R exactly. This header only
# documents what genuinely differs in Yasso20.
#
# KEY DIFFERENCES FROM YASSO15:
#   * Monthly climate input. Yasso20 replaces the four-point Gaussian
#     approximation of intra-annual temperature variation (used in Yasso07/15)
#     with a direct average over all 12 monthly temperatures. The climate
#     object passed through the engine is monthly format (12 rows per year)
#     rather than the annual summary used by Yasso07/15.
#
#   * steady_state_n = STEADY_STATE_YEARS * 12L. Because climate_by_plot
#     holds monthly rows, the engine's seq_len(steady_state_n) slice must
#     index monthly rows to span the intended number of years. With 12 rows
#     per year, STEADY_STATE_YEARS * 12 gives the correct window.
#     compute_xi_mean_yasso20 groups by month via tapply, so receiving
#     240 rows (20 years x 12 months) and returning scalar xi_awe/n/h
#     at mean climate is exactly correct.
#
#   * yasso15.so is reused (no separate Fortran for Yasso20). Both
#     yasso15_wrapper.R and yasso20_wrapper.R must be sourced; the latter
#     defines YASSO20_PARAM_NAMES/DEFAULT_PARAMS (aliased to Yasso15's)
#     and the compute_xi_yasso20 / compute_xi_mean_yasso20 functions.
#
#   * yasso20_run_engine calls yasso15_run directly with the xi_arrays
#     pre-computed by the engine. This avoids the redundant xi recomputation
#     that would occur if yasso20_run() were called instead (it recomputes
#     xi internally from climate_df). The Fortran call is identical either way.
#
#   * Annual precip for yasso15_run's leaching term is aggregated from monthly
#     precip and merged into inputs_by_plot, same approach as Yasso15.
#
# PARAMETER STRUCTURE: identical to Yasso15.
#   35 parameters, same names, same defaults (YASSO20_DEFAULT_PARAMS aliases
#   YASSO15_DEFAULT_PARAMS). Same 28 free / 11 fixed split. Same priors.
#   This is by design: the only structural difference between Yasso15 and
#   Yasso20 in the intercomparison is the intra-annual temperature integration
#   method. All else equal, posterior differences in climate response parameters
#   are attributable to that design choice alone.
#
# OUTPUTS:
#   ./Calibration_real_data/runs/Yasso20_chains_<RUN_ID>.rds
#   ./Calibration_real_data/runs/Yasso20_posterior_<RUN_ID>.rds
#   ./Calibration_real_data/runs/Yasso20_metadata_<RUN_ID>.rds
#   ./Calibration_real_data/model_inputs/Yasso20_inputs_<RUN_ID>.rds
#   ./Calibration_real_data/diagnostics/Yasso20/Yasso20_report_<RUN_ID>.txt
#   ./Calibration_real_data/diagnostics/Yasso20/Yasso20_traces_<RUN_ID>.png
#   ./Calibration_real_data/diagnostics/Yasso20/Yasso20_marginals_<RUN_ID>.png
#   ./Calibration_real_data/diagnostics/Yasso20/Yasso20_forward_at_defaults_<RUN_ID>.png
#
# DOWNSTREAM:
#   Rscript run_Yasso20_predictive.R <RUN_ID>
#   Rscript run_residual_analysis.R Yasso20 <RUN_ID>
# =============================================================================

source("./Calibration_real_data/calibration_engine.R")
source("./Model_functions_real_data/input_compatibility_layer.R")
# Yasso20 reuses yasso15.so; yasso15_wrapper.R must be loaded first because
# yasso20_wrapper.R aliases YASSO15_PARAM_NAMES and calls yasso15_run().
source("./Model_functions_real_data/Decomposition_functions/Yasso/yasso15_wrapper.R")
source("./Model_functions_real_data/Decomposition_functions/Yasso/yasso20_wrapper.R")
dyn.load("./Model_functions_real_data/Decomposition_functions/Yasso/yasso15.so")

library(dplyr)


# =============================================================================
# 0.  Configuration
# =============================================================================

MODEL_NAME <- "Yasso20"
RUN_ID     <- format(Sys.time(), "%Y%m%d_%H%M%S")

# --- Test (Mac) settings ---
N_PLOTS_TEST <- NA
N_CHAINS     <- 4L
N_ITER       <- 10000L
N_BURNIN     <- 2000L
N_LOG        <- 200L

# --- Production (Roihu) settings: uncomment to switch ---
# N_PLOTS_TEST <- NA
# N_CHAINS     <- 4L
# N_ITER       <- 50000L
# N_BURNIN     <- 10000L
# N_LOG        <- 1000L

CORES_PER_CHAIN <- if (grepl("puhti|mahti", Sys.info()["nodename"])) parallelly::availableCores() else parallel::detectCores() - 1L
# CORES_PER_CHAIN <- 90L   # Roihu

# Steady-state climate window in YEARS. Converted to monthly rows below
# when passed to make_likelihood, because climate_by_plot is monthly format.
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
# Identical to Yasso15: 35 parameters, same names, same fixed/free split.
# YASSO20_DEFAULT_PARAMS and YASSO20_PARAM_NAMES are aliases of their Yasso15
# equivalents (defined in yasso20_wrapper.R). Using the _YASSO20_ names here
# for clarity and to make any future divergence between the model families easy
# to introduce without refactoring.
# =============================================================================

p_default        <- YASSO20_DEFAULT_PARAMS
FIXED_RATE_NAMES <- c("alpha_A","alpha_W","alpha_E","alpha_N",
                      "p_H","alpha_H",
                      "w1","w2","w3","w4","w5")
fixed_rates      <- p_default[FIXED_RATE_NAMES]

B <- 1.0 - fixed_rates["p_H"]

FRAC_COL_A <- c("p_AW","p_AE","p_AN")
FRAC_COL_W <- c("p_WA","p_WE","p_WN")
FRAC_COL_E <- c("p_EA","p_EW","p_EN")
FRAC_COL_N <- c("p_NA","p_NW","p_NE")

param_spec <- list(
  list(names = FRAC_COL_A,                         type = "stick_break",   budget = B),
  list(names = FRAC_COL_W,                         type = "stick_break",   budget = B),
  list(names = FRAC_COL_E,                         type = "stick_break",   budget = B),
  list(names = FRAC_COL_N,                         type = "stick_break",   budget = B),
  list(names = c("beta1",  "beta2",  "gamma"),      type = "unconstrained"),  # AWE climate
  list(names = c("betaN1", "betaN2", "gammaN"),     type = "unconstrained"),  # N   climate
  list(names = c("betaH1", "betaH2", "gammaH"),     type = "unconstrained"),  # H   climate
  list(names = c("delta1", "delta2", "r"),          type = "unconstrained"),  # size modifier
  list(names = c("sigma_init"),                     type = "log"),
  list(names = c("sigma_input"),                    type = "log")
)

transforms       <- build_transforms(param_spec)
to_original      <- transforms$to_original
to_unconstrained <- transforms$to_unconstrained
log_jacobian     <- transforms$log_jacobian
FREE_NAMES       <- transforms$param_names
N_FREE           <- transforms$n_params

free_defaults <- c(
  p_default[c(unlist(lapply(param_spec[1:8], `[[`, "names")))],
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
# Identical to Yasso15; uses YASSO20_PARAM_NAMES for the reorder step.
# =============================================================================

assemble_model_params <- function(p_free) {
  yasso_params <- c(
    fixed_rates,
    p_free[c(unlist(lapply(param_spec[1:8], `[[`, "names")))]
  )
  yasso_params <- yasso_params[YASSO20_PARAM_NAMES]
  c(yasso_params, sigma_input = unname(p_free["sigma_input"]))
}


# =============================================================================
# 3.  Data preparation
# =============================================================================
# Yasso20 uses monthly climate (map_climate_monthly: 12 rows per year) and
# Yasso20-specific litter inputs (map_inputs_yasso20).
#
# Two additions relative to Yasso15 data prep:
#   (a) Annual precip is aggregated from monthly precip and merged into
#       inputs_by_plot. yasso15_run() (called by yasso20_run_engine) needs
#       it for the leaching term.
#   (b) precip_mean in litter_means is the annual total from mean monthly
#       values: sum(tapply(monthly_precip, month, mean)) over the steady-state
#       window. This matches how compute_xi_mean_yasso20 constructs the
#       synthetic mean year internally.
# =============================================================================

message("Loading data...")

input_raw <- read.csv("./Data/model_inputs/input_raw_monthly.csv")
site_raw  <- read.csv("./Data/model_inputs/site_raw.csv")

litter_cols <- grep("^C_nwl|^C_fwl|^C_cwl", names(input_raw), value = TRUE)
annual_litter <- input_raw %>%
  group_by(plot_id, year) %>%
  summarise(total_litter = sum(across(all_of(litter_cols)), na.rm = TRUE),
            .groups = "drop")

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
  message(sprintf("WARNING: %d NA values in litter columns. Filtering...", n_lit_na))
  bad_rows <- !complete.cases(input_calib[, litter_cols])
  plots_with_lit_na <- unique(input_calib$plot_id[bad_rows])
  calib_plots <- setdiff(calib_plots, plots_with_lit_na)
  input_calib <- input_calib[input_calib$plot_id %in% calib_plots, ]
  message(sprintf("After litter filter: %d plots", length(calib_plots)))
}

# Monthly climate for Yasso20 (12 rows per year per plot).
Yasso20_monthly <- map_climate_monthly(input_calib)
Yasso20_inputs  <- map_inputs_yasso20(input_calib)

# (a) Aggregate monthly precip to annual and merge into inputs.
#     yasso15_run() (called inside yasso20_run_engine) needs annual precip
#     for the leaching term. Merging here keeps the engine interface unchanged.
annual_precip_df <- Yasso20_monthly %>%
  group_by(plot_id, year) %>%
  summarise(precip = sum(precip, na.rm = TRUE), .groups = "drop")

Yasso20_inputs <- merge(
  Yasso20_inputs,
  annual_precip_df,
  by  = c("plot_id", "year"),
  all.x = TRUE
)

plots_real <- as.character(unique(Yasso20_monthly$plot_id))
message(sprintf("Plots after mapping: %d", length(plots_real)))

# Sanity check: verify that inputs and monthly climate cover identical year sets
# per plot. A mismatch means match() in obs_meta will silently produce NAs,
# causing those plots to be skipped in the likelihood with no warning.
year_check <- vapply(plots_real, function(pid) {
  yrs_inp  <- sort(unique(Yasso20_inputs[Yasso20_inputs$plot_id == pid, "year"]))
  yrs_clim <- sort(unique(Yasso20_monthly[Yasso20_monthly$plot_id == pid, "year"]))
  identical(yrs_inp, yrs_clim)
}, logical(1))

if (!all(year_check)) {
  bad <- plots_real[!year_check]
  stop(sprintf(
    "Year mismatch between inputs and monthly climate for %d plot(s): %s",
    length(bad), paste(bad, collapse = ", ")))
}
message(sprintf("Year coverage check: OK (%d plots)", length(plots_real)))

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

# climate_by_plot: monthly rows (12 per year). The engine's seq_len(n_ss)
# slice is sized in monthly rows (STEADY_STATE_YEARS * 12) so the window
# covers the correct number of years.
climate_by_plot <- split(Yasso20_monthly, Yasso20_monthly$plot_id)
inputs_by_plot  <- split(Yasso20_inputs,  Yasso20_inputs$plot_id)

# litter_means: (b) precip_mean for Yasso20 is the annual total reconstructed
# from mean monthly values over the steady-state window. This matches the
# synthetic mean-year logic inside compute_xi_mean_yasso20:
#   mean_clim$precip = mean monthly precip  ->  sum = annual total.
litter_means <- lapply(plots_real, function(pid) {
  inp  <- Yasso20_inputs[as.character(Yasso20_inputs$plot_id) == pid, ]
  clim <- Yasso20_monthly[as.character(Yasso20_monthly$plot_id) == pid, ]
  n_ss_rows <- min(STEADY_STATE_YEARS * 12L, nrow(clim))
  clim_ss   <- clim[seq_len(n_ss_rows), ]
  # Use years appearing in clim_ss to subset litter inputs consistently.
  ss_years <- unique(clim_ss$year)
  inp_ss   <- inp[inp$year %in% ss_years, ]
  list(
    nwl_mean    = colMeans(inp_ss[, c("nwl_A","nwl_W","nwl_E","nwl_N")]),
    fwl_mean    = colMeans(inp_ss[, c("fwl_A","fwl_W","fwl_E","fwl_N")]),
    cwl_mean    = colMeans(inp_ss[, c("cwl_A","cwl_W","cwl_E","cwl_N")]),
    precip_mean = sum(tapply(clim_ss$precip, clim_ss$month, mean, na.rm = TRUE))
  )
})
names(litter_means) <- plots_real

# obs_meta: index into annual rows. For Yasso20, the model output (from
# yasso15_run) is annual, indexed by year -- same structure as Yasso07/15.
# We match observation years against the annual inputs frame (not monthly).
obs_meta <- lapply(plots_real, function(pid) {
  inp_plot <- Yasso20_inputs[as.character(Yasso20_inputs$plot_id) == pid, ]
  obs_plot <- SOC_obs_all[as.character(SOC_obs_all$plot_id) == pid, ]
  list(
    idx      = match(obs_plot$year, inp_plot$year),
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
# Same four-function pattern as Yasso07/15. Key differences:
#
#   compute_xi_yasso20_engine:
#     Takes monthly climate_df, returns list(xi_awe, xi_n, xi_h) via
#     compute_xi_yasso20. The engine passes the full per-plot monthly clim
#     frame; compute_xi_yasso20 handles the year-by-year loop internally.
#
#   compute_xi_mean_yasso20_engine:
#     Receives the steady-state monthly window (STEADY_STATE_YEARS * 12 rows).
#     Delegates to compute_xi_mean_yasso20, which builds a synthetic mean year
#     (one row per calendar month) and calls compute_xi_yasso20 once.
#
#   steady_state_yasso20_engine:
#     Passes pre-computed xi_ss list and precip_mean (from lm) to the
#     refactored yasso20_steady_state(), which delegates to yasso15_steady_state().
#
#   yasso20_run_engine:
#     Calls yasso15_run directly with the xi_arrays pre-computed by the engine.
#     This avoids the redundant xi recomputation inside yasso20_run() and keeps
#     the sigma_input scaling on litter inputs consistent with Yasso07/15.
# =============================================================================

compute_xi_yasso20_engine <- function(clim, model_params) {
  compute_xi_yasso20(
    climate_df = clim,
    params     = model_params[YASSO20_PARAM_NAMES]
  )
}

compute_xi_mean_yasso20_engine <- function(clim_ss, model_params) {
  compute_xi_mean_yasso20(
    clim_ss = clim_ss,
    params  = model_params[YASSO20_PARAM_NAMES]
  )
}

steady_state_yasso20_engine <- function(model_params, lm, xi_ss) {
  si <- model_params["sigma_input"]
  yasso20_steady_state(
    params      = model_params[YASSO20_PARAM_NAMES],
    nwl_mean    = lm$nwl_mean * si,
    fwl_mean    = lm$fwl_mean * si,
    cwl_mean    = lm$cwl_mean * si,
    xi_ss       = xi_ss,
    precip_mean = lm$precip_mean
  )
}

LITTER_COLS <- c("nwl_A","nwl_W","nwl_E","nwl_N",
                 "fwl_A","fwl_W","fwl_E","fwl_N",
                 "cwl_A","cwl_W","cwl_E","cwl_N")

# Calls yasso15_run directly with the engine's pre-computed xi_arrays and
# annual precip (merged into inputs during data prep). This is equivalent to
# yasso20_run but avoids the xi recomputation that yasso20_run would trigger.
yasso20_run_engine <- function(inputs, model_params, C_init, xi_arrays) {
  si <- model_params["sigma_input"]
  inputs_scaled <- inputs
  inputs_scaled[, LITTER_COLS] <- inputs[, LITTER_COLS] * si
  yasso15_run(
    input_df  = inputs_scaled,
    params    = model_params[YASSO20_PARAM_NAMES],
    C_init    = C_init,
    xi_arrays = xi_arrays,
    precip    = inputs$precip   # annual, aggregated from monthly in Section 3
  )
}


# =============================================================================
# 5.  Prior specification
# =============================================================================
# Identical to Yasso15. See run_Yasso15_calibration.R Section 5 for rationale.
# Keeping priors identical across Yasso15 and Yasso20 is a deliberate choice:
# the only structural difference between the two models is the xi computation
# method, and the prior should not introduce additional asymmetry.
# =============================================================================

sigma_ppm <- setNames(rep(1.0, N_FREE), FREE_NAMES)

sigma_ppm["beta1"]   <- 0.2;  sigma_ppm["beta2"]   <- 0.05
sigma_ppm["betaN1"]  <- 0.2;  sigma_ppm["betaN2"]  <- 0.05
sigma_ppm["betaH1"]  <- 0.2;  sigma_ppm["betaH2"]  <- 0.05

sigma_ppm["sigma_init"]  <- 0.5
sigma_ppm["sigma_input"] <- 0.5

stopifnot(length(sigma_ppm) == N_FREE)
stopifnot(all(sigma_ppm > 0))
stopifnot(all(names(sigma_ppm) == FREE_NAMES))

message("\nPrior widths (unconstrained space):")
print(round(sigma_ppm, 3))


# =============================================================================
# 6.  Pre-MCMC sanity check
# =============================================================================
# Same logic as Yasso15. Note that n_ss for the sanity check is in monthly
# rows (STEADY_STATE_YEARS * 12), matching the make_likelihood call below.
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

  xi_array  <- compute_xi_yasso20_engine(clim, sanity_params)
  n_ss      <- min(STEADY_STATE_YEARS * 12L, nrow(clim))
  xi_for_ss <- compute_xi_mean_yasso20_engine(
    clim[seq_len(n_ss), , drop = FALSE], sanity_params)
  C_init    <- steady_state_yasso20_engine(sanity_params, lm, xi_for_ss)
  run_out   <- yasso20_run_engine(inputs, sanity_params, C_init, xi_array)

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
                        sprintf("%s_forward_at_defaults_%s.png",
                                MODEL_NAME, RUN_ID))
png(sanity_png, width = 10L * PX_PER_IN, height = 8L * PX_PER_IN, res = PX_PER_IN)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for (sr in sanity_results) {
  yr  <- sr$run_out$year
  soc <- sr$run_out$total_soc
  obs_yr  <- yr[sr$meta$idx]
  obs_soc <- sr$meta$soc_obs
  ylim_r  <- range(c(soc, obs_soc), na.rm = TRUE)
  plot(yr, soc, type = "l", col = "steelblue", lwd = 2, ylim = ylim_r,
       xlab = "Year", ylab = "SOC (tC/ha)",
       main = sprintf("Plot %s (defaults)", sr$pid))
  if (length(obs_soc) > 0)
    points(obs_yr, obs_soc, pch = 19, col = "firebrick", cex = 1.3)
  legend("topleft",
         legend = sprintf("xi_awe=%.3f  xi_n=%.3f  xi_h=%.3f",
                          sr$xi_for_ss$xi_awe,
                          sr$xi_for_ss$xi_n,
                          sr$xi_for_ss$xi_h),
         bty = "n", cex = 0.75)
}
mtext(sprintf("%s forward run at published defaults  |  %s", MODEL_NAME, RUN_ID),
      side = 3, outer = TRUE, line = -1.5, cex = 0.95, font = 2)
dev.off()
message(sprintf("Sanity plot saved: %s", sanity_png))
message(sprintf("Forward-run sanity: %s",
                if (forward_run_pass) "PASS" else "WARN -- inspect plot"))


# =============================================================================
# 7.  Save processed inputs
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
# steady_state_n = STEADY_STATE_YEARS * 12L because climate_by_plot holds
# monthly rows. The engine slices clim[seq_len(steady_state_n), ] to get the
# steady-state window; 12 rows per year gives the correct number of years.
# compute_xi_mean_yasso20 handles the monthly grouping internally.
# =============================================================================

ll_fn <- make_likelihood(
  n_cores         = CORES_PER_CHAIN,
  to_original     = to_original,
  log_jacobian    = log_jacobian,
  assemble_params = assemble_model_params,
  compute_xi      = compute_xi_yasso20_engine,
  compute_xi_mean = compute_xi_mean_yasso20_engine,
  steady_state    = steady_state_yasso20_engine,
  run_model       = yasso20_run_engine,
  sigma_obs_fixed = sigma_obs_fixed,
  plots           = plots,
  climate_by_plot = climate_by_plot,
  inputs_by_plot  = inputs_by_plot,
  litter_means    = litter_means,
  obs_meta        = obs_meta,
  steady_state_n  = STEADY_STATE_YEARS * 12L   # monthly rows, not years
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
    for (j in seq_len(N_FREE)) {
      m[, j] <- rnorm(n, mean = best_x[j], sd = sigma_ppm[j] * 0.1)
    }
    m
  },
  lower = rep(-Inf, N_FREE),
  upper = rep( Inf, N_FREE)
)


# =============================================================================
# 9.  Run MCMC
# =============================================================================

t_run <- system.time({
  mcmc_out <- run_mcmc_chains(
    ll_fn         = ll_fn,
    prior         = prior,
    best_x        = best_x,
    param_names   = FREE_NAMES,
    mcmc_settings = mcmc_settings,
    run_config    = run_config)
})["elapsed"]

chain_results <- mcmc_out$chain_results
chain_health  <- mcmc_out$chain_health

message(sprintf("\nAll chains complete. Wallclock: %.1f min\n", t_run / 60))


# =============================================================================
# 10.  Diagnostics, save, and report
# =============================================================================

HIGHLIGHT <- c("beta1","beta2","gamma",
               "betaN1","betaN2","gammaN",
               "betaH1","betaH2","gammaH",
               "delta1","delta2","r",
               "sigma_init","sigma_input")

GROUP1 <- c("p_AW","p_AE","p_AN",
            "p_WA","p_WE","p_WN",
            "p_EA","p_EW","p_EN",
            "p_NA","p_NW","p_NE")
GROUP2 <- HIGHLIGHT

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
message("Next: run_Yasso20_predictive.R with RUN_ID=", RUN_ID)
