# =============================================================================
# run_Yasso15_transient_calibration.R
#
# CALIBRATION STAGE — Yasso15 / HIKET pipeline.
#
# Structure mirrors run_Yasso07_calibration.R exactly. Read that script for
# the full rationale behind every architectural decision. This header only
# documents what genuinely differs in Yasso15.
#
# KEY DIFFERENCES FROM YASSO07:
#   * 35-parameter model (vs 24). Yasso15 adds pool-group-specific climate
#     response: AWE, N, and H each have their own beta1/beta2/gamma triplet
#     (betaN1/betaN2/gammaN and betaH1/betaH2/gammaH). All three triplets
#     are free parameters here (same rationale as beta1/beta2/gamma in 07).
#
#   * Leaching weight parameters w1-w5 exist in the vector but are fixed at
#     defaults. leac = 0.0 throughout, making them unidentifiable in practice.
#
#   * xi_for_ss is a named list (xi_awe, xi_n, xi_h), not a scalar.
#     The engine's is_valid_xi() helper handles both forms.
#
#   * yasso15_run() requires an annual precip vector for the leaching term.
#     This is merged into inputs_by_plot from Yasso07_climate$precip during
#     data preparation (see Section 3), so the engine's run_model interface
#     remains unchanged.
#
#   * litter_means stores precip_mean (mean annual precip over the steady-state
#     window) because yasso15_steady_state() passes it to Fortran for leaching.
#
# FIXED PARAMETERS (11 total):
#   alpha_A, alpha_W, alpha_E, alpha_N -- base decomposition rates.
#   p_H, alpha_H                       -- humus pool fraction and rate.
#   w1, w2, w3, w4, w5                 -- leaching weights (leac = 0.0).
#   Same rationale as Yasso07: rates are well-constrained, humus parameters
#   are not the focus of the intercomparison, leaching is inactive.
#
# FREE PARAMETERS (28 total: 26 model + 2 nuisance):
#   12 inter-pool flow fractions (stick-break, same groups as Yasso07)
#    3 AWE climate response:  beta1,  beta2,  gamma
#    3 N   climate response:  betaN1, betaN2, gammaN
#    3 H   climate response:  betaH1, betaH2, gammaH
#    3 woody size modifier:   delta1, delta2, r
#    2 nuisance:              sigma_init, sigma_input
#
# OUTPUTS (same naming convention as Yasso07):
#   ./Calibration_real_data/runs/Yasso15_chains_<RUN_ID>.rds
#   ./Calibration_real_data/runs/Yasso15_posterior_<RUN_ID>.rds
#   ./Calibration_real_data/runs/Yasso15_metadata_<RUN_ID>.rds
#   ./Calibration_real_data/model_inputs/Yasso15_inputs_<RUN_ID>.rds
#   ./Calibration_real_data/diagnostics/Yasso15/Yasso15_report_<RUN_ID>.txt
#   ./Calibration_real_data/diagnostics/Yasso15/Yasso15_traces_<RUN_ID>.png
#   ./Calibration_real_data/diagnostics/Yasso15/Yasso15_marginals_<RUN_ID>.png
#   ./Calibration_real_data/diagnostics/Yasso15/Yasso15_forward_at_defaults_<RUN_ID>.png
#
# DOWNSTREAM:
#   Rscript run_Yasso15_predictive.R <RUN_ID>
#   Rscript run_residual_analysis.R Yasso15 <RUN_ID>
# =============================================================================

source("./Calibration_real_data_transient/calibration_engine_transient.R")
source("./Model_functions_real_data_transient/input_compatibility_layer.R")
source("./Model_functions_real_data_transient/Decomposition_functions/Yasso/yasso15_wrapper_transient.R")
dyn.load("./Model_functions_real_data/Decomposition_functions/Yasso/yasso15.so")

library(dplyr)


# =============================================================================
# 0.  Configuration
# =============================================================================

MODEL_NAME <- "Yasso15"
RUN_ID     <- format(Sys.time(), "%Y%m%d_%H%M%S")

# read the calibration configuration (N_PLOTS_TEST, N_CHAINS, N_ITER, N_BURNIN, N_LOG)
source("./Calibration_real_data_transient/calib_config.R")


CORES_PER_CHAIN <- if (grepl("puhti|mahti", Sys.info()["nodename"])) parallelly::availableCores() else parallel::detectCores() - 1L
# CORES_PER_CHAIN <- 90L   # Roihu

STEADY_STATE_YEARS <- 20L
PREINIT_YEAR       <- 1917L  # Finland independence -- start of conceptual pre-run
N_PREINIT_SMOOTH   <- 5L     # years averaged for t0 litter endpoint (not used in steady-state approach)

DIR_LOGS   <- "./Calibration_real_data_transient/progress_logs/"
DIR_DIAG   <- file.path("./Calibration_real_data_transient/diagnostics", MODEL_NAME)
DIR_RUNS   <- "./Calibration_real_data_transient/runs/"
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
# Yasso15 has 35 parameters. Fixed (11): the four base decomp rates, humus
# pool parameters, and the five leaching weights (inactive because leac = 0.0).
# Free (26 model + 2 nuisance = 28 total): flow fractions, all three climate
# response triplets, woody size modifier, sigma_init, sigma_input.
#
# The three climate response groups (AWE / N / H) are all unconstrained because
# beta and gamma parameters can take any real value in the original formulation.
# Tighter priors on beta2 terms (quadratic in T) avoid biophysically implausible
# temperature curves; see Section 5.
# =============================================================================

p_default        <- YASSO15_DEFAULT_PARAMS
FIXED_RATE_NAMES <- c("alpha_A","alpha_W","alpha_E","alpha_N",
                      "p_H","alpha_H",
                      "w1","w2","w3","w4","w5")
fixed_rates      <- p_default[FIXED_RATE_NAMES]

# Stick-breaking budget: same logic as Yasso07.
# Each donor column's outflow fractions must sum to <= B = 1 - p_H.
B <- 1.0 - fixed_rates["p_H"]

FRAC_COL_A <- c("p_AW","p_AE","p_AN")
FRAC_COL_W <- c("p_WA","p_WE","p_WN")
FRAC_COL_E <- c("p_EA","p_EW","p_EN")
FRAC_COL_N <- c("p_NA","p_NW","p_NE")

param_spec <- list(
  list(names = FRAC_COL_A,    type = "stick_break",   budget = B),
  list(names = FRAC_COL_W,    type = "stick_break",   budget = B),
  list(names = FRAC_COL_E,    type = "stick_break",   budget = B),
  list(names = FRAC_COL_N,    type = "stick_break",   budget = B),
  # AWE climate -- individual entries for per-parameter transforms.
  # beta1/betaN1/betaH1 > 0 enforced by log (T increases decomp in Finland).
  # gamma/gammaN/gammaH < 0 enforced by tight prior; mathematically unconstrained.
  list(names = "beta1",       type = "log"),           # AWE T linear;    enforce > 0
  list(names = "beta2",       type = "unconstrained"), # AWE T quadratic; can be negative
  list(names = "gamma",       type = "unconstrained"), # AWE precip;      sign via prior
  # N climate
  list(names = "betaN1",      type = "log"),           # N T linear;      enforce > 0
  list(names = "betaN2",      type = "unconstrained"), # N T quadratic
  list(names = "gammaN",      type = "unconstrained"), # N precip
  # H climate
  list(names = "betaH1",      type = "log"),           # H T linear;      enforce > 0
  list(names = "betaH2",      type = "unconstrained"), # H T quadratic
  list(names = "gammaH",      type = "unconstrained"), # H precip
  # Woody size modifier -- individual entries.
  # delta2 > 0: NaN in Fortran if delta2 < 0 for d > ~0.4 cm.
  # r log-transformed: only -ABS(r) used.
  list(names = "delta1",      type = "unconstrained"), # linear;  can be negative
  list(names = "delta2",      type = "log"),           # quadratic; enforce > 0
  list(names = "r",           type = "log"),           # power;   only |r| used
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

# Prior centres in PHYSICAL space.
# Climate/size: posterior means from FMI repository Yasso15.dat (10,000 MCMC
# samples after dropping MAP row; confirmed to match y15par.csv to 4 sig. fig.
# after column remapping -- see Priors_model_matching.R).
# Fractions: Yasso15 published defaults (p_default).
free_defaults <- c(
  p_default[MODEL_FREE_NAMES],
  sigma_init  = 1.00,   # CHANGED: prior centre = 1.0 (same productivity as today)
  sigma_input = 1.00
)

# Override climate/size centres with Yasso15 posterior means
free_defaults["beta1"]   <-   0.09062
free_defaults["beta2"]   <-  -0.000215
free_defaults["gamma"]   <-  -1.80897
free_defaults["betaN1"]  <-   0.04878
free_defaults["betaN2"]  <-  -0.0000792
free_defaults["gammaN"]  <-  -1.17294
free_defaults["betaH1"]  <-   0.03518
free_defaults["betaH2"]  <-  -0.000208
free_defaults["gammaH"]  <- -12.54094
free_defaults["delta1"]  <-  -0.43883
free_defaults["delta2"]  <-   1.26838
free_defaults["r"]       <-   0.25687

best_x <- to_unconstrained(free_defaults)

stopifnot(
  max(abs(to_original(best_x)[FREE_NAMES] - free_defaults[FREE_NAMES])) < 1e-10)
message("Transform round-trip: OK")


# =============================================================================
# 2.  Parameter assembly
# =============================================================================
# Combines fixed_rates (11 elements) and free model parameters (24 elements)
# into a 35-element vector in YASSO15_PARAM_NAMES order, then appends
# sigma_input. MODEL_FREE_NAMES excludes sigma_init and sigma_input.
# =============================================================================

assemble_model_params <- function(p_free) {
  yasso_params <- c(fixed_rates, p_free[MODEL_FREE_NAMES])
  yasso_params <- yasso_params[YASSO15_PARAM_NAMES]
  c(yasso_params,
    sigma_input = unname(p_free["sigma_input"]),
    sigma_init  = unname(p_free["sigma_init"]))   # NEW: needed by transient_init
}


# =============================================================================
# 3.  Data preparation
# =============================================================================
# Yasso15 uses the same annual climate format as Yasso07 (temp_mean,
# temp_amplitude, precip) and the same litter input format (map_inputs_yasso07).
#
# Two additions relative to the Yasso07 run script:
#   (a) Annual precip is merged into inputs_by_plot. yasso15_run() needs it
#       for the leaching term even though leac = 0.0 renders it numerically
#       inactive. Storing it here avoids any change to the engine interface.
#   (b) precip_mean (mean annual precip over the steady-state window) is stored
#       in litter_means. yasso15_steady_state() passes it to Fortran for the
#       same leaching term in the steady-state call.
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

# Yasso15 uses the same annual climate and input format as Yasso07.
Yasso15_climate <- map_climate_yasso07(input_calib)
Yasso15_inputs  <- map_inputs_yasso07(input_calib)

# (a) Merge annual precip into inputs so yasso15_run_engine can pass it to
#     yasso15_run() without touching the engine interface. The precip column
#     in Yasso15_climate aligns by plot_id and year.
Yasso15_inputs <- merge(
  Yasso15_inputs,
  Yasso15_climate[, c("plot_id", "year", "precip")],
  by = c("plot_id", "year"),
  all.x = TRUE
)

plots_real <- as.character(unique(Yasso15_climate$plot_id))
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

# sigma_obs is fixed at the SAME observed CV used for Yasso07.
# This is architecturally critical for the intercomparison: all models are
# evaluated against the same observation error standard so that likelihood
# differences reflect model structure, not error model choices.
obs_cv          <- sd(SOC_obs_all$soc_obs_tCha) / mean(SOC_obs_all$soc_obs_tCha)
sigma_obs_fixed <- obs_cv

message(sprintf("SOC observations: %d  |  obs CV: %.3f  |  sigma_obs fixed: %.3f",
                nrow(SOC_obs_all), obs_cv, sigma_obs_fixed))

climate_by_plot <- split(Yasso15_climate, Yasso15_climate$plot_id)
inputs_by_plot  <- split(Yasso15_inputs,  Yasso15_inputs$plot_id)

# litter_means: (b) adds precip_mean alongside nwl/fwl/cwl means.
# precip_mean = mean annual precip over the steady-state window, passed to
# yasso15_steady_state() for the Fortran leaching term.
litter_means <- lapply(plots_real, function(pid) {
  inp  <- Yasso15_inputs[as.character(Yasso15_inputs$plot_id) == pid, ]
  clim <- Yasso15_climate[as.character(Yasso15_climate$plot_id) == pid, ]
  n_ss    <- min(STEADY_STATE_YEARS, nrow(inp))
  inp_ss  <- inp[seq_len(n_ss), ]
  clim_ss <- clim[seq_len(n_ss), ]
  list(
    nwl_mean    = colMeans(inp_ss[, c("nwl_A","nwl_W","nwl_E","nwl_N")]),
    fwl_mean    = colMeans(inp_ss[, c("fwl_A","fwl_W","fwl_E","fwl_N")]),
    cwl_mean    = colMeans(inp_ss[, c("cwl_A","cwl_W","cwl_E","cwl_N")]),
    precip_mean = mean(clim_ss$precip)   # for Fortran leaching term
  )
})
names(litter_means) <- plots_real

obs_meta <- lapply(plots_real, function(pid) {
  clim     <- Yasso15_climate[as.character(Yasso15_climate$plot_id) == pid, ]
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
# Four engine functions adapting the Yasso15 wrapper to the engine's generic
# interface. Pattern is identical to Yasso07; the differences are:
#
#   compute_xi_yasso15_engine:
#     Returns a named list (xi_awe, xi_n, xi_h) rather than a scalar.
#     The engine's is_valid_xi() handles both forms.
#
#   compute_xi_mean_yasso15_engine:
#     Returns a list (xi_awe, xi_n, xi_h) at mean climate. xi(mean_climate)
#     not mean(xi_annual) -- same Jensen's inequality fix as Yasso07, now
#     applied to all three pool groups.
#
#   steady_state_yasso15_engine:
#     Passes the pre-computed xi_ss list and precip_mean (from lm) to the
#     refactored yasso15_steady_state(). No climate computation inside the
#     wrapper -- that responsibility sits with compute_xi_mean, as in Yasso07.
#
#   yasso15_run_engine:
#     Passes xi_arrays (the list from compute_xi_yasso15_engine) and
#     inputs$precip (merged during data prep) to yasso15_run().
# =============================================================================

compute_xi_yasso15_engine <- function(clim, model_params) {
  compute_xi_yasso15(
    temp_mean = clim$temp_mean,
    temp_amp  = clim$temp_amplitude,
    precip    = clim$precip,
    params    = model_params[YASSO15_PARAM_NAMES]
  )
}

compute_xi_mean_yasso15_engine <- function(clim_ss, model_params) {
  compute_xi_mean_yasso15(
    clim_ss = clim_ss,
    params  = model_params[YASSO15_PARAM_NAMES]
  )
}

# CHANGED: yasso15_transient_init replaces yasso15_steady_state.
# sigma_init scales litter at which steady state is computed (prior centre 1.0).
steady_state_yasso15_engine <- function(model_params, lm, xi_ss) {
  si         <- unname(model_params["sigma_input"])
  sigma_init <- unname(model_params["sigma_init"])
  yasso15_transient_init(model_params, lm, xi_ss)
}

LITTER_COLS <- c("nwl_A","nwl_W","nwl_E","nwl_N",
                 "fwl_A","fwl_W","fwl_E","fwl_N",
                 "cwl_A","cwl_W","cwl_E","cwl_N")

yasso15_run_engine <- function(inputs, model_params, C_init, xi_arrays) {
  si <- model_params["sigma_input"]
  inputs_scaled <- inputs
  inputs_scaled[, LITTER_COLS] <- inputs[, LITTER_COLS] * si
  yasso15_run(
    input_df  = inputs_scaled,
    params    = model_params[YASSO15_PARAM_NAMES],
    C_init    = C_init,
    xi_arrays = xi_arrays,
    precip    = inputs$precip   # merged from climate in Section 3
  )
}


# =============================================================================
# 5.  Prior specification
# =============================================================================
# Default sigma_ppm = 1.0 for all parameters, then tightened for specific
# groups following the same rationale as Yasso07:
#
#   beta2 / betaN2 / betaH2 (= 0.05): quadratic temperature terms. Tighter
#     prior prevents biophysically implausible modifier curves. Applies to
#     all three pool groups since they share the same functional form.
#
#   beta1 / betaN1 / betaH1 (= 0.2): linear temperature terms. Same reasoning
#     as Yasso07. All three groups tightened equally -- the intercomparison
#     treats the pool groups symmetrically at the prior level.
#
#   sigma_init (= 0.5): log-scale. Physical 95% CI ~[0.04, 0.27].
#   sigma_input (= 0.5): log-scale. Physical 95% CI ~[0.37, 2.72].
#
# Note: gamma / gammaN / gammaH left at 1.0. The precipitation response is
# weakly identified on Finnish plots (low inter-annual P variability) and a
# wider prior lets the data speak more freely here.
# =============================================================================

# VALUES (Yasso15 -- from SD of FMI repository posterior, Yasso15.dat):
# Log-transformed params: sigma = SD(log(samples)); physical params: SD(samples).
# gammaH capped at 1.5: near-unidentifiable at Finnish P (~600mm) because
# [1-exp(gammaH*0.6)] ~ 1 for any gammaH < -5; large SD is an artefact.
# Transfer fractions remain at 1.0 (Finnish data informs pool structure).
sigma_ppm <- setNames(rep(1.0, N_FREE), FREE_NAMES)
sigma_ppm["beta1"]       <- 0.04593  # log space
sigma_ppm["beta2"]       <- 0.00014  # physical space
sigma_ppm["gamma"]       <- 0.07022  # physical space
sigma_ppm["betaN1"]      <- 0.10502  # log space
sigma_ppm["betaN2"]      <- 0.00008  # physical space
sigma_ppm["gammaN"]      <- 0.15543  # physical space
sigma_ppm["betaH1"]      <- 0.13477  # log space
sigma_ppm["betaH2"]      <- 0.00016  # physical space
sigma_ppm["gammaH"]      <- 1.50000  # physical space (capped; see above)
sigma_ppm["delta1"]      <- 0.32810  # physical space
sigma_ppm["delta2"]      <- 0.28643  # log space
sigma_ppm["r"]           <- 0.05504  # log space
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
# Same logic as Yasso07. One difference: xi_for_ss is now a list, so the
# legend formats all three pool-group xi values rather than a single scalar.
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

  xi_array  <- compute_xi_yasso15_engine(clim, sanity_params)
  n_ss      <- min(STEADY_STATE_YEARS, nrow(clim))
  xi_for_ss <- compute_xi_mean_yasso15_engine(
    clim[seq_len(n_ss), , drop = FALSE], sanity_params)
  C_init    <- steady_state_yasso15_engine(sanity_params, lm, xi_for_ss)
  run_out   <- yasso15_run_engine(inputs, sanity_params, C_init, xi_array)

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
  # xi_for_ss is a list: show all three pool-group values
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
  compute_xi      = compute_xi_yasso15_engine,
  compute_xi_mean = compute_xi_mean_yasso15_engine,
  steady_state    = steady_state_yasso15_engine,
  run_model       = yasso15_run_engine,
  sigma_obs_fixed = sigma_obs_fixed,
  plots           = plots,
  climate_by_plot = climate_by_plot,
  inputs_by_plot  = inputs_by_plot,
  litter_means    = litter_means,
  obs_meta        = obs_meta,
  steady_state_n  = STEADY_STATE_YEARS,
  transient_init  = TRUE              # NEW: sigma_init enters model, removed from sd_vec
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
# as it is the fastest labile pool. Column N is not constrained -- see
# Yasso07 note. p_EW is structurally fixed at 0 in Yasso20 so the E column
# reduces to p_EA only; kept as sum for Yasso15 where p_EW is free.
check_recycling_fractions <- function(p_phys) {
  (p_phys["p_WA"]                   < 0.5) &&
    (p_phys["p_EA"] + p_phys["p_EW"] < 0.5)
}

prior <- createPrior(
  density = function(x) {
    p_phys <- to_original(x)
    if (!check_recycling_fractions(p_phys)) return(-Inf)
    sum(dnorm(x, mean = best_x, sd = sigma_ppm, log = TRUE))
  },
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
# HIGHLIGHT includes all climate response parameters across the three pool
# groups. These are the parameters most directly comparable across models
# in the intercomparison, and the ones with the clearest policy interpretation
# (temperature sensitivity of decomposition).
# =============================================================================

HIGHLIGHT <- c("beta1","beta2","gamma",
               "betaN1","betaN2","gammaN",
               "betaH1","betaH2","gammaH",
               "delta1","delta2","r",
               "sigma_init","sigma_input")

# group1: inter-pool flow fractions (same 12 as Yasso07; Yasso15 adds
# separate climate parameters per pool group but keeps the same fraction
# naming convention). Plotted in 4-column grid.
# group2: all climate response parameters across the three pool groups,
# plus the two nuisance parameters. Plotted in 4-column grid.
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
message("Next: run_Yasso15_transient_predictive.R with RUN_ID=", RUN_ID)
