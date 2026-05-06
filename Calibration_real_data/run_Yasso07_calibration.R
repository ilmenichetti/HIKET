# =============================================================================
# run_Yasso07_calibration.R
#
# CALIBRATION STAGE of the Yasso07 / HIKET pipeline.
#
# WHAT THIS SCRIPT DOES (in order):
#   1.  Read pre-processed model inputs (input_raw_monthly.csv, site_raw.csv)
#       built by Data_work.R, and filter to calibration-ready plots.
#   2.  Set up the parameter specification: which parameters are free
#       (calibrated), which are fixed (published values), and how each one
#       transforms between unconstrained MCMC space and physical space.
#   3.  Define model interface wrappers: glue functions that adapt the
#       Yasso07 wrapper functions to the engine's likelihood factory.
#   4.  Run a pre-MCMC sanity check: forward simulation at published defaults
#       for 4 sample plots, with a PNG plot showing trajectories vs observations.
#   5.  Run MCMC: 4 DEzs chains via BayesianTools, with per-chain progress
#       logs and health tracking (initial/final ll, -Inf rejection rate, etc.).
#   6.  Compute convergence diagnostics: Gelman-Rubin R-hat (univariate +
#       multivariate), effective sample size, posterior summary tables.
#   7.  Save chain objects, posterior samples, and metadata to disk.
#   8.  Write a single text report summarising everything.
#
# WHAT THIS SCRIPT DOES NOT DO (split into other scripts):
#   - Posterior predictive simulation       -> run_Yasso07_predictive.R
#   - Observed-vs-predicted scatter plots   -> run_Yasso07_predictive.R
#   - Random Forest residual analysis       -> run_residual_analysis.R
#   This split exists because once MCMC has converged, you can iterate on
#   downstream analysis (e.g. tweaking RF covariates, refining plots) without
#   re-running calibration. On Mac the calibration is ~30 minutes; on Roihu
#   with the full ~3700 plots it will be hours.
#
# OUTPUTS:
#   ./Calibration/runs/Yasso07_chains_<RUN_ID>.rds
#       Raw BayesianTools chain objects -- needed for chain extension or
#       re-running diagnostics.
#   ./Calibration/runs/Yasso07_posterior_<RUN_ID>.rds
#       Matrix of post-burnin posterior samples in PHYSICAL space (after
#       to_original transform). Rows = samples, columns = parameters.
#       This is the file used by downstream scripts.
#   ./Calibration/runs/Yasso07_metadata_<RUN_ID>.rds
#       Run configuration, R version, wallclock, summary stats. For
#       reproducibility and bookkeeping.
#   ./Calibration/model_inputs/Yasso07_inputs_<RUN_ID>.rds
#       All pre-built per-plot lookup structures (climate_by_plot,
#       inputs_by_plot, litter_means, obs_meta, plot_info, etc.) and run
#       configuration. The predictive script reads this so it doesn't need
#       to re-source raw data files.
#   ./Calibration/diagnostics/Yasso07/Yasso07_report_<RUN_ID>.txt
#       Single human-readable text report. The predictive and residual
#       analysis scripts append further sections to this same file.
#   ./Calibration/diagnostics/Yasso07/Yasso07_traces_<RUN_ID>.png
#       Trace plots for highlighted parameters (one panel per parameter,
#       one line per chain). Visual convergence check.
#   ./Calibration/diagnostics/Yasso07/Yasso07_marginals_<RUN_ID>.png
#       Posterior vs prior density plots for highlighted parameters.
#       Shows what the data contributed beyond the prior.
#   ./Calibration/diagnostics/Yasso07/Yasso07_forward_at_defaults_<RUN_ID>.png
#       Sanity check: 4 plot panels showing the SOC trajectory under the
#       published default parameters, with observations overlaid. If the
#       trajectories are far from the observations, MCMC will struggle
#       and this is the first thing to look at.
#
# DOWNSTREAM (in this order):
#   Rscript run_Yasso07_predictive.R <RUN_ID>     # ~5 min
#   Rscript run_residual_analysis.R Yasso07 <RUN_ID>   # ~2 min
#
# KEY ARCHITECTURAL DECISIONS (locked in after discussion, Apr 2026):
#   * MULTIPLICATIVE-NORMAL likelihood:
#       SOC_obs ~ Normal(SOC_hat, (SOC_hat * sigma)^2)
#     where sigma = sigma_obs at non-first observations and
#           sigma = sqrt(sigma_obs^2 + sigma_init^2) at the first observation.
#     Reverted from a log-normal workaround that was added when the data had
#     a x10 litter multiplier bug (since fixed in Data_work.R).
#
#   * STEADY-STATE INITIALISATION: xi(mean_climate), NOT mean(xi_annual).
#     Because xi is convex in T and concave in P (Jensen's inequality),
#     averaging xi over years gives a biased estimate of xi at mean climate.
#     The new compute_xi_mean_yasso07 function in the wrapper takes the
#     steady-state climate window, averages T / T_amp / P, then calls
#     compute_xi_yasso07 once. This is wired through the engine via the
#     compute_xi_mean argument to make_likelihood().
#
#   * sigma_obs FIXED across all models (= observed SOC coefficient of
#     variation). NOT calibrated. This ensures structural model differences
#     show up as residuals rather than being absorbed into the error term,
#     which is what we want for the model intercomparison.
#
#   * sigma_init FREE: captures genuine uncertainty about the initial
#     condition (the steady-state assumption is strong; the actual 1985
#     SOC may differ from steady state). Added in quadrature to sigma_obs
#     at the first observation per plot.
#
#   * sigma_input FREE, GLOBAL: a single multiplicative scaling on litter
#     inputs. Captures aleatoric uncertainty in NFI litter estimates
#     (allometric biomass equations have ~20-50% intrinsic uncertainty).
#     Same multiplier applied to all plots, all years -- not stratified by
#     species or region, because that would add identifiability questions
#     that distract from the model intercomparison story.
#
#   * 4 CHAINS in both test and production. Three chains gives noisier
#     R-hat upper bounds; four gives meaningfully tighter convergence
#     diagnostics for negligible extra cost.
#
#   * STEADY_STATE_YEARS = 20: the first 20 years of climate (1985-2004)
#     are used for the steady-state climate scalar. A round number that
#     captures pre-warming-acceleration climate.
# =============================================================================

source("./Calibration_real_data/calibration_engine.R")
source("./Model_functions_real_data/input_compatibility_layer.R")
source("./Model_functions_real_data/Decomposition_functions/Yasso/yasso07_wrapper.R")
dyn.load("./Model_functions_real_data/Decomposition_functions/Yasso/yasso07.so")

library(dplyr)


# =============================================================================
# 0.  Configuration
# =============================================================================
# Two profiles: Mac test (small, fast) and Roihu production (full data, long).
# Switch by uncommenting the relevant block. RUN_ID is generated from the
# current timestamp so multiple runs don't collide on disk.
#
# DIRECTORY STRUCTURE created by this script:
#   ./Calibration/progress_logs/      one log file per chain, written live
#   ./Calibration/diagnostics/Yasso07/   model-specific diagnostic outputs
#   ./Calibration/runs/               saved chain / posterior / metadata
#   ./Calibration/model_inputs/       pre-built input lists for predictive
#   ./Calibration/plots/              (currently unused in calibration; kept
#                                      for legacy compatibility)
# =============================================================================

MODEL_NAME <- "Yasso07"
RUN_ID     <- format(Sys.time(), "%Y%m%d_%H%M%S")

# --- Test (Mac) settings ---
# 20 plots subsampled from the ~400 calib-ready plots. 10000 iterations is
# enough for convergence on this small dataset. Burn-in 2000 = 20%.
# N_LOG = 200: log progress every 200 likelihood evaluations.
N_PLOTS_TEST <- 5L
N_CHAINS     <- 4L
N_ITER       <- 1000L
N_BURNIN     <- 200L
N_LOG        <- 200L

# --- Production (Roihu) settings: uncomment to switch ---
# Full ~3700-plot dataset, longer chains, less frequent logging.
# N_PLOTS_TEST <- NA
# N_CHAINS     <- 4L
# N_ITER       <- 50000L
# N_BURNIN     <- 10000L
# N_LOG        <- 1000L

# CORES_PER_CHAIN: number of cores used by mclapply inside the likelihood for
# the per-plot loop. On Mac, "detectCores() - 1" leaves one core free for the
# OS so the laptop stays responsive. On Roihu with hybrid parallelism (chains
# AND plots in parallel), use 90 to keep total cores at 4*90 = 360 of 384.
CORES_PER_CHAIN <- if (grepl("puhti|mahti", Sys.info()["nodename"])) parallelly::availableCores() else parallel::detectCores() - 1L
# CORES_PER_CHAIN <- 90L   # Roihu

# Steady-state climate window: number of years from the start of the
# transient run used to compute the mean climate for steady-state init.
# 20 years balances "smoothes inter-annual variability" against "still
# pre-recent-warming climate". Sensitivity to this choice should be small
# but is worth checking with 10 / 30 in a sensitivity analysis post-paper.
STEADY_STATE_YEARS <- 20L

# Output directories: model-specific subfolder under diagnostics so multiple
# models can co-exist without filename collisions.
DIR_LOGS   <- "./Calibration_real_data/progress_logs/"
DIR_DIAG   <- file.path("./Calibration_real_data/diagnostics", MODEL_NAME)
DIR_RUNS   <- "./Calibration_real_data/runs/"
DIR_INPUTS <- "./Data/model_inputs/"
DIR_PLOTS  <- "./Calibration_real_data/plots/"
for (d in c(DIR_LOGS, DIR_DIAG, DIR_RUNS, DIR_INPUTS, DIR_PLOTS))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

# run_config and mcmc_settings are passed around between engine functions
# as a single bundle. Keeping them as named lists makes it easy to add new
# fields later without changing every call site.
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
# Yasso07 has 24 parameters total. Of these:
#   FIXED (6):  alpha_A, alpha_W, alpha_E, alpha_N -- decomposition rates
#               for the four litter pools (AWEN). Fixed at published values
#               because rate constants are well-constrained in the literature
#               and calibrating them risks overfitting; structural differences
#               between models should appear in flow fractions and climate
#               response, not in rates.
#               p_H, alpha_H -- humus pool partition fraction and rate.
#               Same rationale.
#
#   FREE (18 from Yasso + 2 nuisance = 20 total):
#               12 inter-pool flow fractions (p_AW, p_AE, ..., p_NE)
#                  -- arranged as four "donor columns" each summing to <= B
#                     where B = 1 - p_H. Use stick-breaking transform to
#                     enforce this compositional constraint while sampling
#                     in unconstrained R^n.
#               beta1, beta2 -- temperature response (linear + quadratic).
#                  Unconstrained; can be any real value.
#               gamma         -- precipitation response.
#               delta1, delta2, r -- woody decomposition response.
#                  Unconstrained.
#               sigma_init    -- initial-condition uncertainty multiplier.
#                  Strictly positive (log transform).
#               sigma_input   -- global litter scaling multiplier.
#                  Strictly positive (log transform); centred at 1.
#
# The param_spec list below is consumed by build_transforms() in the engine,
# which assembles the full transformation machinery (forward, inverse, and
# Jacobian) from these declarations.
# =============================================================================

p_default        <- YASSO07_DEFAULT_PARAMS
FIXED_RATE_NAMES <- c("alpha_A","alpha_W","alpha_E","alpha_N","p_H","alpha_H")
fixed_rates      <- p_default[FIXED_RATE_NAMES]

# Stick-breaking budget B = 1 - p_H. Each donor column's outflow fractions
# must sum to <= B (the rest stays in the pool / is mass-balanced internally).
B <- 1.0 - fixed_rates["p_H"]

# Donor-column groupings: A donates to W, E, N; W to A, E, N; etc.
FRAC_COL_A <- c("p_AW","p_AE","p_AN")
FRAC_COL_W <- c("p_WA","p_WE","p_WN")
FRAC_COL_E <- c("p_EA","p_EW","p_EN")
FRAC_COL_N <- c("p_NA","p_NW","p_NE")

param_spec <- list(
  list(names = FRAC_COL_A,                      type = "stick_break",   budget = B),
  list(names = FRAC_COL_W,                      type = "stick_break",   budget = B),
  list(names = FRAC_COL_E,                      type = "stick_break",   budget = B),
  list(names = FRAC_COL_N,                      type = "stick_break",   budget = B),
  list(names = c("beta1","beta2","gamma"),       type = "unconstrained"),
  list(names = c("delta1","delta2","r"),         type = "unconstrained"),
  list(names = c("sigma_init"),                  type = "log"),
  list(names = c("sigma_input"),                 type = "log")
)

# build_transforms returns four things: forward and inverse transforms, the
# log-Jacobian function (added to the likelihood to keep the posterior
# correct under reparameterisation), and metadata.
transforms       <- build_transforms(param_spec)
to_original      <- transforms$to_original       # unconstrained -> physical
to_unconstrained <- transforms$to_unconstrained   # physical -> unconstrained
log_jacobian     <- transforms$log_jacobian       # log|d_phys / d_uncon|
FREE_NAMES       <- transforms$param_names
N_FREE           <- transforms$n_params

# Default values in PHYSICAL space, used to centre the prior and as the
# starting point for chains. Yasso07 published defaults for the 18 model
# parameters; weakly informed reasonable starts for the 2 nuisance parameters.
free_defaults <- c(
  p_default[c(unlist(lapply(param_spec[1:4], `[[`, "names")))],   # 12 fractions
  p_default[c("beta1","beta2","gamma","delta1","delta2","r")],     # 6 climate/woody
  sigma_init  = 0.10,    # ~10% initial-condition uncertainty
  sigma_input = 1.00     # no scaling adjustment by default
)

# best_x: defaults transformed to UNCONSTRAINED space. This is what the
# sampler sees. The prior is N(best_x, sigma_ppm^2 I) on this scale.
best_x <- to_unconstrained(free_defaults)

# Round-trip sanity: if to_original(to_unconstrained(p)) != p, the transforms
# are wrong (typically a stick-break budget mismatch).
stopifnot(
  max(abs(to_original(best_x)[FREE_NAMES] - free_defaults[FREE_NAMES])) < 1e-10)
message("Transform round-trip: OK")


# =============================================================================
# 2.  Parameter assembly
# =============================================================================
# The MCMC sampler proposes p_free vectors of length 20 (free parameters in
# physical space). The Yasso07 Fortran routine expects a 24-element vector
# in a specific order. This function combines:
#   - fixed_rates (from published values, 6 elements)
#   - p_free (the proposal, 20 elements)
# and reorders to match Yasso07's expected layout.
#
# It also tacks sigma_input onto the end because the engine wrapper functions
# (steady_state_yasso07_engine and yasso07_run_engine, defined in section 4)
# read this from model_params to scale litter inputs.
#
# FUTURE STRATIFICATION: if sigma_input becomes species-stratified, this
# function would receive a stratum argument and pick the appropriate
# multiplier. Currently global so the function is simple.
# =============================================================================

assemble_model_params <- function(p_free) {
  yasso_params <- c(
    fixed_rates,
    p_free[c(unlist(lapply(param_spec[1:6], `[[`, "names")))]
  )
  # Reorder to match Fortran's expected layout (YASSO07_PARAM_NAMES)
  yasso_params <- yasso_params[names(p_default)]
  c(yasso_params, sigma_input = unname(p_free["sigma_input"]))
}


# =============================================================================
# 3.  Data preparation
# =============================================================================
# Reads the data files produced by Data_work.R and builds the per-plot lookup
# structures that the likelihood will read.
#
# Three filtering steps narrow the raw data to calibration-ready plots:
#   (a) site_raw$calib_ready: has SOC obs + has climate + not peatland
#   (b) drop plots with NA in monthly climate (defensive)
#   (c) drop plots with NA in any litter column (defensive)
# After filtering, ~400 plots remain on the test settings.
# =============================================================================

message("Loading data...")

input_raw <- read.csv("./Data/model_inputs/input_raw_monthly.csv")
site_raw  <- read.csv("./Data/model_inputs/site_raw.csv")

# Annual litter total per plot, used in pre-MCMC diagnostic logging.
# (Not used inside MCMC -- the likelihood works with the monthly data.)
litter_cols <- grep("^C_nwl|^C_fwl|^C_cwl", names(input_raw), value = TRUE)
annual_litter <- input_raw %>%
  group_by(plot_id, year) %>%
  summarise(total_litter = sum(across(all_of(litter_cols)), na.rm = TRUE),
            .groups = "drop")

# Filter (a): calib_ready plots only
calib_plots <- site_raw$plot_id[site_raw$calib_ready]
message(sprintf("Calibration-ready plots: %d", length(calib_plots)))
input_calib <- input_raw[input_raw$plot_id %in% calib_plots, ]

# Filter (b): drop plots with NA climate. Should be 0 if Data_work.R worked,
# but defensive in case of edge cases.
n_climate_na <- sum(is.na(input_calib$temp_air))
if (n_climate_na > 0) {
  message(sprintf("WARNING: %d rows with NA climate in calib plots. Filtering...",
                  n_climate_na))
  plots_with_na <- unique(input_calib$plot_id[is.na(input_calib$temp_air)])
  calib_plots   <- setdiff(calib_plots, plots_with_na)
  input_calib   <- input_raw[input_raw$plot_id %in% calib_plots, ]
  message(sprintf("After climate filter: %d plots", length(calib_plots)))
}

# Filter (c): drop plots with NA in any litter column. Same defensive logic.
n_lit_na <- sum(is.na(input_calib[, litter_cols]))
if (n_lit_na > 0) {
  message(sprintf("WARNING: %d NA values in litter columns. Filtering plots...",
                  n_lit_na))
  bad_rows <- !complete.cases(input_calib[, litter_cols])
  plots_with_lit_na <- unique(input_calib$plot_id[bad_rows])
  calib_plots <- setdiff(calib_plots, plots_with_lit_na)
  input_calib <- input_calib[input_calib$plot_id %in% calib_plots, ]
  message(sprintf("After litter filter: %d plots", length(calib_plots)))
}

# Convert monthly inputs to Yasso07's annual timestep format. The compatibility
# layer functions live in input_compatibility_layer.R and aggregate from
# monthly to annual (Yasso) or keep monthly (RothC).
Yasso07_climate <- map_climate_yasso07(input_calib)
Yasso07_inputs  <- map_inputs_yasso07(input_calib)

plots_real <- as.character(unique(Yasso07_climate$plot_id))
message(sprintf("Plots after mapping: %d", length(plots_real)))

# Subsample for Mac testing. Set seed for reproducibility -- this is the same
# 20 plots every time, so successive test runs are directly comparable.
if (!is.na(N_PLOTS_TEST) && N_PLOTS_TEST < length(plots_real)) {
  set.seed(42)   # standardised seed across the project for plot subsetting
  plots <- sort(sample(plots_real, N_PLOTS_TEST))
} else {
  plots <- plots_real
}
plots <- as.character(plots)

# SOC observations: extract the non-NA observations and rank them within plot.
# obs_rank == 1 marks the first observation per plot (1985 in this dataset),
# which gets the inflated sigma (sqrt(sigma_obs^2 + sigma_init^2)) in the
# likelihood to account for steady-state initial condition uncertainty.
SOC_obs_all <- input_calib %>%
  filter(!is.na(soc_obs_tCha)) %>%
  select(plot_id, year, month, soc_obs_tCha) %>%
  arrange(plot_id, year, month) %>%
  group_by(plot_id) %>%
  mutate(obs_rank = row_number()) %>%
  ungroup()

# Fixed observation error: coefficient of variation of SOC across all plots.
# The SAME value is used for all models in the intercomparison. This is
# architecturally critical: if sigma_obs were calibrated, structural model
# differences would be absorbed into the error term and the intercomparison
# would compare error-inflated likelihoods rather than model fits.
obs_cv          <- sd(SOC_obs_all$soc_obs_tCha) / mean(SOC_obs_all$soc_obs_tCha)
sigma_obs_fixed <- obs_cv

message(sprintf("SOC observations: %d  |  obs CV: %.3f  |  sigma_obs fixed: %.3f",
                nrow(SOC_obs_all), obs_cv, sigma_obs_fixed))

# Pre-split the input data by plot. Each plot's data is then accessible as
# climate_by_plot[[pid]] in O(1), avoiding O(n) scans inside the likelihood.
# This is one of the key performance optimisations.
climate_by_plot <- split(Yasso07_climate, Yasso07_climate$plot_id)
inputs_by_plot  <- split(Yasso07_inputs,  Yasso07_inputs$plot_id)

# litter_means: mean litter inputs over the steady-state window per plot.
# Used by yasso07_steady_state to compute C* (the analytical pool distribution
# at equilibrium). Pre-computed once before MCMC because it doesn't depend on
# parameters.
litter_means <- lapply(plots_real, function(pid) {
  inp <- Yasso07_inputs[as.character(Yasso07_inputs$plot_id) == pid, ]
  n_ss   <- min(STEADY_STATE_YEARS, nrow(inp))
  inp_ss <- inp[seq_len(n_ss), ]
  list(
    nwl_mean = colMeans(inp_ss[, c("nwl_A","nwl_W","nwl_E","nwl_N")]),
    fwl_mean = colMeans(inp_ss[, c("fwl_A","fwl_W","fwl_E","fwl_N")]),
    cwl_mean = colMeans(inp_ss[, c("cwl_A","cwl_W","cwl_E","cwl_N")])
  )
})
names(litter_means) <- plots_real

# obs_meta: per-plot list of observation metadata used inside the likelihood.
#   $idx       -- integer indices into climate$year matching each obs year
#   $soc_obs   -- the observed SOC values
#   $is_first  -- TRUE for the first obs per plot (gets inflated sigma)
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

# plot_info: full site_raw row per plot, kept around for the residual analysis
# script. Not consumed by the engine. Stored as a list of one-row lists so
# downstream code can pull out specific covariates without re-reading site_raw.
# The whole row is kept (not a hand-picked subset) so all the new audit-era
# fields (KA, kasvyo_syke, kuvio attributes, derived chemistry, climate
# variability, Köppen, litter quality) automatically flow through.
plot_info <- lapply(plots_real, function(pid) {
  row <- site_raw[as.character(site_raw$plot_id) == pid, ]
  as.list(row[1, , drop = FALSE])
})
names(plot_info) <- plots_real

run_config$n_plots <- length(plots)


# =============================================================================
# 4.  Model interface wrappers
# =============================================================================
# Glue functions that adapt the Yasso07 wrapper functions to the engine's
# generic likelihood factory. Each model in the intercomparison provides
# its own versions of these four functions; the engine treats them as a
# uniform black-box interface.
#
# THE FOUR FUNCTIONS (engine signature in parentheses):
#   compute_xi_yasso07_engine          (compute_xi)
#       Climate modifier time series for the TRANSIENT simulation. One xi
#       per year. Calls compute_xi_yasso07 on the full annual climate.
#
#   compute_xi_mean_yasso07_engine     (compute_xi_mean)
#       SCALAR climate modifier for STEADY-STATE init. Critical: this must
#       compute xi(mean_climate), NOT mean(xi_annual). See the wrapper docs
#       for the full rationale (Jensen's inequality, ~5-15% bias for boreal).
#
#   steady_state_yasso07_engine        (steady_state)
#       Analytical steady-state pool distribution C* given xi_mean and
#       litter_means. Single Fortran call. Multiplies litter inputs by
#       sigma_input before passing to Fortran, so the calibrated multiplier
#       affects both the steady state AND the transient consistently.
#
#   yasso07_run_engine                 (run_model)
#       Full transient simulation 1985 -> 2023 starting from C_init. Single
#       Fortran call returning the full SOC trajectory. Also applies
#       sigma_input to litter inputs so steady state and transient see
#       the same scaling.
# =============================================================================

# Climate modifier for the transient simulation: takes the full per-plot
# climate frame, returns a vector of xi values (one per year).
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

# Steady-state climate modifier: xi(mean_climate), NOT mean(xi_annual).
# This is the architectural fix discussed in the planning thread. The wrapper
# function compute_xi_mean_yasso07 averages T, T_amp, P over the steady-state
# window then calls compute_xi_yasso07 once, which preserves the within-year
# sinusoidal temperature decomposition (T1-T4) while correctly handling the
# inter-annual nonlinearity of xi.
compute_xi_mean_yasso07_engine <- function(clim_ss, model_params) {
  compute_xi_mean_yasso07(
    clim_ss = clim_ss,
    beta1   = model_params["beta1"],
    beta2   = model_params["beta2"],
    gamma   = model_params["gamma"]
  )
}

# Steady-state pool distribution. sigma_input scales the litter inputs so
# that the calibrated multiplier affects steady state and transient runs
# consistently. The wrapper function expects parameters in YASSO07_PARAM_NAMES
# order; assemble_model_params already arranged that.
steady_state_yasso07_engine <- function(model_params, lm, xi_mean) {
  si <- model_params["sigma_input"]
  yasso07_steady_state(
    params   = model_params[YASSO07_PARAM_NAMES],
    nwl_mean = lm$nwl_mean * si,
    fwl_mean = lm$fwl_mean * si,
    cwl_mean = lm$cwl_mean * si,
    xi_mean  = xi_mean
  )
}

# Litter columns in the order Yasso07 expects (NWL, FWL, CWL x AWEN).
LITTER_COLS <- c("nwl_A","nwl_W","nwl_E","nwl_N",
                 "fwl_A","fwl_W","fwl_E","fwl_N",
                 "cwl_A","cwl_W","cwl_E","cwl_N")

# Transient simulation. Same sigma_input scaling as in steady_state for
# consistency. Returns a data frame with year, A/W/E/N/H pools, total_soc,
# and respiration.
yasso07_run_engine <- function(inputs, model_params, C_init, xi_array) {
  si <- model_params["sigma_input"]
  inputs_scaled <- inputs
  inputs_scaled[, LITTER_COLS] <- inputs[, LITTER_COLS] * si
  yasso07_run(inputs_scaled, model_params[YASSO07_PARAM_NAMES], C_init, xi_array)
}


# =============================================================================
# 5.  Prior specification (per-parameter widths)
# =============================================================================
# The prior is N(best_x, diag(sigma_ppm^2)) on UNCONSTRAINED space.
# sigma_ppm is a per-parameter vector, not a scalar. Three widths matter:
#
#   beta1 (=0.2): exp(beta1 * T) with T~10C amplifies 0.2 width to a 95%
#                 multiplicative range of about 0.4x to 7x for the temp
#                 modifier. Wider would let the model produce any modifier
#                 from 0 to 1000x, which is biophysically implausible.
#
#   beta2 (=0.05): the quadratic term. Amplification is even more aggressive
#                 (T^2). Tight prior keeps the modifier on a sensible curve.
#
#   sigma_init (=0.5): on log scale. Physical 95% CI is roughly [0.04, 0.27],
#                 i.e. initial-condition uncertainty between 4% and 27% of
#                 SOC. Reasonable range for boreal forest steady-state error.
#
#   sigma_input (=0.5): on log scale. Physical 95% CI is roughly [0.37, 2.72].
#                 Symmetric on log scale, allows ~2.7x scaling either way to
#                 capture the genuine ~20-50% allometric uncertainty in NFI
#                 litter estimates.
#
#   Stick-break and unconstrained (=1.0): default; weakly informative.
#
# These widths were chosen empirically during script development; they could
# be re-derived using the engine's run_prior_predictive_match() function,
# which finds the smallest sigma giving prior predictive CV >= obs CV.
# =============================================================================

sigma_ppm <- setNames(rep(1.0, N_FREE), FREE_NAMES)
sigma_ppm["beta1"]       <- 0.2
sigma_ppm["beta2"]       <- 0.05
sigma_ppm["sigma_init"]  <- 0.5
sigma_ppm["sigma_input"] <- 0.5    # global multiplier; reflects allometric uncertainty

stopifnot(length(sigma_ppm) == N_FREE)
stopifnot(all(sigma_ppm > 0))
stopifnot(all(names(sigma_ppm) == FREE_NAMES))

message("\nPrior widths (unconstrained space):")
print(round(sigma_ppm, 3))


# =============================================================================
# 6.  Pre-MCMC sanity check
# =============================================================================
# Why this exists:
#   If the published Yasso07 defaults are far from the data, MCMC will
#   struggle to converge and the failure mode is hard to diagnose post-hoc.
#   A 5-second forward run for 4 plots at the prior centre tells us:
#     - is the model numerically stable on this data?
#     - is it within order-of-magnitude of observed SOC?
#     - is steady-state init returning sensible C* values?
#
# Output: <DIR_DIAG>/Yasso07_forward_at_defaults_<RUN_ID>.png
#         A 2x2 grid of plot panels: predicted SOC trajectory (blue line)
#         and observations (red dots) for 4 sample plots, plus the
#         xi_for_ss value used for steady-state init in each panel's legend.
#
# PASS criterion: all four trajectories finite, positive, and within
#                 0.1x to 10x of observations at the obs years.
# =============================================================================

message("\n-------------------------------------------------------------")
message("Pre-MCMC sanity: forward run at defaults...")

set.seed(2025)   # standardised across the project
sanity_pids <- sort(sample(plots, min(4, length(plots))))
sanity_params <- assemble_model_params(free_defaults)

sanity_results <- lapply(sanity_pids, function(pid) {
  clim   <- climate_by_plot[[pid]]
  inputs <- inputs_by_plot[[pid]]
  lm     <- litter_means[[pid]]
  meta   <- obs_meta[[pid]]

  # Step through the same three stages as the likelihood (Stage 1 transient
  # xi, Stage 2 steady state, Stage 3 transient run), at default parameters.
  xi_array  <- compute_xi_yasso07_engine(clim, sanity_params)
  n_ss      <- min(STEADY_STATE_YEARS, nrow(clim))
  xi_for_ss <- compute_xi_mean_yasso07_engine(
    clim[seq_len(n_ss), , drop = FALSE], sanity_params)
  C_init    <- steady_state_yasso07_engine(sanity_params, lm, xi_for_ss)
  run_out   <- yasso07_run_engine(inputs, sanity_params, C_init, xi_array)

  list(pid = pid, run_out = run_out, meta = meta, xi_for_ss = xi_for_ss,
       C_init_total = sum(C_init))
})

# Forward-run sanity status: pass if all four are finite, positive, and within
# 0.1x to 10x of observed SOC at observation years. If this fails, either the
# defaults are off or the data has a unit problem -- inspect the PNG before
# launching MCMC.
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

# Plot: 2x2 panel PNG. xi_for_ss is shown in each legend so any plot with
# a degenerate climate modifier is immediately obvious.
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
  ylim_r <- range(c(soc, obs_soc), na.rm = TRUE)
  plot(yr, soc, type = "l", col = "steelblue", lwd = 2, ylim = ylim_r,
       xlab = "Year", ylab = "SOC (tC/ha)",
       main = sprintf("Plot %s (defaults)", sr$pid))
  if (length(obs_soc) > 0) {
    points(obs_yr, obs_soc, pch = 19, col = "firebrick", cex = 1.3)
  }
  legend("topleft",
         legend = sprintf("xi(mean clim) = %.3f", sr$xi_for_ss),
         bty = "n", cex = 0.8)
}
mtext(sprintf("%s forward run at published defaults  |  %s",
              MODEL_NAME, RUN_ID),
      side = 3, outer = TRUE, line = -1.5, cex = 0.95, font = 2)
dev.off()
message(sprintf("Sanity plot saved: %s", sanity_png))
message(sprintf("Forward-run sanity: %s",
                if (forward_run_pass) "PASS" else "WARN -- inspect plot"))


# =============================================================================
# 7.  Save processed inputs (for predictive script)
# =============================================================================
# Bundle all per-plot lookup structures + run configuration into a single
# .rds file. The predictive script reads this to avoid re-sourcing raw data
# and re-running the data preparation pipeline.
#
# Saving inputs serves a second purpose: provenance. The saved structures
# reflect EXACTLY what the model saw; running the predictive script weeks
# later guarantees consistency with the calibration even if Data_work.R has
# changed in the meantime.
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
# make_likelihood is the engine factory: it returns a byte-compiled function
# of one argument (the unconstrained parameter vector x) that returns the
# log-posterior (= log-likelihood + log-Jacobian for reparameterisation).
#
# All the data structures, model interface functions, and configuration are
# captured in the closure. The sampler only sees a function R^N -> R.
#
# After construction, evaluate the likelihood at the prior centre as a
# sanity check. If it returns -Inf or NaN, something is fundamentally wrong
# (numerical failure in steady-state, pool with negative values, etc.) and
# MCMC cannot proceed.
# =============================================================================

ll_fn <- make_likelihood(
  n_cores         = CORES_PER_CHAIN,
  to_original     = to_original,
  log_jacobian    = log_jacobian,
  assemble_params = assemble_model_params,
  compute_xi      = compute_xi_yasso07_engine,
  compute_xi_mean = compute_xi_mean_yasso07_engine,   # NEW: steady-state fix
  steady_state    = steady_state_yasso07_engine,
  run_model       = yasso07_run_engine,
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

# Pre-MCMC sanity record consumed by write_calibration_report() in section 10.
pre_mcmc_sanity <- list(
  ll_at_defaults     = ll_check,
  n_plots            = length(plots),
  n_obs              = nrow(SOC_obs_all),
  steady_state_years = STEADY_STATE_YEARS,
  forward_run_pass   = forward_run_pass
)

# Prior: per-parameter Gaussian widths. The density function uses element-wise
# vector arithmetic in dnorm. The sampler uses a per-column rnorm because
# rnorm(n*k, mean=vec, sd=vec) recycles element-by-element across the flat
# output, which destroys column-wise structure. The * 0.1 in the sampler
# starts chains in a tight cluster around best_x; DEzs benefits from
# tightly-clustered starts because its proposal distribution is built from
# chain history and aggressive cross-chain proposals early on can be unstable.
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
# DEzs (Differential Evolution with snooker updates) is the sampler used.
# It is gradient-free, population-based, well-suited to correlated 20-30
# dimensional posteriors, and implemented in BayesianTools.
#
# Each chain runs sequentially on Mac (lapply); on Roihu replace the lapply
# inside run_mcmc_chains() with mclapply for outer parallelism.
#
# run_mcmc_chains returns both the chain objects (for diagnostics and saving)
# and a per-chain health summary (initial ll, final ll, -Inf rejection rate,
# wallclock per chain) consumed by the report writer.
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
# run_diagnostics produces:
#   - R-hat (univariate per parameter + multivariate joint)
#   - Effective sample size per parameter
#   - Posterior quantile summary
#   - Trace plots PNG (one panel per HIGHLIGHT parameter)
#   - Marginal plots PNG (posterior vs prior densities)
# All written to DIR_DIAG.
#
# save_results writes three .rds files to DIR_RUNS:
#   - chains_<RUN_ID>.rds      raw BayesianTools chain objects
#   - posterior_<RUN_ID>.rds   physical-space posterior matrix
#   - metadata_<RUN_ID>.rds    run config + summary statistics
#
# write_calibration_report writes one text file to DIR_DIAG with sections:
#   [1] Pre-MCMC sanity
#   [2] Chain health
#   [3] Convergence
#   [4] Posterior shape
# The predictive and residual analysis scripts will append [5] and [6] to
# the same file later.
#
# HIGHLIGHT: the climate response and nuisance parameters are the most
# scientifically interesting and policy-relevant. The 12 stick-break flow
# fractions are also calibrated but produce visually noisy trace/marginal
# plots that obscure the parameters that matter; restricting the highlighted
# set keeps the PNGs concise. All parameters appear in the posterior summary
# table regardless.
# =============================================================================

HIGHLIGHT <- c("beta1","beta2","gamma","delta1","delta2","r",
               "sigma_init","sigma_input")

# group1: the 12 inter-pool flow fractions (stick-break parameterisation).
# Plotted in a 4-column grid (3 rows × 4 cols) so all 12 fit on one page.
# group2: climate response + nuisance parameters — the scientifically most
# interesting set, plotted in a 2-column grid for larger, more readable panels.
# These groupings also define what appears in the KL divergence barplot and
# the marginals PNGs; highlight_params still controls the trace plots.
GROUP1 <- c("p_AW","p_AE","p_AN",
            "p_WA","p_WE","p_WN",
            "p_EA","p_EW","p_EN",
            "p_NA","p_NW","p_NE")
GROUP2 <- HIGHLIGHT   # beta1, beta2, gamma, delta1, delta2, r, sigma_init, sigma_input

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

# Single-file calibration report.
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
message("Next: run_Yasso07_predictive.R with RUN_ID=", RUN_ID)


