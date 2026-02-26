# =============================================================================
# Calibration_Yasso07_production.R
#
# Model-specific calibration script for Yasso07.
# Sources calibration_engine.R for all model-agnostic machinery.
#
# This file contains only what is specific to Yasso07:
#   - Configuration
#   - Parameter specification (param_spec) and defaults
#   - assemble_model_params()
#   - Data preparation
#   - Engine calls
#
# The three model interface functions are the same functions used for
# forward runs and validation:
#   compute_xi_yasso07()    (from yasso07_wrapper.R)
#   yasso07_steady_state()  (from yasso07_wrapper.R)
#   yasso07_run()           (from yasso07_wrapper.R)
#
# ROIHU SWITCH:
#   N_CHAINS        <- 4L
#   CORES_PER_CHAIN <- 90L
#   replace lapply with mclapply(mc.cores = N_CHAINS) in engine
# =============================================================================

source("./Calibration/calibration_engine.R")
source("./Model_functions/input_compatibility_layer.R")
source("./Model_functions/Decomposition_functions/Yasso/yasso07_wrapper.R")
dyn.load("./Model_functions/Decomposition_functions/Yasso/yasso07.so")

library(dplyr)


# =============================================================================
# 0.  Configuration
# =============================================================================

MODEL_NAME  <- "Yasso07"
RUN_ID      <- format(Sys.time(), "%Y%m%d_%H%M%S")

# --- Plot replication (testing) ---
N_PLOTS_TEST <- 25L        # set to NA for production

# --- MCMC ---
N_CHAINS   <- 3L
N_ITER     <- 10000L
N_BURNIN   <- 2000L
N_LOG      <- 200L

# --- Parallelism ---
CORES_PER_CHAIN <- parallel::detectCores() - 1L
# CORES_PER_CHAIN <- 90L   # Roihu

# --- Output directories ---
DIR_LOGS <- "./Calibration/progress_logs/"
DIR_DIAG <- "./Calibration/diagnostics/"
DIR_RUNS <- "./Calibration/runs/"
DIR_INPUTS <- "./Calibration/model_inputs/"
for (d in c(DIR_LOGS, DIR_DIAG, DIR_RUNS))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)


# --- Run configuration (passed to all engine functions) ---
# n_plots is not yet known -- added at end of Section 3
run_config <- list(
  MODEL_NAME      = MODEL_NAME,
  RUN_ID          = RUN_ID,
  DIR_LOGS        = DIR_LOGS,
  DIR_DIAG        = DIR_DIAG,
  DIR_RUNS        = DIR_RUNS,
  DIR_INPUTS      = DIR_INPUTS,
  CORES_PER_CHAIN = CORES_PER_CHAIN
)

# --- MCMC settings (passed to engine functions that need iteration counts) ---
mcmc_settings <- list(
  N_CHAINS = N_CHAINS,
  N_ITER   = N_ITER,
  N_BURNIN = N_BURNIN,
  N_LOG    = N_LOG
)


message("=============================================================")
message(sprintf("  %s Calibration  |  Run: %s", MODEL_NAME, RUN_ID))
message("=============================================================\n")


# =============================================================================
# 1.  Parameter specification
#
#     Declares which parameters are fixed (held at published values throughout
#     MCMC) and which are free (calibrated). For each free parameter, declares
#     its constraint type so the engine can build the correct transformation.
#
#     No computation happens here except building the transform functions and
#     setting the MCMC starting point. The actual transformations happen later:
#     to_unconstrained() is called once here to initialise best_x;
#     to_original() is called millions of times inside the likelihood.
# =============================================================================

# Get the parameter values from each model definition
p_default        <- YASSO07_DEFAULT_PARAMS

# Fixed parameters: held at published values, never sampled.
# Calibrating decomposition rates jointly with transfer fractions is
# not identifiable -- we would be fitting the same signal twice.
FIXED_RATE_NAMES <- c("alpha_A","alpha_W","alpha_E","alpha_N","p_H","alpha_H")
fixed_rates      <- p_default[FIXED_RATE_NAMES]

# Budget for stick-breaking: the total fraction available for inter-pool
# transfers. Whatever doesn't go to humus (p_H) must be distributed among
# the other pools. This is the length of the "stick" to be broken.
B <- 1.0 - fixed_rates["p_H"]

# Free parameter groups: inter-pool transfer fractions.
# Each column (e.g. p_AW, p_AE, p_AN) represents the fractions of pool A's
# outflow going to pools W, E, N respectively. They must be positive and
# sum to <= B. One group per source pool.
FRAC_COL_A <- c("p_AW","p_AE","p_AN")
FRAC_COL_W <- c("p_WA","p_WE","p_WN")
FRAC_COL_E <- c("p_EA","p_EW","p_EN")
FRAC_COL_N <- c("p_NA","p_NW","p_NE")

# Parameter specification table.
# Each entry tells the engine the name(s) and constraint type of a group:
#
#   stick_break : fractions that must be positive and sum to <= budget.
#                 DEzs samples three unconstrained reals; the stick-breaking
#                 transform converts them to valid fractions by sequentially
#                 breaking off pieces of a stick of length B. Guarantees
#                 positivity and sum constraint by construction.
#
#   unconstrained: parameters with no constraint (any real value).
#                 Identity transform -- sampler sees the same values the
#                 model sees. Climate response parameters (beta, gamma) and
#                 size modifiers (delta, r) can legitimately be negative
#                 or zero, so no transformation is needed.
#
#   log         : parameters that must be strictly positive (> 0).
#                 DEzs samples log(sigma); model receives exp(log(sigma)).
#                 Used for observation error terms where negative values
#                 are physically impossible.
param_spec <- list(
  list(names = FRAC_COL_A,                           type = "stick_break",   budget = B),
  list(names = FRAC_COL_W,                           type = "stick_break",   budget = B),
  list(names = FRAC_COL_E,                           type = "stick_break",   budget = B),
  list(names = FRAC_COL_N,                           type = "stick_break",   budget = B),
  list(names = c("beta1","beta2","gamma"),            type = "unconstrained"),
  list(names = c("delta1","delta2","r"),              type = "unconstrained"),
  list(names = c("sigma_init"),                      type = "log")           # sigma_obs fixed; sigma_init free
)

# Build the three transform functions from the spec above.
# After this call:
#   to_original(x)       -- converts a sampler proposal (unconstrained) to
#                           physical parameter values the model can use
#   to_unconstrained(p)  -- converts physical values to unconstrained space;
#                           used once below to initialise best_x
#   log_jacobian(x, p)   -- correction factor for the space warping.
#                           When you change variables, regions of parameter
#                           space stretch or compress. Without this correction
#                           the sampler would over- or under-sample certain
#                           regions. The Jacobian is the volume distortion
#                           factor; adding its log to the likelihood makes
#                           the sampler explore the correct posterior geometry.
transforms       <- build_transforms(param_spec)
to_original      <- transforms$to_original
to_unconstrained <- transforms$to_unconstrained
log_jacobian     <- transforms$log_jacobian
FREE_NAMES       <- transforms$param_names
N_FREE           <- transforms$n_params

# Collect the free parameter defaults in the same order as param_spec.
# sigma_init has no published default -- 10% multiplicative initialisation
# error is a reasonable starting value.
# sigma_obs is NOT included: it is fixed to obs_cv, not sampled.
free_defaults <- c(
  p_default[c(unlist(lapply(param_spec[1:4], `[[`, "names")))],
  p_default[c("beta1","beta2","gamma","delta1","delta2","r")],
  sigma_init = 0.10     # sigma_obs is fixed (= obs_cv), not sampled
)

# Convert published defaults to unconstrained space.
# This is the point where DEzs chains will start, and the centre of
# the Gaussian prior. Starting at the published defaults is a sensible
# choice: they represent the best prior knowledge of Yasso07 parameters
# from the original calibration literature.
best_x <- to_unconstrained(free_defaults)

# Verify the round-trip: physical -> unconstrained -> physical should
# recover the original values to machine precision. If this fails,
# something is wrong with the transform definitions or the budget B.
stopifnot(
  max(abs(to_original(best_x)[FREE_NAMES] - free_defaults[FREE_NAMES])) < 1e-10)
message("Transform round-trip: OK")



# =============================================================================
# 2.  Parameter assembly
#
#     Combines free (calibrated) and fixed parameters into the full
#     named vector expected by the Yasso07 wrapper functions.
# =============================================================================

assemble_model_params <- function(p_free) {
  full <- c(fixed_rates,
            p_free[c(unlist(lapply(param_spec[1:6], `[[`, "names")))])
  full[names(p_default)]
}


# =============================================================================
# 3.  Data preparation
# =============================================================================

message("Loading data...")
input_raw       <- read.csv("synthetic_input_data_template.csv")
site_raw        <- read.csv("synthetic_site_data_template.csv")
Yasso07_climate <- map_climate_yasso07(input_raw)
Yasso07_inputs  <- map_inputs_yasso07(input_raw)
plots_real      <- unique(Yasso07_climate$plot_id)


# -------------------------------------------------------------------------------------
# TRICK FOR VALIDATION ON TEST DATA!!!
# plots is a vector of plot IDs passed to the likelihood -- either the full
# production set or N_PLOTS_TEST replications of plots_real for timing tests.
# -------------------------------------------------------------------------------------

plots <- if (!is.na(N_PLOTS_TEST)) {
  message(sprintf("TEST MODE: %d plots replicated to %d",
                  length(plots_real), N_PLOTS_TEST))
  rep(plots_real, length.out = N_PLOTS_TEST)
} else {
  message(sprintf("PRODUCTION: %d plots", length(plots_real)))
  plots_real
}

SOC_obs_all <- input_raw %>%
  filter(!is.na(soc_obs_tCha)) %>%
  select(plot_id, year, month, soc_obs_tCha) %>%
  arrange(plot_id, year, month) %>%
  group_by(plot_id) %>%
  mutate(obs_rank = row_number()) %>%
  ungroup()

obs_cv <- sd(SOC_obs_all$soc_obs_tCha) / mean(SOC_obs_all$soc_obs_tCha)

# sigma_obs is fixed to obs_cv: the proportional observation error is set equal
# to the observed SOC coefficient of variation. This is the same value used for
# the prior predictive matching width, and is applied identically across all
# models so that structural differences are not absorbed into the error term.
sigma_obs_fixed <- obs_cv

message(sprintf("SOC observations: %d  |  obs CV: %.3f  |  sigma_obs fixed: %.3f\n",
                nrow(SOC_obs_all), obs_cv, sigma_obs_fixed))

# Pre-split for O(1) lookup inside likelihood
climate_by_plot <- split(Yasso07_climate, Yasso07_climate$plot_id)
inputs_by_plot  <- split(Yasso07_inputs,  Yasso07_inputs$plot_id)

litter_means <- lapply(plots_real, function(pid) {
  inp <- Yasso07_inputs[Yasso07_inputs$plot_id == pid, ]
  list(
    nwl_mean = colMeans(inp[, c("nwl_A","nwl_W","nwl_E","nwl_N")]),
    fwl_mean = colMeans(inp[, c("fwl_A","fwl_W","fwl_E","fwl_N")]),
    cwl_mean = colMeans(inp[, c("cwl_A","cwl_W","cwl_E","cwl_N")])
  )
})
names(litter_means) <- plots_real

obs_meta <- lapply(plots_real, function(pid) {
  clim     <- Yasso07_climate[Yasso07_climate$plot_id == pid, ]
  obs_plot <- SOC_obs_all[SOC_obs_all$plot_id == pid, ]
  list(
    idx      = match(obs_plot$year, clim$year),
    soc_obs  = obs_plot$soc_obs_tCha,
    is_first = obs_plot$obs_rank == 1L
  )
})
names(obs_meta) <- plots_real

# Finalise run_config now that n_plots is known
run_config$n_plots <- length(plots)


# =============================================================================
# 4.  Steady state wrapper
#
#     Adapts yasso07_steady_state() to the engine's calling convention:
#       steady_state(model_params, litter_means_list, xi_mean)
# =============================================================================

steady_state_yasso07 <- function(model_params, lm, xi_mean) {
  yasso07_steady_state(
    params   = model_params,
    nwl_mean = lm$nwl_mean,
    fwl_mean = lm$fwl_mean,
    cwl_mean = lm$cwl_mean,
    xi_mean  = xi_mean
  )
}


# =============================================================================
# 5.  Prior predictive matching
# =============================================================================
#
# Goal: find sigma_ppm, the width of the Gaussian prior in unconstrained space.
# We want the prior to be weakly informative -- wide enough that it doesn't
# dominate the posterior, but not so wide that the sampler wastes time in
# physically implausible regions.
#
# The criterion: draw parameters from the prior, run the model, and check
# whether the resulting SOC predictions are as variable as the observed data
# (measured by CV = sd/mean). If the prior predictive CV >= observed CV,
# the prior is broad enough. sigma_ppm is the smallest value that achieves this.
#
# We use a single plot for speed -- we only need the marginal prior predictive
# variance, not the full joint distribution across plots.

run_one_plot_yasso07 <- function(p_free) {
  # Assemble full parameter vector (free + fixed rates) for the model
  model_params <- assemble_model_params(p_free)
  pid          <- plots_real[1]            # representative plot
  clim         <- climate_by_plot[[pid]]
  inputs       <- inputs_by_plot[[pid]]
  lm           <- litter_means[[pid]]
  meta         <- obs_meta[[pid]]
  
  # Stage 1: compute climate modifier time series for this parameter proposal
  xi_array <- compute_xi_yasso07(clim, model_params)
  xi_mean  <- mean(xi_array)
  if (!is.finite(xi_mean) || xi_mean <= 0) return(NA_real_)  # degenerate proposal
  
  # Stage 2: analytical steady-state initialisation (avoids spinup)
  C_init <- steady_state_yasso07(model_params, lm, xi_mean)
  if (any(!is.finite(C_init)) || any(C_init < 0)) return(NA_real_)  # invalid pools
  
  # Stage 3: run transient simulation, return SOC at first observation year
  # (index [1]: we only need one scalar prediction per prior draw for the CV)
  run_out <- yasso07_run(inputs, model_params, C_init, xi_array)
  soc     <- run_out$total_soc[meta$idx[1]]
  if (!is.finite(soc) || soc <= 0) NA_real_ else soc
}

set.seed(42)   # reproducibility: same sigma_ppm across reruns
sigma_ppm <- run_prior_predictive_match(
  best_x       = best_x,          # prior centre = published defaults (unconstrained)
  to_original  = to_original,     # needed inside the engine to back-transform draws
  run_one_plot = run_one_plot_yasso07,
  obs_cv       = obs_cv           # target CV derived from Finnish inventory SOC data
)

# =============================================================================
# 5.1  Adapting the model interfaces to the sampler
# =============================================================================

# Adapter: wraps compute_xi_yasso07 to match the engine's calling convention
# compute_xi(clim_df, params_vector) -> xi_array
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


# =============================================================================
# 5.5  Save processed inputs
# =============================================================================

save_inputs(
  inputs_list = list(
    climate_by_plot = climate_by_plot,
    inputs_by_plot  = inputs_by_plot,
    litter_means    = litter_means,
    obs_meta        = obs_meta,
    plots           = plots,
    plots_real      = plots_real,
    obs_cv          = obs_cv,
    sigma_obs_fixed = sigma_obs_fixed,
    sigma_ppm       = sigma_ppm,
    free_defaults   = free_defaults,
    best_x          = best_x
  ),
  run_config = run_config
)

# =============================================================================
# 6.  Build likelihood and prior
# =============================================================================

ll_fn <- make_likelihood(
  n_cores         = CORES_PER_CHAIN,
  to_original     = to_original,
  log_jacobian    = log_jacobian,
  assemble_params = assemble_model_params,
  compute_xi      = compute_xi_yasso07_engine,
  steady_state    = steady_state_yasso07,
  run_model       = yasso07_run,
  sigma_obs_fixed = sigma_obs_fixed,
  plots           = plots,
  climate_by_plot = climate_by_plot,
  inputs_by_plot  = inputs_by_plot,
  litter_means    = litter_means,
  obs_meta        = obs_meta
)

# Sanity check at defaults
ll_check <- ll_fn(best_x)
message(sprintf("Likelihood at defaults: %.3f", ll_check))
if (!is.finite(ll_check)) stop("FAIL: likelihood not finite at defaults.")
message("Likelihood check: PASS\n")

prior <- createPrior(
  density = function(x) sum(dnorm(x, mean = best_x, sd = sigma_ppm, log = TRUE)),
  sampler = function(n = 1)
    matrix(rnorm(n * N_FREE, mean = best_x, sd = sigma_ppm * 0.1),
           nrow = n, ncol = N_FREE),
  lower = rep(-Inf, N_FREE),
  upper = rep( Inf, N_FREE)
)


# =============================================================================
# 7.  Run MCMC
# =============================================================================

# mcmc_settings <- list(
#   N_CHAINS = N_CHAINS, N_ITER = N_ITER,
#   N_BURNIN = N_BURNIN, N_LOG  = N_LOG)
# 
# run_config <- list(
#   MODEL_NAME      = MODEL_NAME,
#   RUN_ID          = RUN_ID,
#   DIR_LOGS        = DIR_LOGS,
#   DIR_DIAG        = DIR_DIAG,
#   DIR_RUNS        = DIR_RUNS,
#   DIR_INPUTS      = DIR_INPUTS,
#   CORES_PER_CHAIN = CORES_PER_CHAIN,
#   n_plots         = length(plots)
# )

t_run <- system.time({
  chain_results <- run_mcmc_chains(
    ll_fn         = ll_fn,
    prior         = prior,
    best_x        = best_x,
    param_names   = FREE_NAMES,
    mcmc_settings = mcmc_settings,
    run_config    = run_config)
})["elapsed"]

message(sprintf("\nAll chains complete. Wallclock: %.1f min\n", t_run / 60))


# =============================================================================
# 8.  Diagnostics and save
# =============================================================================

# Parameters to highlight in trace/marginal plots
HIGHLIGHT <- c("beta1","beta2","gamma","delta1","delta2","r",
               "sigma_init")   # sigma_obs removed (fixed, not sampled)

diag_out <- run_diagnostics(
  chain_results    = chain_results,
  to_original      = to_original,
  param_names      = FREE_NAMES,
  best_x           = best_x,
  sigma_ppm        = sigma_ppm,
  run_config       = run_config,
  mcmc_settings    = mcmc_settings,    # add this
  highlight_params = HIGHLIGHT
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

message("\n=============================================================")
message(sprintf("  %s complete  |  Run: %s  |  %.1f min",
                MODEL_NAME, RUN_ID, t_run / 60))
message("=============================================================\n")