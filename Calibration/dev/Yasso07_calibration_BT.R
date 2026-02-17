# =============================================================================
# Calibration_Yasso07_BT.R
#
# Bayesian calibration of Yasso07 using BayesianTools (DEzs sampler).
#
# What this script does, in plain language:
#
#   We have ~3700 forest plots with measured soil carbon (SOC) at two or more
#   time points. We want to find parameter values for Yasso07 that make the
#   model's predictions match those measurements as well as possible -- but
#   rather than finding a single "best" set of parameters, we want the full
#   probability distribution over plausible parameter values. That distribution
#   is the "posterior", and MCMC is how we explore it.
#
#   The sampler works in "unconstrained space": all parameters are transformed
#   so they can range freely from -Inf to +Inf. This avoids hard boundaries
#   that would cause the sampler to get stuck. Before calling the model we
#   always back-transform to physical (original) space.
#
# Design principles (unchanged from standalone MH version):
#
#   1. FIXED k: Decomposition rates fixed at published values.
#      These define what each pool IS chemically. We don't touch them.
#
#   2. FREE parameters: Transfer fractions (12), climate sensitivity
#      (beta1, beta2, gamma), woody modifiers (delta1, delta2, r),
#      and two nuisance sigmas that describe observation noise.
#
#   3. TRANSFORMED sampling: Sampler works in unconstrained space.
#      Back-transformation happens in R before any model call.
#
#   4. WHY BayesianTools / DEzs instead of standalone MH:
#      The standalone Metropolis-Hastings sampler used isotropic (round)
#      proposals and needed manual tuning. DEzs (Differential Evolution
#      with snooker updates) runs multiple chains simultaneously and
#      uses the spread of those chains to automatically propose good moves.
#      It handles correlated parameters and different parameter scales
#      much better than simple MH, and requires almost no tuning.
#
# =============================================================================

library(dplyr)
library(BayesianTools)   # provides createPrior, createBayesianSetup, runMCMC

source("./Model_functions/input_compatibility_layer.R")
source("./Model_functions/Decomposition_functions/Yasso/yasso07_wrapper.R")
dyn.load("./Model_functions/Decomposition_functions/Yasso/yasso07.so")


# =============================================================================
# 1.  Load and map inputs
#     (identical to standalone version)
# =============================================================================

input_raw <- read.csv("synthetic_input_data_template.csv")
site_raw  <- read.csv("synthetic_site_data_template.csv")

Yasso07_climate <- map_climate_yasso07(input_raw)
Yasso07_inputs  <- map_inputs_yasso07(input_raw)
plots           <- unique(Yasso07_climate$plot_id)


# =============================================================================
# 2.  Extract SOC observations
# =============================================================================

SOC_obs_all <- input_raw %>%
  filter(!is.na(soc_obs_tCha)) %>%
  select(plot_id, year, month, soc_obs_tCha) %>%
  arrange(plot_id, year, month) %>%
  group_by(plot_id) %>%
  mutate(obs_rank = row_number()) %>%
  ungroup()


# =============================================================================
# 3.  Pre-compute plot-level quantities
#
#     We do this ONCE here, outside the MCMC loop, because these quantities
#     don't depend on parameters -- they're just derived from the input data.
#     Doing it inside the loop would waste time recomputing them millions of
#     times.
# =============================================================================

# For each plot: mean litter inputs by pool (A, W, E, N) across all years.
# These means are used to compute the steady-state initial conditions.
litter_means <- lapply(plots, function(pid) {
  inputs <- Yasso07_inputs[Yasso07_inputs$plot_id == pid, ]
  list(
    nwl_mean = colMeans(inputs[, c("nwl_A","nwl_W","nwl_E","nwl_N")]),
    fwl_mean = colMeans(inputs[, c("fwl_A","fwl_W","fwl_E","fwl_N")]),
    cwl_mean = colMeans(inputs[, c("cwl_A","cwl_W","cwl_E","cwl_N")])
  )
})
names(litter_means) <- plots

# For each plot: which time steps correspond to observations, and which
# observation is the first one (used for a wider uncertainty on initial SOC).
obs_meta <- lapply(plots, function(pid) {
  clim     <- Yasso07_climate[Yasso07_climate$plot_id == pid, ]
  obs_plot <- SOC_obs_all[SOC_obs_all$plot_id == pid, ]
  list(
    idx      = match(obs_plot$year, clim$year),   # row index in model output
    soc_obs  = obs_plot$soc_obs_tCha,             # observed SOC values
    is_first = obs_plot$obs_rank == 1L            # TRUE for first observation
  )
})
names(obs_meta) <- plots


# =============================================================================
# 4.  Parameter structure
# =============================================================================

p_default <- YASSO07_DEFAULT_PARAMS

# Fixed parameters: decomposition rates encode pool chemical identity.
# We never move these during calibration.
FIXED_RATE_NAMES <- c("alpha_A","alpha_W","alpha_E","alpha_N","p_H","alpha_H")
fixed_rates      <- p_default[FIXED_RATE_NAMES]

# Free parameter groups and their names
FRAC_COL_A     <- c("p_AW","p_AE","p_AN")      # transfer fractions out of pool A
FRAC_COL_W     <- c("p_WA","p_WE","p_WN")      # transfer fractions out of pool W
FRAC_COL_E     <- c("p_EA","p_EW","p_EN")      # transfer fractions out of pool E
FRAC_COL_N     <- c("p_NA","p_NW","p_NE")      # transfer fractions out of pool N
FRAC_NAMES     <- c(FRAC_COL_A, FRAC_COL_W, FRAC_COL_E, FRAC_COL_N)

CLIMATE_NAMES  <- c("beta1","beta2","gamma")    # temperature and precip sensitivity
MODIFIER_NAMES <- c("delta1","delta2","r")      # woody litter size effect
NUISANCE_NAMES <- c("sigma_obs","sigma_init")   # observation and init uncertainty

FREE_NAMES <- c(FRAC_NAMES, CLIMATE_NAMES, MODIFIER_NAMES, NUISANCE_NAMES)
N_FREE     <- length(FREE_NAMES)   # 20 free parameters total

# Default values in original (physical) space, used as starting point
free_defaults <- c(
  p_default[FRAC_NAMES],
  p_default[CLIMATE_NAMES],
  p_default[MODIFIER_NAMES],
  sigma_obs  = 0.05,
  sigma_init = 0.10
)

# Index helpers: which positions in the unconstrained vector correspond to what
IDX_FRACS      <- 1:12
IDX_COL_A      <- 1:3
IDX_COL_W      <- 4:6
IDX_COL_E      <- 7:9
IDX_COL_N      <- 10:12
IDX_CLIMATE    <- 13:15
IDX_MODIFIERS  <- 16:18
IDX_SIGMA_OBS  <- 19L
IDX_SIGMA_INIT <- 20L

# "Budget" for each column: the fraction of carbon NOT going to humus.
# Transfer fractions across AWEN pools cannot exceed this.
B <- 1.0 - fixed_rates["p_H"]


# =============================================================================
# 5.  Transformation functions
#
#     WHY TRANSFORM?
#     Transfer fractions must satisfy: p1 + p2 + p3 <= B (the budget).
#     Sigmas must be positive. If we sampled these directly, the sampler
#     would keep hitting the boundaries and most proposals would be rejected.
#
#     STICK-BREAKING for transfer fractions:
#     Imagine breaking a stick of length B into pieces:
#       - First piece:  p1 = logistic(v1) * B
#         (logistic maps any real number to (0,1), so p1 is in (0, B))
#       - Second piece: p2 = logistic(v2) * (B - p1)
#         (takes a fraction of what's left after p1)
#       - Third piece:  p3 = logistic(v3) * (B - p1 - p2)
#         (takes a fraction of what's left after p1 and p2)
#       - Remainder:    respired = B - p1 - p2 - p3  (always >= 0)
#     This guarantees p1+p2+p3 <= B regardless of what v1, v2, v3 are.
#
#     LOG transform for sigmas:
#     log(sigma) can be any real number, but back-transforming with exp()
#     always gives a positive sigma. Simple and effective.
# =============================================================================

logit     <- function(p) log(p / (1 - p))
inv_logit <- function(x) 1 / (1 + exp(-x))

# Forward: unconstrained v (length 3) -> physical fractions (length 3)
stick_break <- function(v, budget) {
  p    <- numeric(3)
  p[1] <- inv_logit(v[1]) * budget
  p[2] <- inv_logit(v[2]) * (budget - p[1])
  p[3] <- inv_logit(v[3]) * (budget - p[1] - p[2])
  p
}

# Inverse: physical fractions (length 3) -> unconstrained v (length 3)
stick_break_inv <- function(p, budget) {
  v    <- numeric(3)
  v[1] <- logit(p[1] / budget)
  v[2] <- logit(p[2] / (budget - p[1]))
  v[3] <- logit(p[3] / (budget - p[1] - p[2]))
  v
}

# Log-Jacobian of the stick-breaking transformation for one column.
#
# WHAT IS THE JACOBIAN?
# When you sample in v-space and transform to p-space, the transformation
# stretches and squishes the probability density unevenly -- fractions near
# 0 or B get compressed, fractions in the middle get spread out. If you
# ignore this, your prior is accidentally non-uniform in a way you didn't
# intend. The Jacobian corrects for this stretching.
#
# Mathematically it's the log-determinant of the matrix of partial derivatives
# d(p1,p2,p3)/d(v1,v2,v3). Because stick-breaking is triangular, the
# determinant simplifies to a product of three terms.
log_jac_stick <- function(p, budget) {
  rem  <- c(budget, budget - p[1], budget - p[1] - p[2])
  frac <- p / rem
  # Each term is log of: logistic(v_i) * (1 - logistic(v_i)) * rem_i
  # which equals log(frac * (1-frac) * rem)
  sum(log(frac * (1 - frac) * rem))
}

# Convert from physical (original) space to unconstrained (sampling) space
to_unconstrained <- function(free_params) {
  x <- numeric(N_FREE)
  names(x) <- FREE_NAMES
  x[IDX_COL_A]      <- stick_break_inv(free_params[FRAC_COL_A], B)
  x[IDX_COL_W]      <- stick_break_inv(free_params[FRAC_COL_W], B)
  x[IDX_COL_E]      <- stick_break_inv(free_params[FRAC_COL_E], B)
  x[IDX_COL_N]      <- stick_break_inv(free_params[FRAC_COL_N], B)
  x[IDX_CLIMATE]    <- free_params[CLIMATE_NAMES]   # unbounded, no transform
  x[IDX_MODIFIERS]  <- free_params[MODIFIER_NAMES]  # unbounded, no transform
  x[IDX_SIGMA_OBS]  <- log(free_params["sigma_obs"])
  x[IDX_SIGMA_INIT] <- log(free_params["sigma_init"])
  x
}

# Convert from unconstrained (sampling) space back to physical (original) space
to_original <- function(x) {
  p <- numeric(N_FREE)
  names(p) <- FREE_NAMES
  p[FRAC_COL_A]     <- stick_break(x[IDX_COL_A], B)
  p[FRAC_COL_W]     <- stick_break(x[IDX_COL_W], B)
  p[FRAC_COL_E]     <- stick_break(x[IDX_COL_E], B)
  p[FRAC_COL_N]     <- stick_break(x[IDX_COL_N], B)
  p[CLIMATE_NAMES]  <- x[IDX_CLIMATE]
  p[MODIFIER_NAMES] <- x[IDX_MODIFIERS]
  p["sigma_obs"]    <- exp(x[IDX_SIGMA_OBS])
  p["sigma_init"]   <- exp(x[IDX_SIGMA_INIT])
  p
}

# Assemble the full model parameter vector (fixed + free, in correct order)
assemble_model_params <- function(p_free) {
  full <- c(fixed_rates,
            p_free[c(FRAC_NAMES, CLIMATE_NAMES, MODIFIER_NAMES)])
  full[names(p_default)]
}

# Starting point: published defaults transformed to unconstrained space
best_x <- to_unconstrained(free_defaults)

# Sanity check: transforming back should recover the originals exactly
stopifnot(max(abs(to_original(best_x)[FREE_NAMES] -
                    free_defaults[FREE_NAMES])) < 1e-10)
message("Transformation round-trip check: OK")


# =============================================================================
# 6.  Prior predictive matching
#
#     HOW WIDE SHOULD THE PRIOR BE?
#     We want the prior to allow plausible variation in SOC predictions, but
#     not be so wide that the sampler wastes time in impossible regions.
#
#     We match the prior width (sigma_ppm) so that the coefficient of
#     variation (CV = sd/mean) of prior-predicted SOC values across a plot
#     matches the observed CV across all plots. This is a data-driven way
#     to set the prior width without making it too tight or too loose.
#
#     We search over a grid of sigma values, draw random parameter sets from
#     each, run the model, and find the sigma that hits the target CV.
# =============================================================================

obs_cv <- sd(SOC_obs_all$soc_obs_tCha) / mean(SOC_obs_all$soc_obs_tCha)
message(sprintf("\nObserved SOC CV: %.3f (%.1f%%)", obs_cv, obs_cv * 100))

run_prior_sample_soc <- function(sigma_val) {
  # Draw random parameters by perturbing best_x in unconstrained space
  x_sample     <- best_x + rnorm(N_FREE, 0, sigma_val)
  p_free       <- to_original(x_sample)
  model_params <- assemble_model_params(p_free)
  
  # Run on the first plot as a representative sample
  pid    <- plots[1]
  clim   <- Yasso07_climate[Yasso07_climate$plot_id == pid, ]
  inputs <- Yasso07_inputs[Yasso07_inputs$plot_id   == pid, ]
  lm     <- litter_means[[pid]]
  meta   <- obs_meta[[pid]]
  
  xi_array <- tryCatch(
    compute_xi_yasso07(
      temp_mean = clim$temp_mean,
      temp_amp  = clim$temp_amplitude,
      precip    = clim$precip,
      beta1     = model_params["beta1"],
      beta2     = model_params["beta2"],
      gamma     = model_params["gamma"]
    ),
    error = function(e) NULL
  )
  if (is.null(xi_array)) return(NA_real_)
  xi_mean <- mean(xi_array)
  if (!is.finite(xi_mean) || xi_mean <= 0) return(NA_real_)
  
  C_init <- tryCatch(
    yasso07_steady_state(
      params   = model_params,
      nwl_mean = lm$nwl_mean,
      fwl_mean = lm$fwl_mean,
      cwl_mean = lm$cwl_mean,
      xi_mean  = xi_mean
    ),
    error = function(e) NULL
  )
  if (is.null(C_init) || any(!is.finite(C_init)) || any(C_init < 0)) return(NA_real_)
  
  run_out <- tryCatch(
    yasso07_run(
      input_df = inputs,
      params   = model_params,
      C_init   = C_init,
      xi_array = xi_array
    ),
    error = function(e) NULL
  )
  if (is.null(run_out)) return(NA_real_)
  
  soc <- run_out$total_soc[meta$idx[1]]
  if (!is.finite(soc) || soc <= 0) return(NA_real_)
  soc
}

message("Running prior predictive matching...")
set.seed(42)
sigma_grid <- seq(0.1, 3.0, by = 0.1)
n_ppm      <- 200L

ppm_cv <- sapply(sigma_grid, function(sigma_val) {
  message(sprintf("  sigma = %.1f ...", sigma_val), appendLF = FALSE)
  preds <- vapply(seq_len(n_ppm), function(i) run_prior_sample_soc(sigma_val),
                  numeric(1))
  valid <- preds[is.finite(preds) & preds > 0]
  cv    <- if (length(valid) < 20L) NA_real_ else sd(valid) / mean(valid)
  message(sprintf(" CV = %.3f (n_valid = %d/%d)", cv, length(valid), n_ppm))
  cv
})

sigma_ppm <- sigma_grid[which(ppm_cv >= obs_cv)[1]]
if (is.na(sigma_ppm)) {
  message("WARNING: CV never reached target -- using sigma_ppm = 1.0")
  sigma_ppm <- 1.0
}
message(sprintf("\nPrior predictive matching complete:"))
message(sprintf("  Target CV    : %.3f", obs_cv))
message(sprintf("  Matched sigma: %.1f", sigma_ppm))
message(sprintf("  Achieved CV  : %.3f", ppm_cv[which(sigma_grid == sigma_ppm)]))


# =============================================================================
# 7.  Log-likelihood function for BayesianTools
#
#     BayesianTools separates the posterior into:
#       log-posterior = log-prior + log-likelihood
#     and handles the prior itself. We only need to supply the log-likelihood.
#
#     IMPORTANT: We also add the log-Jacobian here.
#     The Jacobian correction belongs mathematically with the likelihood
#     (it's a property of our parameterisation, not of the data), and
#     adding it here keeps the prior object clean.
#
#     HOW THE LIKELIHOOD WORKS:
#     For each plot, we:
#       1. Compute climate modifiers (xi) from temperature and precipitation
#       2. Initialise the model at steady state (using mean litter inputs)
#       3. Run the model forward through time
#       4. Compare predicted SOC to observed SOC using a Gaussian error model
#          where the standard deviation scales with predicted SOC (proportional
#          error -- bigger stocks have bigger absolute uncertainty)
#
#     The first observation at each plot uses a wider error (sigma_obs^2 +
#     sigma_init^2) because we are also uncertain about the initial conditions.
#     Later observations only carry sigma_obs because the model has "warmed up".
#
#     We return -Inf (impossible) if anything goes wrong: non-finite xi,
#     negative initial stocks, model crash, non-finite predictions, etc.
# =============================================================================

log_likelihood_yasso07 <- function(x) {
  
  # --- Back-transform from unconstrained to physical space ---
  p_free       <- to_original(x)
  sigma_obs    <- p_free["sigma_obs"]
  sigma_init   <- p_free["sigma_init"]
  model_params <- assemble_model_params(p_free)
  
  # --- Jacobian correction (see Section 5 for explanation) ---
  # This accounts for the uneven stretching caused by our transformations.
  # Adding it here ensures the sampler explores the right distribution.
  log_jac <- log_jac_stick(p_free[FRAC_COL_A], B) +
    log_jac_stick(p_free[FRAC_COL_W], B) +
    log_jac_stick(p_free[FRAC_COL_E], B) +
    log_jac_stick(p_free[FRAC_COL_N], B) +
    log(sigma_obs) +    # Jacobian of log transform: d(sigma)/d(log_sigma) = sigma
    log(sigma_init)
  
  # --- Loop over all plots and accumulate log-likelihood ---
  log_lik <- 0
  
  for (pid in plots) {
    
    clim   <- Yasso07_climate[Yasso07_climate$plot_id == pid, ]
    inputs <- Yasso07_inputs[Yasso07_inputs$plot_id   == pid, ]
    lm     <- litter_means[[pid]]
    meta   <- obs_meta[[pid]]
    
    # Guard: if any observation time step couldn't be matched to a model step
    if (any(is.na(meta$idx))) return(-Inf)
    
    # Step 1: compute annual climate modifiers (xi)
    # xi encodes how temperature and precipitation slow or speed decomposition
    xi_array <- tryCatch(
      compute_xi_yasso07(
        temp_mean = clim$temp_mean,
        temp_amp  = clim$temp_amplitude,
        precip    = clim$precip,
        beta1     = model_params["beta1"],
        beta2     = model_params["beta2"],
        gamma     = model_params["gamma"]
      ),
      error = function(e) NULL
    )
    if (is.null(xi_array))                         return(-Inf)
    xi_mean <- mean(xi_array)
    if (!is.finite(xi_mean) || xi_mean <= 0)       return(-Inf)
    
    # Step 2: initialise at steady state
    # We assume the plot was near equilibrium at the start of the simulation.
    # Steady state = litter inputs balance decomposition losses.
    C_init <- tryCatch(
      yasso07_steady_state(
        params   = model_params,
        nwl_mean = lm$nwl_mean,
        fwl_mean = lm$fwl_mean,
        cwl_mean = lm$cwl_mean,
        xi_mean  = xi_mean
      ),
      error = function(e) NULL
    )
    if (is.null(C_init) ||
        any(!is.finite(C_init)) ||
        any(C_init < 0))                           return(-Inf)
    
    # Step 3: run the model forward through time
    run_out <- tryCatch(
      yasso07_run(
        input_df = inputs,
        params   = model_params,
        C_init   = C_init,
        xi_array = xi_array
      ),
      error = function(e) NULL
    )
    if (is.null(run_out) ||
        any(!is.finite(run_out$total_soc)))        return(-Inf)
    
    # Step 4: extract predictions at observation time steps
    SOC_hat <- run_out$total_soc[meta$idx]
    if (any(!is.finite(SOC_hat)) || any(SOC_hat <= 0)) return(-Inf)
    
    # Step 5: compute observation standard deviation
    # We use a proportional error model: sd = SOC_hat * sigma
    # This means larger carbon stocks have proportionally larger uncertainty,
    # which is realistic (measurement and process errors scale with stock size).
    # The first observation also carries initialisation uncertainty (sigma_init).
    sd_vec <- ifelse(
      meta$is_first,
      SOC_hat * sqrt(sigma_obs^2 + sigma_init^2),  # first obs: obs + init noise
      SOC_hat * sigma_obs                           # later obs: obs noise only
    )
    
    # Step 6: add Gaussian log-likelihood for this plot
    # dnorm(observed, mean=predicted, sd=uncertainty) measures how likely
    # the observations are given these parameter values
    log_lik <- log_lik +
      sum(dnorm(meta$soc_obs, mean = SOC_hat, sd = sd_vec, log = TRUE))
  }
  
  if (!is.finite(log_lik)) return(-Inf)
  
  # Return log-likelihood plus Jacobian correction
  log_lik + log_jac
}


# =============================================================================
# 8.  Prior specification
#
#     We use a multivariate Gaussian (independent dimensions) in unconstrained
#     space, centred at the published defaults (best_x) with width sigma_ppm
#     determined by prior predictive matching above.
#
#     "density"  tells BayesianTools how probable a given point is under the
#                prior -- used to compute the full log-posterior.
#     "sampler"  tells BayesianTools how to draw starting points for chains --
#                this is only used for initialisation, not during sampling.
#     "lower/upper" are set to -Inf/+Inf because our unconstrained parameters
#                have no hard boundaries -- the constraints are handled by
#                the transformations in Section 5.
# =============================================================================

prior_yasso07 <- createPrior(
  
  density = function(x) {
    # Log-probability of x under a Gaussian prior centred at published defaults
    sum(dnorm(x, mean = best_x, sd = sigma_ppm, log = TRUE))
  },
  
  sampler = function(n = 1) {
    # Draw n starting points by perturbing the defaults slightly
    # (used only to initialise the MCMC chains)
    matrix(
      rnorm(n * N_FREE, mean = best_x, sd = sigma_ppm * 0.1),
      nrow = n, ncol = N_FREE
    )
  },
  
  lower = rep(-Inf, N_FREE),
  upper = rep( Inf, N_FREE)
)


# =============================================================================
# 9.  BayesianTools setup
#
#     createBayesianSetup packages the likelihood and prior together into
#     one object. BayesianTools then handles:
#       - Computing log-posterior = log-prior + log-likelihood at each step
#       - Running multiple chains in parallel
#       - Tracking acceptance rates
#       - Storing the chain samples
# =============================================================================

bayes_setup <- createBayesianSetup(
  likelihood = log_likelihood_yasso07,
  prior      = prior_yasso07,
  names      = FREE_NAMES,          # parameter names for diagnostics and plots
  parallel   = FALSE                # set to TRUE if you have parallel backend
)

# Quick sanity check: log-posterior should be finite at the default parameters
lp_check <- bayes_setup$posterior$density(best_x)
message(sprintf("\nLog-posterior at defaults: %.3f", lp_check))
if (!is.finite(lp_check)) stop("Log-posterior is not finite at defaults -- check model setup.")


# =============================================================================
# 10.  Run MCMC with DEzs sampler
#
#      WHAT IS DEzs?
#      Differential Evolution with snooker updates (ter Braak & Vrugt 2008).
#      It runs several chains simultaneously (default: 3). At each step, each
#      chain proposes a move based on the *difference* between two randomly
#      chosen other chains. This means:
#        - Proposals automatically align with the posterior's shape
#        - Correlated parameters are handled well
#        - No manual tuning of proposal width is needed
#        - The chains can escape local optima through crossover moves
#
#      KEY SETTINGS:
#        iterations: total number of samples per chain
#        burnin:     initial samples to discard (chain is still "warming up")
#        nrChains:   number of independent chains (more = better convergence check)
#        eps:        small random perturbation added to DE moves (avoids collapse)
# =============================================================================

mcmc_settings <- list(
  iterations = 50000,    # samples per chain -- increase for production
  burnin     = 10000,    # discard first 10000 as the chain finds the posterior
  nrChains   = 3,        # 3 chains lets us check convergence with R-hat
  eps        = 0.001     # tiny jitter on DE proposals (BayesianTools default)
)

message("\n--- Starting DEzs MCMC ---")
message(sprintf("  %d chains x %d iterations = %d total evaluations",
                mcmc_settings$nrChains,
                mcmc_settings$iterations,
                mcmc_settings$nrChains * mcmc_settings$iterations))
message(sprintf("  Burning in first %d samples per chain\n", mcmc_settings$burnin))

result_bt <- runMCMC(
  bayesianSetup = bayes_setup,
  sampler       = "DEzs",
  settings      = mcmc_settings
)


# =============================================================================
# 11.  Convergence diagnostics
#
#      Before trusting the posterior, check that the chains have converged.
#
#      R-hat (Gelman-Rubin statistic):
#        Compares variance within chains to variance between chains.
#        R-hat close to 1.0 means all chains reached the same region.
#        R-hat > 1.1 is a warning sign -- run longer or check the model.
#
#      Effective sample size (ESS):
#        MCMC samples are correlated (each step depends on the last), so
#        10000 samples might carry the same information as only 500 independent
#        ones. ESS estimates how many "effective" independent samples you have.
#        Aim for ESS > 200 for reliable posterior summaries.
# =============================================================================

message("\n--- Convergence diagnostics ---")

# BayesianTools built-in summary (includes R-hat and ESS)
summary(result_bt)

# Gelman-Rubin R-hat (requires multiple chains)
if (mcmc_settings$nrChains > 1) {
  tryCatch({
    gr <- gelmanDiagnostics(result_bt)
    message("\nGelman-Rubin R-hat (should be < 1.1 for all parameters):")
    print(round(gr$psrf, 3))
  }, error = function(e) {
    message("Could not compute Gelman-Rubin: ", conditionMessage(e))
  })
}


# =============================================================================
# 12.  Extract and summarise posterior
#
#      The chain is stored in unconstrained space. We back-transform each
#      sample to physical space before computing summaries, so the results
#      are in interpretable units (fractions, 1/year, etc.).
# =============================================================================

message("\n--- Posterior summary (physical space, post-burnin) ---")

# Extract the MCMC samples as a matrix (rows = samples, cols = parameters)
# getSample with coda=FALSE returns a plain matrix, burnin already removed
chain_raw <- getSample(result_bt, coda = FALSE)

# Back-transform every row to physical space
chain_physical <- t(apply(chain_raw, 1, to_original))
colnames(chain_physical) <- FREE_NAMES

# Posterior quantiles for each parameter
print(round(
  apply(chain_physical, 2, quantile,
        probs = c(0.025, 0.25, 0.50, 0.75, 0.975)),
  4
))


# =============================================================================
# 13.  Diagnostic plots
#
#      Trace plots: parameter value vs. iteration -- should look like
#        "fuzzy caterpillars", not trends or stuck regions.
#      Marginal densities: posterior distribution for each parameter.
#      Correlation plot: shows which parameters trade off against each other.
# =============================================================================

message("\n--- Generating diagnostic plots ---")

# Trace plots and marginal densities (BayesianTools built-in)
plot(result_bt)

# Correlation structure among posterior samples
# Climate parameters only -- the ones with borderline R-hat
chain_raw <- getSample(result_bt, coda = FALSE)

# Back-transform to physical space
chain_physical <- t(apply(chain_raw, 1, to_original))
colnames(chain_physical) <- FREE_NAMES

# Just the climate + modifier block
climate_cols <- c("beta1", "beta2", "gamma", "delta1", "delta2", "r")
pairs(chain_physical[, climate_cols],
      pch = ".", col = rgb(0, 0, 0, 0.1),
      main = "Climate and modifier parameters")

# One panel per pool -- much more readable than all 12 at once
par(mfrow = c(2, 2))
pairs(chain_physical[, c("p_AW","p_AE","p_AN")], pch = ".", main = "Pool A outflows")
pairs(chain_physical[, c("p_WA","p_WE","p_WN")], pch = ".", main = "Pool W outflows")
pairs(chain_physical[, c("p_EA","p_EW","p_EN")], pch = ".", main = "Pool E outflows")
pairs(chain_physical[, c("p_NA","p_NW","p_NE")], pch = ".", main = "Pool N outflows")
par(mfrow = c(1, 1))

# # Marginal plots for key parameters of interest
# marginalPlot(result_bt,
#              prior = TRUE,    # overlay prior so you can see how much data helped
#              whichParameters = c(IDX_CLIMATE, IDX_SIGMA_OBS, IDX_SIGMA_INIT))

gelmanDiagnostics(result_bt)


# Benchmark a single evaluation
t1 <- system.time(log_likelihood_yasso07(best_x))
message(sprintf("Single evaluation: %.3f seconds", t1["elapsed"]))

# How many plots in your test dataset?
message(sprintf("Number of plots in test data: %d", length(plots)))

# Projected time per evaluation at 3700 plots
# (assuming linear scaling, which is approximately true)
t_per_plot    <- t1["elapsed"] / length(plots)
t_at_3700     <- t_per_plot * 3700
target_evals  <- 150000

message(sprintf("Time per plot:              %.6f seconds", t_per_plot))
message(sprintf("Projected time at 3700 plots: %.3f seconds per eval", t_at_3700))
message(sprintf("Projected total runtime:      %.1f hours",
                t_at_3700 * target_evals / 3600))
# =============================================================================
# 14.  Save results
# =============================================================================

message("\n--- Saving results ---")

# Save the full BayesianTools result object (can be reloaded for further analysis)
saveRDS(result_bt,      file = "Yasso07_BT_result.rds")

# Save the back-transformed chain for downstream analysis
saveRDS(chain_physical, file = "Yasso07_posterior_samples.rds")

message("Saved:")
message("  Yasso07_BT_result.rds        -- full BayesianTools result object")
message("  Yasso07_posterior_samples.rds -- posterior samples in physical space")
message("\nDone.")





# =============================================================================
# 15. test parallelization
# =============================================================================


library(parallel)

# Detect available cores -- leave one free for the OS
n_cores <- parallel::detectCores() - 1

log_likelihood_yasso07_parallel <- function(x) {
  
  p_free       <- to_original(x)
  sigma_obs    <- p_free["sigma_obs"]
  sigma_init   <- p_free["sigma_init"]
  model_params <- assemble_model_params(p_free)
  
  log_jac <- log_jac_stick(p_free[FRAC_COL_A], B) +
    log_jac_stick(p_free[FRAC_COL_W], B) +
    log_jac_stick(p_free[FRAC_COL_E], B) +
    log_jac_stick(p_free[FRAC_COL_N], B) +
    log(sigma_obs) +
    log(sigma_init)
  
  # mclapply: like lapply but forks n_cores processes
  # Each worker gets one plot and returns its log-likelihood contribution
  log_liks <- parallel::mclapply(plots, function(pid) {
    
    clim   <- Yasso07_climate[Yasso07_climate$plot_id == pid, ]
    inputs <- Yasso07_inputs[Yasso07_inputs$plot_id   == pid, ]
    lm     <- litter_means[[pid]]
    meta   <- obs_meta[[pid]]
    
    if (any(is.na(meta$idx))) return(-Inf)
    
    xi_array <- tryCatch(
      compute_xi_yasso07(
        temp_mean = clim$temp_mean,
        temp_amp  = clim$temp_amplitude,
        precip    = clim$precip,
        beta1     = model_params["beta1"],
        beta2     = model_params["beta2"],
        gamma     = model_params["gamma"]
      ),
      error = function(e) NULL
    )
    if (is.null(xi_array))                   return(-Inf)
    xi_mean <- mean(xi_array)
    if (!is.finite(xi_mean) || xi_mean <= 0) return(-Inf)
    
    C_init <- tryCatch(
      yasso07_steady_state(
        params   = model_params,
        nwl_mean = lm$nwl_mean,
        fwl_mean = lm$fwl_mean,
        cwl_mean = lm$cwl_mean,
        xi_mean  = xi_mean
      ),
      error = function(e) NULL
    )
    if (is.null(C_init) ||
        any(!is.finite(C_init)) ||
        any(C_init < 0))                     return(-Inf)
    
    run_out <- tryCatch(
      yasso07_run(
        input_df = inputs,
        params   = model_params,
        C_init   = C_init,
        xi_array = xi_array
      ),
      error = function(e) NULL
    )
    if (is.null(run_out) ||
        any(!is.finite(run_out$total_soc)))  return(-Inf)
    
    SOC_hat <- run_out$total_soc[meta$idx]
    if (any(!is.finite(SOC_hat)) ||
        any(SOC_hat <= 0))                   return(-Inf)
    
    sd_vec <- ifelse(
      meta$is_first,
      SOC_hat * sqrt(sigma_obs^2 + sigma_init^2),
      SOC_hat * sigma_obs
    )
    
    sum(dnorm(meta$soc_obs, mean = SOC_hat, sd = sd_vec, log = TRUE))
    
  }, mc.cores = n_cores)
  
  # If any plot returned -Inf, the whole likelihood is -Inf
  log_liks <- unlist(log_liks)
  if (any(!is.finite(log_liks))) return(-Inf)
  
  sum(log_liks) + log_jac
}



# Compare sequential vs parallel on your test dataset
system.time(log_likelihood_yasso07(best_x))          # sequential
system.time(log_likelihood_yasso07_parallel(best_x)) # parallel



# Simulate increasing plot counts by replicating your existing plots
# (not real data, just to benchmark scaling behaviour)
plots_original <- plots

for (n_plots_target in c(10, 50, 100, 500, 1000, 3700)) {
  
  # Replicate plots to reach target size
  plots <<- rep(plots_original, length.out = n_plots_target)
  
  t_seq <- system.time(log_likelihood_yasso07(best_x))["elapsed"]
  t_par <- system.time(log_likelihood_yasso07_parallel(best_x))["elapsed"]
  
  message(sprintf("  %4d plots | seq: %.3fs | par: %.3fs | speedup: %.2fx",
                  n_plots_target, t_seq, t_par, t_seq / t_par))
}

# Restore original
plots <<- plots_original



#profiling
library(profvis)
profvis({
  for (i in 1:10) log_likelihood_yasso07(best_x)
})



# Trying to refactor the likelihood function to be faster
climate_by_plot <- split(Yasso07_climate, Yasso07_climate$plot_id)
inputs_by_plot  <- split(Yasso07_inputs,  Yasso07_inputs$plot_id)


log_likelihood_yasso07_refactored <- function(x) {
  
  # --- Back-transform from unconstrained to physical space ---
  p_free       <- to_original(x)
  sigma_obs    <- p_free["sigma_obs"]
  sigma_init   <- p_free["sigma_init"]
  model_params <- assemble_model_params(p_free)
  
  # --- Jacobian correction (see Section 5 for explanation) ---
  # This accounts for the uneven stretching caused by our transformations.
  # Adding it here ensures the sampler explores the right distribution.
  log_jac <- log_jac_stick(p_free[FRAC_COL_A], B) +
    log_jac_stick(p_free[FRAC_COL_W], B) +
    log_jac_stick(p_free[FRAC_COL_E], B) +
    log_jac_stick(p_free[FRAC_COL_N], B) +
    log(sigma_obs) +    # Jacobian of log transform: d(sigma)/d(log_sigma) = sigma
    log(sigma_init)
  
  # --- Loop over all plots and accumulate log-likelihood ---
  log_lik <- 0
  
  for (pid in plots) {
    
    clim   <- climate_by_plot[[pid]] #!!!!!!!!!!! Here the changes for speed, avoiding some overheads
    inputs <- inputs_by_plot[[pid]]#!!!!!!!!!!! Here the changes for speed, avoiding some overheads
    lm     <- litter_means[[pid]]
    meta   <- obs_meta[[pid]]
    
    # Guard: if any observation time step couldn't be matched to a model step
    if (any(is.na(meta$idx))) return(-Inf)
    
    # Step 1: compute annual climate modifiers (xi)
    # xi encodes how temperature and precipitation slow or speed decomposition
    xi_array <- tryCatch(
      compute_xi_yasso07(
        temp_mean = clim$temp_mean,
        temp_amp  = clim$temp_amplitude,
        precip    = clim$precip,
        beta1     = model_params["beta1"],
        beta2     = model_params["beta2"],
        gamma     = model_params["gamma"]
      ),
      error = function(e) NULL
    )
    if (is.null(xi_array))                         return(-Inf)
    xi_mean <- mean(xi_array)
    if (!is.finite(xi_mean) || xi_mean <= 0)       return(-Inf)
    
    # Step 2: initialise at steady state
    # We assume the plot was near equilibrium at the start of the simulation.
    # Steady state = litter inputs balance decomposition losses.
    C_init <- tryCatch(
      yasso07_steady_state(
        params   = model_params,
        nwl_mean = lm$nwl_mean,
        fwl_mean = lm$fwl_mean,
        cwl_mean = lm$cwl_mean,
        xi_mean  = xi_mean
      ),
      error = function(e) NULL
    )
    if (is.null(C_init) ||
        any(!is.finite(C_init)) ||
        any(C_init < 0))                           return(-Inf)
    
    # Step 3: run the model forward through time
    run_out <- tryCatch(
      yasso07_run(
        input_df = inputs,
        params   = model_params,
        C_init   = C_init,
        xi_array = xi_array
      ),
      error = function(e) NULL
    )
    if (is.null(run_out) ||
        any(!is.finite(run_out$total_soc)))        return(-Inf)
    
    # Step 4: extract predictions at observation time steps
    SOC_hat <- run_out$total_soc[meta$idx]
    if (any(!is.finite(SOC_hat)) || any(SOC_hat <= 0)) return(-Inf)
    
    # Step 5: compute observation standard deviation
    # We use a proportional error model: sd = SOC_hat * sigma
    # This means larger carbon stocks have proportionally larger uncertainty,
    # which is realistic (measurement and process errors scale with stock size).
    # The first observation also carries initialisation uncertainty (sigma_init).
    sd_vec <- ifelse(
      meta$is_first,
      SOC_hat * sqrt(sigma_obs^2 + sigma_init^2),  # first obs: obs + init noise
      SOC_hat * sigma_obs                           # later obs: obs noise only
    )
    
    # Step 6: add Gaussian log-likelihood for this plot
    # dnorm(observed, mean=predicted, sd=uncertainty) measures how likely
    # the observations are given these parameter values
    log_lik <- log_lik +
      sum(dnorm(meta$soc_obs, mean = SOC_hat, sd = sd_vec, log = TRUE))
  }
  
  if (!is.finite(log_lik)) return(-Inf)
  
  # Return log-likelihood plus Jacobian correction
  log_lik + log_jac
}


#profiling
library(profvis)
profvis({
  for (i in 1:10) log_likelihood_yasso07_refactored(best_x)
})



library(compiler)

# Force immediate byte-compilation -- pays the cost once, not on every call
log_likelihood_yasso07_refactored <- cmpfun(log_likelihood_yasso07_refactored)


system.time(for (i in 1:10) log_likelihood_yasso07_refactored(best_x))









##### refactored parallel likelihood

# This works, refactoring also the parallel version
climate_by_plot <- split(Yasso07_climate, Yasso07_climate$plot_id)
inputs_by_plot  <- split(Yasso07_inputs,  Yasso07_inputs$plot_id)

library(parallel)

# Detect available cores -- leave one free for the OS
n_cores <- parallel::detectCores() - 1

log_likelihood_yasso07_parallel_refactored <- function(x) {
  
  p_free       <- to_original(x)
  sigma_obs    <- p_free["sigma_obs"]
  sigma_init   <- p_free["sigma_init"]
  model_params <- assemble_model_params(p_free)
  
  log_jac <- log_jac_stick(p_free[FRAC_COL_A], B) +
    log_jac_stick(p_free[FRAC_COL_W], B) +
    log_jac_stick(p_free[FRAC_COL_E], B) +
    log_jac_stick(p_free[FRAC_COL_N], B) +
    log(sigma_obs) +
    log(sigma_init)
  
  # mclapply: like lapply but forks n_cores processes
  # Each worker gets one plot and returns its log-likelihood contribution
  log_liks <- parallel::mclapply(plots, function(pid) {
    
    clim   <- climate_by_plot[[pid]] # !!!!!Fix for speed, eliminating some overheads
    inputs <- inputs_by_plot[[pid]]# !!!!!Fix for speed, eliminating some overheads
    lm     <- litter_means[[pid]]
    meta   <- obs_meta[[pid]]
    
    if (any(is.na(meta$idx))) return(-Inf)
    
    xi_array <- tryCatch(
      compute_xi_yasso07(
        temp_mean = clim$temp_mean,
        temp_amp  = clim$temp_amplitude,
        precip    = clim$precip,
        beta1     = model_params["beta1"],
        beta2     = model_params["beta2"],
        gamma     = model_params["gamma"]
      ),
      error = function(e) NULL
    )
    if (is.null(xi_array))                   return(-Inf)
    xi_mean <- mean(xi_array)
    if (!is.finite(xi_mean) || xi_mean <= 0) return(-Inf)
    
    C_init <- tryCatch(
      yasso07_steady_state(
        params   = model_params,
        nwl_mean = lm$nwl_mean,
        fwl_mean = lm$fwl_mean,
        cwl_mean = lm$cwl_mean,
        xi_mean  = xi_mean
      ),
      error = function(e) NULL
    )
    if (is.null(C_init) ||
        any(!is.finite(C_init)) ||
        any(C_init < 0))                     return(-Inf)
    
    run_out <- tryCatch(
      yasso07_run(
        input_df = inputs,
        params   = model_params,
        C_init   = C_init,
        xi_array = xi_array
      ),
      error = function(e) NULL
    )
    if (is.null(run_out) ||
        any(!is.finite(run_out$total_soc)))  return(-Inf)
    
    SOC_hat <- run_out$total_soc[meta$idx]
    if (any(!is.finite(SOC_hat)) ||
        any(SOC_hat <= 0))                   return(-Inf)
    
    sd_vec <- ifelse(
      meta$is_first,
      SOC_hat * sqrt(sigma_obs^2 + sigma_init^2),
      SOC_hat * sigma_obs
    )
    
    sum(dnorm(meta$soc_obs, mean = SOC_hat, sd = sd_vec, log = TRUE))
    
  }, mc.cores = n_cores)
  
  # If any plot returned -Inf, the whole likelihood is -Inf
  log_liks <- unlist(log_liks)
  if (any(!is.finite(log_liks))) return(-Inf)
  
  sum(log_liks) + log_jac
}


#remove some overheads
library(compiler)
log_likelihood_yasso07_parallel_refactored <- cmpfun(log_likelihood_yasso07_parallel_refactored)



# Compare sequential vs parallel on your test dataset
system.time(log_likelihood_yasso07_refactored(best_x))          # sequential
system.time(log_likelihood_yasso07_parallel_refactored(best_x)) # parallel



# Simulate increasing plot counts by replicating your existing plots
# (not real data, just to benchmark scaling behaviour)
plots_original <- plots

for (n_plots_target in c(10, 50, 100, 500, 1000, 3700)) {
  
  # Replicate plots to reach target size
  plots <<- rep(plots_original, length.out = n_plots_target)
  
  t_seq <- system.time(log_likelihood_yasso07_refactored(best_x))["elapsed"]
  t_par <- system.time(log_likelihood_yasso07_parallel_refactored(best_x))["elapsed"]
  
  message(sprintf("  %4d plots | seq: %.3fs | par: %.3fs | speedup: %.2fx",
                  n_plots_target, t_seq, t_par, t_seq / t_par))
}

# Restore original
plots <<- plots_original


t_par_mac      <- 0.359   # seconds per eval, 11 cores on Mac
n_cores_mac    <- 11
target_evals   <- 150000  # 50k iterations x 3 chains


# The sequential floor (overhead that can't be parallelised)
# At 3700 plots: seq = 1.763s, par = 0.359s
# Amdahl: t_par = t_seq * ((1-f) + f/n)
# Solve for f:
n_cores_mac = 11
t_seq_3700 <- 1.763
t_par_11   <- 0.359
f <- (1 - t_par_11/t_seq_3700) / (1 - 1/n_cores_mac)
message(sprintf("Parallel fraction: %.3f (%.1f%% of work is parallelisable)", f, f*100))

# Sequential floor (irreducible minimum per evaluation)
t_floor <- t_seq_3700 * (1 - f)
message(sprintf("Sequential floor:  %.3f seconds", t_floor))

# Project to Roihu core counts
message("\nProjected total runtime at 150k evaluations:")
message(sprintf("  %-10s  %-12s  %-12s  %-10s", 
                "Cores", "Per eval", "Total", "vs Mac"))
for (n_cores in c(11, 16, 32, 64, 128)) {
  t_per_eval  <- t_seq_3700 * ((1 - f) + f / n_cores)
  total_hours <- t_per_eval * target_evals / 3600
  speedup     <- t_par_11 / t_per_eval
  message(sprintf("  %-10d  %-12s  %-12s  %-10s",
                  n_cores,
                  sprintf("%.3fs", t_per_eval),
                  sprintf("%.1f hrs", total_hours),
                  sprintf("%.1fx", speedup)))
}


