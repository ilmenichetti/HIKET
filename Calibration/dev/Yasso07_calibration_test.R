# =============================================================================
# Calibration_Yasso07.R
#
# Design principles:
#
#   1. FIXED k: Decomposition rates (alpha_A/W/E/N, p_H, alpha_H) fixed at
#      published values. These encode pool chemical identity and are
#      independently constrained by experimental data. Fixing them preserves
#      model identity and eliminates the MRT ordering problem.
#
#   2. FREE parameters: Transfer fractions (12), climate sensitivity
#      (beta1, beta2, gamma), woody modifiers (delta1, delta2, r).
#
#   3. TRANSFORMED sampling: Sampler works entirely in unconstrained space.
#      Back-transformation in R before any Fortran call -- Fortran unchanged.
#
#      Transfer fractions: stick-breaking per column, enforcing the physical
#      constraint that outflows from each pool cannot exceed 1 - p_H.
#      Each column uses 3 independent logit-transformed values v1, v2, v3:
#
#        B = 1 - p_H                      (budget after humus formation)
#        p_1 = logistic(v1) * B
#        p_2 = logistic(v2) * (B - p_1)
#        p_3 = logistic(v3) * (B - p_1 - p_2)
#        respiration = B - p_1 - p_2 - p_3  (always >= 0)
#
#      Climate/modifier params (unbounded): identity
#      Nuisance sigma > 0: log transform
#
#   4. PRIOR WIDTH via prior predictive matching: sigma chosen so that the
#      CV of prior-predicted SOC matches observed CV across plots.
#
# Sampler: standalone Metropolis-Hastings (no BayesianTools).
#
# Free parameter vector layout (unconstrained space, length 20):
#   [1:3]   v_AW, v_AE, v_AN   stick-breaking for column A
#   [4:6]   v_WA, v_WE, v_WN   stick-breaking for column W
#   [7:9]   v_EA, v_EW, v_EN   stick-breaking for column E
#   [10:12] v_NA, v_NW, v_NE   stick-breaking for column N
#   [13]    beta1
#   [14]    beta2
#   [15]    gamma
#   [16]    delta1
#   [17]    delta2
#   [18]    r
#   [19]    log(sigma_obs)
#   [20]    log(sigma_init)
# =============================================================================

library(dplyr)

source("./Model_functions/input_compatibility_layer.R")
source("./Model_functions/Decomposition_functions/Yasso/yasso07_wrapper.R")
dyn.load("./Model_functions/Decomposition_functions/Yasso/yasso07.so")


# =============================================================================
# 1.  Load and map inputs
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
# =============================================================================

litter_means <- lapply(plots, function(pid) {
  inputs <- Yasso07_inputs[Yasso07_inputs$plot_id == pid, ]
  list(
    nwl_mean = colMeans(inputs[, c("nwl_A","nwl_W","nwl_E","nwl_N")]),
    fwl_mean = colMeans(inputs[, c("fwl_A","fwl_W","fwl_E","fwl_N")]),
    cwl_mean = colMeans(inputs[, c("cwl_A","cwl_W","cwl_E","cwl_N")])
  )
})
names(litter_means) <- plots

obs_meta <- lapply(plots, function(pid) {
  clim     <- Yasso07_climate[Yasso07_climate$plot_id == pid, ]
  obs_plot <- SOC_obs_all[SOC_obs_all$plot_id == pid, ]
  list(
    idx      = match(obs_plot$year, clim$year),
    soc_obs  = obs_plot$soc_obs_tCha,
    is_first = obs_plot$obs_rank == 1L
  )
})
names(obs_meta) <- plots


# =============================================================================
# 4.  Parameter structure
# =============================================================================

p_default <- YASSO07_DEFAULT_PARAMS

# -- Fixed: decomposition rates (pool chemical identity) --
FIXED_RATE_NAMES <- c("alpha_A","alpha_W","alpha_E","alpha_N","p_H","alpha_H")
fixed_rates      <- p_default[FIXED_RATE_NAMES]

# -- Free: partitioning and climate response --
# Transfer fractions -- named by their destination_source convention
# Column A (outflows from A): p_AW, p_AE, p_AN
# Column W (outflows from W): p_WA, p_WE, p_WN
# Column E (outflows from E): p_EA, p_EW, p_EN
# Column N (outflows from N): p_NA, p_NW, p_NE
FRAC_COL_A     <- c("p_AW","p_AE","p_AN")
FRAC_COL_W     <- c("p_WA","p_WE","p_WN")
FRAC_COL_E     <- c("p_EA","p_EW","p_EN")
FRAC_COL_N     <- c("p_NA","p_NW","p_NE")
FRAC_NAMES     <- c(FRAC_COL_A, FRAC_COL_W, FRAC_COL_E, FRAC_COL_N)

CLIMATE_NAMES  <- c("beta1","beta2","gamma")
MODIFIER_NAMES <- c("delta1","delta2","r")
NUISANCE_NAMES <- c("sigma_obs","sigma_init")

FREE_NAMES <- c(FRAC_NAMES, CLIMATE_NAMES, MODIFIER_NAMES, NUISANCE_NAMES)
N_FREE     <- length(FREE_NAMES)   # 20

# Default values in ORIGINAL space
free_defaults <- c(
  p_default[FRAC_NAMES],
  p_default[CLIMATE_NAMES],
  p_default[MODIFIER_NAMES],
  sigma_obs  = 0.05,
  sigma_init = 0.10
)

# Index helpers in unconstrained vector
IDX_FRACS      <- 1:12
IDX_COL_A      <- 1:3
IDX_COL_W      <- 4:6
IDX_COL_E      <- 7:9
IDX_COL_N      <- 10:12
IDX_CLIMATE    <- 13:15
IDX_MODIFIERS  <- 16:18
IDX_SIGMA_OBS  <- 19L
IDX_SIGMA_INIT <- 20L

# Budget per column: fraction not going to humus
B <- 1.0 - fixed_rates["p_H"]


# =============================================================================
# 5.  Transformation functions
#
#     Stick-breaking for transfer fractions:
#       Given unconstrained v1, v2, v3 for a column with budget B:
#         p1 = logistic(v1) * B
#         p2 = logistic(v2) * (B - p1)
#         p3 = logistic(v3) * (B - p1 - p2)
#       Column sum = p1 + p2 + p3 <= B < 1  (always valid)
#       Remainder  = B - p1 - p2 - p3 is respired
#
#     Inverse (original -> unconstrained):
#       Given p1, p2, p3 and budget B:
#         v1 = logit(p1 / B)
#         v2 = logit(p2 / (B - p1))
#         v3 = logit(p3 / (B - p1 - p2))
# =============================================================================

logit     <- function(p) log(p / (1 - p))
inv_logit <- function(x) 1 / (1 + exp(-x))

# Stick-breaking: unconstrained v (length 3) -> fractions (length 3)
stick_break <- function(v, budget) {
  p <- numeric(3)
  p[1] <- inv_logit(v[1]) * budget
  p[2] <- inv_logit(v[2]) * (budget - p[1])
  p[3] <- inv_logit(v[3]) * (budget - p[1] - p[2])
  p
}

# Inverse stick-breaking: fractions (length 3) -> unconstrained v (length 3)
stick_break_inv <- function(p, budget) {
  v    <- numeric(3)
  v[1] <- logit(p[1] / budget)
  v[2] <- logit(p[2] / (budget - p[1]))
  v[3] <- logit(p[3] / (budget - p[1] - p[2]))
  v
}

# Log-Jacobian of stick-breaking transformation for one column
# d(p1,p2,p3)/d(v1,v2,v3) -- product of logistic derivatives * budget terms
log_jac_stick <- function(p, budget) {
  # logistic derivative: logistic(v) * (1 - logistic(v))
  # = (p/budget) * (1 - p/budget) * budget  for p1
  # etc. -- simplified:
  rem <- c(budget, budget - p[1], budget - p[1] - p[2])
  frac <- p / rem
  sum(log(frac * (1 - frac) * rem))
}

to_unconstrained <- function(free_params) {
  x <- numeric(N_FREE)
  names(x) <- FREE_NAMES
  
  x[IDX_COL_A]  <- stick_break_inv(free_params[FRAC_COL_A], B)
  x[IDX_COL_W]  <- stick_break_inv(free_params[FRAC_COL_W], B)
  x[IDX_COL_E]  <- stick_break_inv(free_params[FRAC_COL_E], B)
  x[IDX_COL_N]  <- stick_break_inv(free_params[FRAC_COL_N], B)
  
  x[IDX_CLIMATE]    <- free_params[CLIMATE_NAMES]
  x[IDX_MODIFIERS]  <- free_params[MODIFIER_NAMES]
  x[IDX_SIGMA_OBS]  <- log(free_params["sigma_obs"])
  x[IDX_SIGMA_INIT] <- log(free_params["sigma_init"])
  x
}

to_original <- function(x) {
  p <- numeric(N_FREE)
  names(p) <- FREE_NAMES
  
  p[FRAC_COL_A] <- stick_break(x[IDX_COL_A], B)
  p[FRAC_COL_W] <- stick_break(x[IDX_COL_W], B)
  p[FRAC_COL_E] <- stick_break(x[IDX_COL_E], B)
  p[FRAC_COL_N] <- stick_break(x[IDX_COL_N], B)
  
  p[CLIMATE_NAMES]  <- x[IDX_CLIMATE]
  p[MODIFIER_NAMES] <- x[IDX_MODIFIERS]
  p["sigma_obs"]    <- exp(x[IDX_SIGMA_OBS])
  p["sigma_init"]   <- exp(x[IDX_SIGMA_INIT])
  p
}

# Assemble full model parameter vector (fixed rates + free, correct order)
assemble_model_params <- function(p_free) {
  full <- c(fixed_rates,
            p_free[c(FRAC_NAMES, CLIMATE_NAMES, MODIFIER_NAMES)])
  full[names(p_default)]
}

# Starting point in unconstrained space
best_x <- to_unconstrained(free_defaults)

# Verify round-trip
stopifnot(max(abs(to_original(best_x)[FREE_NAMES] -
                    free_defaults[FREE_NAMES])) < 1e-10)
message("Transformation round-trip check: OK")

# Verify column sums at defaults
p_check      <- to_original(best_x)
params_check <- assemble_model_params(p_check)
message("Column sum check at defaults (all must be < 1):")
message(sprintf("  From A: %.6f", params_check["p_AW"] + params_check["p_AE"] +
                  params_check["p_AN"] + fixed_rates["p_H"]))
message(sprintf("  From W: %.6f", params_check["p_WA"] + params_check["p_WE"] +
                  params_check["p_WN"] + fixed_rates["p_H"]))
message(sprintf("  From E: %.6f", params_check["p_EA"] + params_check["p_EW"] +
                  params_check["p_EN"] + fixed_rates["p_H"]))
message(sprintf("  From N: %.6f", params_check["p_NA"] + params_check["p_NW"] +
                  params_check["p_NE"] + fixed_rates["p_H"]))


# =============================================================================
# 6.  Prior predictive matching
#
#     Find sigma_ppm such that CV of prior-predicted SOC matches observed CV.
# =============================================================================

obs_cv <- sd(SOC_obs_all$soc_obs_tCha) / mean(SOC_obs_all$soc_obs_tCha)
message(sprintf("\nObserved SOC CV: %.3f (%.1f%%)", obs_cv, obs_cv * 100))

run_prior_sample_soc <- function(sigma_val) {
  x_sample     <- best_x + rnorm(N_FREE, 0, sigma_val)
  p_free       <- to_original(x_sample)
  model_params <- assemble_model_params(p_free)
  
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
message(sprintf("  Achieved CV  : %.3f",
                ppm_cv[which(sigma_grid == sigma_ppm)]))


# =============================================================================
# 7.  Log-posterior
#
#     log p(x | data) = log_prior + log_jacobian + log_likelihood
#
#     Prior:    Gaussian(best_x, sigma_ppm^2) on unconstrained params
#     Jacobian: stick-breaking (fractions) + log (nuisance)
#     Likelihood: Gaussian proportional error on SOC observations
# =============================================================================

make_log_posterior_yasso07 <- function(sigma_prior  = sigma_ppm,
                                       report_every = 500L) {
  
  .plots        <- plots
  .climate      <- Yasso07_climate
  .inputs       <- Yasso07_inputs
  .litter_means <- litter_means
  .obs_meta     <- obs_meta
  .center       <- best_x
  .call_count   <- 0L
  .last_lik     <- NA_real_
  .start_time   <- proc.time()[["elapsed"]]
  
  function(x) {
    
    .call_count <<- .call_count + 1L
    
    # Back-transform
    p_free       <- to_original(x)
    sigma_obs    <- p_free["sigma_obs"]
    sigma_init   <- p_free["sigma_init"]
    model_params <- assemble_model_params(p_free)
    
    # Log-prior: Gaussian on unconstrained params
    log_prior <- sum(dnorm(x, mean = .center, sd = sigma_prior, log = TRUE))
    
    # Log-Jacobian: stick-breaking (4 columns) + log transforms (2 nuisance)
    log_jac <- log_jac_stick(p_free[FRAC_COL_A], B) +
      log_jac_stick(p_free[FRAC_COL_W], B) +
      log_jac_stick(p_free[FRAC_COL_E], B) +
      log_jac_stick(p_free[FRAC_COL_N], B) +
      log(sigma_obs) +
      log(sigma_init)
    
    # Log-likelihood
    log_lik <- 0
    
    for (pid in .plots) {
      
      clim   <- .climate[.climate$plot_id == pid, ]
      inputs <- .inputs[.inputs$plot_id   == pid, ]
      lm     <- .litter_means[[pid]]
      meta   <- .obs_meta[[pid]]
      
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
      if (is.null(xi_array)) return(-Inf)
      
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
      if (is.null(C_init) || any(!is.finite(C_init)) || any(C_init < 0)) return(-Inf)
      
      run_out <- tryCatch(
        yasso07_run(
          input_df = inputs,
          params   = model_params,
          C_init   = C_init,
          xi_array = xi_array
        ),
        error = function(e) NULL
      )
      if (is.null(run_out) || any(!is.finite(run_out$total_soc))) return(-Inf)
      
      SOC_hat <- run_out$total_soc[meta$idx]
      if (any(!is.finite(SOC_hat)) || any(SOC_hat <= 0)) return(-Inf)
      
      sd_vec <- ifelse(
        meta$is_first,
        SOC_hat * sqrt(sigma_obs^2 + sigma_init^2),
        SOC_hat * sigma_obs
      )
      
      log_lik <- log_lik +
        sum(dnorm(meta$soc_obs, mean = SOC_hat, sd = sd_vec, log = TRUE))
    }
    
    if (!is.finite(log_lik)) return(-Inf)
    
    log_post <- log_prior + log_jac + log_lik
    
    .last_lik <<- log_lik
    if (.call_count %% report_every == 0L) {
      elapsed <- proc.time()[["elapsed"]] - .start_time
      rate    <- .call_count / elapsed
      message(sprintf(
        "  [Yasso07] call %6d | log-lik: %9.3f | log-prior: %7.2f | %.1f calls/sec | elapsed: %.0f s",
        .call_count, log_lik, log_prior, rate, elapsed
      ))
    }
    
    log_post
  }
}


# =============================================================================
# 8.  Standalone Metropolis-Hastings sampler
# =============================================================================

run_metropolis <- function(log_post, start, n_iter,
                           proposal_sd, report_every = 500L,
                           seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  n_params <- length(start)
  chain    <- matrix(NA_real_, nrow = n_iter, ncol = n_params)
  colnames(chain) <- names(start)
  
  current    <- start
  current_lp <- log_post(current)
  n_accept   <- 0L
  
  for (i in seq_len(n_iter)) {
    proposal    <- current + rnorm(n_params, 0, proposal_sd)
    proposal_lp <- tryCatch(log_post(proposal), error = function(e) -Inf)
    
    if (is.finite(proposal_lp) &&
        log(runif(1)) < (proposal_lp - current_lp)) {
      current    <- proposal
      current_lp <- proposal_lp
      n_accept   <- n_accept + 1L
    }
    
    chain[i, ] <- current
    
    if (i %% report_every == 0L) {
      message(sprintf(
        "  [MH] iter %6d | log-post: %9.3f | acceptance: %.3f",
        i, current_lp, n_accept / i
      ))
    }
  }
  
  list(chain = chain, acceptance = n_accept / n_iter, final_lp = current_lp)
}


# =============================================================================
# 9.  Run
# =============================================================================

log_post <- make_log_posterior_yasso07(
  sigma_prior  = sigma_ppm,
  report_every = 500L
)

message("\n--- Verification at default parameters ---")
lp_start <- log_post(best_x)
message(sprintf("log-posterior at defaults: %.3f", lp_start))

# Initial proposal width: 10% of prior sigma
proposal_sd <- rep(sigma_ppm * 0.10, N_FREE)
names(proposal_sd) <- FREE_NAMES

# Acceptance rate calibration
message("\n--- Acceptance calibration (500 iterations) ---")
cal <- run_metropolis(
  log_post     = log_post,
  start        = best_x,
  n_iter       = 500L,
  proposal_sd  = proposal_sd,
  report_every = 500L,
  seed         = 42L
)
message(sprintf("Acceptance rate: %.3f (target: 0.20-0.40)", cal$acceptance))

if      (cal$acceptance < 0.10) { proposal_sd <- proposal_sd * 0.3; message("Low  -- scaling down to 30%") } else
  if      (cal$acceptance < 0.20) { proposal_sd <- proposal_sd * 0.7; message("Low  -- scaling down to 70%") } else
    if      (cal$acceptance < 0.40) { message("Good -- no change") } else
    { proposal_sd <- proposal_sd * 1.5; message("High -- scaling up to 150%") }

# Production run
message("\n--- Production run (10000 iterations) ---")
result <- run_metropolis(
  log_post     = log_post,
  start        = best_x,
  n_iter       = 10000L,
  proposal_sd  = proposal_sd,
  report_every = 1000L,
  seed         = 123L
)

message(sprintf("\nFinal acceptance rate: %.3f", result$acceptance))

# Back-transform chain to original parameter space
chain_original <- t(apply(result$chain, 1, to_original))

message("\n--- Posterior summary (original space, burnin = 2000) ---")
burnin <- 2000L
print(round(apply(chain_original[burnin:nrow(chain_original), ], 2,
                  quantile, probs = c(0.025, 0.25, 0.50, 0.75, 0.975)), 4))

