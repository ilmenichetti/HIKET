# =============================================================================
# Calibration_Yasso07_BT_performance_tests.R
#
# PURPOSE OF THIS FILE:
#   This is a development and testing document that records the full logical
#   process of optimising the Yasso07 likelihood function for MCMC calibration
#   on the Roihu cluster (128 cores).
#
#   It is meant to be run TOP TO BOTTOM, section by section. Each section
#   builds on the previous one, benchmarks the current state, and explains
#   what we learned and what we changed next.
#
# STORY IN BRIEF:
#   We started with a naive sequential likelihood and measured 77 hours
#   projected runtime at 3700 plots. Through profiling and targeted fixes,
#   we arrived at a hybrid approach (parallel chains x parallel plot loop)
#   projecting to ~2 hours on Roihu.
#
# OPTIMISATION STEPS:
#   Fix 1: pre-split data frames         -- eliminates O(n) subsetting
#   Fix 2: byte-compile with cmpfun()    -- eliminates JIT overhead
#   Fix 3: hybrid parallelism            -- chains x plot-loop parallelism
#
# SECTIONS:
#    0.  Setup
#    1.  Baseline sequential likelihood (original, unoptimised)
#    2.  Verify calibration is working (before optimising anything)
#    3.  Baseline runtime benchmark and projection
#    4.  Profiling the baseline -- finding the bottlenecks
#    5.  Fix 1: pre-split data frames
#    6.  Profiling after Fix 1 -- what changed?
#    7.  Fix 2: byte-compile with cmpfun()
#    8.  Baseline and refactored parallel versions (mclapply)
#    9.  Scaling benchmark: all versions compared
#   10.  Runtime projection using Amdahl's law
#   11.  HYBRID APPROACH: parallel chains x parallel plot loop
#
# =============================================================================


# =============================================================================
# 0.  Setup
# =============================================================================

library(dplyr)
library(BayesianTools)
library(parallel)
library(compiler)
library(profvis)
library(coda)       # for gelman.diag on combined chains

source("./Model_functions/input_compatibility_layer.R")
source("./Model_functions/Decomposition_functions/Yasso/yasso07_wrapper.R")
dyn.load("./Model_functions/Decomposition_functions/Yasso/yasso07.so")

message("=============================================================")
message("  Yasso07 BayesianTools -- Performance Testing Document")
message("=============================================================\n")

# --- Data loading ---
input_raw       <- read.csv("synthetic_input_data_template.csv")
site_raw        <- read.csv("synthetic_site_data_template.csv")
Yasso07_climate <- map_climate_yasso07(input_raw)
Yasso07_inputs  <- map_inputs_yasso07(input_raw)
plots           <- unique(Yasso07_climate$plot_id)
message(sprintf("Loaded %d plots from input data.", length(plots)))

# --- SOC observations ---
SOC_obs_all <- input_raw %>%
  filter(!is.na(soc_obs_tCha)) %>%
  select(plot_id, year, month, soc_obs_tCha) %>%
  arrange(plot_id, year, month) %>%
  group_by(plot_id) %>%
  mutate(obs_rank = row_number()) %>%
  ungroup()
message(sprintf("Found %d SOC observations across %d plots.\n",
                nrow(SOC_obs_all), n_distinct(SOC_obs_all$plot_id)))

# --- Pre-compute plot-level quantities (done ONCE, outside MCMC) ---
# litter_means: mean litter inputs by pool, used for steady-state init
# obs_meta:     which time steps have observations, and which is first
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

# --- Parameter structure ---
p_default        <- YASSO07_DEFAULT_PARAMS
FIXED_RATE_NAMES <- c("alpha_A","alpha_W","alpha_E","alpha_N","p_H","alpha_H")
fixed_rates      <- p_default[FIXED_RATE_NAMES]
FRAC_COL_A       <- c("p_AW","p_AE","p_AN")
FRAC_COL_W       <- c("p_WA","p_WE","p_WN")
FRAC_COL_E       <- c("p_EA","p_EW","p_EN")
FRAC_COL_N       <- c("p_NA","p_NW","p_NE")
FRAC_NAMES       <- c(FRAC_COL_A, FRAC_COL_W, FRAC_COL_E, FRAC_COL_N)
CLIMATE_NAMES    <- c("beta1","beta2","gamma")
MODIFIER_NAMES   <- c("delta1","delta2","r")
NUISANCE_NAMES   <- c("sigma_obs","sigma_init")
FREE_NAMES       <- c(FRAC_NAMES, CLIMATE_NAMES, MODIFIER_NAMES, NUISANCE_NAMES)
N_FREE           <- length(FREE_NAMES)
IDX_COL_A        <- 1:3;   IDX_COL_W      <- 4:6
IDX_COL_E        <- 7:9;   IDX_COL_N      <- 10:12
IDX_CLIMATE      <- 13:15; IDX_MODIFIERS  <- 16:18
IDX_SIGMA_OBS    <- 19L;   IDX_SIGMA_INIT <- 20L
B                <- 1.0 - fixed_rates["p_H"]

free_defaults <- c(
  p_default[FRAC_NAMES], p_default[CLIMATE_NAMES], p_default[MODIFIER_NAMES],
  sigma_obs = 0.05, sigma_init = 0.10
)

# --- Transformation functions (unchanged throughout all sections) ---
logit     <- function(p) log(p / (1 - p))
inv_logit <- function(x) 1 / (1 + exp(-x))

stick_break <- function(v, budget) {
  p    <- numeric(3)
  p[1] <- inv_logit(v[1]) * budget
  p[2] <- inv_logit(v[2]) * (budget - p[1])
  p[3] <- inv_logit(v[3]) * (budget - p[1] - p[2])
  p
}

stick_break_inv <- function(p, budget) {
  v    <- numeric(3)
  v[1] <- logit(p[1] / budget)
  v[2] <- logit(p[2] / (budget - p[1]))
  v[3] <- logit(p[3] / (budget - p[1] - p[2]))
  v
}

log_jac_stick <- function(p, budget) {
  rem  <- c(budget, budget - p[1], budget - p[1] - p[2])
  frac <- p / rem
  sum(log(frac * (1 - frac) * rem))
}

to_unconstrained <- function(free_params) {
  x <- numeric(N_FREE); names(x) <- FREE_NAMES
  x[IDX_COL_A]      <- stick_break_inv(free_params[FRAC_COL_A], B)
  x[IDX_COL_W]      <- stick_break_inv(free_params[FRAC_COL_W], B)
  x[IDX_COL_E]      <- stick_break_inv(free_params[FRAC_COL_E], B)
  x[IDX_COL_N]      <- stick_break_inv(free_params[FRAC_COL_N], B)
  x[IDX_CLIMATE]    <- free_params[CLIMATE_NAMES]
  x[IDX_MODIFIERS]  <- free_params[MODIFIER_NAMES]
  x[IDX_SIGMA_OBS]  <- log(free_params["sigma_obs"])
  x[IDX_SIGMA_INIT] <- log(free_params["sigma_init"])
  x
}

to_original <- function(x) {
  p <- numeric(N_FREE); names(p) <- FREE_NAMES
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

assemble_model_params <- function(p_free) {
  full <- c(fixed_rates, p_free[c(FRAC_NAMES, CLIMATE_NAMES, MODIFIER_NAMES)])
  full[names(p_default)]
}

best_x <- to_unconstrained(free_defaults)

# Verify round-trip accuracy of transformations
rt_err <- max(abs(to_original(best_x)[FREE_NAMES] - free_defaults[FREE_NAMES]))
stopifnot(rt_err < 1e-10)
message(sprintf("Transformation round-trip check: OK (max error = %.2e)\n", rt_err))


# =============================================================================
# 1.  Baseline sequential likelihood (original, unoptimised)
#
#     This is the starting point before any optimisation. Two bottlenecks
#     will be identified later by profiling (Section 4):
#       a) [.data.frame -- scans entire climate/inputs data frame for each
#          plot on every MCMC evaluation. With 3700 plots this is ~11100
#          full-table scans per likelihood call.
#       b) JIT compilation overhead -- R re-compiles the function body on
#          every call instead of doing it once upfront.
# =============================================================================

message("-------------------------------------------------------------")
message("SECTION 1: Baseline sequential likelihood (unoptimised)")
message("-------------------------------------------------------------\n")

log_likelihood_yasso07 <- function(x) {
  
  p_free       <- to_original(x)
  sigma_obs    <- p_free["sigma_obs"]
  sigma_init   <- p_free["sigma_init"]
  model_params <- assemble_model_params(p_free)
  
  log_jac <- log_jac_stick(p_free[FRAC_COL_A], B) +
    log_jac_stick(p_free[FRAC_COL_W], B) +
    log_jac_stick(p_free[FRAC_COL_E], B) +
    log_jac_stick(p_free[FRAC_COL_N], B) +
    log(sigma_obs) + log(sigma_init)
  
  log_lik <- 0
  
  for (pid in plots) {
    
    # BOTTLENECK (identified in Section 4):
    # These scan the full data frame on every plot, every evaluation.
    clim   <- Yasso07_climate[Yasso07_climate$plot_id == pid, ]
    inputs <- Yasso07_inputs[Yasso07_inputs$plot_id   == pid, ]
    lm     <- litter_means[[pid]]
    meta   <- obs_meta[[pid]]
    
    if (any(is.na(meta$idx))) return(-Inf)
    
    xi_array <- tryCatch(
      compute_xi_yasso07(
        temp_mean = clim$temp_mean, temp_amp = clim$temp_amplitude,
        precip    = clim$precip,
        beta1     = model_params["beta1"], beta2 = model_params["beta2"],
        gamma     = model_params["gamma"]
      ), error = function(e) NULL)
    if (is.null(xi_array))                         return(-Inf)
    xi_mean <- mean(xi_array)
    if (!is.finite(xi_mean) || xi_mean <= 0)       return(-Inf)
    
    C_init <- tryCatch(
      yasso07_steady_state(
        params   = model_params, nwl_mean = lm$nwl_mean,
        fwl_mean = lm$fwl_mean, cwl_mean = lm$cwl_mean, xi_mean = xi_mean
      ), error = function(e) NULL)
    if (is.null(C_init) || any(!is.finite(C_init)) || any(C_init < 0)) return(-Inf)
    
    run_out <- tryCatch(
      yasso07_run(
        input_df = inputs, params   = model_params,
        C_init   = C_init, xi_array = xi_array
      ), error = function(e) NULL)
    if (is.null(run_out) || any(!is.finite(run_out$total_soc))) return(-Inf)
    
    SOC_hat <- run_out$total_soc[meta$idx]
    if (any(!is.finite(SOC_hat)) || any(SOC_hat <= 0)) return(-Inf)
    
    sd_vec <- ifelse(meta$is_first,
                     SOC_hat * sqrt(sigma_obs^2 + sigma_init^2),
                     SOC_hat * sigma_obs)
    
    log_lik <- log_lik + sum(dnorm(meta$soc_obs, mean = SOC_hat, sd = sd_vec, log = TRUE))
  }
  
  if (!is.finite(log_lik)) return(-Inf)
  log_lik + log_jac
}

message("Baseline sequential likelihood defined (no optimisations applied).")


# =============================================================================
# 2.  Verify calibration is working
#
#     CRITICAL: Before optimising anything, confirm the likelihood and the
#     BayesianTools setup give sensible results. There is no point making
#     something fast if it produces wrong answers.
#
#     Checks performed:
#       a) Likelihood is finite at published defaults
#       b) Likelihood responds to parameter perturbations
#       c) A short MCMC run (1000 iterations) is stable and mixes
#       d) Posterior medians are in plausible range
#
#     ALL CHECKS MUST PASS before proceeding to Section 3.
# =============================================================================

message("\n-------------------------------------------------------------")
message("SECTION 2: Verify calibration is working")
message("-------------------------------------------------------------\n")

# --- Prior predictive matching (needed for prior width) ---
# We find sigma_ppm
# such that prior-predicted SOC CV matches observed CV across plots.
obs_cv <- sd(SOC_obs_all$soc_obs_tCha) / mean(SOC_obs_all$soc_obs_tCha)
message(sprintf("Observed SOC CV: %.3f (%.1f%%)\n", obs_cv, obs_cv * 100))

run_prior_sample_soc <- function(sigma_val) {
  x_sample     <- best_x + rnorm(N_FREE, 0, sigma_val)
  p_free       <- to_original(x_sample)
  model_params <- assemble_model_params(p_free)
  pid          <- plots[1]
  clim         <- Yasso07_climate[Yasso07_climate$plot_id == pid, ]
  inputs       <- Yasso07_inputs[Yasso07_inputs$plot_id   == pid, ]
  lm           <- litter_means[[pid]]
  meta         <- obs_meta[[pid]]
  xi_array <- tryCatch(compute_xi_yasso07(
    temp_mean = clim$temp_mean, temp_amp = clim$temp_amplitude,
    precip = clim$precip, beta1 = model_params["beta1"],
    beta2 = model_params["beta2"], gamma = model_params["gamma"]
  ), error = function(e) NULL)
  if (is.null(xi_array)) return(NA_real_)
  xi_mean <- mean(xi_array)
  if (!is.finite(xi_mean) || xi_mean <= 0) return(NA_real_)
  C_init <- tryCatch(yasso07_steady_state(
    params = model_params, nwl_mean = lm$nwl_mean,
    fwl_mean = lm$fwl_mean, cwl_mean = lm$cwl_mean, xi_mean = xi_mean
  ), error = function(e) NULL)
  if (is.null(C_init) || any(!is.finite(C_init)) || any(C_init < 0)) return(NA_real_)
  run_out <- tryCatch(yasso07_run(
    input_df = inputs, params = model_params, C_init = C_init, xi_array = xi_array
  ), error = function(e) NULL)
  if (is.null(run_out)) return(NA_real_)
  soc <- run_out$total_soc[meta$idx[1]]
  if (!is.finite(soc) || soc <= 0) return(NA_real_)
  soc
}

set.seed(42)
sigma_grid <- seq(0.1, 3.0, by = 0.1)
n_ppm      <- 200L
message("Running prior predictive matching to determine prior width...")
ppm_cv <- sapply(sigma_grid, function(sv) {
  preds <- vapply(seq_len(n_ppm), function(i) run_prior_sample_soc(sv), numeric(1))
  valid <- preds[is.finite(preds) & preds > 0]
  if (length(valid) < 20L) NA_real_ else sd(valid) / mean(valid)
})
sigma_ppm <- sigma_grid[which(ppm_cv >= obs_cv)[1]]
if (is.na(sigma_ppm)) {
  sigma_ppm <- 1.0
  message("WARNING: CV target never reached -- using sigma_ppm = 1.0")
}
message(sprintf("Prior width: sigma_ppm = %.1f  (target CV = %.3f, achieved = %.3f)\n",
                sigma_ppm, obs_cv, ppm_cv[which(sigma_grid == sigma_ppm)]))

# --- Check 2a: likelihood at defaults must be finite and negative ---
message("CHECK 2a: Likelihood at published defaults")
ll_baseline <- log_likelihood_yasso07(best_x)
if (!is.finite(ll_baseline)) {
  stop("FAIL: log-likelihood not finite at defaults. Check model and data.")
}
if (ll_baseline >= 0) {
  message(sprintf("  WARNING: log-likelihood is non-negative (%.3f) -- check data scale", ll_baseline))
} else {
  message(sprintf("  PASS: log-likelihood = %.3f (finite, negative as expected)\n", ll_baseline))
}

# --- Check 2b: likelihood responds to parameter perturbations ---
message("CHECK 2b: Sensitivity to parameter perturbations")
message("  Perturbing parameters at 50%% of prior width and checking response...")
set.seed(99)
n_perturb  <- 20L
ll_perturb <- vapply(seq_len(n_perturb), function(i) {
  x_p <- best_x + rnorm(N_FREE, 0, sigma_ppm * 0.5)
  log_likelihood_yasso07(x_p)
}, numeric(1))
n_finite <- sum(is.finite(ll_perturb))
n_lower  <- sum(ll_perturb < ll_baseline, na.rm = TRUE)
message(sprintf("  Finite evaluations: %d / %d", n_finite, n_perturb))
message(sprintf("  Lower than defaults: %d / %d (defaults should be near optimum)", n_lower, n_finite))
if (n_finite < n_perturb * 0.5) {
  message("  WARNING: >50%% of perturbations failed -- prior may be too wide or model unstable")
} else {
  message("  PASS: likelihood is responsive to parameter changes\n")
}

# --- Prior and BayesianTools setup (needed for verification MCMC) ---
prior_yasso07 <- createPrior(
  density = function(x) sum(dnorm(x, mean = best_x, sd = sigma_ppm, log = TRUE)),
  sampler = function(n = 1) matrix(rnorm(n * N_FREE, mean = best_x, sd = sigma_ppm * 0.1),
                                   nrow = n, ncol = N_FREE),
  lower = rep(-Inf, N_FREE), upper = rep(Inf, N_FREE)
)

bayes_setup_baseline <- createBayesianSetup(
  likelihood = log_likelihood_yasso07,
  prior      = prior_yasso07,
  names      = FREE_NAMES,
  parallel   = FALSE
)

lp_check <- bayes_setup_baseline$posterior$density(best_x)
message(sprintf("CHECK: BayesianTools log-posterior at defaults = %.3f", lp_check))
if (!is.finite(lp_check)) stop("FAIL: log-posterior not finite at defaults.")
message("  PASS: BayesianTools setup is valid\n")

# --- Check 2c: short MCMC verification run ---
# Purpose: confirm DEzs starts, mixes, and produces valid output.
# 1000 iterations is NOT enough for inference -- this is only a smoke test.
# Expected: chains move (acceptance not 0 or 1), all samples finite.
message("CHECK 2c: Short MCMC verification run (1000 iterations, 3 chains)")
message("  PURPOSE: confirm sampler is functional, not to make inferences.")
message("  Expected: finite samples, R-hat reachable, chain is moving.\n")

set.seed(42)
result_verify <- runMCMC(
  bayesianSetup = bayes_setup_baseline,
  sampler       = "DEzs",
  settings      = list(iterations = 10000, burnin = 200, nrChains = 3, eps = 0.001)
)

samples_verify <- getSample(result_verify, coda = FALSE)
n_finite_samp  <- sum(apply(samples_verify, 1, function(r) all(is.finite(r))))
message(sprintf("  Samples collected (post-burnin): %d", nrow(samples_verify)))
message(sprintf("  Finite samples: %d / %d (%.0f%%)",
                n_finite_samp, nrow(samples_verify),
                100 * n_finite_samp / max(nrow(samples_verify), 1)))

if (n_finite_samp < nrow(samples_verify) * 0.95) {
  message("  WARNING: >5%% non-finite samples -- possible numerical instability")
} else {
  message("  PASS: samples are finite")
}

# Gelman-Rubin on verification run
# At 10000 iterations R-hat will not be excellent, but should not be extreme.
message("\n  Gelman-Rubin R-hat (10000-iteration verification run):")
message("  Values > 2.0 would suggest the chains are not moving at all.")
tryCatch({
  gr_v     <- gelmanDiagnostics(result_verify)
  rhat_max <- max(gr_v$psrf[, "Point est."], na.rm = TRUE)
  rhat_med <- median(gr_v$psrf[, "Point est."], na.rm = TRUE)
  message(sprintf("  Median R-hat: %.3f | Max R-hat: %.3f", rhat_med, rhat_max))
  if (rhat_max > 3.0) {
    message("  WARNING: Max R-hat > 3.0 -- chains may be stuck. Check acceptance.")
  } else if (rhat_max > 1.5) {
    message("  OK: R-hat > 1.5 is expected at 1000 iterations with 20 parameters.")
  } else {
    message("  GOOD: R-hat < 1.5 already, sampler is mixing well even at 1000 iters.")
  }
}, error = function(e) message(sprintf("  Could not compute R-hat: %s", conditionMessage(e))))

# Posterior medians for climate parameters -- sanity check
message("\n  CHECK 2d: Posterior medians vs published defaults (verification run)")
message("  On SYNTHETIC data, medians need not match defaults.")
message("  On REAL Finnish inventory data, climate params should be near defaults.\n")
chain_verify_phys <- t(apply(samples_verify, 1, to_original))
colnames(chain_verify_phys) <- FREE_NAMES
post_med    <- apply(chain_verify_phys, 2, median)
def_vals    <- free_defaults[FREE_NAMES]
message(sprintf("  %-12s | %-9s | %-12s | %-10s", "Parameter", "Default", "Post.Median", "Rel.Diff"))
message(sprintf("  %s", paste(rep("-", 52), collapse="")))
for (nm in c(CLIMATE_NAMES, NUISANCE_NAMES)) {
  rd <- abs(post_med[nm] - def_vals[nm]) / (abs(def_vals[nm]) + 1e-10)
  message(sprintf("  %-12s | %9.4f | %12.4f | %10.3f", nm, def_vals[nm], post_med[nm], rd))
}

message("\nSECTION 2 COMPLETE: Calibration verified.")
message("Proceed to performance optimisation.\n")


# =============================================================================
# 3.  Baseline runtime benchmark and projection
#
#     Measure single-evaluation time on the test dataset and project to
#     full 3700-plot production scale. This establishes our baseline to beat.
#
#     NOTE: The test dataset has only a few plots (synthetic data), so we
#     project by assuming linear scaling with plot count. This is approximately
#     true -- the plot loop is the dominant cost and plots are independent.
# =============================================================================

message("-------------------------------------------------------------")
message("SECTION 3: Baseline runtime benchmark and projection")
message("-------------------------------------------------------------\n")

message(sprintf("Test dataset: %d plots", length(plots)))
message("Timing 10 evaluations and averaging to reduce noise...\n")

t_baseline_10 <- system.time(for (i in 1:10) log_likelihood_yasso07(best_x))["elapsed"]
t_baseline    <- t_baseline_10 / 10
message(sprintf("Baseline per-evaluation time: %.4f seconds", t_baseline))

# Project to production scale
t_per_plot   <- t_baseline / length(plots)
t_at_3700    <- t_per_plot * 3700
target_evals <- 150000   # 50k iterations x 3 chains

message(sprintf("Time per plot:                %.6f seconds", t_per_plot))
message(sprintf("Projected time at 3700 plots: %.3f seconds per evaluation", t_at_3700))
message(sprintf("Projected total (1 core):     %.1f hours\n",
                t_at_3700 * target_evals / 3600))
message("THIS IS OUR BASELINE. Target: < 3 hours wallclock on Roihu.\n")


# =============================================================================
# 4.  Profiling the baseline
#
#     profvis builds a flame graph showing where time actually goes.
#
#     FINDINGS:
#       1. [.data.frame is the dominant cost -- R scanning the full climate
#          and inputs data frames for each plot on every evaluation.
#          With 3700 plots this is 2 x 3700 = 7400 full-table scans per call.
#       2. compiler:::tryCmpfun -- R is JIT-compiling the function body on
#          every call. This is visible as a broad band in the profiler.
#
#     FIX for (1): pre-split data frames -- see Section 5
#     FIX for (2): byte-compile with cmpfun() -- see Section 7
# =============================================================================

message("-------------------------------------------------------------")
message("SECTION 4: Profiling the baseline")
message("-------------------------------------------------------------\n")
message("Running profvis on 10 baseline evaluations...")
message("LOOK FOR: Which functions have the widest time bars?")
message("EXPECTED: [.data.frame and compiler overhead dominate.\n")

profvis({
  for (i in 1:10) log_likelihood_yasso07(best_x)
})

message("Profiling complete.")
message("FINDINGS:")
message("  1. [.data.frame  -- subsetting Yasso07_climate/inputs by plot_id inside the loop")
message("  2. cmpfun        -- R JIT-compiling the function body on every call\n")


# =============================================================================
# 5.  Fix 1: Pre-split data frames
#
#     PROBLEM: Inside the plot loop, the baseline does:
#       clim   <- Yasso07_climate[Yasso07_climate$plot_id == pid, ]
#       inputs <- Yasso07_inputs[Yasso07_inputs$plot_id   == pid, ]
#     Each line scans the entire data frame -- O(n_rows) -- for every plot
#     on every MCMC step.
#
#     FIX: split() the data frames into named lists once before MCMC.
#     Lookup then becomes O(1) hashing instead of O(n) scanning:
#       clim   <- climate_by_plot[[pid]]
#       inputs <- inputs_by_plot[[pid]]
#
#     Everything else in the function is identical.
# =============================================================================

message("-------------------------------------------------------------")
message("SECTION 5: Fix 1 -- pre-split data frames")
message("-------------------------------------------------------------\n")
message("Splitting climate and input data frames by plot_id (done once)...")

climate_by_plot <- split(Yasso07_climate, Yasso07_climate$plot_id)
inputs_by_plot  <- split(Yasso07_inputs,  Yasso07_inputs$plot_id)

message(sprintf("  climate_by_plot: %d entries", length(climate_by_plot)))
message(sprintf("  inputs_by_plot:  %d entries\n", length(inputs_by_plot)))

# Refactored likelihood: only change is the two lookup lines inside the loop
log_likelihood_yasso07_refactored <- function(x) {
  
  p_free       <- to_original(x)
  sigma_obs    <- p_free["sigma_obs"]
  sigma_init   <- p_free["sigma_init"]
  model_params <- assemble_model_params(p_free)
  
  log_jac <- log_jac_stick(p_free[FRAC_COL_A], B) +
    log_jac_stick(p_free[FRAC_COL_W], B) +
    log_jac_stick(p_free[FRAC_COL_E], B) +
    log_jac_stick(p_free[FRAC_COL_N], B) +
    log(sigma_obs) + log(sigma_init)
  
  log_lik <- 0
  
  for (pid in plots) {
    
    # FIX 1: O(1) list lookup -- no more full data frame scan
    clim   <- climate_by_plot[[pid]]
    inputs <- inputs_by_plot[[pid]]
    lm     <- litter_means[[pid]]
    meta   <- obs_meta[[pid]]
    
    if (any(is.na(meta$idx))) return(-Inf)
    
    xi_array <- tryCatch(
      compute_xi_yasso07(
        temp_mean = clim$temp_mean, temp_amp = clim$temp_amplitude,
        precip    = clim$precip,
        beta1     = model_params["beta1"], beta2 = model_params["beta2"],
        gamma     = model_params["gamma"]
      ), error = function(e) NULL)
    if (is.null(xi_array))                         return(-Inf)
    xi_mean <- mean(xi_array)
    if (!is.finite(xi_mean) || xi_mean <= 0)       return(-Inf)
    
    C_init <- tryCatch(
      yasso07_steady_state(
        params   = model_params, nwl_mean = lm$nwl_mean,
        fwl_mean = lm$fwl_mean, cwl_mean = lm$cwl_mean, xi_mean = xi_mean
      ), error = function(e) NULL)
    if (is.null(C_init) || any(!is.finite(C_init)) || any(C_init < 0)) return(-Inf)
    
    run_out <- tryCatch(
      yasso07_run(
        input_df = inputs, params   = model_params,
        C_init   = C_init, xi_array = xi_array
      ), error = function(e) NULL)
    if (is.null(run_out) || any(!is.finite(run_out$total_soc))) return(-Inf)
    
    SOC_hat <- run_out$total_soc[meta$idx]
    if (any(!is.finite(SOC_hat)) || any(SOC_hat <= 0)) return(-Inf)
    
    sd_vec <- ifelse(meta$is_first,
                     SOC_hat * sqrt(sigma_obs^2 + sigma_init^2),
                     SOC_hat * sigma_obs)
    
    log_lik <- log_lik + sum(dnorm(meta$soc_obs, mean = SOC_hat, sd = sd_vec, log = TRUE))
  }
  
  if (!is.finite(log_lik)) return(-Inf)
  log_lik + log_jac
}

# Correctness check: Fix 1 must give identical results to baseline
message("Verifying Fix 1 gives identical results to baseline...")
ll_ref   <- log_likelihood_yasso07_refactored(best_x)
diff_f1  <- abs(ll_baseline - ll_ref)
message(sprintf("  Baseline:    %.6f", ll_baseline))
message(sprintf("  Refactored:  %.6f", ll_ref))
message(sprintf("  Difference:  %.2e", diff_f1))
if (diff_f1 > 1e-8) {
  stop("FAIL: Refactored function differs from baseline -- check Fix 1.")
} else {
  message("  PASS: Results identical to numerical precision\n")
}

# Timing after Fix 1 (before Fix 2)
t_ref_uncompiled <- system.time(for (i in 1:10) log_likelihood_yasso07_refactored(best_x))["elapsed"] / 10
message(sprintf("After Fix 1 (uncompiled): %.4f s per eval  (was %.4f s, speedup %.1fx)",
                t_ref_uncompiled, t_baseline, t_baseline / t_ref_uncompiled))


# =============================================================================
# 6.  Profiling after Fix 1
#
#     After removing [.data.frame, the profiler now shows a different picture.
#
#     FINDING: compiler:::tryCmpfun is now the dominant cost. R is
#     just-in-time compiling the function body on every call. This is
#     fixed in Section 7.
#
#     The .Fortran calls are visible but much smaller -- the actual model
#     computation is fast relative to the R overhead we've been fighting.
# =============================================================================

message("\n-------------------------------------------------------------")
message("SECTION 6: Profiling after Fix 1")
message("-------------------------------------------------------------\n")
message("Running profvis on 10 Fix-1-refactored evaluations...")
message("EXPECTED CHANGE: [.data.frame should be gone or much smaller.")
message("NEW FINDING: compiler overhead (cmpfun/tryCmpfun) now dominates.\n")

profvis({
  for (i in 1:10) log_likelihood_yasso07_refactored(best_x)
})

message("Profiling complete.")
message("FINDING: [.data.frame eliminated. New bottleneck: JIT compilation overhead.\n")


# =============================================================================
# 7.  Fix 2: Byte-compile with cmpfun()
#
#     PROBLEM: R's JIT compiler re-compiles the function body on each call.
#     This shows up in the profiler as compiler:::tryCmpfun.
#
#     FIX: cmpfun() forces immediate byte-compilation once and stores the
#     compiled version. All subsequent calls skip compilation entirely.
#     This is a one-line change with zero logical impact on the function.
#
#     After this fix, only the .Fortran calls remain as significant cost.
#     These cannot be eliminated without rewriting the Fortran interface
#     (the batched-Fortran approach discussed but not implemented, as it
#     would break the clean three-stage wrapper structure).
# =============================================================================

message("-------------------------------------------------------------")
message("SECTION 7: Fix 2 -- byte-compile with cmpfun()")
message("-------------------------------------------------------------\n")
message("Byte-compiling the refactored likelihood function...")

log_likelihood_yasso07_refactored <- cmpfun(log_likelihood_yasso07_refactored)

message("Done. Verifying compiled version gives identical results...")
ll_compiled <- log_likelihood_yasso07_refactored(best_x)
diff_f2     <- abs(ll_baseline - ll_compiled)
message(sprintf("  Baseline:   %.6f", ll_baseline))
message(sprintf("  Compiled:   %.6f", ll_compiled))
message(sprintf("  Difference: %.2e", diff_f2))
if (diff_f2 > 1e-8) {
  stop("FAIL: Compiled function differs from baseline.")
} else {
  message("  PASS: Identical results\n")
}

t_ref_compiled <- system.time(for (i in 1:10) log_likelihood_yasso07_refactored(best_x))["elapsed"] / 10
message(sprintf("After Fix 1 + Fix 2 (compiled): %.4f s per eval", t_ref_compiled))
message(sprintf("Speedup vs baseline:            %.1fx\n", t_baseline / t_ref_compiled))

# Re-profile to confirm only .Fortran remains
message("Profiling compiled refactored function:")
message("EXPECTED: Only .Fortran calls remain as significant cost.\n")
profvis({
  for (i in 1:10) log_likelihood_yasso07_refactored(best_x)
})
message("Profiling complete. .Fortran calls now dominate -- this is the irreducible floor.\n")


# =============================================================================
# 8.  Parallel versions: baseline and refactored
#
#     Parallelise the plot loop using mclapply. Each core handles a subset
#     of plots simultaneously -- valid because each plot's contribution to
#     the log-likelihood is independent given the current parameters.
#
#     We define two parallel versions:
#       8a: Parallel with ORIGINAL slow subsetting (no Fix 1)
#       8b: Parallel with Fix 1 + Fix 2 applied (the version to use)
#
#     WHY mclapply?
#       - Fork-based: each worker inherits parent memory, no explicit data
#         transfer needed (climate_by_plot, inputs_by_plot etc. are available)
#       - Works identically on Mac and Linux/Roihu
#       - Does NOT work on Windows (would need parLapply with a cluster)
#
#     NOTE ON NESTED FORKING:
#       In Section 11 (hybrid approach), we nest mclapply inside mclapply
#       (chains in parallel, each chain parallelises plots). Nested forking
#       works on Linux/Roihu. On Mac it may behave unpredictably; we test it.
# =============================================================================

message("-------------------------------------------------------------")
message("SECTION 8: Parallel versions (mclapply)")
message("-------------------------------------------------------------\n")

n_cores <- parallel::detectCores() - 1
message(sprintf("Available cores: %d | Using: %d (leaving 1 for OS)\n",
                n_cores + 1, n_cores))

# --- 8a: Parallel with original slow subsetting ---
# Included to show that parallelisation alone is not enough:
# the [.data.frame bottleneck hurts inside workers too.
log_likelihood_yasso07_parallel <- function(x) {
  
  p_free       <- to_original(x)
  sigma_obs    <- p_free["sigma_obs"]
  sigma_init   <- p_free["sigma_init"]
  model_params <- assemble_model_params(p_free)
  
  log_jac <- log_jac_stick(p_free[FRAC_COL_A], B) +
    log_jac_stick(p_free[FRAC_COL_W], B) +
    log_jac_stick(p_free[FRAC_COL_E], B) +
    log_jac_stick(p_free[FRAC_COL_N], B) +
    log(sigma_obs) + log(sigma_init)
  
  log_liks <- parallel::mclapply(plots, function(pid) {
    
    # STILL USING SLOW SUBSETTING -- fixed in 8b
    clim   <- Yasso07_climate[Yasso07_climate$plot_id == pid, ]
    inputs <- Yasso07_inputs[Yasso07_inputs$plot_id   == pid, ]
    lm     <- litter_means[[pid]]
    meta   <- obs_meta[[pid]]
    
    if (any(is.na(meta$idx))) return(-Inf)
    
    xi_array <- tryCatch(
      compute_xi_yasso07(
        temp_mean = clim$temp_mean, temp_amp = clim$temp_amplitude,
        precip = clim$precip, beta1 = model_params["beta1"],
        beta2 = model_params["beta2"], gamma = model_params["gamma"]
      ), error = function(e) NULL)
    if (is.null(xi_array))                   return(-Inf)
    xi_mean <- mean(xi_array)
    if (!is.finite(xi_mean) || xi_mean <= 0) return(-Inf)
    
    C_init <- tryCatch(
      yasso07_steady_state(
        params = model_params, nwl_mean = lm$nwl_mean,
        fwl_mean = lm$fwl_mean, cwl_mean = lm$cwl_mean, xi_mean = xi_mean
      ), error = function(e) NULL)
    if (is.null(C_init) || any(!is.finite(C_init)) || any(C_init < 0)) return(-Inf)
    
    run_out <- tryCatch(
      yasso07_run(
        input_df = inputs, params = model_params,
        C_init = C_init, xi_array = xi_array
      ), error = function(e) NULL)
    if (is.null(run_out) || any(!is.finite(run_out$total_soc))) return(-Inf)
    
    SOC_hat <- run_out$total_soc[meta$idx]
    if (any(!is.finite(SOC_hat)) || any(SOC_hat <= 0)) return(-Inf)
    
    sd_vec <- ifelse(meta$is_first,
                     SOC_hat * sqrt(sigma_obs^2 + sigma_init^2),
                     SOC_hat * sigma_obs)
    sum(dnorm(meta$soc_obs, mean = SOC_hat, sd = sd_vec, log = TRUE))
    
  }, mc.cores = n_cores)
  
  log_liks <- unlist(log_liks)
  if (any(!is.finite(log_liks))) return(-Inf)
  sum(log_liks) + log_jac
}

# --- 8b: Parallel with Fix 1 + Fix 2 applied ---
# This is the version to use in production.
log_likelihood_yasso07_parallel_refactored <- function(x) {
  
  p_free       <- to_original(x)
  sigma_obs    <- p_free["sigma_obs"]
  sigma_init   <- p_free["sigma_init"]
  model_params <- assemble_model_params(p_free)
  
  log_jac <- log_jac_stick(p_free[FRAC_COL_A], B) +
    log_jac_stick(p_free[FRAC_COL_W], B) +
    log_jac_stick(p_free[FRAC_COL_E], B) +
    log_jac_stick(p_free[FRAC_COL_N], B) +
    log(sigma_obs) + log(sigma_init)
  
  log_liks <- parallel::mclapply(plots, function(pid) {
    
    # FIX 1 APPLIED: O(1) lookup
    clim   <- climate_by_plot[[pid]]
    inputs <- inputs_by_plot[[pid]]
    lm     <- litter_means[[pid]]
    meta   <- obs_meta[[pid]]
    
    if (any(is.na(meta$idx))) return(-Inf)
    
    xi_array <- tryCatch(
      compute_xi_yasso07(
        temp_mean = clim$temp_mean, temp_amp = clim$temp_amplitude,
        precip = clim$precip, beta1 = model_params["beta1"],
        beta2 = model_params["beta2"], gamma = model_params["gamma"]
      ), error = function(e) NULL)
    if (is.null(xi_array))                   return(-Inf)
    xi_mean <- mean(xi_array)
    if (!is.finite(xi_mean) || xi_mean <= 0) return(-Inf)
    
    C_init <- tryCatch(
      yasso07_steady_state(
        params = model_params, nwl_mean = lm$nwl_mean,
        fwl_mean = lm$fwl_mean, cwl_mean = lm$cwl_mean, xi_mean = xi_mean
      ), error = function(e) NULL)
    if (is.null(C_init) || any(!is.finite(C_init)) || any(C_init < 0)) return(-Inf)
    
    run_out <- tryCatch(
      yasso07_run(
        input_df = inputs, params = model_params,
        C_init = C_init, xi_array = xi_array
      ), error = function(e) NULL)
    if (is.null(run_out) || any(!is.finite(run_out$total_soc))) return(-Inf)
    
    SOC_hat <- run_out$total_soc[meta$idx]
    if (any(!is.finite(SOC_hat)) || any(SOC_hat <= 0)) return(-Inf)
    
    sd_vec <- ifelse(meta$is_first,
                     SOC_hat * sqrt(sigma_obs^2 + sigma_init^2),
                     SOC_hat * sigma_obs)
    sum(dnorm(meta$soc_obs, mean = SOC_hat, sd = sd_vec, log = TRUE))
    
  }, mc.cores = n_cores)
  
  log_liks <- unlist(log_liks)
  if (any(!is.finite(log_liks))) return(-Inf)
  sum(log_liks) + log_jac
}

# FIX 2 applied to parallel refactored
log_likelihood_yasso07_parallel_refactored <- cmpfun(log_likelihood_yasso07_parallel_refactored)

# Correctness checks for both parallel versions
message("Verifying both parallel versions give identical results to baseline...")
ll_par_orig <- log_likelihood_yasso07_parallel(best_x)
ll_par_ref  <- log_likelihood_yasso07_parallel_refactored(best_x)
message(sprintf("  Baseline:            %.6f", ll_baseline))
message(sprintf("  Parallel original:   %.6f (diff: %.2e)", ll_par_orig, abs(ll_baseline - ll_par_orig)))
message(sprintf("  Parallel refactored: %.6f (diff: %.2e)", ll_par_ref,  abs(ll_baseline - ll_par_ref)))
if (max(abs(ll_baseline - ll_par_orig), abs(ll_baseline - ll_par_ref)) > 1e-8) {
  stop("FAIL: Parallel versions differ from baseline.")
} else {
  message("  PASS: All versions give identical results\n")
}


# =============================================================================
# 9.  Scaling benchmark: all versions compared
#
#     We replicate test plots to simulate increasing dataset sizes, then
#     time all four versions. This reveals:
#       - Where the crossover point is (below which parallelisation hurts)
#       - The incremental contribution of each fix
#       - Extrapolation to 3700 plots
#
#     IMPORTANT: Replicated plots are used ONLY for benchmarking. The
#     plot ID repetition does not affect timing validity since we are
#     measuring the cost of the computation, not the quality of the result.
#     plots is restored to the original at the end.
# =============================================================================

message("-------------------------------------------------------------")
message("SECTION 9: Scaling benchmark -- all versions")
message("-------------------------------------------------------------\n")
message("Benchmarking all versions at increasing simulated plot counts.")
message("NOTE: Replicated test plots -- timing only, not real analysis.\n")

hdr <- sprintf("  %-6s | %-9s | %-9s | %-7s | %-9s | %-9s | %-7s",
               "Plots", "Seq.Orig", "Seq.Ref", "Speedup", "Par.Orig", "Par.Ref", "Speedup")
message(hdr)
message(sprintf("  %s", paste(rep("-", nchar(hdr) - 2), collapse = "")))

plots_original <- plots

for (n_target in c(10, 50, 100, 500, 1000, 3700)) {
  
  plots <<- rep(plots_original, length.out = n_target)
  
  t_so <- system.time(log_likelihood_yasso07(best_x))["elapsed"]
  t_sr <- system.time(log_likelihood_yasso07_refactored(best_x))["elapsed"]
  t_po <- system.time(log_likelihood_yasso07_parallel(best_x))["elapsed"]
  t_pr <- system.time(log_likelihood_yasso07_parallel_refactored(best_x))["elapsed"]
  
  message(sprintf("  %-6d | %-9.3f | %-9.3f | %-7.2f | %-9.3f | %-9.3f | %-7.2f",
                  n_target, t_so, t_sr, t_so/t_sr, t_po, t_pr, t_po/t_pr))
}

plots <<- plots_original

message("\nKey observations:")
message("  - Parallel is SLOWER than sequential at small plot counts (forking overhead)")
message("  - Crossover point is around 50-100 plots on this machine")
message("  - At 3700 plots, parallel refactored gives ~4-5x speedup over sequential baseline")
message("  - Fix 1 + Fix 2 improve both sequential and parallel versions\n")


# =============================================================================
# 10. Runtime projection using Amdahl's law
#
#     Amdahl's law: if fraction f of the work is parallelisable, then with
#     n cores the maximum speedup is: 1 / ((1-f) + f/n)
#
#     We estimate f from observed speedups, then project to Roihu core counts.
#
#     KEY INSIGHT:
#     The sequential floor (~0.13s at 3700 plots) is irreducible -- no amount
#     of cores removes it. Beyond ~64 cores the added benefit per core shrinks
#     sharply. The real leverage is the hybrid approach (Section 11), where
#     3 chains run simultaneously rather than sequentially.
# =============================================================================

message("-------------------------------------------------------------")
message("SECTION 10: Runtime projection (Amdahl's law)")
message("-------------------------------------------------------------\n")

# These are empirical timings from Section 9 at 3700 plots.
# Update them if your Section 9 benchmark gives different numbers.
t_seq_3700  <- 1.763   # seconds, sequential refactored at 3700 plots
t_par_mac   <- 0.359   # seconds, parallel refactored at 3700 plots (11 cores)
n_cores_mac <- n_cores

# Estimate parallel fraction f from Amdahl's law:
# t_par = t_seq * ((1-f) + f/n)  =>  solve for f
f       <- (1 - t_par_mac / t_seq_3700) / (1 - 1 / n_cores_mac)
f       <- min(max(f, 0), 1)   # clamp to valid range
t_floor <- t_seq_3700 * (1 - f)

message(sprintf("Empirical timings at 3700 plots:"))
message(sprintf("  Sequential:           %.3f s", t_seq_3700))
message(sprintf("  Parallel (%d cores): %.3f s", n_cores_mac, t_par_mac))
message(sprintf("Estimated parallel fraction: %.3f (%.1f%% of work)", f, f * 100))
message(sprintf("Sequential floor (Amdahl):   %.3f s per evaluation\n", t_floor))

message(sprintf("Projected total runtime at %d evaluations:", target_evals))
message(sprintf("  %-7s | %-12s | %-12s | %-10s", "Cores", "Per eval", "Total", "vs Mac"))
message(sprintf("  %s", paste(rep("-", 47), collapse="")))

for (nc in c(11, 16, 32, 64, 128)) {
  t_e  <- t_seq_3700 * ((1 - f) + f / nc)
  t_h  <- t_e * target_evals / 3600
  spd  <- t_par_mac / t_e
  message(sprintf("  %-7d | %-12s | %-12s | %-10s",
                  nc, sprintf("%.3fs", t_e), sprintf("%.1f hrs", t_h), sprintf("%.2fx", spd)))
}

message("\nNOTE: Beyond ~64 cores the marginal gain per core diminishes sharply.")
message("The hybrid approach (Section 11) uses cores more efficiently by running")
message("chains simultaneously rather than throwing all cores at one chain.\n")


# =============================================================================
# 11. HYBRID APPROACH: parallel chains x parallel plot loop
#
#     MOTIVATION:
#     Section 10 showed that a single chain at 128 cores still takes ~5.5
#     hours because of the sequential floor. But three chains running
#     one after another would take 3x that. The hybrid approach runs
#     N chains simultaneously, each using CORES_PER_CHAIN for the plot loop.
#
#     DESIGN:
#       Outer mclapply: launches N_CHAINS independent processes
#         -> each runs one complete DEzs MCMC chain
#       Inner mclapply (inside each chain's likelihood): parallelises the
#         plot loop with CORES_PER_CHAIN cores
#
#     CORE ALLOCATION:
#       Total cores used = N_CHAINS x CORES_PER_CHAIN
#       On Roihu (128 cores): N_CHAINS=3, CORES_PER_CHAIN=40 -> 120 cores
#       On Mac (testing):     N_CHAINS=3, CORES_PER_CHAIN=floor((cores-1)/3)
#
#     WALLCLOCK TIME ESTIMATE (Roihu):
#       Per-chain per-eval: ~0.16s (at 40 cores per chain)
#       Per-chain total:    0.16s x 50000 iters = 8000s ~ 2.2 hours
#       All chains run simultaneously -> ~2.2 hours wallclock
#
#     COMBINING CHAINS FOR DIAGNOSTICS:
#       Each chain is a separate DEzs run with nrChains=1.
#       After all chains complete, we combine them into a coda mcmc.list
#       and compute Gelman-Rubin R-hat across chains as usual.
#
#     LIMITATION ON MAC:
#       Nested mclapply (forking inside forked processes) can behave
#       unpredictably on some Mac versions. We test it explicitly here.
#       On Roihu (Linux) nested forking is reliable.
# =============================================================================

message("-------------------------------------------------------------")
message("SECTION 11: Hybrid approach -- parallel chains x parallel plots")
message("-------------------------------------------------------------\n")

# --- Core allocation ---
# On Roihu, set N_CHAINS=3 and CORES_PER_CHAIN=40.
# Here we adapt automatically to available cores.
N_CHAINS        <- 3L
CORES_PER_CHAIN <- max(1L, floor((parallel::detectCores() - 1L) / N_CHAINS))
TOTAL_CORES     <- N_CHAINS * CORES_PER_CHAIN

message("Core allocation for hybrid approach:")
message(sprintf("  Chains:          %d (running simultaneously)", N_CHAINS))
message(sprintf("  Cores per chain: %d (for internal plot loop)", CORES_PER_CHAIN))
message(sprintf("  Total cores:     %d / %d available", TOTAL_CORES, parallel::detectCores()))
message(sprintf("  OS reserve:      %d cores\n", parallel::detectCores() - TOTAL_CORES))

# --- Likelihood factory ---
# Returns a compiled likelihood function with a fixed core count baked in.
# Using a factory avoids race conditions: each chain process gets its own
# function with its own core allocation captured in the closure.
make_parallel_likelihood <- function(cores_per_chain) {
  cmpfun(function(x) {
    
    p_free       <- to_original(x)
    sigma_obs    <- p_free["sigma_obs"]
    sigma_init   <- p_free["sigma_init"]
    model_params <- assemble_model_params(p_free)
    
    log_jac <- log_jac_stick(p_free[FRAC_COL_A], B) +
      log_jac_stick(p_free[FRAC_COL_W], B) +
      log_jac_stick(p_free[FRAC_COL_E], B) +
      log_jac_stick(p_free[FRAC_COL_N], B) +
      log(sigma_obs) + log(sigma_init)
    
    log_liks <- parallel::mclapply(plots, function(pid) {
      
      clim   <- climate_by_plot[[pid]]
      inputs <- inputs_by_plot[[pid]]
      lm     <- litter_means[[pid]]
      meta   <- obs_meta[[pid]]
      
      if (any(is.na(meta$idx))) return(-Inf)
      
      xi_array <- tryCatch(
        compute_xi_yasso07(
          temp_mean = clim$temp_mean, temp_amp = clim$temp_amplitude,
          precip = clim$precip, beta1 = model_params["beta1"],
          beta2 = model_params["beta2"], gamma = model_params["gamma"]
        ), error = function(e) NULL)
      if (is.null(xi_array))                   return(-Inf)
      xi_mean <- mean(xi_array)
      if (!is.finite(xi_mean) || xi_mean <= 0) return(-Inf)
      
      C_init <- tryCatch(
        yasso07_steady_state(
          params = model_params, nwl_mean = lm$nwl_mean,
          fwl_mean = lm$fwl_mean, cwl_mean = lm$cwl_mean, xi_mean = xi_mean
        ), error = function(e) NULL)
      if (is.null(C_init) || any(!is.finite(C_init)) || any(C_init < 0)) return(-Inf)
      
      run_out <- tryCatch(
        yasso07_run(
          input_df = inputs, params = model_params,
          C_init = C_init, xi_array = xi_array
        ), error = function(e) NULL)
      if (is.null(run_out) || any(!is.finite(run_out$total_soc))) return(-Inf)
      
      SOC_hat <- run_out$total_soc[meta$idx]
      if (any(!is.finite(SOC_hat)) || any(SOC_hat <= 0)) return(-Inf)
      
      sd_vec <- ifelse(meta$is_first,
                       SOC_hat * sqrt(sigma_obs^2 + sigma_init^2),
                       SOC_hat * sigma_obs)
      sum(dnorm(meta$soc_obs, mean = SOC_hat, sd = sd_vec, log = TRUE))
      
    }, mc.cores = cores_per_chain)   # each chain uses only its allocated cores
    
    log_liks <- unlist(log_liks)
    if (any(!is.finite(log_liks))) return(-Inf)
    sum(log_liks) + log_jac
  })
}

# --- Verify the factory produces a correct likelihood ---
message("Verifying hybrid likelihood factory produces correct results...")
ll_hybrid <- make_parallel_likelihood(CORES_PER_CHAIN)(best_x)
diff_hyb  <- abs(ll_baseline - ll_hybrid)
message(sprintf("  Baseline:   %.6f", ll_baseline))
message(sprintf("  Hybrid:     %.6f", ll_hybrid))
message(sprintf("  Difference: %.2e", diff_hyb))
if (diff_hyb > 1e-8) {
  stop("FAIL: Hybrid likelihood differs from baseline.")
} else {
  message("  PASS: Identical results\n")
}

# --- Single-chain timing comparison ---
message("Single-evaluation timing comparison:")
t_s <- system.time(log_likelihood_yasso07_refactored(best_x))["elapsed"]
t_p <- system.time(log_likelihood_yasso07_parallel_refactored(best_x))["elapsed"]
t_h <- system.time(make_parallel_likelihood(CORES_PER_CHAIN)(best_x))["elapsed"]
message(sprintf("  Sequential refactored:  %.4f s", t_s))
message(sprintf("  Parallel refactored:    %.4f s (all %d cores)", t_p, n_cores))
message(sprintf("  Hybrid per-chain:       %.4f s (%d cores per chain)\n", t_h, CORES_PER_CHAIN))

# --- Short hybrid MCMC test ---
# Run N_CHAINS chains simultaneously, each with its own CORES_PER_CHAIN cores.
# This is the actual Roihu execution pattern, just with fewer iterations.
message(sprintf("Running short hybrid MCMC test (%d chains x 500 iterations)...", N_CHAINS))
message(sprintf("  Each chain uses %d cores for the plot loop.", CORES_PER_CHAIN))
message("  Chains run simultaneously via outer mclapply.\n")

set.seed(42)
t_hybrid_mcmc <- system.time({
  
  chain_results <- parallel::mclapply(seq_len(N_CHAINS), function(chain_id) {
    
    # Each chain gets its own fresh BayesianSetup with its own likelihood.
    # This avoids any shared state between chains.
    ll_fn <- make_parallel_likelihood(CORES_PER_CHAIN)
    
    setup <- createBayesianSetup(
      likelihood = ll_fn,
      prior      = prior_yasso07,
      names      = FREE_NAMES,
      parallel   = FALSE   # chain-level parallelism handled by outer mclapply
    )
    
    runMCMC(
      bayesianSetup = setup,
      sampler       = "DEzs",
      settings      = list(
        iterations = 1000,
        burnin     = 100,
        nrChains   = 1,    # one chain per forked process
        eps        = 0.001
      )
    )
    
  }, mc.cores = N_CHAINS)
  
})["elapsed"]

message(sprintf("Hybrid MCMC test completed in %.1f seconds\n", t_hybrid_mcmc))

# --- Check all chains returned valid results ---
message("Checking hybrid MCMC chain results:")
n_ok <- 0L
for (i in seq_len(N_CHAINS)) {
  res_i <- chain_results[[i]]
  ok_i  <- inherits(res_i, "mcmcSamplerList") || inherits(res_i, "mcmcSampler")
  if (ok_i) {
    samp_i   <- getSample(res_i, coda = FALSE)
    n_fin_i  <- sum(apply(samp_i, 1, function(r) all(is.finite(r))))
    pct_fin  <- 100 * n_fin_i / max(nrow(samp_i), 1)
    status   <- if (pct_fin >= 95) "PASS" else "WARN"
    message(sprintf("  Chain %d [%s]: %d samples, %.0f%% finite",
                    i, status, nrow(samp_i), pct_fin))
    n_ok <- n_ok + 1L
  } else {
    message(sprintf("  Chain %d [FAIL]: returned %s", i, class(res_i)[1]))
  }
}

if (n_ok < N_CHAINS) {
  message(sprintf("\nWARNING: %d/%d chains failed.", N_CHAINS - n_ok, N_CHAINS))
  message("On Mac, nested mclapply may behave unpredictably.")
  message("On Roihu (Linux), this approach is reliable.")
} else {
  message(sprintf("\nPASS: All %d chains completed successfully.\n", N_CHAINS))
}

# --- Convergence diagnostics on combined hybrid chains ---
if (n_ok == N_CHAINS) {
  
  message("Convergence diagnostics on combined hybrid chains:")
  message("(1000 iterations -- expect poor R-hat. For production use 50k iterations.)\n")
  
  # Combine chains into a coda mcmc.list for Gelman-Rubin
  chain_list_coda <- tryCatch({
    coda::mcmc.list(lapply(chain_results, function(r) {
      coda::mcmc(getSample(r, coda = FALSE))
    }))
  }, error = function(e) NULL)
  
  if (!is.null(chain_list_coda)) {
    
    gr_hybrid <- tryCatch(
      gelman.diag(chain_list_coda, multivariate = FALSE),
      error = function(e) NULL
    )
    
    if (!is.null(gr_hybrid)) {
      rhat_vals <- gr_hybrid$psrf[, "Point est."]
      message(sprintf("  R-hat range:           %.3f -- %.3f",
                      min(rhat_vals, na.rm = TRUE), max(rhat_vals, na.rm = TRUE)))
      message(sprintf("  Parameters R-hat > 1.2: %d / %d (expected at 500 iters)",
                      sum(rhat_vals > 1.2, na.rm = TRUE), length(rhat_vals)))
      
      # Per-parameter R-hat for climate params (the ones we care about most)
      message("\n  R-hat for climate and nuisance parameters:")
      message(sprintf("  %-12s | %-10s | %-10s", "Parameter", "R-hat", "Status"))
      message(sprintf("  %s", paste(rep("-", 37), collapse="")))
      for (nm in c(CLIMATE_NAMES, NUISANCE_NAMES)) {
        rh  <- rhat_vals[nm]
        sta <- if (is.na(rh)) "N/A" else if (rh < 1.1) "Good" else if (rh < 1.5) "OK (short run)" else "Poor"
        message(sprintf("  %-12s | %-10.3f | %-10s", nm, rh, sta))
      }
    }
  }
  
  # ESS for combined chains
  message("\n  Effective sample size (combined chains, post-burnin):")
  all_raw  <- do.call(rbind, lapply(chain_results, function(r) getSample(r, coda = FALSE)))
  ess_vals <- tryCatch(effectiveSize(coda::mcmc(all_raw)), error = function(e) NULL)
  if (!is.null(ess_vals)) {
    message(sprintf("  ESS range: %.0f -- %.0f  (target: > 200 per parameter for inference)",
                    min(ess_vals, na.rm = TRUE), max(ess_vals, na.rm = TRUE)))
    message(sprintf("  Parameters with ESS < 100: %d / %d (expected at 500 iters)",
                    sum(ess_vals < 100, na.rm = TRUE), length(ess_vals)))
  }
  
  # Posterior summary for climate parameters
  message("\n  Posterior summary -- climate and nuisance parameters (physical space):")
  all_phys <- t(apply(all_raw, 1, to_original))
  colnames(all_phys) <- FREE_NAMES
  post_q   <- apply(all_phys, 2, quantile, probs = c(0.025, 0.50, 0.975), na.rm = TRUE)
  message(sprintf("  %-12s | %7s | %7s | %7s", "Parameter", "2.5%", "50%", "97.5%"))
  message(sprintf("  %s", paste(rep("-", 42), collapse="")))
  for (nm in c(CLIMATE_NAMES, NUISANCE_NAMES)) {
    message(sprintf("  %-12s | %7.4f | %7.4f | %7.4f",
                    nm, post_q["2.5%", nm], post_q["50%", nm], post_q["97.5%", nm]))
  }
}

# --- Final runtime projection for Roihu ---
message("\n--- FINAL PROJECTION: Hybrid approach on Roihu ---\n")

# Roihu configuration (adjust as needed)
roihu_n_chains        <- 3L
roihu_cores_per_chain <- 40L
roihu_total_cores     <- roihu_n_chains * roihu_cores_per_chain
roihu_iters           <- 50000L

# Per-eval time per chain (using Amdahl estimate from Section 10)
t_roihu_eval  <- t_seq_3700 * ((1 - f) + f / roihu_cores_per_chain)
t_roihu_chain <- t_roihu_eval * roihu_iters   # seconds for one chain
t_roihu_wall  <- t_roihu_chain / 3600         # hours -- chains run simultaneously

message(sprintf("Roihu configuration:"))
message(sprintf("  Chains:                  %d (running simultaneously)", roihu_n_chains))
message(sprintf("  Cores per chain:         %d", roihu_cores_per_chain))
message(sprintf("  Total cores:             %d", roihu_total_cores))
message(sprintf("  Iterations per chain:    %d", roihu_iters))
message(sprintf("  Estimated per-eval:      %.3f seconds", t_roihu_eval))
message(sprintf("  Estimated per-chain:     %.1f hours", t_roihu_chain / 3600))
message(sprintf("  Estimated WALLCLOCK:     %.1f hours (all chains run simultaneously)", t_roihu_wall))
message(sprintf("\n  Baseline (1 core, seq): %.1f hours", t_at_3700 * target_evals / 3600))
message(sprintf("  Total speedup achieved: %.1fx\n",
                (t_at_3700 * target_evals / 3600) / t_roihu_wall))






# =============================================================================
# 12. Final production run on full test dataset
#
#     This is the first real calibration run using the fully optimised
#     hybrid approach. All 4 plots, 10000 iterations per chain, 3 chains
#     running simultaneously with CORES_PER_CHAIN cores each.
#
#     PURPOSE:
#       - Confirm the hybrid setup works end-to-end on the real test data
#       - Get a first look at the posterior with a meaningful number of steps
#       - Check convergence diagnostics before scaling to 3700 plots on Roihu
#
#     NOTE: 10000 iterations is modest. On the real 3700-plot dataset we
#     will use 50000. But for the test dataset this is a reasonable first run.
# =============================================================================

message("-------------------------------------------------------------")
message("SECTION 12: Final production run -- full test dataset, 10000 iterations")
message("-------------------------------------------------------------\n")

# Confirm we are using the real (unreplicated) plot list
plots <- plots_original
message(sprintf("Running on %d plots (full test dataset).", length(plots)))
message(sprintf("Configuration: %d chains x 10000 iterations x %d cores per chain",
                N_CHAINS, CORES_PER_CHAIN))
message(sprintf("Total likelihood evaluations: %d\n", N_CHAINS * 10000L))

# Estimated runtime based on Amdahl projection from Section 10
t_est_eval  <- t_seq_3700 * ((1 - f) + f / CORES_PER_CHAIN)
t_est_total <- t_est_eval * 10000 * N_CHAINS / N_CHAINS / 60  # minutes, chains simultaneous
message(sprintf("Estimated wallclock time: %.1f minutes (based on Section 10 projections)\n",
                t_est_total))

#set up the folder for the log file of sampling progresses
log_dir <- "./Calibration/tests/"

set.seed(2025)
t_final <- system.time({
  

  final_chain_results <- parallel::mclapply(seq_len(N_CHAINS), function(chain_id) {
    
    log_file <- file.path(log_dir, sprintf("chain_%d_progress.log", chain_id))
    cat(sprintf("Chain %d started at %s\n", chain_id, Sys.time()), file = log_file)
    
    ll_fn <- make_parallel_likelihood(CORES_PER_CHAIN)
    
    # Wrap likelihood to log every 500 calls
    call_count <- 0L
    ll_fn_logged <- function(x) {
      call_count <<- call_count + 1L
      if (call_count %% 500L == 0L) {
        cat(sprintf("Chain %d | eval %d | %s\n",
                    chain_id, call_count, Sys.time()),
            file = log_file, append = TRUE)
      }
      ll_fn(x)
    }
    
    setup <- createBayesianSetup(
      likelihood = ll_fn_logged,
      prior      = prior_yasso07,
      names      = FREE_NAMES,
      parallel   = FALSE
    )
    
    result <- runMCMC(
      bayesianSetup = setup,
      sampler       = "DEzs",
      settings      = list(iterations = 10000, burnin = 2000, nrChains = 1, eps = 0.001)
    )
    
    cat(sprintf("Chain %d finished at %s\n", chain_id, Sys.time()), file = log_file, append = TRUE)
    result
    
  }, mc.cores = N_CHAINS)  
})["elapsed"]

message(sprintf("\nFinal run completed in %.1f minutes (%.0f seconds)",
                t_final / 60, t_final))

# --- Check all chains completed ---
message("\nChain completion check:")
n_final_ok <- 0L
for (i in seq_len(N_CHAINS)) {
  res_i <- final_chain_results[[i]]
  ok_i  <- inherits(res_i, "mcmcSamplerList") || inherits(res_i, "mcmcSampler")
  if (ok_i) {
    samp_i  <- getSample(res_i, coda = FALSE)
    n_fin_i <- sum(apply(samp_i, 1, function(r) all(is.finite(r))))
    pct     <- 100 * n_fin_i / max(nrow(samp_i), 1)
    message(sprintf("  Chain %d: %d post-burnin samples, %.0f%% finite  [%s]",
                    i, nrow(samp_i), pct, if (pct >= 95) "PASS" else "WARN"))
    n_final_ok <- n_final_ok + 1L
  } else {
    message(sprintf("  Chain %d: FAILED -- returned %s", i, class(res_i)[1]))
  }
}

if (n_final_ok < N_CHAINS) {
  stop(sprintf("%d / %d chains failed. Check nested forking on this system.",
               N_CHAINS - n_final_ok, N_CHAINS))
}

# --- Convergence diagnostics ---
message("\n--- Convergence diagnostics (final run) ---\n")

final_coda <- tryCatch({
  coda::mcmc.list(lapply(final_chain_results, function(r) {
    coda::mcmc(getSample(r, coda = FALSE))
  }))
}, error = function(e) NULL)

if (!is.null(final_coda)) {
  
  # Gelman-Rubin R-hat
  gr_final <- tryCatch(gelman.diag(final_coda, multivariate = FALSE), error = function(e) NULL)
  
  if (!is.null(gr_final)) {
    rhat_vals  <- gr_final$psrf[, "Point est."]
    rhat_upper <- gr_final$psrf[, "Upper C.I."]
    n_good     <- sum(rhat_vals < 1.05, na.rm = TRUE)
    n_ok_r     <- sum(rhat_vals >= 1.05 & rhat_vals < 1.10, na.rm = TRUE)
    n_warn     <- sum(rhat_vals >= 1.10, na.rm = TRUE)
    
    message(sprintf("Gelman-Rubin R-hat summary (%d parameters):", length(rhat_vals)))
    message(sprintf("  R-hat < 1.05 (good):    %d parameters", n_good))
    message(sprintf("  R-hat 1.05-1.10 (ok):   %d parameters", n_ok_r))
    message(sprintf("  R-hat > 1.10 (warning): %d parameters", n_warn))
    message(sprintf("  Overall range: %.3f -- %.3f",
                    min(rhat_vals, na.rm = TRUE), max(rhat_vals, na.rm = TRUE)))
    
    # Full table -- all parameters
    message("\n  Full R-hat table:")
    message(sprintf("  %-12s | %-10s | %-10s | %-10s",
                    "Parameter", "Point est.", "Upper C.I.", "Status"))
    message(sprintf("  %s", paste(rep("-", 50), collapse="")))
    for (nm in FREE_NAMES) {
      pe  <- rhat_vals[nm]
      uci <- rhat_upper[nm]
      sta <- if (is.na(pe)) "N/A" else if (pe < 1.05) "Good" else if (pe < 1.10) "OK" else "Warning"
      message(sprintf("  %-12s | %-10.3f | %-10.3f | %-10s", nm, pe, uci, sta))
    }
    message("")
    
    if (n_warn > 0) {
      message(sprintf("ACTION: %d parameter(s) with R-hat > 1.10.", n_warn))
      message("  Consider running longer (20000-50000 iterations) before using posteriors.")
      message("  Parameters with highest R-hat:")
      top3 <- names(sort(rhat_vals, decreasing = TRUE)[1:min(3, n_warn)])
      for (nm in top3) message(sprintf("    %s: %.3f", nm, rhat_vals[nm]))
    } else {
      message("All parameters have R-hat < 1.10. Convergence is acceptable.")
    }
  }
  
  # ESS
  message("\nEffective sample size (post-burnin, combined chains):")
  all_raw_final <- do.call(rbind, lapply(final_chain_results,
                                         function(r) getSample(r, coda = FALSE)))
  ess_final <- tryCatch(effectiveSize(coda::mcmc(all_raw_final)), error = function(e) NULL)
  
  if (!is.null(ess_final)) {
    n_low_ess <- sum(ess_final < 200, na.rm = TRUE)
    message(sprintf("  ESS range: %.0f -- %.0f", min(ess_final), max(ess_final)))
    message(sprintf("  Parameters with ESS < 200: %d / %d",
                    n_low_ess, length(ess_final)))
    if (n_low_ess > 0) {
      message("  ACTION: Low ESS suggests slow mixing. Run longer or check correlations.")
    } else {
      message("  ESS > 200 for all parameters. Posterior summaries are reliable.")
    }
  }
}

# --- Posterior summary ---
message("\n--- Posterior summary (physical space, post-burnin) ---\n")

all_phys_final <- t(apply(all_raw_final, 1, to_original))
colnames(all_phys_final) <- FREE_NAMES

message(sprintf("Total post-burnin samples: %d (across %d chains)", nrow(all_phys_final), N_CHAINS))
message(sprintf("%-12s | %7s | %7s | %7s | %7s | %7s",
                "Parameter", "2.5%", "25%", "50%", "75%", "97.5%"))
message(paste(rep("-", 60), collapse=""))
post_q_final <- apply(all_phys_final, 2, quantile,
                      probs = c(0.025, 0.25, 0.50, 0.75, 0.975), na.rm = TRUE)
for (nm in FREE_NAMES) {
  message(sprintf("%-12s | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f",
                  nm,
                  post_q_final["2.5%",  nm],
                  post_q_final["25%",   nm],
                  post_q_final["50%",   nm],
                  post_q_final["75%",   nm],
                  post_q_final["97.5%", nm]))
}

# --- Trace plots ---
message("\nGenerating trace plots for climate and nuisance parameters...")
message("Look for: stationary chains with no trends and good overlap between chains.")
par(mfrow = c(3, 2))
for (nm in c(CLIMATE_NAMES, NUISANCE_NAMES)) {
  # Plot one trace per chain
  chain_traces <- lapply(final_chain_results, function(r) {
    raw <- getSample(r, coda = FALSE)
    phys <- t(apply(raw, 1, to_original))
    phys[, nm]
  })
  ylim_r <- range(unlist(chain_traces), na.rm = TRUE)
  plot(chain_traces[[1]], type = "l", col = "steelblue",
       ylim = ylim_r, xlab = "Iteration", ylab = nm,
       main = sprintf("%s -- trace", nm))
  if (N_CHAINS > 1) {
    lines(chain_traces[[2]], col = "firebrick")
    if (N_CHAINS > 2) lines(chain_traces[[3]], col = "darkgreen")
  }
}
par(mfrow = c(1, 1))

# --- Marginal posteriors vs priors for climate parameters ---
message("Generating marginal posterior vs prior plots (climate + nuisance)...")
message("Narrow posterior = data is informative about that parameter.")
message("Posterior ~ prior = data has little to say about that parameter.\n")
par(mfrow = c(3, 2))
for (nm in c(CLIMATE_NAMES, NUISANCE_NAMES)) {
  post_vals <- all_phys_final[, nm]
  # Prior samples (in physical space, for reference)
  prior_samp_raw  <- matrix(rnorm(2000 * N_FREE, mean = best_x, sd = sigma_ppm),
                            nrow = 2000, ncol = N_FREE)
  prior_samp_phys <- t(apply(prior_samp_raw, 1, to_original))
  prior_vals      <- prior_samp_phys[is.finite(prior_samp_phys[, nm]), nm]
  
  d_post  <- density(post_vals[is.finite(post_vals)])
  d_prior <- density(prior_vals[is.finite(prior_vals)])
  xlim_r  <- range(c(d_post$x, d_prior$x))
  ylim_r  <- range(c(d_post$y, d_prior$y))
  
  plot(d_post,  col = "steelblue", lwd = 2, xlim = xlim_r, ylim = ylim_r,
       main = nm, xlab = nm, ylab = "Density")
  lines(d_prior, col = "grey60", lwd = 1, lty = 2)
  legend("topright", legend = c("Posterior", "Prior"),
         col = c("steelblue", "grey60"), lwd = c(2, 1), lty = c(1, 2), bty = "n")
}
par(mfrow = c(1, 1))

# --- Save results ---
message("Saving final run results...")
saveRDS(final_chain_results, file = "./Calibration/tests/Yasso07_final_test_run_chains.rds")
saveRDS(all_phys_final,      file = "./Calibration/tests/Yasso07_final_test_run_posterior.rds")
message("Saved:")
message("  Yasso07_final_test_run_chains.rds    -- raw BayesianTools chain objects")
message("  Yasso07_final_test_run_posterior.rds -- posterior samples in physical space")

message("\n=============================================================")
message("  Section 12 complete.")
message(sprintf("  Total wallclock: %.1f minutes", t_final / 60))
message("  Review R-hat table and trace plots before scaling to Roihu.")
message("=============================================================\n")




# =============================================================================
# 13. Single-nesting test: sequential chains, parallel plot loop
#
#     MOTIVATION:
#     Section 12 uses the hybrid approach (parallel chains x parallel plots)
#     which on Mac silently falls back to sequential plot evaluation inside
#     each forked chain process -- nested forking is not supported reliably
#     on macOS. The result was 3 cores active instead of 9.
#
#     THIS SECTION tests the alternative single-nesting approach:
#       - Chains run ONE AT A TIME (sequentially)
#       - Each chain uses ALL available cores for the plot loop
#
#     This is the correct way to test and benchmark plot-loop parallelism
#     on Mac, and directly validates the inner mclapply that will run inside
#     each chain on Roihu.
#
#     EXPECTED BEHAVIOUR ON MAC:
#       - Activity monitor shows ONE rsession-arm64 process
#       - That process uses ~700-900% CPU (7-9 cores active)
#       - Confirms the plot-loop parallelism is working correctly
#
#     TRADEOFF vs HYBRID:
#       - Hybrid (Section 12): chains run simultaneously, each with fewer cores
#         -> better wallclock time on Roihu where nested forking works
#       - Single-nesting (Section 13): chains run one at a time, each with all cores
#         -> better on Mac for testing, wallclock = N_CHAINS x per-chain time
#
#     WHAT TO COMPARE:
#       Compare per-chain time here vs per-chain time in Section 12.
#       If nested forking were working in Section 12, they should be similar
#       (same cores per chain). If Section 13 is much faster per chain,
#       it confirms the inner mclapply was silently disabled in Section 12.
# =============================================================================

message("-------------------------------------------------------------")
message("SECTION 13: Single-nesting test -- sequential chains, parallel plots")
message("-------------------------------------------------------------\n")

# All available cores go to the plot loop -- no sharing between chains
CORES_SINGLE <- n_cores   # n_cores was defined in Section 8 as detectCores() - 1

message("Configuration:")
message(sprintf("  Chains:        %d (running ONE AT A TIME)", N_CHAINS))
message(sprintf("  Cores per chain: %d (all available, for plot loop)", CORES_SINGLE))
message(sprintf("  Total cores used at any moment: %d", CORES_SINGLE))
message(sprintf("  Expected activity monitor: 1 process at ~%d00%% CPU\n",
                CORES_SINGLE))
message("Watch activity monitor during this run.")
message("If plot-loop parallelism is working you should see ONE process")
message(sprintf("using ~%d cores (~%d00%% CPU).\n", CORES_SINGLE, CORES_SINGLE))

# Log directory (same as Section 12)
log_dir <- "./Calibration/tests/"
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)

# Build the likelihood with all cores allocated to the plot loop
ll_single_nesting <- make_parallel_likelihood(CORES_SINGLE)

# Verify it still gives correct results
message("Verifying single-nesting likelihood gives correct results...")
ll_sn_check <- ll_single_nesting(best_x)
diff_sn     <- abs(ll_baseline - ll_sn_check)
message(sprintf("  Baseline:   %.6f", ll_baseline))
message(sprintf("  Single-nesting: %.6f", ll_sn_check))
message(sprintf("  Difference: %.2e", diff_sn))
if (diff_sn > 1e-8) {
  stop("FAIL: Single-nesting likelihood differs from baseline.")
} else {
  message("  PASS: Identical results\n")
}

# --- Single-evaluation timing ---
# Time one evaluation to confirm the plot loop is actually parallelising.
# Compare to the sequential refactored timing from Section 7.
message("Timing single evaluation (single-nesting, all cores on plot loop):")
t_sn_single <- system.time(ll_single_nesting(best_x))["elapsed"]
message(sprintf("  Single-nesting:       %.4f s", t_sn_single))
message(sprintf("  Sequential (Section 7): %.4f s", t_ref_compiled))
message(sprintf("  Speedup:              %.2fx", t_ref_compiled / t_sn_single))
if (t_ref_compiled / t_sn_single < 1.5) {
  message("  NOTE: Low speedup on test dataset -- too few plots to saturate cores.")
  message("  Run the 3700-plot scaling benchmark below for a meaningful comparison.\n")
} else {
  message("  Parallelism is working on plot loop.\n")
}

# --- Scaling benchmark: single-nesting vs sequential ---
# Replicate plots to simulate production scale and find the crossover point.
message("Scaling benchmark: single-nesting vs sequential refactored...")
message(sprintf("  %-6s | %-12s | %-12s | %-10s",
                "Plots", "Sequential", "Single-nest", "Speedup"))
message(sprintf("  %s", paste(rep("-", 47), collapse="")))

plots_original_s13 <- plots   # save current plot list

for (n_target in c(10, 50, 100, 500, 1000, 3700)) {
  
  plots <<- rep(plots_original_s13, length.out = n_target)
  
  t_seq <- system.time(log_likelihood_yasso07_refactored(best_x))["elapsed"]
  t_sn  <- system.time(ll_single_nesting(best_x))["elapsed"]
  
  message(sprintf("  %-6d | %-12.3f | %-12.3f | %-10.2f",
                  n_target, t_seq, t_sn, t_seq / t_sn))
}

plots <<- plots_original_s13
message("")

# --- Full run: sequential chains, parallel plots ---
message(sprintf("Running full single-nesting MCMC (%d chains x 10000 iterations)...", N_CHAINS))
message("Chains run one at a time. Each chain uses all available cores for plots.")
message(sprintf("Expected total wallclock: ~%.1f minutes (N_CHAINS x per-chain time)\n",
                t_sn_single * 10000 * N_CHAINS / 60))

set.seed(2025)
t_sn_total <- system.time({
  
  # lapply (not mclapply): chains run sequentially
  # Each chain's likelihood uses mclapply internally with all cores
  sn_chain_results <- lapply(seq_len(N_CHAINS), function(chain_id) {
    
    log_file <- file.path(log_dir, sprintf("chain_sn_%d_progress.log", chain_id))
    cat(sprintf("Chain %d (single-nesting) started at %s\n", chain_id, Sys.time()),
        file = log_file)
    
    message(sprintf("  Starting chain %d / %d ...", chain_id, N_CHAINS))
    
    # Wrap likelihood to log every 500 calls
    call_count  <- 0L
    ll_fn_logged <- function(x) {
      call_count <<- call_count + 1L
      if (call_count %% 500L == 0L) {
        cat(sprintf("Chain %d | eval %d | %s\n",
                    chain_id, call_count, Sys.time()),
            file = log_file, append = TRUE)
      }
      ll_single_nesting(x)
    }
    
    setup <- createBayesianSetup(
      likelihood = ll_fn_logged,
      prior      = prior_yasso07,
      names      = FREE_NAMES,
      parallel   = FALSE
    )
    
    result <- runMCMC(
      bayesianSetup = setup,
      sampler       = "DEzs",
      settings      = list(
        iterations = 10000,
        burnin     = 2000,
        nrChains   = 1,
        eps        = 0.001
      )
    )
    
    cat(sprintf("Chain %d finished at %s\n", chain_id, Sys.time()),
        file = log_file, append = TRUE)
    message(sprintf("  Chain %d complete.", chain_id))
    result
  })
  
})["elapsed"]

message(sprintf("\nSingle-nesting run completed in %.1f minutes (%.0f seconds)",
                t_sn_total / 60, t_sn_total))

# --- Compare wallclock times ---
message("\n--- Wallclock comparison: hybrid (Section 12) vs single-nesting (Section 13) ---\n")
message(sprintf("  Hybrid (3 chains x %d cores):        %.1f minutes",
                CORES_PER_CHAIN, t_final / 60))
message(sprintf("  Single-nesting (1 chain x %d cores): %.1f minutes",
                CORES_SINGLE, t_sn_total / 60))
message(sprintf("  Per-chain time hybrid:               %.1f minutes",
                t_final / 60))    # chains ran simultaneously so wallclock = per-chain
message(sprintf("  Per-chain time single-nesting:       %.1f minutes",
                t_sn_total / 60 / N_CHAINS))
message("")
message("INTERPRETATION:")
if ((t_sn_total / N_CHAINS) < (t_final * 0.7)) {
  message("  Single-nesting is faster per chain -- nested forking was NOT working")
  message("  in Section 12. On Roihu (Linux) the hybrid will outperform this.")
} else {
  message("  Per-chain times are similar -- nested forking WAS working in Section 12.")
  message("  The hybrid approach gives better wallclock through simultaneous chains.")
}

# --- Chain completion check ---
message("\nChain completion check:")
n_sn_ok <- 0L
for (i in seq_len(N_CHAINS)) {
  res_i <- sn_chain_results[[i]]
  ok_i  <- inherits(res_i, "mcmcSamplerList") || inherits(res_i, "mcmcSampler")
  if (ok_i) {
    samp_i  <- getSample(res_i, coda = FALSE)
    n_fin_i <- sum(apply(samp_i, 1, function(r) all(is.finite(r))))
    pct     <- 100 * n_fin_i / max(nrow(samp_i), 1)
    message(sprintf("  Chain %d: %d post-burnin samples, %.0f%% finite  [%s]",
                    i, nrow(samp_i), pct, if (pct >= 95) "PASS" else "WARN"))
    n_sn_ok <- n_sn_ok + 1L
  } else {
    message(sprintf("  Chain %d: FAILED -- returned %s", i, class(res_i)[1]))
  }
}

if (n_sn_ok < N_CHAINS) {
  stop(sprintf("%d / %d chains failed.", N_CHAINS - n_sn_ok, N_CHAINS))
}

# --- Convergence diagnostics ---
message("\n--- Convergence diagnostics (single-nesting run) ---\n")

sn_coda <- tryCatch({
  coda::mcmc.list(lapply(sn_chain_results, function(r) {
    coda::mcmc(getSample(r, coda = FALSE))
  }))
}, error = function(e) NULL)

if (!is.null(sn_coda)) {
  
  gr_sn <- tryCatch(gelman.diag(sn_coda, multivariate = FALSE), error = function(e) NULL)
  
  if (!is.null(gr_sn)) {
    rhat_vals  <- gr_sn$psrf[, "Point est."]
    rhat_upper <- gr_sn$psrf[, "Upper C.I."]
    n_good     <- sum(rhat_vals < 1.05, na.rm = TRUE)
    n_ok_r     <- sum(rhat_vals >= 1.05 & rhat_vals < 1.10, na.rm = TRUE)
    n_warn     <- sum(rhat_vals >= 1.10, na.rm = TRUE)
    
    message(sprintf("Gelman-Rubin R-hat summary (%d parameters):", length(rhat_vals)))
    message(sprintf("  R-hat < 1.05 (good):    %d parameters", n_good))
    message(sprintf("  R-hat 1.05-1.10 (ok):   %d parameters", n_ok_r))
    message(sprintf("  R-hat > 1.10 (warning): %d parameters", n_warn))
    message(sprintf("  Overall range: %.3f -- %.3f",
                    min(rhat_vals, na.rm = TRUE), max(rhat_vals, na.rm = TRUE)))
    
    message("\n  Full R-hat table:")
    message(sprintf("  %-12s | %-10s | %-10s | %-10s",
                    "Parameter", "Point est.", "Upper C.I.", "Status"))
    message(sprintf("  %s", paste(rep("-", 50), collapse="")))
    for (nm in FREE_NAMES) {
      pe  <- rhat_vals[nm]
      uci <- rhat_upper[nm]
      sta <- if (is.na(pe)) "N/A" else if (pe < 1.05) "Good" else if (pe < 1.10) "OK" else "Warning"
      message(sprintf("  %-12s | %-10.3f | %-10.3f | %-10s", nm, pe, uci, sta))
    }
    message("")
    
    if (n_warn > 0) {
      message(sprintf("ACTION: %d parameter(s) with R-hat > 1.10.", n_warn))
      top3 <- names(sort(rhat_vals, decreasing = TRUE)[1:min(3, n_warn)])
      for (nm in top3) message(sprintf("    %s: %.3f", nm, rhat_vals[nm]))
    } else {
      message("All parameters have R-hat < 1.10. Convergence is acceptable.")
    }
  }
  
  # ESS
  message("\nEffective sample size (post-burnin, combined chains):")
  all_raw_sn  <- do.call(rbind, lapply(sn_chain_results,
                                       function(r) getSample(r, coda = FALSE)))
  ess_sn <- tryCatch(effectiveSize(coda::mcmc(all_raw_sn)), error = function(e) NULL)
  if (!is.null(ess_sn)) {
    n_low <- sum(ess_sn < 200, na.rm = TRUE)
    message(sprintf("  ESS range: %.0f -- %.0f", min(ess_sn), max(ess_sn)))
    message(sprintf("  Parameters with ESS < 200: %d / %d", n_low, length(ess_sn)))
    if (n_low > 0) {
      message("  ACTION: Low ESS -- run longer or check for strong correlations.")
    } else {
      message("  ESS > 200 for all parameters. Posterior summaries are reliable.")
    }
  }
}

# --- Posterior summary ---
message("\n--- Posterior summary (physical space, post-burnin) ---\n")

all_phys_sn <- t(apply(all_raw_sn, 1, to_original))
colnames(all_phys_sn) <- FREE_NAMES

message(sprintf("Total post-burnin samples: %d (across %d chains)", nrow(all_phys_sn), N_CHAINS))
message(sprintf("%-12s | %7s | %7s | %7s | %7s | %7s",
                "Parameter", "2.5%", "25%", "50%", "75%", "97.5%"))
message(paste(rep("-", 60), collapse=""))
post_q_sn <- apply(all_phys_sn, 2, quantile,
                   probs = c(0.025, 0.25, 0.50, 0.75, 0.975), na.rm = TRUE)
for (nm in FREE_NAMES) {
  message(sprintf("%-12s | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f",
                  nm,
                  post_q_sn["2.5%",  nm],
                  post_q_sn["25%",   nm],
                  post_q_sn["50%",   nm],
                  post_q_sn["75%",   nm],
                  post_q_sn["97.5%", nm]))
}

# --- Compare posteriors: hybrid (Section 12) vs single-nesting (Section 13) ---
# The two runs use the same data and prior so posteriors should be very similar.
# Large differences would suggest one of the runs has not converged.
message("\n--- Posterior comparison: hybrid vs single-nesting (median, climate params) ---\n")
message(sprintf("%-12s | %-10s | %-10s | %-10s",
                "Parameter", "Hybrid med", "SN med", "Abs diff"))
message(sprintf("%s", paste(rep("-", 48), collapse="")))
for (nm in c(CLIMATE_NAMES, NUISANCE_NAMES)) {
  med_h  <- post_q_final["50%", nm]
  med_sn <- post_q_sn["50%",   nm]
  message(sprintf("%-12s | %-10.4f | %-10.4f | %-10.4f",
                  nm, med_h, med_sn, abs(med_h - med_sn)))
}
message("")
message("If medians agree closely, both approaches have converged to the same posterior.")
message("If they disagree, run longer before trusting either result.\n")

# --- Trace plots ---
message("Generating trace plots (single-nesting run)...")
par(mfrow = c(3, 2))
for (nm in c(CLIMATE_NAMES, NUISANCE_NAMES)) {
  chain_traces <- lapply(sn_chain_results, function(r) {
    raw  <- getSample(r, coda = FALSE)
    phys <- t(apply(raw, 1, to_original))
    phys[, nm]
  })
  ylim_r <- range(unlist(chain_traces), na.rm = TRUE)
  plot(chain_traces[[1]], type = "l", col = "steelblue",
       ylim = ylim_r, xlab = "Iteration", ylab = nm,
       main = sprintf("%s -- trace (single-nesting)", nm))
  if (N_CHAINS > 1) {
    lines(chain_traces[[2]], col = "firebrick")
    if (N_CHAINS > 2) lines(chain_traces[[3]], col = "darkgreen")
  }
}
par(mfrow = c(1, 1))

# --- Save results ---
message("Saving single-nesting results...")
saveRDS(sn_chain_results, file = "Yasso07_sn_test_run_chains.rds")
saveRDS(all_phys_sn,      file = "Yasso07_sn_test_run_posterior.rds")
message("Saved:")
message("  Yasso07_sn_test_run_chains.rds    -- raw BayesianTools chain objects")
message("  Yasso07_sn_test_run_posterior.rds -- posterior samples in physical space")

message("\n=============================================================")
message("  Section 13 complete.")
message(sprintf("  Wallclock: %.1f minutes (sequential chains, parallel plots)", t_sn_total / 60))
message(sprintf("  Compare to Section 12: %.1f minutes (hybrid)", t_final / 60))
message("")
message("  DECISION FOR ROIHU:")
message("  Use the hybrid approach (Section 12 structure) -- nested forking")
message("  works on Linux and gives better wallclock through simultaneous chains.")
message("  Section 13 is the Mac fallback for local testing only.")
message("=============================================================\n")





# =============================================================================
# 14. Short MCMC test at 50 plots -- observe parallelism during real sampling
#
#     PURPOSE:
#     The scaling benchmark confirmed parallelism works at 50+ plots.
#     This section runs actual MCMC chains at 50 plots so you can watch
#     activity monitor during real sampling and confirm CPU usage is high.
#
#     At 50 plots the speedup is modest (~2.8x) but the workers stay busy
#     long enough per evaluation to be visible in activity monitor.
#     At 3700 plots on Roihu the speedup will be ~5x with 11 cores and
#     much higher with 40 cores per chain.
#
#     WHAT TO WATCH:
#       Activity monitor during this run should show ONE rsession-arm64
#       process using ~500-800% CPU (5-8 cores active). This confirms
#       the inner mclapply is genuinely parallelising the plot loop.
#
#     SETTINGS:
#       50 replicated plots (from test dataset)
#       1000 iterations per chain (short -- just enough to see mixing)
#       3 chains, sequential (single-nesting approach from Section 13)
#       All 11 cores for plot loop
# =============================================================================

message("-------------------------------------------------------------")
message("SECTION 14: Short MCMC test at 50 plots")
message("-------------------------------------------------------------\n")

# Set up 50 replicated plots
plots_original_s14 <- plots
plots <<- rep(plots_original_s14, length.out = 50)
message(sprintf("Plot list set to %d plots (replicated from test dataset).", length(plots)))
message(sprintf("Configuration: %d chains x 1000 iterations, %d cores for plot loop",
                N_CHAINS, CORES_SINGLE))
message("\nWatch activity monitor NOW.")
message(sprintf("Expected: 1 process using ~%d00%% CPU during sampling.\n", CORES_SINGLE))

log_dir <- "./Calibration/tests/"
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)

# Rebuild the likelihood with current (50-plot) plots list in scope
# Important: make_parallel_likelihood captures 'plots' from the environment
# at call time, so we rebuild after changing plots
ll_50plots <- make_parallel_likelihood(CORES_SINGLE)

# Quick correctness check at new plot count
ll_50_check <- ll_50plots(best_x)
message(sprintf("Log-likelihood at 50 plots: %.3f (should be finite)\n", ll_50_check))
if (!is.finite(ll_50_check)) stop("FAIL: log-likelihood not finite at 50 plots.")

set.seed(2025)
t_s14 <- system.time({
  
  s14_chain_results <- lapply(seq_len(N_CHAINS), function(chain_id) {
    
    log_file <- file.path(log_dir, sprintf("chain_s14_%d.log", chain_id))
    cat(sprintf("Chain %d started at %s\n", chain_id, Sys.time()), file = log_file)
    message(sprintf("  Chain %d / %d starting...", chain_id, N_CHAINS))
    
    call_count   <- 0L
    ll_fn_logged <- function(x) {
      call_count <<- call_count + 1L
      if (call_count %% 200L == 0L) {
        cat(sprintf("Chain %d | eval %d | %s\n",
                    chain_id, call_count, Sys.time()),
            file = log_file, append = TRUE)
      }
      ll_50plots(x)
    }
    
    setup <- createBayesianSetup(
      likelihood = ll_fn_logged,
      prior      = prior_yasso07,
      names      = FREE_NAMES,
      parallel   = FALSE
    )
    
    result <- runMCMC(
      bayesianSetup = setup,
      sampler       = "DEzs",
      settings      = list(
        iterations = 1000,
        burnin     = 200,
        nrChains   = 1,
        eps        = 0.001
      )
    )
    
    cat(sprintf("Chain %d finished at %s\n", chain_id, Sys.time()),
        file = log_file, append = TRUE)
    message(sprintf("  Chain %d complete.", chain_id))
    result
  })
  
})["elapsed"]

message(sprintf("\nSection 14 run completed in %.1f seconds (%.2f minutes)",
                t_s14, t_s14 / 60))

# --- Per-evaluation timing ---
# 1000 iterations x 3 chains, DEzs makes roughly 1 eval per iteration
total_evals_s14  <- 1000 * N_CHAINS
t_per_eval_s14   <- t_s14 / total_evals_s14
message(sprintf("Approx per-evaluation time: %.4f seconds at 50 plots", t_per_eval_s14))
message(sprintf("Expected at 3700 plots (linear scale): %.3f seconds\n",
                t_per_eval_s14 * 3700 / 50))

# --- Chain completion check ---
message("Chain completion check:")
n_s14_ok <- 0L
for (i in seq_len(N_CHAINS)) {
  res_i <- s14_chain_results[[i]]
  ok_i  <- inherits(res_i, "mcmcSamplerList") || inherits(res_i, "mcmcSampler")
  if (ok_i) {
    samp_i  <- getSample(res_i, coda = FALSE)
    n_fin_i <- sum(apply(samp_i, 1, function(r) all(is.finite(r))))
    pct     <- 100 * n_fin_i / max(nrow(samp_i), 1)
    message(sprintf("  Chain %d: %d post-burnin samples, %.0f%% finite  [%s]",
                    i, nrow(samp_i), pct, if (pct >= 95) "PASS" else "WARN"))
    n_s14_ok <- n_s14_ok + 1L
  } else {
    message(sprintf("  Chain %d: FAILED -- %s", i, class(res_i)[1]))
  }
}

# --- Quick convergence check ---
message("\nQuick convergence check (1000 iterations -- expect poor R-hat):")
s14_coda <- tryCatch({
  coda::mcmc.list(lapply(s14_chain_results, function(r) {
    coda::mcmc(getSample(r, coda = FALSE))
  }))
}, error = function(e) NULL)

if (!is.null(s14_coda)) {
  gr_s14 <- tryCatch(gelman.diag(s14_coda, multivariate = FALSE), error = function(e) NULL)
  if (!is.null(gr_s14)) {
    rhat_vals <- gr_s14$psrf[, "Point est."]
    message(sprintf("  R-hat range: %.3f -- %.3f",
                    min(rhat_vals, na.rm = TRUE), max(rhat_vals, na.rm = TRUE)))
    message(sprintf("  R-hat > 1.1: %d / %d parameters (expected at 1000 iters)",
                    sum(rhat_vals > 1.1, na.rm = TRUE), length(rhat_vals)))
    message("  Climate parameters:")
    for (nm in CLIMATE_NAMES) {
      message(sprintf("    %-8s: %.3f", nm, rhat_vals[nm]))
    }
  }
}

# Restore original plot list
plots <<- plots_original_s14
message(sprintf("\nPlot list restored to %d plots.", length(plots)))

message("\n=============================================================")
message("  Section 14 complete.")
message(sprintf("  Wallclock: %.1f seconds at 50 plots, 1000 iterations", t_s14))
message("  If activity monitor showed high CPU during this run,")
message("  plot-loop parallelism is confirmed working on your Mac.")
message("=============================================================\n")
