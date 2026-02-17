# =============================================================================
# Calibration_Yasso07_production.R
#
# PURPOSE:
#   Production calibration script for Yasso07 using the fully optimised
#   single-nesting approach:
#     - Chains run sequentially (lapply)
#     - Each chain uses all available cores for the plot loop (mclapply)
#
#   This is the Mac-compatible implementation. On Roihu (Linux, 384 cores
#   per node) switch to the hybrid approach:
#     N_CHAINS        <- 4L
#     CORES_PER_CHAIN <- 90L
#   and replace lapply with mclapply for the outer chain loop.
#
# OUTPUT STRUCTURE:
#   ./Calibration/progress_logs/   -- per-chain log files written during sampling
#   ./Calibration/diagnostics/     -- R-hat tables, ESS, posterior plots
#   ./Calibration/runs/            -- saved chain objects and posterior samples
#
#   All output files are prefixed with "Yasso07_" to distinguish from
#   other model runs (Yasso15, Yasso20, RothC, Q-model) in the same folders.
#
# TEST CONFIGURATION (current):
#   100 replicated plots -- enough to stress the parallelism visibly
#   10000 iterations per chain, 3 chains
#   All available cores for plot loop
#
# ROIHU CONFIGURATION (April 2026):
#   3700 real plots
#   50000 iterations per chain, 4 chains
#   N_CHAINS=4, CORES_PER_CHAIN=90 (hybrid approach, 360 cores)
#
# =============================================================================

library(dplyr)
library(BayesianTools)
library(parallel)
library(compiler)
library(coda)

source("./Model_functions/input_compatibility_layer.R")
source("./Model_functions/Decomposition_functions/Yasso/yasso07_wrapper.R")
dyn.load("./Model_functions/Decomposition_functions/Yasso/yasso07.so")

# =============================================================================
# 0.  Configuration
#
#     Change these settings when moving to Roihu.
#     Everything else in the script adapts automatically.
# =============================================================================

# --- Run identity ---
MODEL_NAME  <- "Yasso07"          # used as prefix for all output files
RUN_ID      <- format(Sys.time(), "%Y%m%d_%H%M%S")  # unique timestamp per run

# --- Plot replication (testing only -- remove for real data) ---
# Set N_PLOTS_TEST to NA to use the real unreplicated dataset
N_PLOTS_TEST <- 100L              # replicate to this many plots for Mac testing

# --- MCMC settings ---
N_CHAINS   <- 3L                  # number of chains
N_ITER     <- 10000L              # iterations per chain
N_BURNIN   <- 2000L              # burn-in per chain (discarded)
N_LOG      <- 200L               # log progress every N_LOG evaluations

# --- Parallelism ---
# Mac (single-nesting): chains sequential, all cores to plot loop
# Roihu (hybrid):       replace lapply with mclapply, set CORES_PER_CHAIN=90
CORES_PER_CHAIN <- parallel::detectCores() - 1L   # all cores to plot loop on Mac
# CORES_PER_CHAIN <- 90L                           # Roihu setting

# --- Output directories ---
DIR_LOGS  <- "./Calibration/progress_logs/"
DIR_DIAG  <- "./Calibration/diagnostics/"
DIR_RUNS  <- "./Calibration/runs/"

# Create output directories if they don't exist
for (d in c(DIR_LOGS, DIR_DIAG, DIR_RUNS)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

message("=============================================================")
message(sprintf("  %s Calibration -- Production Run", MODEL_NAME))
message(sprintf("  Run ID: %s", RUN_ID))
message(sprintf("  %s", Sys.time()))
message("=============================================================\n")
message(sprintf("Output directories:"))
message(sprintf("  Logs:        %s", DIR_LOGS))
message(sprintf("  Diagnostics: %s", DIR_DIAG))
message(sprintf("  Runs:        %s\n", DIR_RUNS))


# =============================================================================
# 1.  Load and prepare data
# =============================================================================

message("Loading and preparing data...")

input_raw       <- read.csv("synthetic_input_data_template.csv")
site_raw        <- read.csv("synthetic_site_data_template.csv")
Yasso07_climate <- map_climate_yasso07(input_raw)
Yasso07_inputs  <- map_inputs_yasso07(input_raw)
plots_real      <- unique(Yasso07_climate$plot_id)

# Replicate plots for Mac testing if requested
if (!is.na(N_PLOTS_TEST)) {
  plots <- rep(plots_real, length.out = N_PLOTS_TEST)
  message(sprintf("TEST MODE: %d real plots replicated to %d plots.",
                  length(plots_real), N_PLOTS_TEST))
} else {
  plots <- plots_real
  message(sprintf("PRODUCTION MODE: %d real plots.", length(plots)))
}

# SOC observations
SOC_obs_all <- input_raw %>%
  filter(!is.na(soc_obs_tCha)) %>%
  select(plot_id, year, month, soc_obs_tCha) %>%
  arrange(plot_id, year, month) %>%
  group_by(plot_id) %>%
  mutate(obs_rank = row_number()) %>%
  ungroup()

message(sprintf("SOC observations: %d across %d unique plots.\n",
                nrow(SOC_obs_all), n_distinct(SOC_obs_all$plot_id)))

# Pre-split data frames for O(1) lookup inside likelihood
# (avoids O(n) data frame scanning on every MCMC evaluation)
climate_by_plot <- split(Yasso07_climate, Yasso07_climate$plot_id)
inputs_by_plot  <- split(Yasso07_inputs,  Yasso07_inputs$plot_id)

# Pre-compute plot-level quantities (done once, outside MCMC loop)
litter_means <- lapply(plots_real, function(pid) {
  inputs <- Yasso07_inputs[Yasso07_inputs$plot_id == pid, ]
  list(
    nwl_mean = colMeans(inputs[, c("nwl_A","nwl_W","nwl_E","nwl_N")]),
    fwl_mean = colMeans(inputs[, c("fwl_A","fwl_W","fwl_E","fwl_N")]),
    cwl_mean = colMeans(inputs[, c("cwl_A","cwl_W","cwl_E","cwl_N")])
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


# =============================================================================
# 2.  Parameter structure and transformations
# =============================================================================

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
stopifnot(max(abs(to_original(best_x)[FREE_NAMES] - free_defaults[FREE_NAMES])) < 1e-10)
message("Transformation round-trip check: OK")


# =============================================================================
# 3.  Prior predictive matching
# =============================================================================

message("\nRunning prior predictive matching...")
obs_cv <- sd(SOC_obs_all$soc_obs_tCha) / mean(SOC_obs_all$soc_obs_tCha)

run_prior_sample_soc <- function(sigma_val) {
  x_sample     <- best_x + rnorm(N_FREE, 0, sigma_val)
  p_free       <- to_original(x_sample)
  model_params <- assemble_model_params(p_free)
  pid          <- plots_real[1]
  clim         <- climate_by_plot[[pid]]
  inputs       <- inputs_by_plot[[pid]]
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
ppm_cv <- sapply(seq_along(sigma_grid), function(j) {
  sv <- sigma_grid[j]
  cat(sprintf("\r  Prior predictive matching: sigma = %.1f (%d/%d)...",
              sv, j, length(sigma_grid)))
  flush.console()
  preds <- vapply(seq_len(200L), function(i) run_prior_sample_soc(sv), numeric(1))
  valid <- preds[is.finite(preds) & preds > 0]
  if (length(valid) < 20L) NA_real_ else sd(valid) / mean(valid)
})
cat("\n")  # move to new line after the \r progress bar

sigma_ppm <- sigma_grid[which(ppm_cv >= obs_cv)[1]]
if (is.na(sigma_ppm)) { sigma_ppm <- 1.0 }
message(sprintf("Prior width: sigma_ppm = %.1f  (obs CV = %.3f)\n", sigma_ppm, obs_cv))


# =============================================================================
# 4.  Likelihood function (parallel plot loop, byte-compiled)
#
#     make_likelihood() returns a compiled likelihood function with a fixed
#     core count baked in via closure. Using a factory means each chain gets
#     its own independent function with no shared state.
#
#     The Jacobian correction is included here (not in the prior) to keep
#     the prior object clean. See performance tests document for explanation.
# =============================================================================

make_likelihood <- function(n_cores) {
  # Capture n_cores explicitly before compiling -- cmpfun can break closures
  nc <- n_cores
  cmpfun(function(x) {
    
    p_free       <- to_original(x)
    sigma_obs    <- p_free["sigma_obs"]
    sigma_init   <- p_free["sigma_init"]
    model_params <- assemble_model_params(p_free)
    
    # Jacobian correction for stick-breaking and log transforms
    log_jac <- log_jac_stick(p_free[FRAC_COL_A], B) +
      log_jac_stick(p_free[FRAC_COL_W], B) +
      log_jac_stick(p_free[FRAC_COL_E], B) +
      log_jac_stick(p_free[FRAC_COL_N], B) +
      log(sigma_obs) + log(sigma_init)
    
    # Parallel plot loop: each worker evaluates one plot's log-likelihood
    log_liks <- parallel::mclapply(plots, function(pid) {
      
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
      if (is.null(xi_array))                   return(-Inf)
      xi_mean <- mean(xi_array)
      if (!is.finite(xi_mean) || xi_mean <= 0) return(-Inf)
      
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
      
      sum(dnorm(meta$soc_obs, mean = SOC_hat, sd = sd_vec, log = TRUE))
      
    }, mc.cores = nc)
    
    log_liks <- unlist(log_liks)
    if (any(!is.finite(log_liks))) return(-Inf)
    sum(log_liks) + log_jac
  })
}

# Verify likelihood at defaults
ll_check <- make_likelihood(CORES_PER_CHAIN)(best_x)
message(sprintf("Likelihood check at defaults: %.3f", ll_check))
if (!is.finite(ll_check)) stop("FAIL: likelihood not finite at defaults.")
message("Likelihood check: PASS\n")


# =============================================================================
# 5.  Prior
# =============================================================================

prior_yasso07 <- createPrior(
  density = function(x) sum(dnorm(x, mean = best_x, sd = sigma_ppm, log = TRUE)),
  sampler = function(n = 1) matrix(rnorm(n * N_FREE, mean = best_x, sd = sigma_ppm * 0.1),
                                   nrow = n, ncol = N_FREE),
  lower = rep(-Inf, N_FREE),
  upper = rep( Inf, N_FREE)
)


# =============================================================================
# 6.  Run MCMC chains
#
#     Single-nesting approach:
#       - Outer lapply:   chains run sequentially (Mac compatible)
#       - Inner mclapply: plot loop parallelised with CORES_PER_CHAIN cores
#
#     For Roihu, replace lapply with mclapply(mc.cores = N_CHAINS) and
#     set CORES_PER_CHAIN = 90. Everything else is identical.
#
#     Each chain writes progress to its own log file in DIR_LOGS.
#     Log files are named: Yasso07_chain_01_RUNID.log
# =============================================================================

message("-------------------------------------------------------------")
message(sprintf("Starting MCMC: %d chains x %d iterations", N_CHAINS, N_ITER))
message(sprintf("  Plot loop: %d cores per chain", CORES_PER_CHAIN))
message(sprintf("  Burn-in:   %d iterations per chain", N_BURNIN))
message(sprintf("  Logging:   every %d evaluations", N_LOG))
message(sprintf("  Log files: %s%s_chain_XX_%s.log", DIR_LOGS, MODEL_NAME, RUN_ID))
message("-------------------------------------------------------------\n")

# Rough runtime estimate based on benchmark timings
t_est_eval  <- 0.052   # seconds per eval at 100 plots (from Section 9 benchmark)
t_est_chain <- t_est_eval * N_ITER / 60
message(sprintf("Estimated wallclock: ~%.0f minutes (sequential chains at %d plots)\n",
                t_est_chain * N_CHAINS, length(plots)))

set.seed(2025)
t_run <- system.time({
  
  # lapply: chains run one at a time
  # On Roihu: replace with mclapply(mc.cores = N_CHAINS)
  chain_results <- lapply(seq_len(N_CHAINS), function(chain_id) {
    
    # Per-chain log file -- named with model, chain number, and run ID
    log_file <- file.path(
      DIR_LOGS,
      sprintf("%s_chain_%02d_%s.log", MODEL_NAME, chain_id, RUN_ID)
    )
    
    cat(sprintf("============================================================\n"), file = log_file)
    cat(sprintf("  %s Calibration -- Chain %d of %d\n", MODEL_NAME, chain_id, N_CHAINS), file = log_file, append = TRUE)
    cat(sprintf("  Run ID:    %s\n", RUN_ID), file = log_file, append = TRUE)
    cat(sprintf("  Started:   %s\n", Sys.time()), file = log_file, append = TRUE)
    cat(sprintf("  Plots:     %d\n", length(plots)), file = log_file, append = TRUE)
    cat(sprintf("  Iterations:%d\n", N_ITER), file = log_file, append = TRUE)
    cat(sprintf("  Cores:     %d\n", CORES_PER_CHAIN), file = log_file, append = TRUE)
    cat(sprintf("============================================================\n\n"), file = log_file, append = TRUE)
    
    message(sprintf("  Chain %d / %d starting  -->  %s",
                    chain_id, N_CHAINS, log_file))
    
    # Build likelihood for this chain
    ll_fn <- make_likelihood(CORES_PER_CHAIN)
    
    # Wrap with logging: writes to log file every N_LOG evaluations
    call_count   <- 0L
    t_chain_start <- proc.time()[["elapsed"]]
    
    ll_fn_logged <- function(x) {
      call_count <<- call_count + 1L
      result <- ll_fn(x)
      if (call_count %% N_LOG == 0L) {
        elapsed <- proc.time()[["elapsed"]] - t_chain_start
        rate    <- call_count / elapsed
        cat(sprintf("[%s]  eval %5d  |  log-lik: %10.3f  |  %.1f eval/s  |  elapsed: %.1f min\n",
                    format(Sys.time(), "%H:%M:%S"),
                    call_count, result, rate, elapsed / 60),
            file = log_file, append = TRUE)
      }
      result
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
        iterations = N_ITER,
        burnin     = N_BURNIN,
        nrChains   = 1,
        eps        = 0.001
      )
    )
    
    elapsed_chain <- proc.time()[["elapsed"]] - t_chain_start
    cat(sprintf("\nChain %d finished at %s  (%.1f minutes)\n",
                chain_id, Sys.time(), elapsed_chain / 60),
        file = log_file, append = TRUE)
    
    message(sprintf("  Chain %d / %d complete  (%.1f minutes)",
                    chain_id, N_CHAINS, elapsed_chain / 60))
    result
  })
  
})["elapsed"]

message(sprintf("\nAll chains complete. Total wallclock: %.1f minutes\n", t_run / 60))


# =============================================================================
# 7.  Convergence diagnostics
#
#     Results written to DIR_DIAG, prefixed with MODEL_NAME and RUN_ID.
#     Files produced:
#       Yasso07_rhat_RUNID.txt          -- R-hat table for all parameters
#       Yasso07_ess_RUNID.txt           -- ESS for all parameters
#       Yasso07_posterior_summary_RUNID.txt  -- posterior quantile table
#       Yasso07_traces_RUNID.pdf        -- trace plots (climate + nuisance)
#       Yasso07_marginals_RUNID.pdf     -- posterior vs prior marginal plots
# =============================================================================

message("-------------------------------------------------------------")
message("Computing convergence diagnostics...")
message(sprintf("  Output: %s%s_*_%s.*", DIR_DIAG, MODEL_NAME, RUN_ID))
message("-------------------------------------------------------------\n")

# Check chain completion
n_ok <- 0L
for (i in seq_len(N_CHAINS)) {
  res_i <- chain_results[[i]]
  ok_i  <- inherits(res_i, "mcmcSamplerList") || inherits(res_i, "mcmcSampler")
  if (ok_i) {
    samp_i  <- getSample(res_i, coda = FALSE)
    n_fin_i <- sum(apply(samp_i, 1, function(r) all(is.finite(r))))
    pct     <- 100 * n_fin_i / max(nrow(samp_i), 1)
    message(sprintf("  Chain %d: %d post-burnin samples, %.0f%% finite  [%s]",
                    i, nrow(samp_i), pct, if (pct >= 95) "PASS" else "WARN"))
    n_ok <- n_ok + 1L
  } else {
    message(sprintf("  Chain %d: FAILED", i))
  }
}
if (n_ok < N_CHAINS) stop(sprintf("%d / %d chains failed.", N_CHAINS - n_ok, N_CHAINS))

# Combine chains
all_raw  <- do.call(rbind, lapply(chain_results, function(r) getSample(r, coda = FALSE)))
all_phys <- t(apply(all_raw, 1, to_original))
colnames(all_phys) <- FREE_NAMES

chain_coda <- coda::mcmc.list(lapply(chain_results, function(r) {
  coda::mcmc(getSample(r, coda = FALSE))
}))

# --- R-hat ---
message("Computing Gelman-Rubin R-hat...")
gr <- tryCatch(gelman.diag(chain_coda, multivariate = FALSE), error = function(e) NULL)

rhat_file <- file.path(DIR_DIAG, sprintf("%s_rhat_%s.txt", MODEL_NAME, RUN_ID))
sink(rhat_file)
cat(sprintf("%s Calibration -- Gelman-Rubin R-hat\n", MODEL_NAME))
cat(sprintf("Run ID: %s\n", RUN_ID))
cat(sprintf("Date:   %s\n", Sys.time()))
cat(sprintf("Plots:  %d | Chains: %d | Iterations: %d | Burnin: %d\n\n",
            length(plots), N_CHAINS, N_ITER, N_BURNIN))

if (!is.null(gr)) {
  rhat_vals  <- gr$psrf[, "Point est."]
  rhat_upper <- gr$psrf[, "Upper C.I."]
  n_good <- sum(rhat_vals < 1.05, na.rm = TRUE)
  n_okr  <- sum(rhat_vals >= 1.05 & rhat_vals < 1.10, na.rm = TRUE)
  n_warn <- sum(rhat_vals >= 1.10, na.rm = TRUE)
  cat(sprintf("Summary: %d good (<1.05) | %d ok (1.05-1.10) | %d warning (>1.10)\n\n",
              n_good, n_okr, n_warn))
  cat(sprintf("%-12s  %-10s  %-10s  %-10s\n", "Parameter", "Point est.", "Upper C.I.", "Status"))
  cat(paste(rep("-", 48), collapse=""), "\n")
  for (nm in FREE_NAMES) {
    pe  <- rhat_vals[nm]
    uci <- rhat_upper[nm]
    sta <- if (is.na(pe)) "N/A" else if (pe < 1.05) "Good" else if (pe < 1.10) "OK" else "Warning"
    cat(sprintf("%-12s  %-10.4f  %-10.4f  %-10s\n", nm, pe, uci, sta))
  }
  # Print to console summary
  cat(sprintf("\nMultivariate psrf: %.4f\n", gr$mpsrf))
}
sink()

# Console summary
if (!is.null(gr)) {
  message(sprintf("R-hat range: %.3f -- %.3f  (%d warnings)",
                  min(gr$psrf[,1], na.rm=TRUE), max(gr$psrf[,1], na.rm=TRUE), n_warn))
  message(sprintf("R-hat file: %s", rhat_file))
  if (n_warn > 0) {
    message(sprintf("WARNING: %d parameter(s) with R-hat > 1.10 -- consider longer run", n_warn))
  }
}

# --- ESS ---
message("\nComputing effective sample size...")
ess_vals <- tryCatch(effectiveSize(coda::mcmc(all_raw)), error = function(e) NULL)

ess_file <- file.path(DIR_DIAG, sprintf("%s_ess_%s.txt", MODEL_NAME, RUN_ID))
sink(ess_file)
cat(sprintf("%s Calibration -- Effective Sample Size\n", MODEL_NAME))
cat(sprintf("Run ID: %s\n", RUN_ID))
cat(sprintf("Date:   %s\n\n", Sys.time()))
if (!is.null(ess_vals)) {
  cat(sprintf("Total post-burnin samples: %d\n", nrow(all_raw)))
  cat(sprintf("ESS range: %.0f -- %.0f\n", min(ess_vals), max(ess_vals)))
  cat(sprintf("Parameters with ESS < 200: %d / %d\n\n",
              sum(ess_vals < 200), length(ess_vals)))
  cat(sprintf("%-12s  %-10s\n", "Parameter", "ESS"))
  cat(paste(rep("-", 25), collapse=""), "\n")
  for (nm in FREE_NAMES) {
    cat(sprintf("%-12s  %-10.0f\n", nm, ess_vals[nm]))
  }
}
sink()

if (!is.null(ess_vals)) {
  n_low <- sum(ess_vals < 200, na.rm = TRUE)
  message(sprintf("ESS range: %.0f -- %.0f  (%d parameters < 200)",
                  min(ess_vals), max(ess_vals), n_low))
  message(sprintf("ESS file: %s", ess_file))
  if (n_low > 0) message("WARNING: Low ESS -- consider running longer")
}

# --- Posterior summary ---
message("\nComputing posterior summary...")
post_q <- apply(all_phys, 2, quantile, probs = c(0.025, 0.25, 0.50, 0.75, 0.975), na.rm = TRUE)

post_file <- file.path(DIR_DIAG, sprintf("%s_posterior_summary_%s.txt", MODEL_NAME, RUN_ID))
sink(post_file)
cat(sprintf("%s Calibration -- Posterior Summary (physical space)\n", MODEL_NAME))
cat(sprintf("Run ID: %s\n", RUN_ID))
cat(sprintf("Date:   %s\n", Sys.time()))
cat(sprintf("Total post-burnin samples: %d\n\n", nrow(all_phys)))
cat(sprintf("%-12s  %8s  %8s  %8s  %8s  %8s\n",
            "Parameter", "2.5%", "25%", "50%", "75%", "97.5%"))
cat(paste(rep("-", 60), collapse=""), "\n")
for (nm in FREE_NAMES) {
  cat(sprintf("%-12s  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f\n",
              nm,
              post_q["2.5%",  nm], post_q["25%",   nm], post_q["50%",   nm],
              post_q["75%",   nm], post_q["97.5%", nm]))
}
sink()
message(sprintf("Posterior summary file: %s", post_file))

# --- Trace plots (PDF) ---
message("\nGenerating trace plots...")
trace_file <- file.path(DIR_DIAG, sprintf("%s_traces_%s.pdf", MODEL_NAME, RUN_ID))
pdf(trace_file, width = 10, height = 8)
par(mfrow = c(3, 2), mar = c(4, 4, 3, 1))
for (nm in c(CLIMATE_NAMES, NUISANCE_NAMES)) {
  chain_traces <- lapply(chain_results, function(r) {
    raw  <- getSample(r, coda = FALSE)
    phys <- t(apply(raw, 1, to_original))
    phys[, nm]
  })
  ylim_r <- range(unlist(chain_traces), na.rm = TRUE)
  cols   <- c("steelblue", "firebrick", "darkgreen", "darkorange")[seq_len(N_CHAINS)]
  plot(chain_traces[[1]], type = "l", col = cols[1],
       ylim = ylim_r, xlab = "Post-burnin iteration", ylab = nm,
       main = sprintf("%s | %s | %s", MODEL_NAME, nm, RUN_ID))
  for (ci in seq_len(N_CHAINS)[-1]) lines(chain_traces[[ci]], col = cols[ci])
  legend("topright", legend = sprintf("Chain %d", seq_len(N_CHAINS)),
         col = cols, lwd = 1, bty = "n", cex = 0.8)
}
dev.off()
message(sprintf("Trace plots: %s", trace_file))

# --- Marginal posterior vs prior plots (PDF) ---
message("Generating marginal plots...")
marg_file <- file.path(DIR_DIAG, sprintf("%s_marginals_%s.pdf", MODEL_NAME, RUN_ID))
set.seed(99)
prior_raw  <- matrix(rnorm(3000 * N_FREE, mean = best_x, sd = sigma_ppm),
                     nrow = 3000, ncol = N_FREE)
prior_phys <- t(apply(prior_raw, 1, to_original))

pdf(marg_file, width = 10, height = 8)
par(mfrow = c(3, 2), mar = c(4, 4, 3, 1))
for (nm in c(CLIMATE_NAMES, NUISANCE_NAMES)) {
  post_v  <- all_phys[is.finite(all_phys[, nm]), nm]
  prior_v <- prior_phys[is.finite(prior_phys[, nm]), nm]
  d_post  <- density(post_v)
  d_prior <- density(prior_v)
  xlim_r  <- range(c(d_post$x, d_prior$x))
  ylim_r  <- range(c(d_post$y, d_prior$y))
  plot(d_post, col = "steelblue", lwd = 2, xlim = xlim_r, ylim = ylim_r,
       main = sprintf("%s | %s", MODEL_NAME, nm),
       xlab = nm, ylab = "Density")
  lines(d_prior, col = "grey60", lwd = 1, lty = 2)
  abline(v = free_defaults[nm], col = "firebrick", lwd = 1, lty = 3)
  legend("topright",
         legend = c("Posterior", "Prior", "Published default"),
         col = c("steelblue", "grey60", "firebrick"),
         lwd = c(2, 1, 1), lty = c(1, 2, 3), bty = "n", cex = 0.8)
}
dev.off()
message(sprintf("Marginal plots: %s", marg_file))


# =============================================================================
# 8.  Save results
#
#     Files saved to DIR_RUNS:
#       Yasso07_chains_RUNID.rds      -- raw BayesianTools result objects
#       Yasso07_posterior_RUNID.rds   -- posterior samples in physical space
#       Yasso07_metadata_RUNID.rds    -- run configuration and timing
# =============================================================================

message("\n-------------------------------------------------------------")
message("Saving results...")
message(sprintf("  Output: %s%s_*_%s.rds", DIR_RUNS, MODEL_NAME, RUN_ID))
message("-------------------------------------------------------------\n")

# Metadata: everything needed to reproduce or extend this run
metadata <- list(
  model        = MODEL_NAME,
  run_id       = RUN_ID,
  timestamp    = Sys.time(),
  n_plots      = length(plots),
  n_plots_real = length(plots_real),
  n_plots_test = N_PLOTS_TEST,
  n_chains     = N_CHAINS,
  n_iter       = N_ITER,
  n_burnin     = N_BURNIN,
  cores_per_chain = CORES_PER_CHAIN,
  sigma_ppm    = sigma_ppm,
  wallclock_min = t_run / 60,
  r_version    = R.version$version.string,
  platform     = .Platform$OS.type,
  rhat_range   = if (!is.null(gr)) range(gr$psrf[,1], na.rm=TRUE) else NA,
  ess_range    = if (!is.null(ess_vals)) range(ess_vals, na.rm=TRUE) else NA,
  n_samples    = nrow(all_phys)
)

saveRDS(chain_results, file.path(DIR_RUNS, sprintf("%s_chains_%s.rds",    MODEL_NAME, RUN_ID)))
saveRDS(all_phys,      file.path(DIR_RUNS, sprintf("%s_posterior_%s.rds", MODEL_NAME, RUN_ID)))
saveRDS(metadata,      file.path(DIR_RUNS, sprintf("%s_metadata_%s.rds",  MODEL_NAME, RUN_ID)))

message(sprintf("Saved: %s_chains_%s.rds",    MODEL_NAME, RUN_ID))
message(sprintf("Saved: %s_posterior_%s.rds", MODEL_NAME, RUN_ID))
message(sprintf("Saved: %s_metadata_%s.rds",  MODEL_NAME, RUN_ID))


# =============================================================================
# 9.  Final summary
# =============================================================================

message("\n=============================================================")
message(sprintf("  %s Calibration Complete", MODEL_NAME))
message(sprintf("  Run ID:    %s", RUN_ID))
message(sprintf("  Wallclock: %.1f minutes", t_run / 60))
message(sprintf("  Samples:   %d post-burnin (across %d chains)", nrow(all_phys), N_CHAINS))
if (!is.null(gr)) {
  message(sprintf("  R-hat:     %.3f -- %.3f  (%s)",
                  min(gr$psrf[,1], na.rm=TRUE), max(gr$psrf[,1], na.rm=TRUE),
                  if (n_warn == 0) "all good" else sprintf("%d warnings", n_warn)))
}
if (!is.null(ess_vals)) {
  message(sprintf("  ESS:       %.0f -- %.0f  (%s)",
                  min(ess_vals), max(ess_vals),
                  if (sum(ess_vals < 200) == 0) "all > 200" else sprintf("%d < 200", sum(ess_vals < 200))))
}
message("")
message("  To switch to Roihu production mode:")
message("    1. Set N_PLOTS_TEST <- NA")
message("    2. Set N_ITER       <- 50000L")
message("    3. Set N_CHAINS     <- 4L")
message("    4. Set CORES_PER_CHAIN <- 90L")
message("    5. Replace lapply with mclapply(mc.cores = N_CHAINS)")
message("       in Section 6 (nested forking works on Linux)")
message("=============================================================\n")