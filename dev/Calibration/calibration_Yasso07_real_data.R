# =============================================================================
# Calibration_Yasso07_real_data.R
#
# Model-specific calibration script for Yasso07 using real Finnish NFI data.
# Sources calibration_engine.R for model-agnostic MCMC machinery.
#
# KEY FEATURES:
#   - Reads input_raw_monthly.csv and site_raw.csv (from data_preparation.R)
#   - Filters to calibration-ready plots (has_climate & has_soc)
#   - sigma_input: global multiplicative litter scaling (free parameter)
#   - STEADY_STATE_YEARS: configurable initialization period (default 20)
#   - plot_info lookup carries stratification metadata into the likelihood
#   - assemble_model_params designed for future stratification extension
#   - Posterior predictive output: 100 draws, long format + summary
#   - sigma_ppm is a PER-PARAMETER vector (not scalar)
#
# REQUIRED ENGINE MODIFICATION: steady_state_n parameter in make_likelihood()
# (already applied in the latest calibration_engine.R)
#
# MAC TEST:  N_PLOTS_TEST <- 20L, N_ITER <- 10000L, 3 chains
# ROIHU:    N_PLOTS_TEST <- NA,   N_ITER <- 50000L, 4 chains, CORES_PER_CHAIN <- 90L
# =============================================================================

source("./Calibration_real_data/calibration_engine.R")
source("./Model_functions/input_compatibility_layer.R")
source("./Model_functions/Decomposition_functions/Yasso/yasso07_wrapper.R")
dyn.load("./Model_functions/Decomposition_functions/Yasso/yasso07.so")

library(dplyr)


# =============================================================================
# 0.  Configuration
# =============================================================================

MODEL_NAME <- "Yasso07"
RUN_ID     <- format(Sys.time(), "%Y%m%d_%H%M%S")

N_PLOTS_TEST <- 20L

N_CHAINS <- 3L
N_ITER   <- 10000L
N_BURNIN <- 2000L
N_LOG    <- 200L

CORES_PER_CHAIN <- parallel::detectCores() - 1L
# CORES_PER_CHAIN <- 90L   # Roihu

STEADY_STATE_YEARS <- 20L
N_PP_DRAWS         <- 100L

DIR_LOGS   <- "./Calibration_real_data/progress_logs/"
DIR_DIAG   <- "./Calibration_real_data/diagnostics/"
DIR_RUNS   <- "./Calibration_real_data/runs/"
DIR_INPUTS <- "./Calibration/model_inputs/"
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
message(sprintf("  %s Calibration (real data)  |  Run: %s", MODEL_NAME, RUN_ID))
message("=============================================================\n")


# =============================================================================
# 1.  Parameter specification
# =============================================================================

p_default        <- YASSO07_DEFAULT_PARAMS
FIXED_RATE_NAMES <- c("alpha_A","alpha_W","alpha_E","alpha_N","p_H","alpha_H")
fixed_rates      <- p_default[FIXED_RATE_NAMES]

B <- 1.0 - fixed_rates["p_H"]

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

transforms       <- build_transforms(param_spec)
to_original      <- transforms$to_original
to_unconstrained <- transforms$to_unconstrained
log_jacobian     <- transforms$log_jacobian
FREE_NAMES       <- transforms$param_names
N_FREE           <- transforms$n_params

free_defaults <- c(
  p_default[c(unlist(lapply(param_spec[1:4], `[[`, "names")))],
  p_default[c("beta1","beta2","gamma","delta1","delta2","r")],
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
# FUTURE STRATIFICATION: single-point extension would be
#   assemble_model_params(p_free, stratum = NULL)

assemble_model_params <- function(p_free) {
  yasso_params <- c(
    fixed_rates,
    p_free[c(unlist(lapply(param_spec[1:6], `[[`, "names")))]
  )
  yasso_params <- yasso_params[names(p_default)]
  # unname() strips the inner name so outer c() naming works
  c(yasso_params, sigma_input = unname(p_free["sigma_input"]))
}


# =============================================================================
# 3.  Data preparation
# =============================================================================

message("Loading real data...")

input_raw <- read.csv("./Data/model_inputs/input_raw_monthly.csv")
site_raw  <- read.csv("./Data/model_inputs/site_raw.csv")


# Annual litter total per plot (all fractions summed)
litter_cols <- grep("^C_nwl|^C_fwl|^C_cwl", names(input_raw), value = TRUE)

annual_litter <- input_raw %>%
  group_by(plot_id, year) %>%
  summarise(total_litter = sum(across(all_of(litter_cols)), na.rm = TRUE),
            .groups = "drop")

summary(annual_litter$total_litter)


# Filter to calibration-ready plots
calib_plots <- site_raw$plot_id[site_raw$calib_ready]
message(sprintf("Calibration-ready plots: %d", length(calib_plots)))

input_calib <- input_raw[input_raw$plot_id %in% calib_plots, ]

# Drop any plots with climate NAs
n_climate_na <- sum(is.na(input_calib$temp_air))
if (n_climate_na > 0) {
  message(sprintf("WARNING: %d rows with NA climate in calib plots. Filtering...",
                  n_climate_na))
  plots_with_na <- unique(input_calib$plot_id[is.na(input_calib$temp_air)])
  calib_plots   <- setdiff(calib_plots, plots_with_na)
  input_calib   <- input_raw[input_raw$plot_id %in% calib_plots, ]
  message(sprintf("After climate filter: %d plots", length(calib_plots)))
}

# Map through compatibility layer (monthly → annual per model)
Yasso07_climate <- map_climate_yasso07(input_calib)
Yasso07_inputs  <- map_inputs_yasso07(input_calib)

# Force character plot IDs for consistent keying
plots_real <- as.character(unique(Yasso07_climate$plot_id))
message(sprintf("Plots after mapping: %d", length(plots_real)))

# Subsample for Mac testing
if (!is.na(N_PLOTS_TEST) && N_PLOTS_TEST < length(plots_real)) {
  set.seed(42)
  plots <- sort(sample(plots_real, N_PLOTS_TEST))
} else {
  plots <- plots_real
}
plots <- as.character(plots)

# SOC observations
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

# Pre-split data frames for O(1) lookup (split() keys become character)
climate_by_plot <- split(Yasso07_climate, Yasso07_climate$plot_id)
inputs_by_plot  <- split(Yasso07_inputs,  Yasso07_inputs$plot_id)

# Litter means (first STEADY_STATE_YEARS only, for steady-state init)
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

# Observation metadata (FIXED: filter clim on Yasso07_climate, not Yasso07_inputs)
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

# Plot-level metadata for stratification (not consumed by engine yet)
plot_info <- lapply(plots_real, function(pid) {
  row <- site_raw[as.character(site_raw$plot_id) == pid, ]
  list(
    region       = row$region,
    species_code = row$species_code,
    species_name = row$species_name,
    soil_code    = row$soil_code,
    temp_zone    = row$temp_zone,
    shallow      = row$shallow,
    lat          = row$lat_WGS84,
    mean_temp    = row$mean_temp
  )
})
names(plot_info) <- plots_real

run_config$n_plots <- length(plots)


# =============================================================================
# 4.  Model interface wrappers
#
#     sigma_input scales litter inputs in both steady_state and transient run.
#     compute_xi is unaffected (climate modifier does not depend on litter).
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

steady_state_yasso07_engine <- function(model_params, lm, xi_mean) {
  si <- model_params["sigma_input"]                            # scaling factor
  yasso07_steady_state(
    params   = model_params[YASSO07_PARAM_NAMES],
    nwl_mean = lm$nwl_mean * si,                                # scale inputs
    fwl_mean = lm$fwl_mean * si,
    cwl_mean = lm$cwl_mean * si,
    xi_mean  = xi_mean
  )
}

LITTER_COLS <- c("nwl_A","nwl_W","nwl_E","nwl_N",
                 "fwl_A","fwl_W","fwl_E","fwl_N",
                 "cwl_A","cwl_W","cwl_E","cwl_N")

yasso07_run_engine <- function(inputs, model_params, C_init, xi_array) {
  si <- model_params["sigma_input"]
  inputs_scaled <- inputs
  inputs_scaled[, LITTER_COLS] <- inputs[, LITTER_COLS] * si    # scale inputs
  yasso07_run(inputs_scaled, model_params[YASSO07_PARAM_NAMES], C_init, xi_array)
}


# =============================================================================
# 5.  Prior specification (per-parameter widths)
# =============================================================================
# Rationale:
#   beta1, beta2:  exp(beta1*T + beta2*T^2) amplifies width dramatically.
#                  sigma=1 allows T=10°C modifier to span ~17 orders of magnitude.
#   sigma_init, sigma_input: log transform → tight enough for reasonable CI.
#   Stick-break and unconstrained: 1.0 is already weakly informative.

sigma_ppm <- setNames(rep(1.0, N_FREE), FREE_NAMES)
sigma_ppm["beta1"]       <- 0.2    # temp-mod 95% CI at 10°C: ~0.4× to ~7×
sigma_ppm["beta2"]       <- 0.05   # keeps quadratic term biologically plausible
sigma_ppm["sigma_init"]  <- 0.5    # physical 95% CI: ~0.04 to ~0.27
sigma_ppm["sigma_input"] <- 0.5    # physical 95% CI: ~0.37 to ~2.72

# Sanity check
stopifnot(length(sigma_ppm) == N_FREE)
stopifnot(all(sigma_ppm > 0))
stopifnot(all(names(sigma_ppm) == FREE_NAMES))

message("\nPrior widths (unconstrained space):")
print(round(sigma_ppm, 3))


# =============================================================================
# 5.5  Save processed inputs
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
# 6.  Build likelihood and prior
# =============================================================================

ll_fn <- make_likelihood(
  n_cores         = CORES_PER_CHAIN,
  to_original     = to_original,
  log_jacobian    = log_jacobian,
  assemble_params = assemble_model_params,
  compute_xi      = compute_xi_yasso07_engine,
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
message(sprintf("Likelihood at defaults: %.3f", ll_check))
if (!is.finite(ll_check)) stop("FAIL: likelihood not finite at defaults.")
message("Likelihood check: PASS\n")

# Prior: density is correct with dnorm + vector sd (element-wise on x vs best_x).
# Sampler needs per-column loop because rnorm(n*k, mean=vec, sd=vec) recycles
# element-by-element across the flat output — destroying column-wise structure.
prior <- createPrior(
  density = function(x) sum(dnorm(x, mean = best_x, sd = sigma_ppm, log = TRUE)),
  sampler = function(n = 1) {
    m <- matrix(0, nrow = n, ncol = N_FREE)               # per-column sampling
    for (j in seq_len(N_FREE)) {
      m[, j] <- rnorm(n, mean = best_x[j], sd = sigma_ppm[j] * 0.1)
    }
    m
  },
  lower = rep(-Inf, N_FREE),
  upper = rep( Inf, N_FREE)
)


# =============================================================================
# 7.  Run MCMC
# =============================================================================

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

HIGHLIGHT <- c("beta1","beta2","gamma","delta1","delta2","r",
               "sigma_init","sigma_input")

diag_out <- run_diagnostics(
  chain_results    = chain_results,
  to_original      = to_original,
  param_names      = FREE_NAMES,
  best_x           = best_x,
  sigma_ppm        = sigma_ppm,
  run_config       = run_config,
  mcmc_settings    = mcmc_settings,
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


# =============================================================================
# 9.  Posterior predictive output
# =============================================================================

message("-------------------------------------------------------------")
message("Generating posterior predictions...")
message(sprintf("  %d draws × %d plots × ~%d years",
                N_PP_DRAWS, length(plots_real),
                nrow(Yasso07_inputs) / length(plots_real)))
message("-------------------------------------------------------------\n")

set.seed(2025)
draw_idx <- sample(nrow(diag_out$all_phys), N_PP_DRAWS)
draws    <- diag_out$all_phys[draw_idx, ]

posterior_predictions <- do.call(rbind, lapply(seq_len(N_PP_DRAWS), function(d) {
  if (d %% 10 == 0) message(sprintf("  Draw %d / %d", d, N_PP_DRAWS))
  
  p_free       <- draws[d, ]
  model_params <- assemble_model_params(p_free)
  
  do.call(rbind, mclapply(plots_real, function(pid) {
    clim   <- climate_by_plot[[pid]]
    inputs <- inputs_by_plot[[pid]]
    lm     <- litter_means[[pid]]
    
    xi_array <- tryCatch(
      compute_xi_yasso07_engine(clim, model_params),
      error = function(e) NULL)
    if (is.null(xi_array)) return(NULL)
    
    n_ss    <- min(STEADY_STATE_YEARS, length(xi_array))
    xi_mean <- mean(xi_array[seq_len(n_ss)])
    if (!is.finite(xi_mean) || xi_mean <= 0) return(NULL)
    
    C_init <- tryCatch(
      steady_state_yasso07_engine(model_params, lm, xi_mean),
      error = function(e) NULL)
    if (is.null(C_init) || any(!is.finite(C_init)) || any(C_init < 0)) return(NULL)
    
    run_out <- tryCatch(
      yasso07_run_engine(inputs, model_params, C_init, xi_array),
      error = function(e) NULL)
    if (is.null(run_out)) return(NULL)
    
    data.frame(
      plot_id     = pid,
      year        = run_out$year,
      draw        = d,
      A           = run_out$A,
      W           = run_out$W,
      E           = run_out$E,
      N           = run_out$N,
      H           = run_out$H,
      total_soc   = run_out$total_soc,
      respiration = run_out$respiration
    )
  }, mc.cores = CORES_PER_CHAIN))
}))

message(sprintf("Posterior predictions: %d rows", nrow(posterior_predictions)))

# Summary table (defensive type coercion on plot_id for the join)
soc_obs_lookup <- SOC_obs_all[, c("plot_id", "year", "soc_obs_tCha")]
soc_obs_lookup$plot_id <- as.character(soc_obs_lookup$plot_id)

posterior_summary <- posterior_predictions %>%
  group_by(plot_id, year) %>%
  summarise(
    soc_mean   = mean(total_soc,  na.rm = TRUE),
    soc_median = median(total_soc, na.rm = TRUE),
    soc_sd     = sd(total_soc,    na.rm = TRUE),
    soc_q025   = quantile(total_soc, 0.025, na.rm = TRUE),
    soc_q975   = quantile(total_soc, 0.975, na.rm = TRUE),
    resp_mean  = mean(respiration, na.rm = TRUE),
    n_draws    = n(),
    .groups    = "drop"
  ) %>%
  left_join(soc_obs_lookup, by = c("plot_id", "year"))

plot_metadata <- site_raw[as.character(site_raw$plot_id) %in% plots_real,
                          c("plot_id", "region", "species_code", "species_name",
                            "soil_code", "temp_zone", "shallow",
                            "mean_litter", "mean_soc_Mgha", "litter_cv",
                            "lat_WGS84", "mean_temp", "mean_precip_annual",
                            "clay", "soil_depth")]

pp_output <- list(
  posterior_predictions = posterior_predictions,
  posterior_summary     = as.data.frame(posterior_summary),
  plot_metadata         = plot_metadata,
  config = list(
    model        = MODEL_NAME,
    run_id       = RUN_ID,
    n_draws      = N_PP_DRAWS,
    n_plots      = length(plots_real),
    steady_state_years = STEADY_STATE_YEARS,
    sigma_obs_fixed    = sigma_obs_fixed,
    sigma_ppm          = sigma_ppm,
    timestamp          = Sys.time()
  )
)

pp_file <- file.path(DIR_RUNS,
                     sprintf("%s_posterior_predictive_%s.rds", MODEL_NAME, RUN_ID))
saveRDS(pp_output, pp_file)
message(sprintf("Posterior predictive saved: [%s]", pp_file))


# =============================================================================
# 10.  Final summary + R-hat / ESS table
# =============================================================================

message("\n=============================================================")
message(sprintf("  %s Calibration Complete (real data)", MODEL_NAME))
message(sprintf("  Run ID:             %s", RUN_ID))
message(sprintf("  Wallclock:          %.1f minutes", t_run / 60))
message(sprintf("  Plots:              %d (of %d calib-ready)",
                length(plots), length(calib_plots)))
message(sprintf("  Free parameters:    %d", N_FREE))
message(sprintf("  sigma_obs (fixed):  %.3f", sigma_obs_fixed))
message(sprintf("  sigma_ppm (prior): %.2f – %.2f  (mean %.2f)",
                min(sigma_ppm), max(sigma_ppm), mean(sigma_ppm)))
message(sprintf("  Steady-state years: %d", STEADY_STATE_YEARS))
message(sprintf("  Posterior draws:    %d", N_PP_DRAWS))
if (!is.null(diag_out$gr)) {
  message(sprintf("  R-hat:              %.3f -- %.3f",
                  min(diag_out$gr$psrf[,1], na.rm=TRUE),
                  max(diag_out$gr$psrf[,1], na.rm=TRUE)))
}
if (!is.null(diag_out$ess_vals)) {
  message(sprintf("  ESS:                %.0f -- %.0f",
                  min(diag_out$ess_vals), max(diag_out$ess_vals)))
}
message("=============================================================\n")

# R-hat / ESS table, sorted by worst R-hat
rhat_df <- data.frame(
  param = rownames(diag_out$gr$psrf),
  rhat  = diag_out$gr$psrf[, 1],
  ess   = diag_out$ess_vals
)
rhat_df <- rhat_df[order(-rhat_df$rhat), ]
print(rhat_df)


# =============================================================================
# 11.  Diagnostic plots (prior/posterior, delta SOC trajectory)
#
#      Uses objects from the calibration run in memory:
#        diag_out$all_phys, best_x, sigma_ppm, to_original,
#        FREE_NAMES, posterior_predictions, posterior_summary
#
#      Outputs PNGs to ./Calibration/plots/ at 150 dpi.
# =============================================================================

plot_metadata$plot_id <- as.character(plot_metadata$plot_id)

posterior <- diag_out$all_phys

message("-------------------------------------------------------------")
message("Generating diagnostic plots...")
message(sprintf("  Output: %s%s_*_%s.png", DIR_PLOTS, MODEL_NAME, RUN_ID))
message("-------------------------------------------------------------")


# ---------------------------------------------------------------------------
# Plot 1: Prior vs posterior densities (highlighted parameters)
# ---------------------------------------------------------------------------
# Per-column prior sampling correctly applies the sigma_ppm vector.
set.seed(42)
n_prior   <- min(5000L, nrow(posterior))
prior_raw <- matrix(0, nrow = n_prior, ncol = length(best_x))
for (j in seq_along(best_x)) {
  prior_raw[, j] <- rnorm(n_prior, mean = best_x[j], sd = sigma_ppm[j])
}
prior_phys <- t(apply(prior_raw, 1, to_original))
colnames(prior_phys) <- FREE_NAMES

plot_params <- c("beta1", "beta2", "gamma",
                 "delta1", "delta2", "r",
                 "sigma_init", "sigma_input")

png(file.path(DIR_PLOTS,
              sprintf("%s_prior_posterior_%s.png", MODEL_NAME, RUN_ID)),
    width = 1800, height = 1200, res = 150)

par(mfrow = c(3, 3), mar = c(4, 4, 3, 1), oma = c(0, 0, 3, 0))

for (nm in plot_params) {
  post_v  <- posterior[is.finite(posterior[, nm]), nm]
  prior_v <- prior_phys[is.finite(prior_phys[, nm]), nm]
  
  # Clip prior tails to keep x-axis readable
  q_prior <- quantile(prior_v, c(0.005, 0.995), na.rm = TRUE)
  prior_v <- prior_v[prior_v >= q_prior[1] & prior_v <= q_prior[2]]
  
  d_post  <- density(post_v)
  d_prior <- density(prior_v)
  xlim_r  <- range(c(d_post$x, d_prior$x))
  ylim_r  <- range(c(d_post$y, d_prior$y))
  
  plot(NA, xlim = xlim_r, ylim = ylim_r,
       xlab = nm, ylab = "Density", main = nm)
  
  # Prior: transparent grey fill + line
  polygon(c(d_prior$x, rev(d_prior$x)),
          c(d_prior$y, rep(0, length(d_prior$y))),
          col = adjustcolor("grey60", alpha.f = 0.3), border = NA)
  lines(d_prior, col = "grey40", lwd = 1)
  
  # Posterior: transparent steelblue fill + line
  polygon(c(d_post$x, rev(d_post$x)),
          c(d_post$y, rep(0, length(d_post$y))),
          col = adjustcolor("steelblue", alpha.f = 0.4), border = NA)
  lines(d_post, col = "steelblue", lwd = 2)
  
  # Median (solid) / mean (dashed) for both prior and posterior
  abline(v = median(post_v),  col = "steelblue", lwd = 1.5, lty = 1)
  abline(v = mean(post_v),    col = "steelblue", lwd = 1.5, lty = 2)
  abline(v = median(prior_v), col = "grey40",    lwd = 1,   lty = 1)
  abline(v = mean(prior_v),   col = "grey40",    lwd = 1,   lty = 2)
}

# Legend panel
plot.new()
legend("center",
       legend = c("Prior density", "Posterior density",
                  "Median (solid)", "Mean (dashed)"),
       col    = c("grey40", "steelblue", "black", "black"),
       lwd    = c(1, 2, 1.5, 1.5),
       lty    = c(1, 1, 1, 2),
       fill   = c(adjustcolor("grey60", 0.3),
                  adjustcolor("steelblue", 0.4), NA, NA),
       border = NA, bty = "n", cex = 1.1)

mtext(sprintf("%s prior vs posterior  |  Run: %s  |  n_post = %d, n_prior = %d",
              MODEL_NAME, RUN_ID, nrow(posterior), n_prior),
      side = 3, outer = TRUE, line = 0.5, cex = 1.1, font = 2)

dev.off()
message(sprintf("  Saved: %s_prior_posterior_%s.png", MODEL_NAME, RUN_ID))

# ---------------------------------------------------------------------------
# Plot 1b: Prior vs posterior densities (flow fraction parameters)
# ---------------------------------------------------------------------------

frac_params <- list(
  A = c("p_AW","p_AE","p_AN"),
  W = c("p_WA","p_WE","p_WN"),
  E = c("p_EA","p_EW","p_EN"),
  N = c("p_NA","p_NW","p_NE")
)
frac_all <- unlist(frac_params, use.names = FALSE)

png(file.path(DIR_PLOTS,
              sprintf("%s_prior_posterior_fractions_%s.png", MODEL_NAME, RUN_ID)),
    width = 1800, height = 1600, res = 150)

par(mfrow = c(4, 3), mar = c(4, 4, 3, 1), oma = c(0, 0, 4, 0))

for (grp in names(frac_params)) {
  for (nm in frac_params[[grp]]) {
    post_v  <- posterior[is.finite(posterior[, nm]), nm]
    prior_v <- prior_phys[is.finite(prior_phys[, nm]), nm]
    
    # Clip extreme prior tails
    q_prior <- quantile(prior_v, c(0.005, 0.995), na.rm = TRUE)
    prior_v <- prior_v[prior_v >= q_prior[1] & prior_v <= q_prior[2]]
    
    d_post  <- density(post_v,  from = 0, to = 1)
    d_prior <- density(prior_v, from = 0, to = 1)
    xlim_r  <- c(0, 1)
    ylim_r  <- range(c(d_post$y, d_prior$y))
    
    plot(NA, xlim = xlim_r, ylim = ylim_r,
         xlab = nm, ylab = "Density",
         main = sprintf("%s  [donor: %s]", nm, grp))
    
    polygon(c(d_prior$x, rev(d_prior$x)),
            c(d_prior$y, rep(0, length(d_prior$y))),
            col = adjustcolor("grey60", alpha.f = 0.3), border = NA)
    lines(d_prior, col = "grey40", lwd = 1)
    
    polygon(c(d_post$x, rev(d_post$x)),
            c(d_post$y, rep(0, length(d_post$y))),
            col = adjustcolor("steelblue", alpha.f = 0.4), border = NA)
    lines(d_post, col = "steelblue", lwd = 2)
    
    abline(v = median(post_v),  col = "steelblue", lwd = 1.5, lty = 1)
    abline(v = mean(post_v),    col = "steelblue", lwd = 1.5, lty = 2)
    abline(v = median(prior_v), col = "grey40",    lwd = 1,   lty = 1)
    abline(v = mean(prior_v),   col = "grey40",    lwd = 1,   lty = 2)
  }
}

mtext(sprintf("%s flow fractions — prior vs posterior  |  Run: %s  |  n_post = %d, n_prior = %d",
              MODEL_NAME, RUN_ID, nrow(posterior), n_prior),
      side = 3, outer = TRUE, line = 1, cex = 1.1, font = 2)

dev.off()
message(sprintf("  Saved: %s_prior_posterior_fractions_%s.png", MODEL_NAME, RUN_ID))


# ---------------------------------------------------------------------------
# Plot 2: Mean ΔSOC trajectory across plots with posterior CI
# ---------------------------------------------------------------------------

pp <- posterior_predictions

# ΔSOC(t) per plot per draw, relative to year 1
pp <- pp %>%
  group_by(plot_id, draw) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(delta_soc = total_soc - first(total_soc)) %>%
  ungroup()

# Mean across plots per draw × year
mean_delta <- pp %>%
  group_by(draw, year) %>%
  summarise(mean_delta_soc = mean(delta_soc, na.rm = TRUE), .groups = "drop")

# Aggregate across draws
traj_summary <- mean_delta %>%
  group_by(year) %>%
  summarise(
    median_delta = median(mean_delta_soc, na.rm = TRUE),
    mean_delta   = mean(mean_delta_soc,   na.rm = TRUE),
    q025 = quantile(mean_delta_soc, 0.025, na.rm = TRUE),
    q250 = quantile(mean_delta_soc, 0.250, na.rm = TRUE),
    q750 = quantile(mean_delta_soc, 0.750, na.rm = TRUE),
    q975 = quantile(mean_delta_soc, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>% arrange(year)

# Observed ΔSOC (mean across plots for each observation year)
obs_baseline <- posterior_summary %>%
  filter(!is.na(soc_obs_tCha)) %>%
  group_by(plot_id) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(obs_delta = soc_obs_tCha - first(soc_obs_tCha)) %>%
  ungroup()

obs_mean_delta <- obs_baseline %>%
  group_by(year) %>%
  summarise(obs_mean = mean(obs_delta, na.rm = TRUE),
            obs_sd   = sd(obs_delta,   na.rm = TRUE),
            obs_n    = n(),
            obs_se   = obs_sd / sqrt(obs_n),
            .groups  = "drop")

png(file.path(DIR_PLOTS,
              sprintf("%s_delta_soc_trajectory_%s.png", MODEL_NAME, RUN_ID)),
    width = 1800, height = 1000, res = 150)

par(mar = c(4, 4, 3, 1))

ylim_r <- range(c(traj_summary$q025, traj_summary$q975,
                  obs_mean_delta$obs_mean + 1.96 * obs_mean_delta$obs_se,
                  obs_mean_delta$obs_mean - 1.96 * obs_mean_delta$obs_se,
                  0), na.rm = TRUE)
ylim_r <- ylim_r + c(-0.05, 0.05) * diff(ylim_r)

plot(NA, xlim = range(traj_summary$year), ylim = ylim_r,
     xlab = "Year",
     ylab = expression(paste(Delta, "SOC relative to year 1 (tC/ha)")),
     main = sprintf("Mean ΔSOC across plots (n = %d)  |  %s  |  Run: %s",
                    length(unique(pp$plot_id)), MODEL_NAME, RUN_ID))

# 95% CI (outer)
polygon(c(traj_summary$year, rev(traj_summary$year)),
        c(traj_summary$q025, rev(traj_summary$q975)),
        col = adjustcolor("steelblue", alpha.f = 0.15), border = NA)

# 50% CI (inner)
polygon(c(traj_summary$year, rev(traj_summary$year)),
        c(traj_summary$q250, rev(traj_summary$q750)),
        col = adjustcolor("steelblue", alpha.f = 0.35), border = NA)

# Median (solid) / mean (dashed)
lines(traj_summary$year, traj_summary$median_delta,
      col = "steelblue", lwd = 2, lty = 1)
lines(traj_summary$year, traj_summary$mean_delta,
      col = "steelblue", lwd = 1.5, lty = 2)

abline(h = 0, lty = 3, col = "grey40")

if (nrow(obs_mean_delta) > 0) {
  points(obs_mean_delta$year, obs_mean_delta$obs_mean,
         pch = 19, col = "firebrick", cex = 1.3)
  segments(obs_mean_delta$year,
           obs_mean_delta$obs_mean - 1.96 * obs_mean_delta$obs_se,
           obs_mean_delta$year,
           obs_mean_delta$obs_mean + 1.96 * obs_mean_delta$obs_se,
           col = "firebrick", lwd = 1.5)
}

legend("topleft",
       legend = c("Posterior median", "Posterior mean",
                  "50% CI", "95% CI",
                  "Observed mean (± 95% CI)"),
       col    = c("steelblue", "steelblue",
                  adjustcolor("steelblue", 0.35),
                  adjustcolor("steelblue", 0.15),
                  "firebrick"),
       lwd    = c(2, 1.5, NA, NA, NA),
       lty    = c(1, 2, NA, NA, NA),
       pch    = c(NA, NA, 15, 15, 19),
       pt.cex = c(NA, NA, 2, 2, 1.3),
       bty    = "n", cex = 0.9)

dev.off()
message(sprintf("  Saved: %s_delta_soc_trajectory_%s.png", MODEL_NAME, RUN_ID))

# ---------------------------------------------------------------------------
# Plot 3: Observed vs predicted SOC (1:1 scatterplot, colored by region)
# ---------------------------------------------------------------------------

obs_pred <- posterior_summary %>%
  filter(!is.na(soc_obs_tCha)) %>%
  left_join(
    plot_metadata[, c("plot_id", "region")],
    by = "plot_id"
  )

region_cols <- c("1" = "#d73027", "2" = "#4575b4")   # 1=South, 2=North
region_labs <- c("1" = "Southern Finland", "2" = "Northern Finland")

png(file.path(DIR_PLOTS,
              sprintf("%s_obs_vs_pred_%s.png", MODEL_NAME, RUN_ID)),
    width = 1200, height = 1200, res = 150)

par(mar = c(4.5, 4.5, 3, 1))

ax_lim <- range(c(obs_pred$soc_obs_tCha, obs_pred$soc_median,
                  obs_pred$soc_q025, obs_pred$soc_q975), na.rm = TRUE)
ax_lim <- c(0, ax_lim[2] * 1.05)

plot(NA, xlim = ax_lim, ylim = ax_lim,
     xlab = "Observed SOC (tC/ha)",
     ylab = "Predicted SOC — posterior median (tC/ha)",
     main = sprintf("Observed vs predicted  |  %s  |  Run: %s",
                    MODEL_NAME, RUN_ID))

abline(0, 1, lty = 2, col = "grey40", lwd = 1.5)

# 95% CI error bars (vertical)
segments(obs_pred$soc_obs_tCha, obs_pred$soc_q025,
         obs_pred$soc_obs_tCha, obs_pred$soc_q975,
         col = adjustcolor(region_cols[as.character(obs_pred$region)], 0.3),
         lwd = 0.8)

# Points
points(obs_pred$soc_obs_tCha, obs_pred$soc_median,
       pch = 16, cex = 0.9,
       col = adjustcolor(region_cols[as.character(obs_pred$region)], 0.7))

# Fit statistics
r2   <- cor(obs_pred$soc_obs_tCha, obs_pred$soc_median, use = "complete.obs")^2
rmse <- sqrt(mean((obs_pred$soc_obs_tCha - obs_pred$soc_median)^2, na.rm = TRUE))
bias <- mean(obs_pred$soc_median - obs_pred$soc_obs_tCha, na.rm = TRUE)

legend("topleft",
       legend = c(region_labs,
                  "",
                  sprintf("R² = %.3f", r2),
                  sprintf("RMSE = %.1f tC/ha", rmse),
                  sprintf("Bias = %+.1f tC/ha", bias),
                  "1:1 line"),
       col    = c(region_cols, NA, NA, NA, NA, "grey40"),
       pch    = c(16, 16, NA, NA, NA, NA, NA),
       lty    = c(NA, NA, NA, NA, NA, NA, 2),
       lwd    = c(NA, NA, NA, NA, NA, NA, 1.5),
       bty    = "n", cex = 0.85)

dev.off()
message(sprintf("  Saved: %s_obs_vs_pred_%s.png", MODEL_NAME, RUN_ID))


# ---------------------------------------------------------------------------
# Printout: predicted vs observed ΔSOC
# ---------------------------------------------------------------------------

message("\n=== ΔSOC predictions vs observations ===")
for (y in obs_mean_delta$year) {
  pred <- traj_summary[traj_summary$year == y, ]
  obs  <- obs_mean_delta[obs_mean_delta$year == y, ]
  message(sprintf("Year %d: observed mean ΔSOC = %+.2f tC/ha (n=%d, SE=%.2f)",
                  y, obs$obs_mean, obs$obs_n, obs$obs_se))
  message(sprintf("          predicted median  = %+.2f  [95%% CI: %+.2f, %+.2f]",
                  pred$median_delta, pred$q025, pred$q975))
}

final <- tail(traj_summary, 1)
message(sprintf("\nFinal year (%d): median ΔSOC = %+.2f tC/ha", final$year, final$median_delta))
message(sprintf("  50%% CI: [%+.2f, %+.2f]", final$q250, final$q750))
message(sprintf("  95%% CI: [%+.2f, %+.2f]", final$q025, final$q975))

message("-------------------------------------------------------------\n")








########### TESTING SPACE
# Pick one plot that has observations
pid <- plots[10]
obs <- SOC_obs_all[SOC_obs_all$plot_id == pid, "soc_obs_tCha"]
obs
pid_pp <- posterior_predictions[posterior_predictions$plot_id == pid, ]
obs_plot <- SOC_obs_all[SOC_obs_all$plot_id == pid, c("year", "soc_obs_tCha")]

# Summarise across draws
pid_summary <- pid_pp %>%
  group_by(year) %>%
  summarise(
    median_soc = median(total_soc),
    q025 = quantile(total_soc, 0.025),
    q975 = quantile(total_soc, 0.975),
    median_H = median(H),
    .groups = "drop"
  )

ylim_r <- range(c(pid_summary$q025, pid_summary$q975, 
                  obs_plot$soc_obs_tCha), na.rm = TRUE)

plot(pid_summary$year, pid_summary$median_soc, type = "l",
     col = "steelblue", lwd = 2, ylim = ylim_r,
     xlab = "Year", ylab = "SOC (tC/ha)",
     main = sprintf("Plot %s — posterior predictions vs observed", pid))

polygon(c(pid_summary$year, rev(pid_summary$year)),
        c(pid_summary$q025, rev(pid_summary$q975)),
        col = adjustcolor("steelblue", 0.2), border = NA)

# H pool separately
lines(pid_summary$year, pid_summary$median_H,
      col = "darkorange", lwd = 1.5, lty = 2)

points(obs_plot$year, obs_plot$soc_obs_tCha,
       pch = 19, col = "firebrick", cex = 1.5)

legend("topleft",
       legend = c("Predicted total (median)", "95% CI", "H pool (median)", "Observed"),
       col    = c("steelblue", adjustcolor("steelblue", 0.2), "darkorange", "firebrick"),
       lwd    = c(2, NA, 1.5, NA), lty = c(1, NA, 2, NA),
       pch    = c(NA, 15, NA, 19), pt.cex = c(NA, 2, NA, 1.5),
       bty    = "n")


