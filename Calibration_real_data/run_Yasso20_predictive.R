# =============================================================================
# run_Yasso20_predictive.R
#
# Posterior predictive stage — Yasso20 / HIKET pipeline.
# Mirrors run_Yasso15_predictive.R. Differences from Yasso15 are marked YASSO20.
#
# DIFFERENCES FROM YASSO15:
#   * Sources yasso15_wrapper.R + yasso20_wrapper.R; loads yasso15.so.
#   * Engine wrappers use Yasso20 functions (monthly climate, yasso15_run direct).
#   * n_ss = STEADY_STATE_YEARS * 12L (monthly rows, not annual).
#   * SOC obs year lookup uses inputs_by_plot[[pid]]$year, not
#     climate_by_plot[[pid]]$year, because climate is monthly (repeated years).
#
# USAGE:
#   Rscript run_Yasso20_predictive.R <RUN_ID>
# =============================================================================

source("./Calibration_real_data/calibration_engine.R")
source("./Model_functions_real_data/input_compatibility_layer.R")
source("./Model_functions_real_data/Decomposition_functions/Yasso/yasso15_wrapper.R")
source("./Model_functions_real_data/Decomposition_functions/Yasso/yasso20_wrapper.R")
dyn.load("./Model_functions_real_data/Decomposition_functions/Yasso/yasso15.so")

library(dplyr)
library(parallel)


# =============================================================================
# 0.  Configuration
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) {
  RUN_ID <- args[1]
} else {
  RUN_ID <- NULL
  if (is.null(RUN_ID)) {
    rns <- list.files(
      "./Calibration_real_data/runs/",
      pattern = "^Yasso20_posterior_[0-9]{8}_[0-9]{6}\\.rds$"
    )
    if (length(rns) == 0)
      stop("No Yasso20 posterior found in ./Calibration_real_data/runs/")
    rns    <- sort(rns, decreasing = TRUE)
    RUN_ID <- sub("^Yasso20_posterior_(.+)\\.rds$", "\\1", rns[1])
    message(sprintf("Auto-detected RUN_ID: %s", RUN_ID))
  }
}

MODEL_NAME <- "Yasso20"
N_PP_DRAWS <- 100L

CORES_PER_CHAIN <- parallel::detectCores() - 1L
# CORES_PER_CHAIN <- 90L   # Roihu

DIR_DIAG   <- file.path("./Calibration_real_data/diagnostics", MODEL_NAME)
DIR_RUNS   <- "./Calibration_real_data/runs/"
DIR_INPUTS <- "./Data/model_inputs/"
for (d in c(DIR_DIAG, DIR_RUNS, DIR_INPUTS))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

run_config <- list(MODEL_NAME = MODEL_NAME, RUN_ID = RUN_ID, DIR_DIAG = DIR_DIAG)

message("=============================================================")
message(sprintf("  %s Posterior Predictive  |  Run: %s", MODEL_NAME, RUN_ID))
message("=============================================================\n")


# =============================================================================
# 1.  Load calibration outputs
# =============================================================================

posterior_phys <- readRDS(file.path(DIR_RUNS,
                                    sprintf("Yasso20_posterior_%s.rds", RUN_ID)))
inputs_pkg     <- readRDS(file.path(DIR_INPUTS,
                                    sprintf("Yasso20_inputs_%s.rds", RUN_ID)))

climate_by_plot <- inputs_pkg$climate_by_plot   # monthly format
inputs_by_plot  <- inputs_pkg$inputs_by_plot    # annual format
litter_means    <- inputs_pkg$litter_means
obs_meta        <- inputs_pkg$obs_meta
plot_info       <- inputs_pkg$plot_info
plots_real      <- inputs_pkg$plots_real
sigma_obs_fixed <- inputs_pkg$sigma_obs_fixed
STEADY_STATE_YEARS <- inputs_pkg$STEADY_STATE_YEARS

message(sprintf("Loaded posterior:   %d samples x %d params",
                nrow(posterior_phys), ncol(posterior_phys)))
message(sprintf("Loaded plots:       %d", length(plots_real)))


# =============================================================================
# 2.  Re-derive parameter assembly
# =============================================================================

p_default        <- YASSO20_DEFAULT_PARAMS
FIXED_RATE_NAMES <- c("alpha_A","alpha_W","alpha_E","alpha_N",
                      "p_H","alpha_H","w1","w2","w3","w4","w5")
fixed_rates      <- p_default[FIXED_RATE_NAMES]

FREE_FRAC_NAMES <- c("p_AW","p_AE","p_AN",
                     "p_WA","p_WE","p_WN",
                     "p_EA","p_EW","p_EN",
                     "p_NA","p_NW","p_NE")

assemble_model_params <- function(p_free) {
  yasso_params <- c(
    fixed_rates,
    p_free[FREE_FRAC_NAMES],
    p_free[c("beta1","beta2","gamma",
             "betaN1","betaN2","gammaN",
             "betaH1","betaH2","gammaH",
             "delta1","delta2","r")]
  )
  yasso_params <- yasso_params[YASSO20_PARAM_NAMES]
  c(yasso_params, sigma_input = unname(p_free["sigma_input"]))
}


# Diagnostic: test first plot to see which step fails
test_pid <- plots_real[1]
test_clim <- climate_by_plot[[test_pid]]
test_inputs <- inputs_by_plot[[test_pid]]
test_lm <- litter_means[[test_pid]]

test_xi <- tryCatch(
  compute_xi_yasso20_engine(test_clim, model_params),
  error = function(e) { message("XI_ERROR: ", e$message); NULL })

if (!is.null(test_xi)) {
  test_xi_ss <- tryCatch(
    compute_xi_mean_yasso20_engine(
      test_clim[seq_len(min(STEADY_STATE_YEARS * 12L, nrow(test_clim))), , drop = FALSE],
      model_params),
    error = function(e) { message("XI_SS_ERROR: ", e$message); NULL })
  
  if (!is.null(test_xi_ss)) {
    test_cinit <- tryCatch(
      steady_state_yasso20_engine(model_params, test_lm, test_xi_ss),
      error = function(e) { message("CINIT_ERROR: ", e$message); NULL })
    
    if (!is.null(test_cinit)) {
      test_run <- tryCatch(
        yasso20_run_engine(test_inputs, model_params, test_cinit, test_xi),
        error = function(e) { message("RUN_ERROR: ", e$message); NULL })
    }
  }
}

# =============================================================================
# 3.  Model interface wrappers  [YASSO20]
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

yasso20_run_engine <- function(inputs, model_params, C_init, xi_arrays) {
  si <- model_params["sigma_input"]
  inputs_scaled <- inputs
  inputs_scaled[, LITTER_COLS] <- inputs[, LITTER_COLS] * si
  yasso15_run(
    input_df  = inputs_scaled,
    params    = model_params[YASSO20_PARAM_NAMES],
    C_init    = C_init,
    xi_arrays = xi_arrays,
    precip    = inputs$precip
  )
}


# =============================================================================
# 4.  Posterior predictive simulation  [YASSO20: n_ss in monthly rows]
# =============================================================================

set.seed(2025)
draw_idx <- sample(nrow(posterior_phys), N_PP_DRAWS)
draws    <- posterior_phys[draw_idx, ]

# Initialize failure counters
xi_fail_count <- 0
xi_ss_fail_count <- 0
cinit_fail_count <- 0
run_fail_count <- 0

message(sprintf("\nGenerating %d posterior predictive draws across %d plots...",
                N_PP_DRAWS, length(plots_real)))
t_pp <- system.time({
  
  # DEBUG: Force write diagnostics
  debug_file <- file.path(DIR_DIAG, sprintf("Yasso20_debug_%s.txt", RUN_ID))
  
  # Check 1: inputs structure before draws
  sink(debug_file, append = FALSE)
  cat("=== YASSO20 DEBUG START ===\n")
  cat("posterior_phys dims:", nrow(posterior_phys), "x", ncol(posterior_phys), "\n")
  cat("plots_real length:", length(plots_real), "\n")
  cat("climate_by_plot is NULL?", is.null(climate_by_plot), "\n")
  cat("First plot climate structure:\n")
  if (length(plots_real) > 0) {
    first_pid <- plots_real[1]
    cat("  First plot ID:", first_pid, "\n")
    if (!is.null(climate_by_plot[[first_pid]])) {
      cat("  Climate rows:", nrow(climate_by_plot[[first_pid]]), "\n")
      cat("  Climate cols:", ncol(climate_by_plot[[first_pid]]), "\n")
      cat("  Climate col names:", paste(colnames(climate_by_plot[[first_pid]]), collapse = ", "), "\n")
    } else {
      cat("  Climate for first plot is NULL!\n")
    }
  }
  sink()
  
  posterior_predictions <- do.call(rbind, lapply(seq_len(N_PP_DRAWS), function(d) {
    if (d %% 10 == 0) message(sprintf("  Draw %d / %d", d, N_PP_DRAWS))

    p_free       <- draws[d, ]
    model_params <- assemble_model_params(p_free)

    do.call(rbind, mclapply(plots_real, function(pid) {
      clim   <- climate_by_plot[[pid]]
      inputs <- inputs_by_plot[[pid]]
      lm     <- litter_means[[pid]]
      
      xi_array <- tryCatch(
        compute_xi_yasso20_engine(clim, model_params),
        error = function(e) NULL)
      if (is.null(xi_array)) {
        assign("xi_fail_count", get("xi_fail_count", envir = .GlobalEnv) + 1, envir = .GlobalEnv)
        return(NULL)
      }
      
      n_ss      <- min(STEADY_STATE_YEARS * 12L, nrow(clim))
      xi_for_ss <- tryCatch(
        compute_xi_mean_yasso20_engine(
          clim[seq_len(n_ss), , drop = FALSE], model_params),
        error = function(e) NULL)
      if (is.null(xi_for_ss) || !is_valid_xi(xi_for_ss)) {
        assign("xi_ss_fail_count", get("xi_ss_fail_count", envir = .GlobalEnv) + 1, envir = .GlobalEnv)
        return(NULL)
      }
      
      C_init <- tryCatch(
        steady_state_yasso20_engine(model_params, lm, xi_for_ss),
        error = function(e) NULL)
      if (is.null(C_init) || any(!is.finite(C_init)) || any(C_init < 0)) {
        assign("cinit_fail_count", get("cinit_fail_count", envir = .GlobalEnv) + 1, envir = .GlobalEnv)
        return(NULL)
      }
      
      run_out <- tryCatch(
        yasso20_run_engine(inputs, model_params, C_init, xi_array),
        error = function(e) NULL)
      if (is.null(run_out)) {
        assign("run_fail_count", get("run_fail_count", envir = .GlobalEnv) + 1, envir = .GlobalEnv)
        return(NULL)
      }
      
      result_df <- data.frame(
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
      assign("success_count", get("success_count", envir = .GlobalEnv) + nrow(result_df), envir = .GlobalEnv)
      result_df
    }, mc.cores = CORES_PER_CHAIN))
  }))
  message(sprintf("\nFailure counts across all plots/draws:"))
  message(sprintf("  XI failures: %d", xi_fail_count))
  message(sprintf("  XI_SS failures: %d", xi_ss_fail_count))
  message(sprintf("  C_init failures: %d", cinit_fail_count))
  message(sprintf("  Run failures: %d", run_fail_count))
  
  write(sprintf(
    "=== FAILURE COUNTERS ===\nXI: %d\nXI_SS: %d\nC_init: %d\nRun: %d",
    xi_fail_count, xi_ss_fail_count, cinit_fail_count, run_fail_count),
    file = file.path(DIR_DIAG, "YASSO20_FAILURES.txt"))
})["elapsed"]
message(sprintf("Posterior predictive complete: %.1f min  (%d rows)",
                t_pp / 60, nrow(posterior_predictions)))

# Check 2: posterior_predictions structure after draws
sink(debug_file, append = TRUE)
cat("\n=== After posterior predictive generation ===\n")
cat("posterior_predictions class:", class(posterior_predictions), "\n")
cat("posterior_predictions is NULL?", is.null(posterior_predictions), "\n")
if (!is.null(posterior_predictions)) {
  cat("posterior_predictions dims:", nrow(posterior_predictions), "x", ncol(posterior_predictions), "\n")
  cat("posterior_predictions colnames:", paste(colnames(posterior_predictions), collapse = ", "), "\n")
  cat("First few rows of posterior_predictions:\n")
  print(head(posterior_predictions))
} else {
  cat("posterior_predictions is NULL - THIS IS THE PROBLEM\n")
}
sink()

message(sprintf("Debug output written to: %s", debug_file))


# =============================================================================
# 5.  Posterior summary table
# =============================================================================
# YASSO20: obs year lookup uses inputs_by_plot[[pid]]$year (annual), not
# climate_by_plot[[pid]]$year (monthly, repeated). obs_meta$idx was built
# against annual inputs in the calibration script.

SOC_obs_all <- do.call(rbind, lapply(names(obs_meta), function(pid) {
  m <- obs_meta[[pid]]
  if (length(m$soc_obs) == 0) return(NULL)
  inp_yrs <- inputs_by_plot[[pid]]$year   # annual -- YASSO20 difference
  data.frame(plot_id = pid,
             year    = inp_yrs[m$idx],
             soc_obs_tCha = m$soc_obs,
             is_first = m$is_first)
}))

posterior_summary <- posterior_predictions %>%
  group_by(plot_id, year) %>%
  summarise(
    soc_mean   = mean(total_soc,   na.rm = TRUE),
    soc_median = median(total_soc, na.rm = TRUE),
    soc_sd     = sd(total_soc,     na.rm = TRUE),
    soc_q025   = quantile(total_soc, 0.025, na.rm = TRUE),
    soc_q975   = quantile(total_soc, 0.975, na.rm = TRUE),
    resp_mean  = mean(respiration, na.rm = TRUE),
    n_draws    = n(),
    .groups    = "drop"
  ) %>%
  left_join(SOC_obs_all, by = c("plot_id", "year"))


# =============================================================================
# 6.  Residuals data frame
# =============================================================================

site_raw <- read.csv("./Data/model_inputs/site_raw.csv")
site_raw$plot_id <- as.character(site_raw$plot_id)

residuals_df <- posterior_summary %>%
  filter(!is.na(soc_obs_tCha)) %>%
  mutate(
    plot_id      = as.character(plot_id),
    log_obs      = log(soc_obs_tCha),
    log_hat_mean = log(soc_mean),
    residual_log = log_obs - log_hat_mean,
    residual_abs = soc_obs_tCha - soc_mean
  ) %>%
  select(plot_id, year, soc_obs_tCha, soc_mean, soc_median,
         soc_q025, soc_q975, log_obs, log_hat_mean,
         residual_log, residual_abs, is_first) %>%
  left_join(site_raw, by = "plot_id")

residuals_csv <- file.path(DIR_DIAG,
                           sprintf("%s_residuals_%s.csv", MODEL_NAME, RUN_ID))
write.csv(residuals_df, residuals_csv, row.names = FALSE)
message(sprintf("Residuals: [%s]  (%d rows)", residuals_csv, nrow(residuals_df)))


# =============================================================================
# 7.  Predictive metrics
# =============================================================================

obs     <- residuals_df$soc_obs_tCha
hat     <- residuals_df$soc_mean
hat_med <- residuals_df$soc_median

R2   <- cor(obs, hat, use = "complete.obs")^2
RMSE <- sqrt(mean((obs - hat)^2,     na.rm = TRUE))
bias <- mean(hat - obs,              na.rm = TRUE)

in_ci       <- residuals_df$soc_obs_tCha >= residuals_df$soc_q025 &
               residuals_df$soc_obs_tCha <= residuals_df$soc_q975
coverage_95 <- mean(in_ci, na.rm = TRUE)

RMSE_med <- sqrt(mean((obs - hat_med)^2, na.rm = TRUE))
bias_med <- mean(hat_med - obs,          na.rm = TRUE)

message("\nPredictive metrics:")
message(sprintf("  R²:            %.3f", R2))
message(sprintf("  RMSE (mean):   %.2f tC/ha", RMSE))
message(sprintf("  Bias (mean):   %+.2f tC/ha", bias))
message(sprintf("  RMSE (median): %.2f tC/ha", RMSE_med))
message(sprintf("  Bias (median): %+.2f tC/ha", bias_med))
message(sprintf("  95%% coverage:  %.3f", coverage_95))


# =============================================================================
# 8.  Plots (identical to Yasso15/07)
# =============================================================================

PX_PER_IN <- 150L

obs_pred_png <- file.path(DIR_DIAG,
                          sprintf("%s_obs_vs_pred_%s.png", MODEL_NAME, RUN_ID))
ktp_labels <- c(
  "1" = "1 Lehdot (OMaT)",       "2" = "2 Lehtomainen (OMT)",
  "3" = "3 Tuorekangas (MT)",    "4" = "4 Kuivahko (VT)",
  "5" = "5 Kuiva (CT)",          "6" = "6 Karukkokangas (ClT)",
  "7" = "7 Kalliot/hietikot",    "8" = "8 Lakimetsät/tunturit"
)
ktp_palette <- c(
  "1" = "#000000", "2" = "#E69F00", "3" = "#56B4E9", "4" = "#009E73",
  "5" = "#F0E442", "6" = "#0072B2", "7" = "#D55E00", "8" = "#CC79A7"
)
ktp_vals <- as.character(residuals_df$kasvup_tyyppi)
ktp_col  <- ktp_palette[ktp_vals]
ktp_col[is.na(ktp_col)] <- "grey60"

png(obs_pred_png, width = 7L * PX_PER_IN, height = 7L * PX_PER_IN, res = PX_PER_IN)
par(mar = c(4.5, 4.5, 3, 1))
ax_lim <- c(0, max(c(residuals_df$soc_obs_tCha, residuals_df$soc_q975),
                   na.rm = TRUE) * 1.05)
plot(NA, xlim = ax_lim, ylim = ax_lim,
     xlab = "Observed SOC (tC/ha)",
     ylab = "Predicted SOC -- posterior median (tC/ha)",
     main = sprintf("%s  |  %s", MODEL_NAME, RUN_ID))
abline(0, 1, lty = 2, col = "grey40", lwd = 1.5)
segments(residuals_df$soc_obs_tCha, residuals_df$soc_q025,
         residuals_df$soc_obs_tCha, residuals_df$soc_q975,
         col = adjustcolor(ktp_col, 0.3), lwd = 0.8)
points(residuals_df$soc_obs_tCha, residuals_df$soc_median,
       pch = 16, cex = 0.9, col = adjustcolor(ktp_col, 0.7))
legend("topleft",
       legend = c(sprintf("R² = %.3f", R2),
                  sprintf("RMSE = %.1f tC/ha", RMSE_med),
                  sprintf("Bias = %+.1f tC/ha", bias_med),
                  sprintf("95%% coverage = %.2f", coverage_95)),
       bty = "n", cex = 0.85)
ktp_present <- sort(unique(na.omit(ktp_vals)))
ktp_present <- ktp_present[ktp_present %in% names(ktp_palette)]
if (length(ktp_present) > 0)
  legend("bottomright", legend = ktp_labels[ktp_present],
         col = ktp_palette[ktp_present],
         pch = 16, bty = "n", cex = 0.75, title = "kasvup_tyyppi")
dev.off()
message(sprintf("Plot saved: %s", obs_pred_png))

# --- Plot B: mean ΔSOC trajectory ---
pp <- posterior_predictions %>%
  group_by(plot_id, draw) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(
    first_soc  = first(total_soc),
    first_year = first(year),
    annual_rate = (total_soc - first_soc) / (year - first_year)
  ) %>%
  filter(year != first_year) %>%   # drop t=0 (0/0 undefined)
  ungroup()

traj_summary <- pp %>%
  group_by(draw, year) %>%
  summarise(mean_annual_rate = mean(annual_rate, na.rm = TRUE), .groups = "drop") %>%
  group_by(year) %>%
  summarise(
    median_delta = median(mean_annual_rate, na.rm = TRUE),
    q025 = quantile(mean_annual_rate, 0.025, na.rm = TRUE),
    q250 = quantile(mean_annual_rate, 0.250, na.rm = TRUE),
    q750 = quantile(mean_annual_rate, 0.750, na.rm = TRUE),
    q975 = quantile(mean_annual_rate, 0.975, na.rm = TRUE),
    .groups = "drop") %>%
  arrange(year)

obs_delta_per_plot <- posterior_summary %>%
  filter(!is.na(soc_obs_tCha)) %>%
  group_by(plot_id) %>%
  arrange(year, .by_group = TRUE) %>%
  summarise(
    delta_obs  = (last(soc_obs_tCha) - first(soc_obs_tCha)) /
      (last(year) - first(year)),
    year_first = first(year),
    year_last  = last(year),
    .groups = "drop") %>%
  filter(year_first != year_last)

obs_mean_delta <- mean(obs_delta_per_plot$delta_obs, na.rm = TRUE)
obs_se_delta   <- sd(obs_delta_per_plot$delta_obs, na.rm = TRUE) /
  sqrt(nrow(obs_delta_per_plot))

traj_png <- file.path(DIR_DIAG,
                      sprintf("%s_delta_soc_trajectory_%s.png", MODEL_NAME, RUN_ID))
png(traj_png, width = 12L * PX_PER_IN, height = 7L * PX_PER_IN, res = PX_PER_IN)
par(mar = c(4, 4, 3, 1))
ylim_r <- range(c(traj_summary$q025, traj_summary$q975,
                  obs_mean_delta + 1.96 * obs_se_delta,
                  obs_mean_delta - 1.96 * obs_se_delta,
                  0), na.rm = TRUE)
ylim_r <- ylim_r + c(-0.05, 0.05) * diff(ylim_r)
plot(NA, xlim = range(traj_summary$year), ylim = ylim_r,
     xlab = "Year",
     ylab = "Mean annual ΔSOC since first year (tC/ha/yr)",
     main = sprintf("Mean annual ΔSOC across plots  |  %s  |  %s",
                    MODEL_NAME, RUN_ID))
polygon(c(traj_summary$year, rev(traj_summary$year)),
        c(traj_summary$q025, rev(traj_summary$q975)),
        col = adjustcolor("steelblue", 0.15), border = NA)
polygon(c(traj_summary$year, rev(traj_summary$year)),
        c(traj_summary$q250, rev(traj_summary$q750)),
        col = adjustcolor("steelblue", 0.35), border = NA)
lines(traj_summary$year, traj_summary$median_delta, col = "steelblue", lwd = 2)
abline(h = 0, lty = 3, col = "grey40")
abline(h = obs_mean_delta, col = "firebrick", lwd = 2.5, lty = 2)
rect(xleft   = min(obs_delta_per_plot$year_first),
     xright  = max(obs_delta_per_plot$year_last),
     ybottom = obs_mean_delta - 1.96 * obs_se_delta,
     ytop    = obs_mean_delta + 1.96 * obs_se_delta,
     col = adjustcolor("firebrick", 0.12), border = NA)
legend("topleft",
       legend = c("Posterior median", "50% CI", "95% CI",
                  sprintf("Observed mean annual ΔSOC ± 95%% CI  (n = %d)",
                          nrow(obs_delta_per_plot))),
       col    = c("steelblue", adjustcolor("steelblue", 0.35),
                  adjustcolor("steelblue", 0.15), "firebrick"),
       lwd = c(2, NA, NA, 2.5), pch = c(NA, 15, 15, NA), lty = c(1, NA, NA, 2),
       pt.cex = c(NA, 2, 2, NA), bty = "n", cex = 0.9)
dev.off()
message(sprintf("Plot saved: %s", traj_png))


# =============================================================================
# 9.  Save
# =============================================================================

pp_output <- list(
  posterior_predictions = posterior_predictions,
  posterior_summary     = as.data.frame(posterior_summary),
  residuals_df          = residuals_df,
  metrics = list(
    R2          = R2,
    RMSE_mean   = RMSE,
    bias_mean   = bias,
    RMSE_median = RMSE_med,
    bias_median = bias_med,
    coverage_95 = coverage_95
  ),
  config = list(
    model              = MODEL_NAME,
    run_id             = RUN_ID,
    n_draws            = N_PP_DRAWS,
    n_plots            = length(plots_real),
    steady_state_years = STEADY_STATE_YEARS,
    sigma_obs_fixed    = sigma_obs_fixed,
    timestamp          = Sys.time()
  )
)

pp_file <- file.path(DIR_RUNS,
                     sprintf("%s_posterior_predictive_%s.rds", MODEL_NAME, RUN_ID))
saveRDS(pp_output, pp_file)
message(sprintf("Posterior predictive saved: %s", pp_file))


# =============================================================================
# 10.  Append [5] to report
# =============================================================================

append_to_report(run_config, paste(c(
  "\n",
  "[5] Predictive performance\n",
  sprintf("    Posterior draws used:   %d\n", N_PP_DRAWS),
  sprintf("    Observations:           %d\n", nrow(residuals_df)),
  sprintf("    R²:                     %.3f\n", R2),
  sprintf("    RMSE (mean predictor):  %.2f tC/ha\n", RMSE),
  sprintf("    Bias (mean):            %+.2f tC/ha\n", bias),
  sprintf("    RMSE (median):          %.2f tC/ha\n", RMSE_med),
  sprintf("    Bias (median):          %+.2f tC/ha\n", bias_med),
  sprintf("    95%% CI coverage:        %.3f   [target ~0.95]\n", coverage_95),
  "\n"
), collapse = ""))

message("\nDone. Next: run_residual_analysis.R Yasso20 ", RUN_ID)
