# =============================================================================
# run_TP2_transient_predictive.R
#
# Posterior predictive stage for TP2 transient calibration.
# Differs from run_TP2_predictive.R only in:
#   - Sources calibration_engine_transient.R and tp2_wrapper_transient.R
#   - DIR_RUNS / DIR_DIAG point to Calibration_real_data_transient/
#   - steady_state engine wrapper calls tp2_transient_init
#
# USAGE:
#   Rscript run_TP2_transient_predictive.R <RUN_ID>
# =============================================================================

source("./Calibration_real_data_transient/calibration_engine_transient.R")
source("./Model_functions_real_data_transient/Decomposition_functions/Yasso/yasso07_wrapper_transient.R")
# No dyn.load: only the pure-R xi functions are used (compute_xi_yasso07,
# compute_xi_mean_yasso07), not yasso07_run which requires the Fortran .so
source("./Model_functions_real_data_transient/Decomposition_functions/SimpleModels/tp2_wrapper_transient.R")
source("./Model_functions_real_data_transient/input_compatibility_layer.R")

library(dplyr)
library(parallel)


# =============================================================================
# 0.  Configuration
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) {
  RUN_ID <- args[1]
} else {
  rns <- list.files(
    "./Calibration_real_data_transient/runs/",
    pattern = "^TP2_posterior_[0-9]{8}_[0-9]{6}\\.rds$"
  )
  if (length(rns) == 0)
    stop("No TP2 transient posterior found in Calibration_real_data_transient/runs/")
  rns    <- sort(rns, decreasing = TRUE)
  RUN_ID <- sub("^TP2_posterior_(.+)\\.rds$", "\\1", rns[1])
  message(sprintf("Auto-detected RUN_ID: %s", RUN_ID))
}

MODEL_NAME      <- "TP2"
N_PP_DRAWS      <- 100L
CORES_PER_CHAIN <- parallel::detectCores() - 1L

# Forward projection settings (Section 4b)
PROJ_YEARS    <- 60L   # years to project beyond last historical year
RECYCLE_YEARS <- 20L   # length of climate window to recycle
POOL_COLS     <- c("A", "H")   # TP2 has two carbon pools

DIR_DIAG   <- file.path("./Calibration_real_data_transient/diagnostics", MODEL_NAME)
DIR_RUNS   <- "./Calibration_real_data_transient/runs/"
DIR_INPUTS <- "./Data/model_inputs/"
for (d in c(DIR_DIAG, DIR_RUNS)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

run_config <- list(MODEL_NAME = MODEL_NAME, RUN_ID = RUN_ID, DIR_DIAG = DIR_DIAG)

message("=============================================================")
message(sprintf("  %s Transient Posterior Predictive  |  Run: %s",
                MODEL_NAME, RUN_ID))
message("=============================================================\n")


# =============================================================================
# 1.  Load calibration outputs
# =============================================================================

posterior_phys <- readRDS(file.path(DIR_RUNS,
                                    sprintf("TP2_posterior_%s.rds", RUN_ID)))
inputs_pkg     <- readRDS(file.path(DIR_INPUTS,
                                    sprintf("TP2_inputs_%s.rds",    RUN_ID)))

climate_by_plot    <- inputs_pkg$climate_by_plot
inputs_by_plot     <- inputs_pkg$inputs_by_plot
litter_means       <- inputs_pkg$litter_means
obs_meta           <- inputs_pkg$obs_meta
plot_info          <- inputs_pkg$plot_info
plots_real         <- inputs_pkg$plots_real
sigma_obs_fixed    <- inputs_pkg$sigma_obs_fixed
STEADY_STATE_YEARS <- inputs_pkg$STEADY_STATE_YEARS

message(sprintf("Loaded posterior:   %d samples x %d params",
                nrow(posterior_phys), ncol(posterior_phys)))
message(sprintf("Loaded plots:       %d", length(plots_real)))


# =============================================================================
# 2.  Parameter assembly (pass-through: all TP2 params free)
# =============================================================================

assemble_model_params <- function(p_free) p_free


# =============================================================================
# 3.  Model interface wrappers
# =============================================================================

compute_xi_tp2_engine <- function(clim, model_params) {
  compute_xi_yasso07(
    temp_mean = clim$temp_mean,
    temp_amp  = clim$temp_amplitude,
    precip    = clim$precip,
    beta1     = model_params["beta1"],
    beta2     = model_params["beta2"],
    gamma     = model_params["gamma"]
  )
}

compute_xi_mean_tp2_engine <- function(clim_ss, model_params) {
  compute_xi_mean_yasso07(
    clim_ss = clim_ss,
    beta1   = model_params["beta1"],
    beta2   = model_params["beta2"],
    gamma   = model_params["gamma"]
  )
}

# TRANSIENT: tp2_transient_init replaces tp2_steady_state
steady_state_tp2_engine <- function(model_params, lm, xi_mean) {
  tp2_transient_init(model_params, lm, xi_mean)
}

tp2_run_engine <- function(inputs, model_params, C_init, xi_array) {
  tp2_run(inputs, model_params, C_init, xi_array)
}


# =============================================================================
# 4.  Posterior predictive simulation (historical period)
# =============================================================================

set.seed(2025)
draw_idx <- sample(nrow(posterior_phys), N_PP_DRAWS)
draws    <- posterior_phys[draw_idx, ]

message(sprintf("\nGenerating %d posterior predictive draws across %d plots...",
                N_PP_DRAWS, length(plots_real)))
t_pp <- system.time({
  posterior_predictions <- do.call(rbind, lapply(seq_len(N_PP_DRAWS), function(d) {
    if (d %% 10 == 0) message(sprintf("  Draw %d / %d", d, N_PP_DRAWS))
    p_free       <- draws[d, ]
    model_params <- assemble_model_params(p_free)

    do.call(rbind, mclapply(plots_real, function(pid) {
      clim   <- climate_by_plot[[pid]]
      inputs <- inputs_by_plot[[pid]]
      lm     <- litter_means[[pid]]

      xi_array <- tryCatch(compute_xi_tp2_engine(clim, model_params),
                           error = function(e) NULL)
      if (is.null(xi_array)) return(NULL)

      n_ss      <- min(STEADY_STATE_YEARS, nrow(clim))
      xi_for_ss <- tryCatch(
        compute_xi_mean_tp2_engine(
          clim[seq_len(n_ss), , drop = FALSE], model_params),
        error = function(e) NULL)
      if (is.null(xi_for_ss) || !is.finite(xi_for_ss) || xi_for_ss <= 0)
        return(NULL)

      C_init <- tryCatch(steady_state_tp2_engine(model_params, lm, xi_for_ss),
                         error = function(e) NULL)
      if (is.null(C_init) || any(!is.finite(C_init)) || any(C_init < 0))
        return(NULL)

      run_out <- tryCatch(tp2_run_engine(inputs, model_params, C_init, xi_array),
                          error = function(e) NULL)
      if (is.null(run_out)) return(NULL)

      data.frame(plot_id   = pid,
                 year      = run_out$year,
                 draw      = d,
                 A         = run_out$A,
                 H         = run_out$H,
                 total_soc = run_out$total_soc)
    }, mc.cores = CORES_PER_CHAIN))
  }))
})[["elapsed"]]
message(sprintf("Posterior predictive complete: %.1f min  (%d rows)",
                t_pp / 60, nrow(posterior_predictions)))


# =============================================================================
# 4b.  60-year forward projection
# =============================================================================
# Extends each posterior draw PROJ_YEARS beyond the last historical year.
#
# The historical run is repeated here to recover the terminal pool state.
# This duplicates Section 4 but keeps both sections self-contained and
# Sections 5-9 completely unchanged.
#
# Climate: last RECYCLE_YEARS of observed climate recycled cyclically.
# Litter:  constant at the last available annual value (last inputs row).
# Initial state: terminal pool values from the repeated historical run.
# =============================================================================

message(sprintf("\nGenerating %d-year forward projections (%d draws x %d plots)...",
                PROJ_YEARS, N_PP_DRAWS, length(plots_real)))
t_proj <- system.time({
  projection_predictions <- do.call(rbind, lapply(seq_len(N_PP_DRAWS), function(d) {
    if (d %% 10 == 0) message(sprintf("  Projection draw %d / %d", d, N_PP_DRAWS))
    p_free       <- draws[d, ]
    model_params <- assemble_model_params(p_free)

    do.call(rbind, mclapply(plots_real, function(pid) {
      clim   <- climate_by_plot[[pid]]
      inputs <- inputs_by_plot[[pid]]
      lm     <- litter_means[[pid]]

      # -- Repeat historical run to recover terminal pool state ---------------
      xi_array <- tryCatch(compute_xi_tp2_engine(clim, model_params),
                           error = function(e) NULL)
      if (is.null(xi_array)) return(NULL)

      n_ss      <- min(STEADY_STATE_YEARS, nrow(clim))
      xi_for_ss <- tryCatch(
        compute_xi_mean_tp2_engine(
          clim[seq_len(n_ss), , drop = FALSE], model_params),
        error = function(e) NULL)
      if (is.null(xi_for_ss) || !is.finite(xi_for_ss) || xi_for_ss <= 0)
        return(NULL)

      C_init <- tryCatch(steady_state_tp2_engine(model_params, lm, xi_for_ss),
                         error = function(e) NULL)
      if (is.null(C_init) || any(!is.finite(C_init)) || any(C_init < 0))
        return(NULL)

      run_hist <- tryCatch(tp2_run_engine(inputs, model_params, C_init, xi_array),
                           error = function(e) NULL)
      if (is.null(run_hist)) return(NULL)

      # -- Build projection climate and litter --------------------------------
      last_year   <- max(run_hist$year)
      proj_yrs    <- seq(last_year + 1L, last_year + PROJ_YEARS)

      clim_recycle <- tail(clim, RECYCLE_YEARS)
      clim_proj    <- clim_recycle[
        ((seq_len(PROJ_YEARS) - 1L) %% RECYCLE_YEARS) + 1L, , drop = FALSE]
      clim_proj$year <- proj_yrs

      inputs_proj      <- tail(inputs, 1L)[rep(1L, PROJ_YEARS), , drop = FALSE]
      inputs_proj$year <- proj_yrs

      # -- Initial state: terminal pool values from historical run ------------
      C_proj_init <- unlist(run_hist[nrow(run_hist), POOL_COLS])

      # -- Forward run --------------------------------------------------------
      xi_proj <- tryCatch(compute_xi_tp2_engine(clim_proj, model_params),
                          error = function(e) NULL)
      if (is.null(xi_proj)) return(NULL)

      run_proj <- tryCatch(
        tp2_run_engine(inputs_proj, model_params, C_proj_init, xi_proj),
        error = function(e) NULL)
      if (is.null(run_proj)) return(NULL)

      data.frame(
        plot_id   = pid,
        year      = run_proj$year,
        draw      = d,
        total_soc = run_proj$total_soc
      )
    }, mc.cores = CORES_PER_CHAIN))
  }))
})[["elapsed"]]
message(sprintf("Projection complete: %.1f min  (%d rows)",
                t_proj / 60, nrow(projection_predictions)))


# =============================================================================
# 5.  Posterior summary
# =============================================================================

SOC_obs_all <- do.call(rbind, lapply(names(obs_meta), function(pid) {
  m <- obs_meta[[pid]]
  if (length(m$soc_obs) == 0) return(NULL)
  clim_yrs <- climate_by_plot[[pid]]$year
  data.frame(plot_id = pid, year = clim_yrs[m$idx],
             soc_obs_tCha = m$soc_obs, is_first = m$is_first)
}))

posterior_summary <- posterior_predictions %>%
  group_by(plot_id, year) %>%
  summarise(
    soc_mean   = mean(total_soc,   na.rm = TRUE),
    soc_median = median(total_soc, na.rm = TRUE),
    soc_sd     = sd(total_soc,     na.rm = TRUE),
    soc_q025   = quantile(total_soc, 0.025, na.rm = TRUE),
    soc_q975   = quantile(total_soc, 0.975, na.rm = TRUE),
    n_draws    = n(), .groups = "drop") %>%
  left_join(SOC_obs_all, by = c("plot_id", "year"))


# =============================================================================
# 6.  Residuals CSV
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

R2          <- cor(obs, hat, use = "complete.obs")^2
RMSE        <- sqrt(mean((obs - hat)^2,     na.rm = TRUE))
bias        <- mean(hat - obs,              na.rm = TRUE)
RMSE_med    <- sqrt(mean((obs - hat_med)^2, na.rm = TRUE))
bias_med    <- mean(hat_med - obs,          na.rm = TRUE)
coverage_95 <- mean(obs >= residuals_df$soc_q025 &
                    obs <= residuals_df$soc_q975, na.rm = TRUE)

message("\nPredictive metrics:")
message(sprintf("  R²:            %.3f", R2))
message(sprintf("  RMSE (mean):   %.2f tC/ha", RMSE))
message(sprintf("  Bias (mean):   %+.2f tC/ha", bias))
message(sprintf("  RMSE (median): %.2f tC/ha", RMSE_med))
message(sprintf("  Bias (median): %+.2f tC/ha", bias_med))
message(sprintf("  95%% coverage:  %.3f", coverage_95))


# =============================================================================
# 8.  Plots
# =============================================================================

PX_PER_IN <- 150L

# --- Obs vs predicted ---
obs_pred_png <- file.path(DIR_DIAG,
                          sprintf("%s_obs_vs_pred_%s.png", MODEL_NAME, RUN_ID))
ktp_palette <- c("1"="#000000","2"="#E69F00","3"="#56B4E9","4"="#009E73",
                 "5"="#F0E442","6"="#0072B2","7"="#D55E00","8"="#CC79A7")
ktp_labels  <- c("1"="1 Lehdot (OMaT)","2"="2 Lehtomainen (OMT)",
                 "3"="3 Tuorekangas (MT)","4"="4 Kuivahko (VT)",
                 "5"="5 Kuiva (CT)","6"="6 Karukkokangas (ClT)",
                 "7"="7 Kalliot/hietikot","8"="8 Lakimetsät/tunturit")
ktp_vals <- as.character(residuals_df$kasvup_tyyppi)
ktp_col  <- ktp_palette[ktp_vals]; ktp_col[is.na(ktp_col)] <- "grey60"

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
         col = ktp_palette[ktp_present], pch = 16, bty = "n", cex = 0.75,
         title = "kasvup_tyyppi")
dev.off()
message(sprintf("Plot saved: %s", obs_pred_png))

# --- Annual ΔSOC trajectory ---
pp <- posterior_predictions %>%
  group_by(plot_id, draw) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(first_soc  = first(total_soc),
         first_year = first(year),
         annual_rate = (total_soc - first_soc) / (year - first_year)) %>%
  filter(year != first_year) %>%
  ungroup()

traj_summary <- pp %>%
  group_by(draw, year) %>%
  summarise(mean_annual_rate = mean(annual_rate, na.rm = TRUE), .groups = "drop") %>%
  group_by(year) %>%
  summarise(median_delta = median(mean_annual_rate, na.rm = TRUE),
            q025 = quantile(mean_annual_rate, 0.025, na.rm = TRUE),
            q250 = quantile(mean_annual_rate, 0.250, na.rm = TRUE),
            q750 = quantile(mean_annual_rate, 0.750, na.rm = TRUE),
            q975 = quantile(mean_annual_rate, 0.975, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(year)

obs_delta_per_plot <- posterior_summary %>%
  filter(!is.na(soc_obs_tCha)) %>%
  group_by(plot_id) %>% arrange(year, .by_group = TRUE) %>%
  summarise(delta_obs  = (last(soc_obs_tCha) - first(soc_obs_tCha)) /
                           (last(year) - first(year)),
            year_first = first(year), year_last = last(year), .groups = "drop") %>%
  filter(year_first != year_last)
obs_mean_delta <- mean(obs_delta_per_plot$delta_obs, na.rm = TRUE)
obs_se_delta   <- sd(obs_delta_per_plot$delta_obs, na.rm = TRUE) /
                  sqrt(nrow(obs_delta_per_plot))

traj_png <- file.path(DIR_DIAG,
                      sprintf("%s_delta_soc_trajectory_%s.png", MODEL_NAME, RUN_ID))
png(traj_png, width = 12L * PX_PER_IN, height = 7L * PX_PER_IN, res = PX_PER_IN)
par(mar = c(4, 4, 3, 1))
ylim_r <- range(c(traj_summary$q025, traj_summary$q975,
                  obs_mean_delta + 1.96*obs_se_delta,
                  obs_mean_delta - 1.96*obs_se_delta, 0), na.rm = TRUE)
ylim_r <- ylim_r + c(-0.05, 0.05) * diff(ylim_r)
plot(NA, xlim = range(traj_summary$year), ylim = ylim_r,
     xlab = "Year",
     ylab = "Mean annual ΔSOC since first year (tC/ha/yr)",
     main = sprintf("Mean annual ΔSOC  |  %s  |  %s", MODEL_NAME, RUN_ID))
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
     ybottom = obs_mean_delta - 1.96*obs_se_delta,
     ytop    = obs_mean_delta + 1.96*obs_se_delta,
     col = adjustcolor("firebrick", 0.12), border = NA)
legend("topleft",
       legend = c("Posterior median","50% CI","95% CI",
                  sprintf("Observed mean annual ΔSOC ± 95%% CI  (n = %d)",
                          nrow(obs_delta_per_plot))),
       col = c("steelblue", adjustcolor("steelblue",0.35),
               adjustcolor("steelblue",0.15), "firebrick"),
       lwd = c(2,NA,NA,2.5), pch = c(NA,15,15,NA), lty = c(1,NA,NA,2),
       pt.cex = c(NA,2,2,NA), bty = "n", cex = 0.9)
dev.off()
message(sprintf("Plot saved: %s", traj_png))


# =============================================================================
# 9.  Save
# =============================================================================

pp_output <- list(
  posterior_predictions  = posterior_predictions,
  projection_predictions = projection_predictions,
  posterior_summary      = as.data.frame(posterior_summary),
  residuals_df           = residuals_df,
  metrics = list(R2 = R2, RMSE_mean = RMSE, bias_mean = bias,
                 RMSE_median = RMSE_med, bias_median = bias_med,
                 coverage_95 = coverage_95),
  config  = list(model = MODEL_NAME, run_id = RUN_ID,
                 n_draws = N_PP_DRAWS, n_plots = length(plots_real),
                 proj_years = PROJ_YEARS, recycle_years = RECYCLE_YEARS,
                 transient_init = TRUE, timestamp = Sys.time())
)
pp_file <- file.path(DIR_RUNS,
                     sprintf("%s_posterior_predictive_%s.rds", MODEL_NAME, RUN_ID))
saveRDS(pp_output, pp_file)
message(sprintf("Saved: %s", pp_file))

append_to_report(run_config, paste(c(
  "\n[5] Predictive performance (transient init)\n",
  sprintf("    Posterior draws:        %d\n",       N_PP_DRAWS),
  sprintf("    Observations:           %d\n",       nrow(residuals_df)),
  sprintf("    R²:                     %.3f\n",     R2),
  sprintf("    RMSE (mean):            %.2f tC/ha\n", RMSE),
  sprintf("    Bias (mean):            %+.2f tC/ha\n", bias),
  sprintf("    RMSE (median):          %.2f tC/ha\n", RMSE_med),
  sprintf("    Bias (median):          %+.2f tC/ha\n", bias_med),
  sprintf("    95%% CI coverage:        %.3f\n",    coverage_95),
  sprintf("    Projection years:       %d\n",       PROJ_YEARS),
  sprintf("    Climate recycle window: %d years\n", RECYCLE_YEARS)), collapse=""))

message("\nDone. Next: run_residual_analysis.R TP2 ", RUN_ID)
