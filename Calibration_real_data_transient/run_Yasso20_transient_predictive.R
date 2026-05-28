# =============================================================================
# run_Yasso20_transient_predictive.R
#
# Posterior predictive stage for Yasso20 transient calibration.
# Key differences from run_Yasso20_predictive.R:
#   - Sources yasso20_wrapper_transient.R (standalone)
#   - DIR_RUNS / DIR_DIAG point to Calibration_real_data_transient/
#   - steady_state engine calls yasso20_transient_init
#   - assemble_model_params includes sigma_init
#   - yasso20_run_engine: yasso20_run takes climate_df + inputs_df (not
#     xi_array separately); it computes xi and precip internally.
# =============================================================================

source("./Calibration_real_data_transient/calibration_engine_transient.R")
source("./Model_functions_real_data/input_compatibility_layer.R")
source("./Model_functions_real_data_transient/Decomposition_functions/Yasso/yasso15_wrapper_transient.R")
source("./Model_functions_real_data_transient/Decomposition_functions/Yasso/yasso20_wrapper_transient.R")

dyn.load("./Model_functions_real_data_transient/Decomposition_functions/Yasso/yasso15.so")

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
    pattern = "^Yasso20_posterior_[0-9]{8}_[0-9]{6}\\.rds$"
  )
  if (length(rns) == 0)
    stop("No Yasso20 transient posterior found.")
  rns    <- sort(rns, decreasing = TRUE)
  RUN_ID <- sub("^Yasso20_posterior_(.+)\\.rds$", "\\1", rns[1])
  message(sprintf("Auto-detected RUN_ID: %s", RUN_ID))
}

MODEL_NAME      <- "Yasso20"
N_PP_DRAWS      <- 100L
CORES_PER_CHAIN <- parallel::detectCores() - 1L

# Forward projection settings (Section 4b)
PROJ_YEARS    <- 60L                        # years to project beyond last historical year
RECYCLE_YEARS <- 20L                        # length of climate window to recycle
POOL_COLS     <- c("A","W","E","N","H")     # Yasso20 five-pool state vector
# Note: climate recycling uses RECYCLE_YEARS * 12L monthly rows (see Section 4b)

LITTER_COLS <- c("nwl_A","nwl_W","nwl_E","nwl_N",
                 "fwl_A","fwl_W","fwl_E","fwl_N",
                 "cwl_A","cwl_W","cwl_E","cwl_N")

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
                                    sprintf("Yasso20_posterior_%s.rds", RUN_ID)))
inputs_pkg     <- readRDS(file.path(DIR_INPUTS,
                                    sprintf("Yasso20_inputs_%s.rds",    RUN_ID)))

climate_by_plot    <- inputs_pkg$climate_by_plot
inputs_by_plot     <- inputs_pkg$inputs_by_plot
litter_means       <- inputs_pkg$litter_means
obs_meta           <- inputs_pkg$obs_meta
plot_info          <- inputs_pkg$plot_info
plots_real         <- inputs_pkg$plots_real
holdout_plots      <- inputs_pkg$holdout_plots   # loaded for validation subsetting
sigma_obs_fixed    <- inputs_pkg$sigma_obs_fixed
STEADY_STATE_YEARS <- inputs_pkg$STEADY_STATE_YEARS

message(sprintf("Loaded posterior:   %d samples x %d params",
                nrow(posterior_phys), ncol(posterior_phys)))
message(sprintf("Loaded plots:       %d", length(plots_real)))


# =============================================================================
# 2.  Parameter assembly
# =============================================================================

FIXED_RATE_NAMES <- setdiff(YASSO20_PARAM_NAMES,
                             setdiff(names(inputs_pkg$free_defaults),
                                     c("sigma_init","sigma_input")))
fixed_rates      <- YASSO20_DEFAULT_PARAMS[FIXED_RATE_NAMES]
MODEL_FREE_NAMES <- setdiff(colnames(posterior_phys), c("sigma_init","sigma_input"))

assemble_model_params <- function(p_free) {
  yasso_params <- c(fixed_rates, p_free[MODEL_FREE_NAMES])
  yasso_params <- yasso_params[YASSO20_PARAM_NAMES]
  c(yasso_params,
    sigma_input = unname(p_free["sigma_input"]),
    sigma_init  = unname(p_free["sigma_init"]))
}


# =============================================================================
# 3.  Model interface wrappers
# =============================================================================
# yasso20_run() computes xi and precip internally from climate_df.
# compute_xi_mean_yasso20 still needed for yasso20_transient_init (steady state).

compute_xi_mean_yasso20_engine <- function(clim_ss, model_params) {
  compute_xi_mean_yasso20(
    clim_ss = clim_ss,
    params  = model_params[YASSO20_PARAM_NAMES]
  )
}

# TRANSIENT: yasso20_transient_init replaces yasso20_steady_state
steady_state_yasso20_engine <- function(model_params, lm, xi_ss) {
  yasso20_transient_init(model_params, lm, xi_ss)
}

# yasso20_run computes xi and precip from climate_df internally.
# Pass climate_df directly; sigma_input scales inputs_df beforehand.
yasso20_run_engine <- function(inputs, model_params, C_init, clim) {
  si <- model_params["sigma_input"]
  inputs_scaled <- inputs
  inputs_scaled[, LITTER_COLS] <- inputs[, LITTER_COLS] * si
  yasso20_run(
    climate_df = clim,
    inputs_df  = inputs_scaled,
    params     = model_params[YASSO20_PARAM_NAMES],
    C_init     = C_init
  )
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

      n_ss      <- min(STEADY_STATE_YEARS, nrow(clim))
      xi_for_ss <- tryCatch(
        compute_xi_mean_yasso20_engine(
          clim[seq_len(n_ss), , drop = FALSE], model_params),
        error = function(e) NULL)
      if (is.null(xi_for_ss)) return(NULL)
      xi_vals <- unlist(xi_for_ss)
      if (any(!is.finite(xi_vals)) || any(xi_vals <= 0)) return(NULL)

      C_init <- tryCatch(
        steady_state_yasso20_engine(model_params, lm, xi_for_ss),
        error = function(e) NULL)
      if (is.null(C_init) || any(!is.finite(C_init)) || any(C_init < 0))
        return(NULL)

      # yasso20_run takes clim directly (computes xi and precip internally)
      run_out <- tryCatch(
        yasso20_run_engine(inputs, model_params, C_init, clim),
        error = function(e) NULL)
      if (is.null(run_out)) return(NULL)

      data.frame(plot_id   = pid,
                 year      = run_out$year,
                 draw      = d,
                 A         = run_out$A,
                 W         = run_out$W,
                 E         = run_out$E,
                 N         = run_out$N,
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
# YASSO20-SPECIFIC: climate_by_plot holds monthly rows (12 per year).
#   Climate recycling operates on RECYCLE_YEARS * 12 monthly rows and
#   produces PROJ_YEARS * 12 projection rows. The year and month columns
#   are overwritten to reflect the projection period.
#   yasso20_run() computes xi internally from clim_proj -- no separate
#   xi computation is needed before the forward run call.
#
# Litter:  constant at the last available annual value (last inputs row).
# Initial state: terminal pool values from the repeated historical run.
#
# No SOC ceiling guard is applied: structurally unstable draws (Yasso20's
# N-pool recycling architecture) produce large but finite total_soc and are
# visible in the projection envelope. This is the intended behaviour for
# the cross-model instability comparison.
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
      n_ss      <- min(STEADY_STATE_YEARS, nrow(clim))
      xi_for_ss <- tryCatch(
        compute_xi_mean_yasso20_engine(
          clim[seq_len(n_ss), , drop = FALSE], model_params),
        error = function(e) NULL)
      if (is.null(xi_for_ss)) return(NULL)
      xi_vals <- unlist(xi_for_ss)
      if (any(!is.finite(xi_vals)) || any(xi_vals <= 0)) return(NULL)

      C_init <- tryCatch(
        steady_state_yasso20_engine(model_params, lm, xi_for_ss),
        error = function(e) NULL)
      if (is.null(C_init) || any(!is.finite(C_init)) || any(C_init < 0))
        return(NULL)

      run_hist <- tryCatch(
        yasso20_run_engine(inputs, model_params, C_init, clim),
        error = function(e) NULL)
      if (is.null(run_hist)) return(NULL)

      # -- Build projection climate (monthly) and litter (annual) ------------
      last_year <- max(run_hist$year)
      proj_yrs  <- seq(last_year + 1L, last_year + PROJ_YEARS)

      # Monthly climate: recycle last RECYCLE_YEARS * 12 rows cyclically
      clim_recycle <- tail(clim, RECYCLE_YEARS * 12L)
      clim_proj    <- clim_recycle[
        ((seq_len(PROJ_YEARS * 12L) - 1L) %% (RECYCLE_YEARS * 12L)) + 1L,
        , drop = FALSE]
      clim_proj$year  <- rep(proj_yrs, each = 12L)
      clim_proj$month <- rep(1:12, PROJ_YEARS)

      # Annual litter: constant at last available year
      inputs_proj      <- tail(inputs, 1L)[rep(1L, PROJ_YEARS), , drop = FALSE]
      inputs_proj$year <- proj_yrs

      # Belt-and-suspenders: tail() + index subset carries original row names,
      # which can confuse downstream rbind / dplyr joins. Reset both frames here.
      rownames(inputs_proj) <- NULL
      rownames(clim_proj)   <- NULL
      
      # -- Initial state: exact terminal per-cohort state from historical run --
      # PATCH (2025-05-25): replaced steady-state workaround with C_final.
      #
      # Previous code re-ran yasso20_steady_state() under "terminal climate"
      # (clim_recycle xi) to approximate C_proj_init. This was wrong in two
      # ways:
      #   (a) steady state != transient end-state after 21 years of forcing;
      #   (b) the result was a 5-element summed vector, but Fortran yasso15_run_r
      #       declares C_init(15) = [C_nwl(1:5) | C_fwl(6:10) | C_cwl(11:15)].
      #       Passing 5 elements caused Fortran to read garbage memory for
      #       slots 6-15, producing the ~10 tC/ha discontinuity at 2006 and
      #       NaN / Inf explosions in ~20% of draw x plot pairs.
      #
      # Fix chain (all already in place before this script):
      #   yasso15.f90  yasso15_run / yasso15_run_r: C_final(15) INTENT(OUT),
      #                populated as [C_nwl | C_fwl | C_cwl] after the time loop.
      #   yasso15_wrapper_transient.R  yasso15_run(): .Fortran() receives
      #                C_final = single(15); attached as attr(out, "C_final")
      #                converted back to double.
      #   yasso20_wrapper_transient.R  yasso20_run(): delegates directly to
      #                yasso15_run() with no intermediate data frame, so the
      #                attribute survives the passthrough automatically.
      #   This script: consume attr(run_hist, "C_final") here.
      #
      # No separate yasso20.f90 exists -- Yasso20 reuses yasso15.so throughout.
      C_proj_init <- attr(run_hist, "C_final")

      # -- Forward run --------------------------------------------------------
      # yasso20_run computes xi internally from clim_proj; no xi call needed
      run_proj <- tryCatch(
        yasso20_run_engine(inputs_proj, model_params, C_proj_init, clim_proj),
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
# 5-9.  Summary, residuals, metrics, plots, save
# =============================================================================

SOC_obs_all <- do.call(rbind, lapply(names(obs_meta), function(pid) {
  m <- obs_meta[[pid]]
  if (length(m$soc_obs) == 0) return(NULL)
  inp_yrs  <- inputs_by_plot[[pid]]$year
  data.frame(plot_id = pid, year = inp_yrs[m$idx],
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

# --- Tag holdout rows (used by residual analysis and multimodel) ---
residuals_df$is_holdout <- residuals_df$plot_id %in% holdout_plots

residuals_csv <- file.path(DIR_DIAG,
                           sprintf("%s_residuals_%s.csv", MODEL_NAME, RUN_ID))
write.csv(residuals_df, residuals_csv, row.names = FALSE)
message(sprintf("Residuals: [%s]  (%d rows)", residuals_csv, nrow(residuals_df)))

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
message(sprintf("  R²: %.3f  |  RMSE: %.2f  |  Bias: %+.2f  |  Cov95: %.3f",
                R2, RMSE_med, bias_med, coverage_95))

# --- Split calibration / holdout and compute metrics for each ---
res_calib   <- residuals_df[!residuals_df$is_holdout, ]
res_holdout <- residuals_df[ residuals_df$is_holdout, ]

compute_metrics <- function(df) {
  obs <- df$soc_obs_tCha; hat <- df$soc_mean; hat_med <- df$soc_median
  list(R2          = cor(obs, hat, use = "complete.obs")^2,
       RMSE_mean   = sqrt(mean((obs - hat)^2,     na.rm = TRUE)),
       bias_mean   = mean(hat - obs,              na.rm = TRUE),
       RMSE_median = sqrt(mean((obs - hat_med)^2, na.rm = TRUE)),
       bias_median = mean(hat_med - obs,          na.rm = TRUE),
       coverage_95 = mean(obs >= df$soc_q025 & obs <= df$soc_q975, na.rm = TRUE))
}

metrics_calib   <- compute_metrics(res_calib)
metrics_holdout <- compute_metrics(res_holdout)

message("\nCalibration metrics:")
message(sprintf("  R²: %.3f  RMSE: %.2f  Bias: %+.2f  Cov: %.3f",
                metrics_calib$R2, metrics_calib$RMSE_median,
                metrics_calib$bias_median, metrics_calib$coverage_95))
message("Holdout (independent validation) metrics:")
message(sprintf("  R²: %.3f  RMSE: %.2f  Bias: %+.2f  Cov: %.3f",
                metrics_holdout$R2, metrics_holdout$RMSE_median,
                metrics_holdout$bias_median, metrics_holdout$coverage_95))

PX_PER_IN <- 150L

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
     main = sprintf("%s transient  |  %s", MODEL_NAME, RUN_ID))
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

pp <- posterior_predictions %>%
  group_by(plot_id, draw) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(first_soc   = first(total_soc),
         first_year  = first(year),
         annual_rate = (total_soc - first_soc) / (year - first_year)) %>%
  filter(year != first_year) %>% ungroup()

traj_summary <- pp %>%
  group_by(draw, year) %>%
  summarise(mean_annual_rate = mean(annual_rate, na.rm = TRUE), .groups = "drop") %>%
  group_by(year) %>%
  summarise(median_delta = median(mean_annual_rate, na.rm = TRUE),
            q025 = quantile(mean_annual_rate, 0.025, na.rm = TRUE),
            q250 = quantile(mean_annual_rate, 0.250, na.rm = TRUE),
            q750 = quantile(mean_annual_rate, 0.750, na.rm = TRUE),
            q975 = quantile(mean_annual_rate, 0.975, na.rm = TRUE),
            .groups = "drop") %>% arrange(year)

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
     xlab = "Year", ylab = "Mean annual ΔSOC since first year (tC/ha/yr)",
     main = sprintf("Mean annual ΔSOC  |  %s transient  |  %s", MODEL_NAME, RUN_ID))
polygon(c(traj_summary$year, rev(traj_summary$year)),
        c(traj_summary$q025, rev(traj_summary$q975)),
        col = adjustcolor("steelblue", 0.15), border = NA)
polygon(c(traj_summary$year, rev(traj_summary$year)),
        c(traj_summary$q250, rev(traj_summary$q750)),
        col = adjustcolor("steelblue", 0.35), border = NA)
lines(traj_summary$year, traj_summary$median_delta, col = "steelblue", lwd = 2)
abline(h = 0, lty = 3, col = "grey40")
abline(h = obs_mean_delta, col = "firebrick", lwd = 2.5, lty = 2)
rect(xleft = min(obs_delta_per_plot$year_first),
     xright = max(obs_delta_per_plot$year_last),
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

# --- Holdout: obs vs predicted (mirrors Section 8 scatter, holdout subset) ---
holdout_pred_png <- file.path(DIR_DIAG,
                              sprintf("%s_obs_vs_pred_holdout_%s.png", MODEL_NAME, RUN_ID))
ktp_vals_h <- as.character(res_holdout$kasvup_tyyppi)
ktp_col_h  <- ktp_palette[ktp_vals_h]; ktp_col_h[is.na(ktp_col_h)] <- "grey60"
png(holdout_pred_png, width = 7L * PX_PER_IN, height = 7L * PX_PER_IN, res = PX_PER_IN)
par(mar = c(4.5, 4.5, 3, 1))
ax_lim_h <- c(0, max(c(res_holdout$soc_obs_tCha, res_holdout$soc_q975),
                     na.rm = TRUE) * 1.05)
plot(NA, xlim = ax_lim_h, ylim = ax_lim_h,
     xlab = "Observed SOC (tC/ha)",
     ylab = "Predicted SOC -- posterior median (tC/ha)",
     main = sprintf("%s  |  Independent validation  |  %s", MODEL_NAME, RUN_ID))
abline(0, 1, lty = 2, col = "grey40", lwd = 1.5)
segments(res_holdout$soc_obs_tCha, res_holdout$soc_q025,
         res_holdout$soc_obs_tCha, res_holdout$soc_q975,
         col = adjustcolor(ktp_col_h, 0.3), lwd = 0.8)
points(res_holdout$soc_obs_tCha, res_holdout$soc_median,
       pch = 16, cex = 0.9, col = adjustcolor(ktp_col_h, 0.7))
legend("topleft",
       legend = c(sprintf("R² = %.3f",           metrics_holdout$R2),
                  sprintf("RMSE = %.1f tC/ha",   metrics_holdout$RMSE_median),
                  sprintf("Bias = %+.1f tC/ha",  metrics_holdout$bias_median),
                  sprintf("95%% coverage = %.2f", metrics_holdout$coverage_95),
                  sprintf("n holdout = %d",       length(holdout_plots))),
       bty = "n", cex = 0.85)
ktp_present_h <- sort(unique(na.omit(ktp_vals_h)))
ktp_present_h <- ktp_present_h[ktp_present_h %in% names(ktp_palette)]
if (length(ktp_present_h) > 0)
  legend("bottomright", legend = ktp_labels[ktp_present_h],
         col = ktp_palette[ktp_present_h], pch = 16, bty = "n", cex = 0.75,
         title = "kasvup_tyyppi")
dev.off()
message(sprintf("Holdout obs vs pred: %s", holdout_pred_png))

# --- Holdout: ΔSOC trajectory (mirrors Section 8 trajectory, holdout subset) ---
pp_holdout <- posterior_predictions[
  posterior_predictions$plot_id %in% holdout_plots, ]
pp_h <- pp_holdout %>%
  group_by(plot_id, draw) %>% arrange(year, .by_group = TRUE) %>%
  mutate(first_soc   = first(total_soc), first_year = first(year),
         annual_rate = (total_soc - first_soc) / (year - first_year)) %>%
  filter(year != first_year) %>% ungroup()
traj_summary_h <- pp_h %>%
  group_by(draw, year) %>%
  summarise(mean_annual_rate = mean(annual_rate, na.rm = TRUE), .groups = "drop") %>%
  group_by(year) %>%
  summarise(median_delta = median(mean_annual_rate, na.rm = TRUE),
            q025 = quantile(mean_annual_rate, 0.025, na.rm = TRUE),
            q250 = quantile(mean_annual_rate, 0.250, na.rm = TRUE),
            q750 = quantile(mean_annual_rate, 0.750, na.rm = TRUE),
            q975 = quantile(mean_annual_rate, 0.975, na.rm = TRUE),
            .groups = "drop") %>% arrange(year)
obs_delta_h <- posterior_summary %>%
  filter(!is.na(soc_obs_tCha), plot_id %in% holdout_plots) %>%
  group_by(plot_id) %>% arrange(year, .by_group = TRUE) %>%
  summarise(delta_obs  = (last(soc_obs_tCha) - first(soc_obs_tCha)) /
              (last(year) - first(year)),
            year_first = first(year), year_last = last(year), .groups = "drop") %>%
  filter(year_first != year_last)
obs_mean_delta_h <- mean(obs_delta_h$delta_obs, na.rm = TRUE)
obs_se_delta_h   <- sd(obs_delta_h$delta_obs,   na.rm = TRUE) / sqrt(nrow(obs_delta_h))
traj_holdout_png <- file.path(DIR_DIAG,
                              sprintf("%s_delta_soc_trajectory_holdout_%s.png", MODEL_NAME, RUN_ID))
png(traj_holdout_png, width = 12L * PX_PER_IN, height = 7L * PX_PER_IN, res = PX_PER_IN)
par(mar = c(4, 4, 3, 1))
ylim_h <- range(c(traj_summary_h$q025, traj_summary_h$q975,
                  obs_mean_delta_h + 1.96*obs_se_delta_h,
                  obs_mean_delta_h - 1.96*obs_se_delta_h, 0), na.rm = TRUE)
ylim_h <- ylim_h + c(-0.05, 0.05) * diff(ylim_h)
plot(NA, xlim = range(traj_summary_h$year), ylim = ylim_h,
     xlab = "Year", ylab = "Mean annual ΔSOC since first year (tC/ha/yr)",
     main = sprintf("Mean annual ΔSOC — independent validation  |  %s  |  %s",
                    MODEL_NAME, RUN_ID))
polygon(c(traj_summary_h$year, rev(traj_summary_h$year)),
        c(traj_summary_h$q025, rev(traj_summary_h$q975)),
        col = adjustcolor("steelblue", 0.15), border = NA)
polygon(c(traj_summary_h$year, rev(traj_summary_h$year)),
        c(traj_summary_h$q250, rev(traj_summary_h$q750)),
        col = adjustcolor("steelblue", 0.35), border = NA)
lines(traj_summary_h$year, traj_summary_h$median_delta, col = "steelblue", lwd = 2)
abline(h = 0, lty = 3, col = "grey40")
abline(h = obs_mean_delta_h, col = "firebrick", lwd = 2.5, lty = 2)
rect(xleft   = min(obs_delta_h$year_first), xright  = max(obs_delta_h$year_last),
     ybottom = obs_mean_delta_h - 1.96*obs_se_delta_h,
     ytop    = obs_mean_delta_h + 1.96*obs_se_delta_h,
     col = adjustcolor("firebrick", 0.12), border = NA)
legend("topleft",
       legend = c("Posterior median", "50% CI", "95% CI",
                  sprintf("Holdout observed ΔSOC ± 95%% CI  (n = %d)",
                          nrow(obs_delta_h))),
       col    = c("steelblue", adjustcolor("steelblue",0.35),
                  adjustcolor("steelblue",0.15), "firebrick"),
       lwd    = c(2,NA,NA,2.5), pch = c(NA,15,15,NA), lty = c(1,NA,NA,2),
       pt.cex = c(NA,2,2,NA), bty = "n", cex = 0.9)
dev.off()
message(sprintf("Holdout ΔSOC trajectory: %s", traj_holdout_png))

# --- Updated bundle: adds holdout slots alongside existing calibration metrics ---
pp_output <- list(
  posterior_predictions  = posterior_predictions,
  projection_predictions = projection_predictions,
  posterior_summary      = as.data.frame(posterior_summary),
  residuals_df           = residuals_df,       # now includes is_holdout column
  holdout_plots          = holdout_plots,
  metrics         = list(                      # calibration — kept for backward compat
    R2 = R2, RMSE_mean = RMSE, bias_mean = bias,
    RMSE_median = RMSE_med, bias_median = bias_med, coverage_95 = coverage_95),
  metrics_calib   = metrics_calib,
  metrics_holdout = metrics_holdout,
  config  = list(model = MODEL_NAME, run_id = RUN_ID,
                 n_draws = N_PP_DRAWS, n_plots_total = length(plots_real),
                 n_calib_plots   = length(plots_real) - length(holdout_plots),
                 n_holdout_plots = length(holdout_plots),
                 proj_years = PROJ_YEARS, recycle_years = RECYCLE_YEARS,
                 transient_init = TRUE, timestamp = Sys.time())
)
saveRDS(pp_output,
        file.path(DIR_RUNS,
                  sprintf("%s_posterior_predictive_%s.rds", MODEL_NAME, RUN_ID)))
message("Saved posterior predictive bundle.")

append_to_report(run_config, paste(c(
  "\n[5] Predictive performance (transient init)\n",
  sprintf("    R²: %.3f  RMSE: %.2f  Bias: %+.2f  Cov95: %.3f\n",
          R2, RMSE_med, bias_med, coverage_95),
  sprintf("    Projection years:       %d\n",       PROJ_YEARS),
  sprintf("    Climate recycle window: %d years\n", RECYCLE_YEARS)), collapse=""))

message("\nDone. Next: run_residual_analysis.R Yasso20 ", RUN_ID)
