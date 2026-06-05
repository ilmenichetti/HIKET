# =============================================================================
# run_yasso20_baseline.R
#
# Deterministic forward run of Yasso20 at published default parameters.
# Reference baseline for comparison with calibrated posterior predictive outputs.
#
# Key differences from run_yasso15_baseline.R:
#   - Monthly climate input (map_climate_monthly, 12 rows per year).
#   - compute_xi_yasso20() averages over 12 monthly temperatures directly
#     rather than the four-point Gaussian approximation used in Yasso15.
#   - yasso20_run(climate_df, inputs_df, params, C_init) takes separate monthly
#     climate and annual inputs; computes xi and annual precip internally.
#   - precip_mean in litter_means is the annual total from mean monthly values:
#     sum(tapply(precip, month, mean)).
#   - Parameter names and values are identical to Yasso15 (same Fortran).
#
# Initialization: steady-state (sigma_init = 1.0 has no effect).
#
# OUTPUTS (diagnostics/Yasso20_baseline/):
#   Yasso20_baseline_predictions.rds
#   Yasso20_baseline_residuals.csv
#   Yasso20_baseline_plot_ranking.csv
#   Yasso20_baseline_mean_soc_trajectory.png
#   Yasso20_baseline_delta_soc_trajectory.png
#   Yasso20_baseline_obs_vs_pred.png
#   Yasso20_baseline_residuals_diagnostic.png
#   Yasso20_baseline_rf_importance.png
#   Yasso20_baseline_sampled_plots.png
#   Yasso20_baseline_worst25_timeseries.png
#   Yasso20_baseline_residuals_histogram.png
#   Yasso20_baseline_report.txt
# =============================================================================

source("./Model_functions_real_data_transient/input_compatibility_layer.R")
source("./Model_functions_real_data_transient/Decomposition_functions/Yasso/yasso15_wrapper_transient.R")
source("./Model_functions_real_data_transient/Decomposition_functions/Yasso/yasso20_wrapper_transient.R")
source("./Prior_specs/Yasso20_priors.R")
dyn.load("./Model_functions_real_data_transient/Decomposition_functions/Yasso/yasso15.so")

suppressPackageStartupMessages({
  library(dplyr)
  library(ranger)
  library(parallel)
})


# =============================================================================
# 0.  Configuration
# =============================================================================

MODEL_NAME <- "Yasso20"
DIR_OUT    <- "./Calibration_real_data_transient/diagnostics/Yasso20_baseline"
dir.create(DIR_OUT, showWarnings = FALSE, recursive = TRUE)

STEADY_STATE_YEARS <- 20L
N_CORES   <- parallel::detectCores() - 1L
PX_PER_IN <- 150L

cat("=============================================================\n")
cat(sprintf("  %s baseline (published defaults, no calibration)\n", MODEL_NAME))
cat(sprintf("  Output: %s\n", DIR_OUT))
cat("=============================================================\n\n")


# =============================================================================
# 1.  Parameters at published defaults
# =============================================================================
# YASSO20_PARAM_NAMES and YASSO20_DEFAULT_PARAMS are identical to Yasso15
# (defined in yasso20_wrapper_transient.R as aliases).

FIXED_RATE_NAMES <- c("alpha_A","alpha_W","alpha_E","alpha_N",
                      "p_H","alpha_H","w1","w2","w3","w4","w5")
fixed_rates      <- YASSO20_DEFAULT_PARAMS[FIXED_RATE_NAMES]
MODEL_FREE_NAMES <- setdiff(names(YASSO20_FREE_DEFAULTS), c("sigma_init","sigma_input"))
free_defaults    <- YASSO20_FREE_DEFAULTS   # sigma_init=1, sigma_input=1 by construction

assemble_model_params <- function(p_free) {
  yasso_params <- c(fixed_rates, p_free[MODEL_FREE_NAMES])
  yasso_params <- yasso_params[YASSO20_PARAM_NAMES]
  c(yasso_params, sigma_input = unname(p_free["sigma_input"]))
}

model_params <- assemble_model_params(free_defaults)

cat("Parameters at published defaults:\n")
cat(sprintf("  sigma_input: %.2f\n", model_params["sigma_input"]))
cat(sprintf("  AWE: beta1=%.5f  beta2=%.7f  gamma=%.5f\n",
            model_params["beta1"], model_params["beta2"], model_params["gamma"]))
cat(sprintf("  N:   betaN1=%.5f  betaN2=%.7f  gammaN=%.5f\n",
            model_params["betaN1"], model_params["betaN2"], model_params["gammaN"]))
cat(sprintf("  H:   betaH1=%.5f  betaH2=%.7f  gammaH=%.5f\n\n",
            model_params["betaH1"], model_params["betaH2"], model_params["gammaH"]))


# =============================================================================
# 2.  Load and filter data
# =============================================================================

cat("Loading data...\n")
input_raw <- read.csv("./Data/model_inputs/input_raw_monthly.csv")
site_raw  <- read.csv("./Data/model_inputs/site_raw.csv")
site_raw$plot_id <- as.character(site_raw$plot_id)

calib_plots <- as.character(site_raw$plot_id[site_raw$calib_ready])
input_calib <- input_raw[as.character(input_raw$plot_id) %in% calib_plots, ]

litter_cols <- grep("^C_(nwl|fwl|cwl)_[AWEN]$", names(input_calib), value = TRUE)
bad_clim    <- unique(input_calib$plot_id[is.na(input_calib$temp_air)])
bad_litter  <- unique(input_calib$plot_id[!complete.cases(input_calib[, litter_cols])])
calib_plots <- setdiff(calib_plots, c(bad_clim, bad_litter))
input_calib <- input_calib[as.character(input_calib$plot_id) %in% calib_plots, ]

# Yasso20 uses monthly climate (12 rows per year) and its own litter mapping
Yasso20_monthly <- map_climate_monthly(input_calib)
Yasso20_inputs  <- map_inputs_yasso20(input_calib)

plots_real <- as.character(Reduce(intersect, list(
  unique(as.character(Yasso20_monthly$plot_id)),
  unique(as.character(Yasso20_inputs$plot_id))
)))
cat(sprintf("Plots with valid data: %d\n", length(plots_real)))

climate_by_plot <- split(Yasso20_monthly, Yasso20_monthly$plot_id)
inputs_by_plot  <- split(Yasso20_inputs,  Yasso20_inputs$plot_id)

# precip_mean: annual total from mean monthly values over the SS window.
# Passed to yasso20_steady_state() for the Fortran leaching term.
litter_means <- lapply(plots_real, function(pid) {
  inp  <- Yasso20_inputs[as.character(Yasso20_inputs$plot_id) == pid, ]
  clim <- Yasso20_monthly[as.character(Yasso20_monthly$plot_id) == pid, ]
  n_ss_inp   <- min(STEADY_STATE_YEARS, nrow(inp))
  n_ss_months <- min(STEADY_STATE_YEARS * 12L, nrow(clim))
  inp_ss  <- inp[seq_len(n_ss_inp), ]
  clim_ss <- clim[seq_len(n_ss_months), ]
  list(
    nwl_mean    = colMeans(inp_ss[, c("nwl_A","nwl_W","nwl_E","nwl_N")]),
    fwl_mean    = colMeans(inp_ss[, c("fwl_A","fwl_W","fwl_E","fwl_N")]),
    cwl_mean    = colMeans(inp_ss[, c("cwl_A","cwl_W","cwl_E","cwl_N")]),
    precip_mean = sum(tapply(clim_ss$precip, clim_ss$month, mean, na.rm = TRUE))
  )
})
names(litter_means) <- plots_real

SOC_obs_all <- input_calib %>%
  filter(!is.na(soc_obs_tCha)) %>%
  select(plot_id, year, soc_obs_tCha) %>%
  mutate(plot_id = as.character(plot_id))


# =============================================================================
# 3.  Deterministic forward run per plot
# =============================================================================

LITTER_COLS <- c("nwl_A","nwl_W","nwl_E","nwl_N",
                 "fwl_A","fwl_W","fwl_E","fwl_N",
                 "cwl_A","cwl_W","cwl_E","cwl_N")

cat(sprintf("\nRunning %s at defaults across %d plots...\n", MODEL_NAME, length(plots_real)))
t0 <- Sys.time()

predictions <- do.call(rbind, mclapply(plots_real, function(pid) {
  clim   <- climate_by_plot[[pid]]   # monthly: year, month, temp_air, precip
  inputs <- inputs_by_plot[[pid]]    # annual:  year, nwl/fwl/cwl AWEN
  lm     <- litter_means[[pid]]

  n_ss_months <- min(STEADY_STATE_YEARS * 12L, nrow(clim))
  xi_for_ss <- tryCatch(
    compute_xi_mean_yasso20(clim[seq_len(n_ss_months), , drop = FALSE],
                            model_params[YASSO20_PARAM_NAMES]),
    error = function(e) NULL)
  if (is.null(xi_for_ss) ||
      !all(is.finite(c(xi_for_ss$xi_awe, xi_for_ss$xi_n, xi_for_ss$xi_h))) ||
      any(c(xi_for_ss$xi_awe, xi_for_ss$xi_n, xi_for_ss$xi_h) <= 0)) return(NULL)

  C_init <- tryCatch(
    yasso20_steady_state(
      params      = model_params[YASSO20_PARAM_NAMES],
      nwl_mean    = lm$nwl_mean,
      fwl_mean    = lm$fwl_mean,
      cwl_mean    = lm$cwl_mean,
      xi_ss       = xi_for_ss,
      precip_mean = lm$precip_mean),
    error = function(e) NULL)
  if (is.null(C_init) || any(!is.finite(C_init)) || any(C_init < 0)) return(NULL)

  si <- model_params["sigma_input"]   # 1.0 at defaults
  inputs_scaled <- inputs
  inputs_scaled[, LITTER_COLS] <- inputs[, LITTER_COLS] * si

  # yasso20_run computes xi and annual precip from monthly climate internally
  run_out <- tryCatch(
    yasso20_run(clim, inputs_scaled, model_params[YASSO20_PARAM_NAMES], C_init),
    error = function(e) NULL)
  if (is.null(run_out)) return(NULL)

  data.frame(
    plot_id = pid, year = run_out$year,
    A = run_out$A, W = run_out$W, E = run_out$E, N = run_out$N, H = run_out$H,
    total_soc = run_out$total_soc, respiration = run_out$respiration
  )
}, mc.cores = N_CORES))

cat(sprintf("Run complete: %.1f sec  (%d plot-years)\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs")), nrow(predictions)))
saveRDS(predictions, file.path(DIR_OUT, "Yasso20_baseline_predictions.rds"))


# =============================================================================
# 4.  Residuals and metrics
# =============================================================================

residuals_df <- predictions %>%
  inner_join(SOC_obs_all, by = c("plot_id", "year")) %>%
  mutate(
    log_obs      = log(soc_obs_tCha),
    log_hat      = log(total_soc),
    residual_log = log_obs - log_hat,
    residual_abs = soc_obs_tCha - total_soc
  ) %>%
  select(plot_id, year, soc_obs_tCha, soc_pred = total_soc,
         log_obs, log_hat, residual_log, residual_abs) %>%
  left_join(site_raw, by = "plot_id")

write.csv(residuals_df, file.path(DIR_OUT, "Yasso20_baseline_residuals.csv"),
          row.names = FALSE)
cat(sprintf("Residuals: %d obs\n", nrow(residuals_df)))

obs    <- residuals_df$soc_obs_tCha
hat    <- residuals_df$soc_pred
finite <- is.finite(obs) & is.finite(hat)
R2   <- cor(obs[finite], hat[finite])^2
RMSE <- sqrt(mean((obs[finite] - hat[finite])^2))
bias <- mean(hat[finite] - obs[finite])
mape <- mean(abs(hat[finite] - obs[finite]) / obs[finite]) * 100

cat(sprintf("\nPredictive metrics at published defaults:\n"))
cat(sprintf("  R²: %.3f  RMSE: %.2f  Bias: %+.2f  MAPE: %.1f%%\n", R2, RMSE, bias, mape))


# =============================================================================
# 5.  Plot: Mean absolute SOC trajectory
# =============================================================================

soc_traj <- predictions %>%
  group_by(year) %>%
  summarise(mean_soc = mean(total_soc, na.rm = TRUE),
            se_soc   = sd(total_soc,   na.rm = TRUE) / sqrt(sum(!is.na(total_soc))),
            .groups  = "drop") %>%
  arrange(year)

obs_by_year <- SOC_obs_all %>%
  group_by(year) %>%
  summarise(mean_obs = mean(soc_obs_tCha, na.rm = TRUE),
            se_obs   = sd(soc_obs_tCha,  na.rm = TRUE) / sqrt(n()),
            n_obs    = n(), .groups = "drop")

ylim_soc <- range(c(soc_traj$mean_soc - soc_traj$se_soc,
                    soc_traj$mean_soc + soc_traj$se_soc,
                    obs_by_year$mean_obs - 1.96 * obs_by_year$se_obs,
                    obs_by_year$mean_obs + 1.96 * obs_by_year$se_obs), na.rm = TRUE)
ylim_soc <- ylim_soc + c(-0.05, 0.05) * diff(ylim_soc)

png(file.path(DIR_OUT, "Yasso20_baseline_mean_soc_trajectory.png"),
    width = 12L * PX_PER_IN, height = 6L * PX_PER_IN, res = PX_PER_IN)
par(mar = c(4, 4, 3, 1))
plot(NA, xlim = range(soc_traj$year), ylim = ylim_soc,
     xlab = "Year", ylab = "Mean SOC across plots (tC/ha)",
     main = sprintf("%s baseline — mean absolute SOC  |  %d plots",
                    MODEL_NAME, length(plots_real)))
polygon(c(soc_traj$year, rev(soc_traj$year)),
        c(soc_traj$mean_soc - soc_traj$se_soc, rev(soc_traj$mean_soc + soc_traj$se_soc)),
        col = adjustcolor("steelblue", 0.20), border = NA)
lines(soc_traj$year, soc_traj$mean_soc, col = "steelblue", lwd = 2.5)
arrows(obs_by_year$year, obs_by_year$mean_obs - 1.96 * obs_by_year$se_obs,
       obs_by_year$year, obs_by_year$mean_obs + 1.96 * obs_by_year$se_obs,
       code = 3, angle = 90, length = 0.05, col = "firebrick", lwd = 2)
points(obs_by_year$year, obs_by_year$mean_obs, pch = 19, cex = 1.3, col = "firebrick")
legend("topleft",
       legend = c("Predicted mean SOC ± 1 SE (cross-plot variability)",
                  "Observed campaign mean ± 95% CI"),
       col = c("steelblue","firebrick"), lwd = c(2.5, 2), pch = c(NA, 19),
       lty = c(1, NA), bty = "n", cex = 0.9)
dev.off()


# =============================================================================
# 6.  Plot: Mean annual ΔSOC trajectory
# =============================================================================

pp <- predictions %>%
  group_by(plot_id) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(first_soc   = first(total_soc),
         first_year  = first(year),
         annual_rate = (total_soc - first_soc) / (year - first_year)) %>%
  filter(year != first_year) %>%
  ungroup()

traj <- pp %>%
  group_by(year) %>%
  summarise(mean_delta = mean(annual_rate, na.rm = TRUE), .groups = "drop") %>%
  arrange(year)

obs_delta_per_plot <- residuals_df %>%
  group_by(plot_id) %>%
  arrange(year, .by_group = TRUE) %>%
  summarise(delta_obs  = (last(soc_obs_tCha) - first(soc_obs_tCha)) /
                         (last(year) - first(year)),
            year_first = first(year), year_last = last(year), .groups = "drop") %>%
  filter(year_first != year_last)

obs_mean_delta <- mean(obs_delta_per_plot$delta_obs, na.rm = TRUE)
obs_se_delta   <- sd(obs_delta_per_plot$delta_obs, na.rm = TRUE) /
                  sqrt(nrow(obs_delta_per_plot))

ylim_d <- range(c(traj$mean_delta,
                  obs_mean_delta + 1.96 * obs_se_delta,
                  obs_mean_delta - 1.96 * obs_se_delta, 0), na.rm = TRUE)
ylim_d <- ylim_d + c(-0.05, 0.05) * diff(ylim_d)

png(file.path(DIR_OUT, "Yasso20_baseline_delta_soc_trajectory.png"),
    width = 12L * PX_PER_IN, height = 6L * PX_PER_IN, res = PX_PER_IN)
par(mar = c(4, 4, 3, 1))
plot(traj$year, traj$mean_delta, type = "l", col = "steelblue", lwd = 2.5,
     ylim = ylim_d, xlab = "Year",
     ylab = "Mean annual ΔSOC since first year (tC/ha/yr)",
     main = sprintf("%s baseline — mean annual ΔSOC  |  %d plots",
                    MODEL_NAME, length(plots_real)))
abline(h = 0, lty = 3, col = "grey40")
abline(h = obs_mean_delta, col = "firebrick", lwd = 2.5, lty = 2)
rect(xleft  = min(obs_delta_per_plot$year_first),
     xright = max(obs_delta_per_plot$year_last),
     ybottom = obs_mean_delta - 1.96 * obs_se_delta,
     ytop    = obs_mean_delta + 1.96 * obs_se_delta,
     col = adjustcolor("firebrick", 0.12), border = NA)
legend("topleft",
       legend = c("Predicted mean annual ΔSOC",
                  sprintf("Observed mean annual ΔSOC ± 95%% CI  (n = %d plots)",
                          nrow(obs_delta_per_plot))),
       col = c("steelblue","firebrick"), lwd = c(2.5, 2.5), lty = c(1, 2),
       bty = "n", cex = 0.9)
dev.off()


# =============================================================================
# 7.  Plot: Observed vs predicted
# =============================================================================

ktp_palette <- c("1"="#000000","2"="#E69F00","3"="#56B4E9","4"="#009E73",
                 "5"="#F0E442","6"="#0072B2","7"="#D55E00","8"="#CC79A7")
ktp_labels  <- c("1"="1 Lehdot (OMaT)","2"="2 Lehtomainen (OMT)",
                 "3"="3 Tuorekangas (MT)","4"="4 Kuivahko (VT)",
                 "5"="5 Kuiva (CT)","6"="6 Karukkokangas (ClT)",
                 "7"="7 Kalliot/hietikot","8"="8 Lakimetsät/tunturit")
ktp_vals <- as.character(residuals_df$kasvup_tyyppi)
ktp_col  <- ktp_palette[ktp_vals]; ktp_col[is.na(ktp_col)] <- "grey60"

png(file.path(DIR_OUT, "Yasso20_baseline_obs_vs_pred.png"),
    width = 7L * PX_PER_IN, height = 7L * PX_PER_IN, res = PX_PER_IN)
par(mar = c(4.5, 4.5, 3, 1))
ax_lim <- c(0, max(c(obs, hat), na.rm = TRUE) * 1.05)
plot(NA, xlim = ax_lim, ylim = ax_lim,
     xlab = "Observed SOC (tC/ha)", ylab = "Predicted SOC (tC/ha)",
     main = sprintf("%s at published defaults", MODEL_NAME))
abline(0, 1, lty = 2, col = "grey40", lwd = 1.5)
points(obs, hat, pch = 16, cex = 0.9, col = adjustcolor(ktp_col, 0.7))
legend("topleft",
       legend = c(sprintf("R² = %.3f", R2), sprintf("RMSE = %.1f tC/ha", RMSE),
                  sprintf("Bias = %+.1f tC/ha", bias), sprintf("MAPE = %.0f%%", mape)),
       bty = "n", cex = 0.85)
ktp_present <- sort(unique(na.omit(ktp_vals)))
ktp_present <- ktp_present[ktp_present %in% names(ktp_palette)]
if (length(ktp_present) > 0)
  legend("bottomright", legend = ktp_labels[ktp_present],
         col = ktp_palette[ktp_present], pch = 16, bty = "n", cex = 0.8,
         title = "kasvup_tyyppi")
dev.off()


# =============================================================================
# 8.  Plot: Residuals diagnostic (6 panels)
# =============================================================================

resid      <- residuals_df$residual_log
fitted_log <- residuals_df$log_hat
mask       <- is.finite(resid) & is.finite(fitted_log)
resid_f    <- resid[mask]; fitted_f <- fitted_log[mask]
resid_df_f <- residuals_df[mask, , drop = FALSE]

png(file.path(DIR_OUT, "Yasso20_baseline_residuals_diagnostic.png"),
    width = 12L * PX_PER_IN, height = 8L * PX_PER_IN, res = PX_PER_IN)
par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))

plot(fitted_f, resid_f, pch = 19, col = adjustcolor("steelblue", 0.5),
     xlab = "log(SOC predicted)", ylab = "Residual (log)", main = "Residuals vs fitted")
abline(h = 0, lty = 2, col = "grey40")
if (length(fitted_f) > 10) {
  lo <- loess(resid_f ~ fitted_f); ord <- order(fitted_f)
  lines(fitted_f[ord], predict(lo)[ord], col = "firebrick", lwd = 2)
}

boxplot(resid_f ~ resid_df_f$year, col = "lightsteelblue",
        xlab = "Year", ylab = "Residual (log)", main = "Residuals by year")
abline(h = 0, lty = 2, col = "grey40")

if ("KA" %in% names(resid_df_f)) {
  boxplot(resid_f ~ as.factor(resid_df_f$KA), col = "lightsteelblue",
          xlab = "Cajander fertility class (KA)", ylab = "Residual (log)",
          main = "Residuals by Cajander class")
  abline(h = 0, lty = 2, col = "grey40")
} else { plot.new(); title("KA not available") }

if ("CoarseFragments" %in% names(resid_df_f)) {
  stone <- resid_df_f$CoarseFragments; stone_mask <- is.finite(stone)
  plot(stone[stone_mask], resid_f[stone_mask], pch = 19,
       col = adjustcolor("steelblue", 0.45),
       xlab = "Coarse fragments (%)", ylab = "Residual (log)",
       main = sprintf("Residuals vs stoniness  (n = %d)", sum(stone_mask)))
  abline(h = 0, lty = 2, col = "grey40")
  if (sum(stone_mask) > 20) {
    lo_s <- loess(resid_f[stone_mask] ~ stone[stone_mask], span = 0.75)
    ord_s <- order(stone[stone_mask])
    lines(stone[stone_mask][ord_s], predict(lo_s)[ord_s], col = "firebrick", lwd = 2)
  }
  rug(stone[stone_mask], col = adjustcolor("steelblue", 0.3), ticksize = 0.02)
} else { plot.new(); title("CoarseFragments not available") }

if ("species_code" %in% names(resid_df_f)) {
  sp_labels <- c("1"="Scots pine","2"="Nor. spruce","3"="Birch","4"="Aspen",
                 "5"="Grey alder","6"="Blk. alder","7"="Other brdlv.")
  sp_factor <- as.factor(resid_df_f$species_code)
  levels(sp_factor) <- sp_labels[levels(sp_factor)]
  boxplot(resid_f ~ sp_factor, col = "lightsteelblue", las = 2,
          xlab = "", ylab = "Residual (log)",
          main = "Residuals by dominant species", cex.axis = 0.78)
  abline(h = 0, lty = 2, col = "grey40")
} else { plot.new(); title("species_code not available") }

qqnorm(resid_f, pch = 19, col = adjustcolor("steelblue", 0.5),
       main = "Q-Q plot of log residuals")
qqline(resid_f, col = "firebrick", lwd = 2)
mtext(sprintf("%s baseline — residual diagnostics", MODEL_NAME),
      side = 3, outer = TRUE, line = -1.5, cex = 1.0, font = 2)
dev.off()


# =============================================================================
# 9.  Random Forest on residuals
# =============================================================================

target <- "residual_log"
candidate_vars <- c(
  "KA","kasvyo_syke","kasvyo_ahti","species_code","soil_code",
  "TexturalClass","temp_zone","koppen_class","ojitustilanne",
  "dev_class_85","kasvup_tyyppi","alaryhma",
  "clay","SiltContent","SandContent","CoarseFragments",
  "EstimatedBulkDensity","profile_depth_cm","ofh_lower_cm","ofh_weight_kgm2",
  "pH.CaCl2.","pH.H2O.","TotalNitrogen","CN_ratio","CEC","base_saturation",
  "ExchangeableAl","ExchangeableFe","ExchangeableCa",
  "mean_temp","mean_precip_annual","GDD5","coldest_month_T",
  "warmest_month_T","T_seasonality","P_seasonality","aridity_index",
  "temp_sum_NFI","elevation_m","lat_WGS84","lon_WGS84",
  "stand_age_85","stand_age_2006_est","basal_area_85","mean_height_85_dm",
  "any_cut_85_95","n_cuts_85_95","any_trt_85_95","soil_prep_pre85",
  "litter_A_frac","litter_W_frac","litter_E_frac","litter_N_frac",
  "woody_share","conifer_share","shallow","n_soc_obs"
)
available <- intersect(candidate_vars, names(residuals_df))
rf_df <- residuals_df[, c(target, available), drop = FALSE]
for (v in available) {
  x <- rf_df[[v]]
  if (is.character(x) || is.logical(x)) rf_df[[v]] <- as.factor(x)
  if (is.integer(x) && length(unique(na.omit(x))) <= 12)
    rf_df[[v]] <- as.factor(as.character(x))
}
keep_vars <- vapply(available, function(v) {
  vec <- rf_df[[v]]
  if (all(is.na(vec))) return(FALSE)
  if (is.factor(vec)) return(nlevels(droplevels(vec)) > 1)
  length(unique(vec[!is.na(vec)])) > 1
}, logical(1))
available <- available[keep_vars]
rf_df <- rf_df[, c(target, available), drop = FALSE]
rf_df <- rf_df[complete.cases(rf_df) & is.finite(rf_df[[target]]), , drop = FALSE]
cat(sprintf("\nRF data: %d obs x %d covariates\n", nrow(rf_df), length(available)))

set.seed(2025)
rf_fit <- ranger(formula = as.formula(paste(target, "~ .")), data = rf_df,
                 num.trees = 500L, importance = "permutation", num.threads = N_CORES)
oob_r2 <- 1 - rf_fit$prediction.error / var(rf_df[[target]])
imp    <- sort(rf_fit$variable.importance, decreasing = TRUE)
imp_df <- data.frame(variable = names(imp), importance = unname(imp))
cat(sprintf("RF OOB R²: %.3f\n", oob_r2))

n_show  <- min(20L, nrow(imp_df))
imp_top <- imp_df[seq_len(n_show), ]; imp_top <- imp_top[order(imp_top$importance), ]
png(file.path(DIR_OUT, "Yasso20_baseline_rf_importance.png"),
    width = 10L * PX_PER_IN, height = max(5L, n_show * 0.3) * PX_PER_IN, res = PX_PER_IN)
par(mar = c(4, 12, 3, 1))
barplot(imp_top$importance, horiz = TRUE, names.arg = imp_top$variable,
        las = 1, cex.names = 0.85, col = "steelblue", border = "white",
        xlab = "Permutation importance",
        main = sprintf("%s baseline — top %d covariates  |  OOB R² = %.3f",
                       MODEL_NAME, n_show, oob_r2))
dev.off()


# =============================================================================
# 10.  16 sampled plots
# =============================================================================

set.seed(2025)
SAMPLE_IDs <- sample(unique(predictions$plot_id), min(16L, length(unique(predictions$plot_id))))
plot_range <- range(c(
  predictions[predictions$plot_id %in% SAMPLE_IDs, "total_soc"],
  SOC_obs_all[as.character(SOC_obs_all$plot_id) %in% SAMPLE_IDs, "soc_obs_tCha"]),
  na.rm = TRUE)

png(file.path(DIR_OUT, "Yasso20_baseline_sampled_plots.png"),
    width = 12L * PX_PER_IN, height = 9L * PX_PER_IN, res = PX_PER_IN)
par(mfrow = c(4, 4), mar = c(3, 3, 2, 1))
for (pid in SAMPLE_IDs) {
  pred_p <- predictions[predictions$plot_id == pid, ]
  obs_p  <- SOC_obs_all[as.character(SOC_obs_all$plot_id) == pid, ]
  plot(pred_p$year, pred_p$total_soc, type = "l", col = "steelblue", lwd = 2,
       ylim = plot_range, xlab = "", ylab = "SOC (tC/ha)",
       main = sprintf("Plot %s", pid), cex.main = 0.85)
  points(obs_p$year, obs_p$soc_obs_tCha, pch = 16, col = "firebrick", cex = 1.3)
}
dev.off()


# =============================================================================
# 11.  Plot ranking and worst-25 panel
# =============================================================================

plot_resid_summary <- residuals_df %>%
  filter(is.finite(residual_log)) %>%
  group_by(plot_id) %>%
  summarise(n_obs          = n(),
            mean_resid     = mean(residual_log),
            abs_resid      = mean(abs(residual_log)),
            mean_obs       = mean(soc_obs_tCha),
            mean_pred      = mean(soc_pred),
            ratio_obs_pred = mean(soc_obs_tCha / soc_pred),
            .groups        = "drop") %>%
  arrange(desc(abs_resid))

write.csv(plot_resid_summary,
          file.path(DIR_OUT, "Yasso20_baseline_plot_ranking.csv"), row.names = FALSE)

worst_ids  <- head(plot_resid_summary$plot_id, 25)
worst_pred <- predictions[predictions$plot_id %in% worst_ids, ]
worst_obs  <- SOC_obs_all[as.character(SOC_obs_all$plot_id) %in% worst_ids, ]
y_range    <- range(c(worst_pred$total_soc, worst_obs$soc_obs_tCha), na.rm = TRUE)

png(file.path(DIR_OUT, "Yasso20_baseline_worst25_timeseries.png"),
    width = 15L * PX_PER_IN, height = 15L * PX_PER_IN, res = PX_PER_IN)
par(mfrow = c(5, 5), mar = c(3, 3, 2, 1), oma = c(0, 0, 3, 0))
for (pid in worst_ids) {
  pred_p <- worst_pred[worst_pred$plot_id == pid, ]
  obs_p  <- worst_obs[as.character(worst_obs$plot_id) == pid, ]
  sr     <- plot_resid_summary[plot_resid_summary$plot_id == pid, ]
  plot(pred_p$year, pred_p$total_soc, type = "l", lwd = 2, col = "steelblue",
       ylim = y_range, xlab = "", ylab = "SOC (tC/ha)",
       main = sprintf("Plot %s\nbias=%.2f  ratio=%.2f",
                      pid, sr$mean_resid, sr$ratio_obs_pred), cex.main = 0.75)
  points(obs_p$year, obs_p$soc_obs_tCha, pch = 16, col = "firebrick", cex = 1.4)
}
mtext(sprintf("%s defaults — 25 worst plots (by mean |log-residual|)", MODEL_NAME),
      side = 3, outer = TRUE, cex = 1.1, font = 2)
dev.off()


# =============================================================================
# 12.  Residual histogram
# =============================================================================

resid_finite <- residuals_df$residual_log[is.finite(residuals_df$residual_log)]
png(file.path(DIR_OUT, "Yasso20_baseline_residuals_histogram.png"),
    width = 8L * PX_PER_IN, height = 6L * PX_PER_IN, res = PX_PER_IN)
par(mar = c(4.5, 4.5, 3, 1))
hist(resid_finite, breaks = 40, col = "steelblue", border = "white",
     xlab = "Residual  log(obs) - log(pred)", ylab = "Number of observations",
     main = sprintf("%s defaults — residual distribution  (n = %d)",
                    MODEL_NAME, length(resid_finite)))
abline(v = 0,                    lty = 2, lwd = 2, col = "grey30")
abline(v = mean(resid_finite),   lty = 1, lwd = 2, col = "firebrick")
abline(v = median(resid_finite), lty = 1, lwd = 2, col = "darkorange")
legend("topleft",
       legend = c("zero", sprintf("mean = %.3f", mean(resid_finite)),
                  sprintf("median = %.3f", median(resid_finite))),
       col = c("grey30","firebrick","darkorange"), lty = c(2,1,1), lwd = 2, bty = "n")
dev.off()


# =============================================================================
# 13.  Text report
# =============================================================================

report_path <- file.path(DIR_OUT, "Yasso20_baseline_report.txt")
sink(report_path)
cat("=========================================================\n")
cat(sprintf("  %s Published-Defaults Baseline Report\n", MODEL_NAME))
cat(sprintf("  Generated: %s\n", format(Sys.time())))
cat("=========================================================\n\n")
cat("[1] Configuration\n")
cat(sprintf("    Plots simulated:           %d\n", length(plots_real)))
cat(sprintf("    Plot-year predictions:     %d\n", nrow(predictions)))
cat(sprintf("    Observations matched:      %d\n", nrow(residuals_df)))
cat(sprintf("    Steady-state window (yr):  %d\n", STEADY_STATE_YEARS))
cat(sprintf("    Climate:                   monthly (12 rows/yr)\n"))
cat(sprintf("    sigma_input:               1.00 (no scaling)\n\n"))
cat("[2] Predictive performance at defaults\n")
cat(sprintf("    R²:    %.3f\n", R2))
cat(sprintf("    RMSE:  %.2f tC/ha\n", RMSE))
cat(sprintf("    Bias:  %+.2f tC/ha\n", bias))
cat(sprintf("    MAPE:  %.1f%%\n\n", mape))
cat("[3] Bias by stratum (median residual_log; positive => obs > pred)\n")
for (strat in c("kasvup_tyyppi","KA","species_code","temp_zone")) {
  if (strat %in% names(residuals_df)) {
    bs <- residuals_df %>%
      filter(is.finite(residual_log)) %>%
      group_by(.data[[strat]]) %>%
      summarise(n = n(), med_resid = median(residual_log), .groups = "drop")
    cat(sprintf("    %s:\n", strat))
    for (i in seq_len(nrow(bs)))
      cat(sprintf("        %-10s  n=%-4d  median_resid_log = %+.3f\n",
                  as.character(bs[[strat]][i]), bs$n[i], bs$med_resid[i]))
  }
}
cat("\n[4] Residual analysis (Random Forest)\n")
cat(sprintf("    Observations:          %d\n", nrow(rf_df)))
cat(sprintf("    Covariates:            %d\n", length(available)))
cat(sprintf("    RF OOB R²:             %.3f\n", oob_r2))
cat("    Top 10 covariates:\n")
for (i in seq_len(min(10, nrow(imp_df))))
  cat(sprintf("        %-25s  %.4f\n", imp_df$variable[i], imp_df$importance[i]))
sink()

cat(sprintf("\nReport: %s\n", report_path))
cat("=============================================================\n")
cat("  Done\n")
cat("=============================================================\n")
