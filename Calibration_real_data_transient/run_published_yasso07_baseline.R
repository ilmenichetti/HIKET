# =============================================================================
# run_published_yasso07_baseline.R
#
# Baseline diagnostic: runs Yasso07 in single mode at the published default parameters
# (free_defaults; sigma_input = 1) across all calibration-ready plots, with
# no calibration  Produces the same plot types as the calibrated
# pipeline (delta_soc_trajectory, obs_vs_pred, residuals_diagnostic,
# rf_importance) but with single-point predictions instead of posterior bands.
#
# PURPOSE:
#   - Reference baseline: how well does published Yasso07 predict observed
#     SOC before any calibration?
#   - Quantify systematic biases by region / zone / species at defaults.
#   - RF on residuals identifies covariates the published model is missing.
#   - If calibration drives sigma_input far from 1 (or other parameters far
#     from defaults), the residual pattern here vs. post-calibration tells
#     us what the data is informing.
#
# OUTPUTS (./Calibration_real_data/diagnostics/Yasso07_published/):
#   Yasso07_published_predictions.rds       (one row per plot-year)
#   Yasso07_published_residuals.csv         (one row per observation)
#   Yasso07_published_delta_soc_trajectory.png
#   Yasso07_published_obs_vs_pred.png
#   Yasso07_published_residuals_diagnostic.png
#   Yasso07_published_rf_importance.png
#   Yasso07_published_report.txt
# =============================================================================

source("./Calibration_real_data/calibration_engine.R")
source("./Model_functions_real_data/input_compatibility_layer.R")
source("./Model_functions_real_data/Decomposition_functions/Yasso/yasso07_wrapper.R")
dyn.load("./Model_functions_real_data/Decomposition_functions/Yasso/yasso07.so")

suppressPackageStartupMessages({
  library(dplyr)
  library(ranger)
  library(parallel)
})


# =============================================================================
# 0.  Configuration
# =============================================================================

MODEL_NAME   <- "Yasso07_published"
DIR_OUT      <- "./Calibration_real_data/diagnostics/Yasso07_published"
dir.create(DIR_OUT, showWarnings = FALSE, recursive = TRUE)

STEADY_STATE_YEARS <- 20L
N_CORES <- parallel::detectCores() - 1L
PX_PER_IN <- 150L

cat("=============================================================\n")
cat("  Published Yasso07 baseline (no calibration)\n")
cat(sprintf("  Output: %s\n", DIR_OUT))
cat("=============================================================\n\n")


# =============================================================================
# 1.  Parameters at published defaults
# =============================================================================
# free_defaults uses the same construction as the calibration script, just
# without the to_unconstrained transform — we run the model directly.

p_default <- YASSO07_DEFAULT_PARAMS

free_defaults <- c(
  p_default[c("p_AW","p_AE","p_AN",
              "p_WA","p_WE","p_WN",
              "p_EA","p_EW","p_EN",
              "p_NA","p_NW","p_NE",
              "beta1","beta2","gamma",
              "delta1","delta2","r")],
  sigma_init  = 0.10,
  sigma_input = 1.00       # by definition, no input scaling at baseline
)

FIXED_RATE_NAMES <- c("alpha_A","alpha_W","alpha_E","alpha_N","p_H","alpha_H")
fixed_rates      <- p_default[FIXED_RATE_NAMES]

assemble_model_params <- function(p_free) {
  yasso_params <- c(fixed_rates,
                    p_free[c("p_AW","p_AE","p_AN",
                             "p_WA","p_WE","p_WN",
                             "p_EA","p_EW","p_EN",
                             "p_NA","p_NW","p_NE",
                             "beta1","beta2","gamma",
                             "delta1","delta2","r")])
  yasso_params <- yasso_params[names(p_default)]
  c(yasso_params, sigma_input = unname(p_free["sigma_input"]))
}

model_params <- assemble_model_params(free_defaults)


# =============================================================================
# 2.  Load and filter data
# =============================================================================

cat("Loading data...\n")
input_raw <- read.csv("./Data/model_inputs/input_raw_monthly.csv")
site_raw  <- read.csv("./Data/model_inputs/site_raw.csv")
site_raw$plot_id <- as.character(site_raw$plot_id)

calib_plots <- as.character(site_raw$plot_id[site_raw$calib_ready])
input_calib <- input_raw[as.character(input_raw$plot_id) %in% calib_plots, ]

# Drop plots with NA climate or NA litter (defensive)
litter_cols <- grep("^C_(nwl|fwl|cwl)_[AWEN]$",
                    names(input_calib), value = TRUE)
bad_clim   <- unique(input_calib$plot_id[is.na(input_calib$temp_air)])
bad_litter <- unique(input_calib$plot_id[!complete.cases(input_calib[, litter_cols])])
calib_plots <- setdiff(calib_plots, c(bad_clim, bad_litter))
input_calib <- input_calib[as.character(input_calib$plot_id) %in% calib_plots, ]

Yasso07_climate <- map_climate_yasso07(input_calib)
Yasso07_inputs  <- map_inputs_yasso07(input_calib)
plots_real      <- as.character(unique(Yasso07_climate$plot_id))

cat(sprintf("Plots with valid data: %d\n", length(plots_real)))

# Pre-split for O(1) lookup
climate_by_plot <- split(Yasso07_climate, Yasso07_climate$plot_id)
inputs_by_plot  <- split(Yasso07_inputs,  Yasso07_inputs$plot_id)

# Litter means for steady-state
litter_means <- lapply(plots_real, function(pid) {
  inp <- Yasso07_inputs[as.character(Yasso07_inputs$plot_id) == pid, ]
  n_ss <- min(STEADY_STATE_YEARS, nrow(inp))
  inp_ss <- inp[seq_len(n_ss), ]
  list(
    nwl_mean = colMeans(inp_ss[, c("nwl_A","nwl_W","nwl_E","nwl_N")]),
    fwl_mean = colMeans(inp_ss[, c("fwl_A","fwl_W","fwl_E","fwl_N")]),
    cwl_mean = colMeans(inp_ss[, c("cwl_A","cwl_W","cwl_E","cwl_N")])
  )
})
names(litter_means) <- plots_real

# SOC observations
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

compute_xi_engine <- function(clim, params) {
  compute_xi_yasso07(
    temp_mean = clim$temp_mean, temp_amp = clim$temp_amplitude,
    precip    = clim$precip,
    beta1 = params["beta1"], beta2 = params["beta2"], gamma = params["gamma"])
}
compute_xi_mean_engine <- function(clim_ss, params) {
  compute_xi_mean_yasso07(
    clim_ss = clim_ss,
    beta1 = params["beta1"], beta2 = params["beta2"], gamma = params["gamma"])
}
steady_state_engine <- function(params, lm, xi_mean) {
  si <- params["sigma_input"]
  yasso07_steady_state(
    params   = params[YASSO07_PARAM_NAMES],
    nwl_mean = lm$nwl_mean * si,
    fwl_mean = lm$fwl_mean * si,
    cwl_mean = lm$cwl_mean * si,
    xi_mean  = xi_mean)
}
run_engine <- function(inputs, params, C_init, xi_array) {
  si <- params["sigma_input"]
  inputs_scaled <- inputs
  inputs_scaled[, LITTER_COLS] <- inputs[, LITTER_COLS] * si
  yasso07_run(inputs_scaled, params[YASSO07_PARAM_NAMES], C_init, xi_array)
}

cat(sprintf("\nRunning Yasso07 at defaults across %d plots...\n", length(plots_real)))
t0 <- Sys.time()

predictions <- do.call(rbind, mclapply(plots_real, function(pid) {
  clim   <- climate_by_plot[[pid]]
  inputs <- inputs_by_plot[[pid]]
  lm     <- litter_means[[pid]]

  xi_array <- tryCatch(compute_xi_engine(clim, model_params),
                       error = function(e) NULL)
  if (is.null(xi_array)) return(NULL)

  n_ss      <- min(STEADY_STATE_YEARS, nrow(clim))
  xi_for_ss <- tryCatch(
    compute_xi_mean_engine(clim[seq_len(n_ss), , drop = FALSE], model_params),
    error = function(e) NULL)
  if (is.null(xi_for_ss) || !is.finite(xi_for_ss) || xi_for_ss <= 0) return(NULL)

  C_init <- tryCatch(steady_state_engine(model_params, lm, xi_for_ss),
                     error = function(e) NULL)
  if (is.null(C_init) || any(!is.finite(C_init)) || any(C_init < 0)) return(NULL)

  run_out <- tryCatch(run_engine(inputs, model_params, C_init, xi_array),
                      error = function(e) NULL)
  if (is.null(run_out)) return(NULL)

  data.frame(
    plot_id = pid,
    year    = run_out$year,
    A = run_out$A, W = run_out$W, E = run_out$E,
    N = run_out$N, H = run_out$H,
    total_soc   = run_out$total_soc,
    respiration = run_out$respiration
  )
}, mc.cores = N_CORES))

cat(sprintf("Run complete: %.1f sec  (%d plot-years)\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs")),
            nrow(predictions)))

saveRDS(predictions, file.path(DIR_OUT, "Yasso07_published_predictions.rds"))




# =============================================================================
# 3.5.  Testing a few plots
# =============================================================================

ALL_IDs <- unique(predictions$plot_id)

SAMPLED_IDs <- sample(ALL_IDs, 16)

plot_range <- range(c(predictions[predictions$plot_id %in% SAMPLED_IDs,]$total_soc,
                      SOC_obs_all[SOC_obs_all$plot_id %in% SAMPLED_IDs,]$soc_obs_tCha))
png(file.path(DIR_OUT, "Yasso07_published_sampled_plots.png"),
    width = 12L * PX_PER_IN, height = 9L * PX_PER_IN, res = PX_PER_IN)
par(mar = c(4, 4, 3, 1))
par(mfrow=c(4,4))
for(i in 1:length(SAMPLED_IDs)){
ID_plot <- SAMPLED_IDs[i]
plot(predictions[predictions$plot_id==ID_plot,]$year,
     predictions[predictions$plot_id==ID_plot,]$total_soc, type="l", ylim=(plot_range),
     main=paste("Plot ID:",SAMPLED_IDs[i] ),
     ylab="SOC stocks t/ha", xlab="year")
points(SOC_obs_all[SOC_obs_all$plot_id==ID_plot,]$year,
       SOC_obs_all[SOC_obs_all$plot_id==ID_plot,]$soc_obs_tCha, col="red", pch=16)
}

dev.off()

# =============================================================================
# 4.  Residuals data frame
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

residuals_csv <- file.path(DIR_OUT, "Yasso07_published_residuals.csv")
write.csv(residuals_df, residuals_csv, row.names = FALSE)
cat(sprintf("Residuals: %d obs  -->  %s\n", nrow(residuals_df), residuals_csv))


# =============================================================================
# 5.  Predictive metrics
# =============================================================================

obs <- residuals_df$soc_obs_tCha
hat <- residuals_df$soc_pred
finite <- is.finite(obs) & is.finite(hat)

R2   <- cor(obs[finite], hat[finite])^2
RMSE <- sqrt(mean((obs[finite] - hat[finite])^2))
bias <- mean(hat[finite] - obs[finite])
mape <- mean(abs(hat[finite] - obs[finite]) / obs[finite]) * 100

cat("\nPredictive metrics at published defaults:\n")
cat(sprintf("  R²:    %.3f\n", R2))
cat(sprintf("  RMSE:  %.2f tC/ha\n", RMSE))
cat(sprintf("  Bias:  %+.2f tC/ha\n", bias))
cat(sprintf("  MAPE:  %.1f%%\n", mape))


# =============================================================================
# 6.  Plot A: Mean ΔSOC trajectory (model continuous, observed single interval)
# =============================================================================

# Model: annual delta from first predicted year per plot
pp <- predictions %>%
  group_by(plot_id) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(delta_soc = total_soc - first(total_soc)) %>%
  ungroup()

traj <- pp %>%
  group_by(year) %>%
  summarise(mean_delta = mean(delta_soc, na.rm = TRUE), .groups = "drop") %>%
  arrange(year)

# Observed: one delta per plot = SOC_last - SOC_first, anchored at midpoint year
obs_delta_per_plot <- residuals_df %>%
  group_by(plot_id) %>%
  arrange(year, .by_group = TRUE) %>%
  summarise(
    delta_obs  = (last(soc_obs_tCha) - first(soc_obs_tCha)) /
      (last(year) - first(year)),   # tC/ha/yr
    year_first = first(year),
    year_last  = last(year),
    .groups = "drop"
  ) %>%
  filter(year_first != year_last)

obs_mean_delta <- mean(obs_delta_per_plot$delta_obs, na.rm = TRUE)
obs_se_delta   <- sd(obs_delta_per_plot$delta_obs, na.rm = TRUE) /
  sqrt(nrow(obs_delta_per_plot))

# x-span for the horizontal line: earliest first-obs to latest last-obs
x_left  <- min(obs_delta_per_plot$year_first)
x_right <- max(obs_delta_per_plot$year_last)

png(file.path(DIR_OUT, "Yasso07_published_delta_soc_trajectory.png"),
    width = 12L * PX_PER_IN, height = 6L * PX_PER_IN, res = PX_PER_IN)
par(mar = c(4, 4, 3, 1))

ylim_r <- range(c(traj$mean_delta,
                  obs_mean_delta + 1.96 * obs_se_delta,
                  obs_mean_delta - 1.96 * obs_se_delta,
                  0), na.rm = TRUE)

plot(traj$year, traj$mean_delta,
     type = "l", col = "steelblue", lwd = 2.5,
     ylim = ylim_r,
     xlab = "Year",
     ylab = expression(paste("Mean ", Delta, "SOC across plots (tC/ha)")),
     main = sprintf("Published Yasso07 baseline -- %d plots, no calibration",
                    length(plots_real)))
abline(h = 0, lty = 3, col = "grey40")

# Observed mean delta as horizontal line spanning the obs window
abline(h = obs_mean_delta,
         col = "firebrick", lwd = 2.5, lty = 2)

# 95% CI as shaded rectangle
rect(xleft = x_left, xright = x_right,
     ybottom = obs_mean_delta - 1.96 * obs_se_delta,
     ytop    = obs_mean_delta + 1.96 * obs_se_delta,
     col = adjustcolor("firebrick", 0.12), border = NA)

legend("topleft",
       legend = c("Predicted mean ΔSOC (annual)",
                  sprintf("Observed mean ΔSOC ± 95%% CI  (n = %d plots)",
                          nrow(obs_delta_per_plot))),
       col  = c("steelblue", "firebrick"),
       lwd  = c(2.5, 2.5),
       lty  = c(1, 2),
       bty  = "n", cex = 0.9)

dev.off()

# =============================================================================
# 7.  Plot B: Observed vs predicted (no error bars)
# =============================================================================

ka_vals <- residuals_df$KA
ka_palette <- c("1" = "#2166ac",   # OMT  rich
                "2" = "#74add1",   # MT
                "3" = "#f4a582",   # VT
                "4" = "#d6604d")   # CT  poor
ka_col <- ka_palette[as.character(ka_vals)]
ka_col[is.na(ka_col)] <- "grey60"

png(file.path(DIR_OUT, "Yasso07_published_obs_vs_pred.png"),
    width = 7L * PX_PER_IN, height = 7L * PX_PER_IN, res = PX_PER_IN)
par(mar = c(4.5, 4.5, 3, 1))
ax_lim <- range(c(obs, hat), na.rm = TRUE)
ax_lim <- c(0, ax_lim[2] * 1.05)
plot(NA, xlim = ax_lim, ylim = ax_lim,
     xlab = "Observed SOC (tC/ha)", ylab = "Predicted SOC (tC/ha)",
     main = "Yasso07 at published defaults")
abline(0, 1, lty = 2, col = "grey40", lwd = 1.5)
points(obs, hat, pch = 16, cex = 0.9, col = adjustcolor(ka_col, 0.7))
legend("topleft",
       legend = c(sprintf("R² = %.3f", R2),
                  sprintf("RMSE = %.1f tC/ha", RMSE),
                  sprintf("Bias = %+.1f tC/ha", bias),
                  sprintf("MAPE = %.0f%%", mape)),
       bty = "n", cex = 0.85)
ka_present <- sort(unique(na.omit(ka_vals)))
ka_labels  <- c("1" = "KA 1: OMT (rich)", "2" = "KA 2: MT",
                "3" = "KA 3: VT",         "4" = "KA 4: CT (poor)")
if (length(ka_present) > 0) {
  legend("bottomright",
         legend = ka_labels[as.character(ka_present)],
         col    = ka_palette[as.character(ka_present)],
         pch = 16, bty = "n", cex = 0.8)
}
dev.off()

# =============================================================================
# 8.  Plot C: Residuals diagnostic (6 panels)
# =============================================================================

resid <- residuals_df$residual_log
fitted_log <- residuals_df$log_hat
mask <- is.finite(resid) & is.finite(fitted_log)
resid_f <- resid[mask]
fitted_f <- fitted_log[mask]
resid_df_f <- residuals_df[mask, , drop = FALSE]
png(file.path(DIR_OUT, "Yasso07_published_residuals_diagnostic.png"),
    width = 12L * PX_PER_IN, height = 8L * PX_PER_IN, res = PX_PER_IN)
par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))

# Panel 1: residuals vs fitted
plot(fitted_f, resid_f, pch = 19, col = adjustcolor("steelblue", 0.5),
     xlab = "log(SOC predicted)", ylab = "Residual (log)",
     main = "Residuals vs fitted")
abline(h = 0, lty = 2, col = "grey40")
if (length(fitted_f) > 10) {
  lo <- loess(resid_f ~ fitted_f)
  ord <- order(fitted_f)
  lines(fitted_f[ord], predict(lo)[ord], col = "firebrick", lwd = 2)
}

# Panel 2: by year
boxplot(resid_f ~ resid_df_f$year, col = "lightsteelblue",
        xlab = "Year", ylab = "Residual (log)", main = "Residuals by year")
abline(h = 0, lty = 2, col = "grey40")

# Panel 3: by Cajander fertility class (KA)
if ("KA" %in% names(resid_df_f)) {
  boxplot(resid_f ~ as.factor(resid_df_f$KA), col = "lightsteelblue",
          xlab = "Cajander fertility class (KA)", ylab = "Residual (log)",
          main = "Residuals by Cajander class")
  abline(h = 0, lty = 2, col = "grey40")
} else {
  plot.new(); title("KA not available")
}

# Panel 4: residuals vs CoarseFragments (stoniness)
if ("CoarseFragments" %in% names(resid_df_f)) {
  stone      <- resid_df_f$CoarseFragments
  stone_mask <- is.finite(stone)
  n_stone    <- sum(stone_mask)
  plot(stone[stone_mask], resid_f[stone_mask],
       pch = 19, col = adjustcolor("steelblue", 0.45),
       xlab = "Coarse fragments / stoniness (%)",
       ylab = "Residual (log)",
       main = sprintf("Residuals vs stoniness  (n = %d)", n_stone))
  abline(h = 0, lty = 2, col = "grey40")
  if (n_stone > 20) {
    lo_s  <- loess(resid_f[stone_mask] ~ stone[stone_mask], span = 0.75)
    ord_s <- order(stone[stone_mask])
    lines(stone[stone_mask][ord_s], predict(lo_s)[ord_s],
          col = "firebrick", lwd = 2)
  }
  rug(stone[stone_mask], col = adjustcolor("steelblue", 0.3), ticksize = 0.02)
} else {
  plot.new(); title("CoarseFragments not available")
}

# Panel 5: by species_code
if ("species_code" %in% names(resid_df_f)) {
  sp_labels <- c("1" = "Scots pine", "2" = "Nor. spruce", "3" = "Birch",
                 "4" = "Aspen",      "5" = "Grey alder",  "6" = "Blk. alder",
                 "7" = "Other brdlv.")
  sp_factor <- as.factor(resid_df_f$species_code)
  levels(sp_factor) <- sp_labels[levels(sp_factor)]
  boxplot(resid_f ~ sp_factor, col = "lightsteelblue", las = 2,
          xlab = "", ylab = "Residual (log)",
          main = "Residuals by dominant species",
          cex.axis = 0.78)
  abline(h = 0, lty = 2, col = "grey40")
} else {
  plot.new(); title("species_code not available")
}

# Panel 6: QQ
qqnorm(resid_f, pch = 19, col = adjustcolor("steelblue", 0.5),
       main = "Q-Q plot of log residuals")
qqline(resid_f, col = "firebrick", lwd = 2)

mtext("Yasso07 published-defaults residual diagnostics",
      side = 3, outer = TRUE, line = -1.5, cex = 1.0, font = 2)
dev.off()

# =============================================================================
# 9.  Random Forest on residuals
# =============================================================================

target <- "residual_log"
candidate_vars <- c(
  "KA", "kasvyo_syke", "kasvyo_ahti", "species_code", "soil_code",
  "TexturalClass", "temp_zone", "koppen_class", "ojitustilanne",
  "dev_class_85", "kasvup_tyyppi", "alaryhma",
  "clay", "SiltContent", "SandContent", "CoarseFragments",
  "EstimatedBulkDensity", "profile_depth_cm",
  "ofh_lower_cm", "ofh_weight_kgm2",
  "pH.CaCl2.", "pH.H2O.", "TotalNitrogen", "CN_ratio", "CEC", "base_saturation",
  "ExchangeableAl", "ExchangeableFe", "ExchangeableCa",
  "mean_temp", "mean_precip_annual", "GDD5", "coldest_month_T",
  "warmest_month_T", "T_seasonality", "P_seasonality", "aridity_index",
  "temp_sum_NFI", "elevation_m",
  "lat_WGS84", "lon_WGS84",
  "stand_age_85", "stand_age_2006_est", "basal_area_85", "mean_height_85_dm",
  "any_cut_85_95", "n_cuts_85_95", "any_trt_85_95", "soil_prep_pre85",
  "litter_A_frac", "litter_W_frac", "litter_E_frac", "litter_N_frac",
  "woody_share", "conifer_share",
  "shallow", "n_soc_obs"
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
rf_df     <- rf_df[, c(target, available), drop = FALSE]

rf_df <- rf_df[complete.cases(rf_df) & is.finite(rf_df[[target]]), ,
               drop = FALSE]
cat(sprintf("\nRF data: %d obs x %d covariates\n",
            nrow(rf_df), length(available)))

set.seed(2025)
rf_fit <- ranger(
  formula     = as.formula(paste(target, "~ .")),
  data        = rf_df,
  num.trees   = 500L,
  importance  = "permutation",
  num.threads = N_CORES)

oob_r2 <- 1 - rf_fit$prediction.error / var(rf_df[[target]])
imp <- sort(rf_fit$variable.importance, decreasing = TRUE)
imp_df <- data.frame(variable = names(imp), importance = unname(imp))

cat(sprintf("RF OOB R²: %.3f\n", oob_r2))

# RF importance plot
n_show <- min(20L, nrow(imp_df))
imp_top <- imp_df[seq_len(n_show), ]
imp_top <- imp_top[order(imp_top$importance), ]

png(file.path(DIR_OUT, "Yasso07_published_rf_importance.png"),
    width  = 10L * PX_PER_IN,
    height = max(5L, n_show * 0.3) * PX_PER_IN,
    res    = PX_PER_IN)
par(mar = c(4, 12, 3, 1))
barplot(imp_top$importance, horiz = TRUE,
        names.arg = imp_top$variable, las = 1, cex.names = 0.85,
        col = "steelblue", border = "white",
        xlab = "Permutation importance",
        main = sprintf("Yasso07 baseline -- top %d covariates  |  OOB R² = %.3f",
                       n_show, oob_r2))
dev.off()


# =============================================================================
# 10.  Text report
# =============================================================================

report_path <- file.path(DIR_OUT, "Yasso07_published_report.txt")
sink(report_path)
cat("=========================================================\n")
cat("  Yasso07 Published-Defaults Baseline Report\n")
cat(sprintf("  Generated: %s\n", format(Sys.time())))
cat("=========================================================\n\n")

cat("[1] Configuration\n")
cat(sprintf("    Plots simulated:           %d\n", length(plots_real)))
cat(sprintf("    Plot-year predictions:     %d\n", nrow(predictions)))
cat(sprintf("    Observations matched:      %d\n", nrow(residuals_df)))
cat(sprintf("    Steady-state window (yr):  %d\n", STEADY_STATE_YEARS))
cat(sprintf("    sigma_input:               1.00 (no scaling)\n"))
cat("\n")

cat("[2] Predictive performance at defaults\n")
cat(sprintf("    R²:    %.3f\n", R2))
cat(sprintf("    RMSE:  %.2f tC/ha\n", RMSE))
cat(sprintf("    Bias:  %+.2f tC/ha\n", bias))
cat(sprintf("    MAPE:  %.1f%%\n", mape))
cat("\n")

cat("[3] Bias by stratum (median residual_log; positive => obs > pred)\n")
for (strat in c("kasvyo_syke", "KA", "species_code", "temp_zone")) {
  if (strat %in% names(residuals_df)) {
    bs <- residuals_df %>%
      filter(is.finite(residual_log)) %>%
      group_by(.data[[strat]]) %>%
      summarise(n = n(),
                med_resid = median(residual_log),
                .groups = "drop")
    cat(sprintf("    %s:\n", strat))
    for (i in seq_len(nrow(bs))) {
      cat(sprintf("        %-10s  n=%-4d  median_resid_log = %+.3f\n",
                  as.character(bs[[strat]][i]),
                  bs$n[i], bs$med_resid[i]))
    }
  }
}
cat("\n")

cat("[4] Residual analysis (Random Forest)\n")
cat(sprintf("    Observations:              %d\n", nrow(rf_df)))
cat(sprintf("    Covariates considered:     %d\n", length(available)))
cat(sprintf("    RF Out-Of-Bag R²:          %.3f\n", oob_r2))
cat("    Top covariates by permutation importance:\n")
for (i in seq_len(min(10, nrow(imp_df)))) {
  cat(sprintf("        %-25s  %.4f\n",
              imp_df$variable[i], imp_df$importance[i]))
}
cat("\n")

cat(sprintf("=========================================================\n"))
sink()

cat(sprintf("\nReport: %s\n", report_path))
cat("=============================================================\n")
cat("  Done\n")
cat("=============================================================\n")





# =============================================================================
# 11.  Rank plots by residual magnitude and inspect worst cases
# =============================================================================

# One summary row per plot: mean log-residual and mean absolute log-residual.
# We use mean rather than sum so plots with 2 obs are comparable to plots with 1.
plot_resid_summary <- residuals_df %>%
  filter(is.finite(residual_log)) %>%
  group_by(plot_id) %>%
  summarise(
    n_obs         = n(),
    mean_resid    = mean(residual_log),      # signed: positive => obs > pred (underprediction)
    abs_resid     = mean(abs(residual_log)), # magnitude regardless of direction
    mean_obs      = mean(soc_obs_tCha),
    mean_pred     = mean(soc_pred),
    ratio_obs_pred = mean(soc_obs_tCha / soc_pred),
    .groups = "drop"
  ) %>%
  arrange(desc(abs_resid))

write.csv(plot_resid_summary,
          file.path(DIR_OUT, "Yasso07_published_plot_ranking.csv"),
          row.names = FALSE)

cat("\n=== WORST 20 PLOTS BY MEAN |log-RESIDUAL| ===\n")
print(head(plot_resid_summary, 20), digits = 3)



# =============================================================================
# 12.  Time series panel: 25 worst plots
# =============================================================================

worst_ids <- head(plot_resid_summary$plot_id, 25)

worst_pred <- predictions[predictions$plot_id %in% worst_ids, ]
worst_obs  <- SOC_obs_all[SOC_obs_all$plot_id %in% worst_ids, ]

y_range <- range(c(worst_pred$total_soc, worst_obs$soc_obs_tCha), na.rm = TRUE)

png(file.path(DIR_OUT, "Yasso07_published_worst25_timeseries.png"),
    width = 15L * PX_PER_IN, height = 15L * PX_PER_IN, res = PX_PER_IN)
par(mfrow = c(5, 5), mar = c(3, 3, 2, 1), oma = c(0, 0, 3, 0))

for (pid in worst_ids) {
  pred_p <- worst_pred[worst_pred$plot_id == pid, ]
  obs_p  <- worst_obs[worst_obs$plot_id  == pid, ]
  
  # Pull summary stats for subtitle
  sr <- plot_resid_summary[plot_resid_summary$plot_id == pid, ]
  
  plot(pred_p$year, pred_p$total_soc,
       type = "l", lwd = 2, col = "steelblue",
       ylim = y_range,
       xlab = "", ylab = "SOC (tC/ha)",
       main = sprintf("Plot %s\nbias=%.2f  ratio=%.2f",
                      pid, sr$mean_resid, sr$ratio_obs_pred),
       cex.main = 0.75)
  points(obs_p$year, obs_p$soc_obs_tCha,
         pch = 16, col = "firebrick", cex = 1.4)
}

mtext("Yasso07 defaults — 25 worst plots (by mean |log-residual|)",
      side = 3, outer = TRUE, cex = 1.1, font = 2)
dev.off()

cat(sprintf("\nWorst-25 panel: %s\n",
            file.path(DIR_OUT, "Yasso07_published_worst25_timeseries.png")))



# =============================================================================
# 13.  Histogram of log-residuals
# =============================================================================

resid_finite <- residuals_df$residual_log[is.finite(residuals_df$residual_log)]

png(file.path(DIR_OUT, "Yasso07_published_residuals_histogram.png"),
    width = 8L * PX_PER_IN, height = 6L * PX_PER_IN, res = PX_PER_IN)
par(mar = c(4.5, 4.5, 3, 1))

hist(resid_finite,
     breaks = 40,
     col    = "steelblue", border = "white",
     xlab   = "Residual  log(obs) - log(pred)",
     ylab   = "Number of observations",
     main   = sprintf("Yasso07 defaults — residual distribution  (n = %d)", 
                      length(resid_finite)))
abline(v = 0,                      lty = 2, lwd = 2, col = "grey30")
abline(v = mean(resid_finite),     lty = 1, lwd = 2, col = "firebrick")
abline(v = median(resid_finite),   lty = 1, lwd = 2, col = "darkorange")
legend("topleft",
       legend = c("zero",
                  sprintf("mean   = %.3f", mean(resid_finite)),
                  sprintf("median = %.3f", median(resid_finite))),
       col    = c("grey30", "firebrick", "darkorange"),
       lty    = c(2, 1, 1), lwd = 2, bty = "n")
dev.off()

cat(sprintf("\nResidual histogram: %s\n",
            file.path(DIR_OUT, "Yasso07_published_residuals_histogram.png")))
