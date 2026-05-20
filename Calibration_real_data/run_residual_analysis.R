# =============================================================================
# run_residual_analysis.R
#
# Model-agnostic residual analysis stage of the HIKET pipeline.
# Reads the residuals CSV produced by run_<MODEL>_predictive.R, runs Random
# Forest variable importance and basic by-stratum residual visualisations,
# and appends a [6] Residual analysis section to the calibration report.
#
# USAGE:
#   Rscript run_residual_analysis.R <MODEL_NAME> <RUN_ID>
#   e.g.  Rscript run_residual_analysis.R Yasso07 20260427_213045
#
# Same script will be reused for Yasso15, Yasso20, RothC once they are
# calibrated.
#
# OUTPUTS:
#   ./Calibration_real_data/diagnostics/<MODEL>/
#     <MODEL>_residuals_diagnostic_<RUN_ID>.png   -- 3x3 residual panel
#     <MODEL>_litter_posterior_<RUN_ID>.png        -- posterior-corrected litter by region
#     <MODEL>_residuals_histogram_<RUN_ID>.png     -- log-residual distribution
#     <MODEL>_rf_importance_<RUN_ID>.png           -- RF permutation importance
#     <MODEL>_rf_importance_<RUN_ID>.csv
#     <MODEL>_rf_summary_<RUN_ID>.rds
#   Appends [6] Residual analysis to the calibration report.
# =============================================================================

library(ranger)
library(dplyr)

source("./Calibration_real_data/calibration_engine.R")   # for append_to_report

# =============================================================================
# 0.  Configuration / argument parsing
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  # Auto-detect: most recent residuals CSV across all model dirs.
  diag_root <- "./Calibration_real_data/diagnostics"
  candidates <- list.files(
    diag_root,
    pattern = "_residuals_[0-9]{8}_[0-9]{6}\\.csv$",
    recursive = TRUE, full.names = TRUE)
  if (length(candidates) == 0)
    stop("No residual CSVs found. Run a predictive script first, ",
         "or pass MODEL RUN_ID.")
  candidates <- candidates[order(file.info(candidates)$mtime, decreasing = TRUE)]
  latest     <- candidates[1]
  MODEL_NAME <- basename(dirname(latest))
  RUN_ID     <- sub(sprintf("^%s_residuals_(.+)\\.csv$", MODEL_NAME),
                    "\\1", basename(latest))
  message(sprintf("Auto-detected: %s / %s", MODEL_NAME, RUN_ID))
} else {
  MODEL_NAME <- args[1]
  RUN_ID     <- args[2]
}

DIR_DIAG <- file.path("./Calibration_real_data/diagnostics", MODEL_NAME)
run_config <- list(MODEL_NAME = MODEL_NAME, RUN_ID = RUN_ID, DIR_DIAG = DIR_DIAG)

PX_PER_IN <- 150L

message("=============================================================")
message(sprintf("  %s Residual Analysis  |  Run: %s", MODEL_NAME, RUN_ID))
message("=============================================================\n")


# =============================================================================
# 1.  Load residuals (already merged with site_raw in predictive script)
# =============================================================================

resid_csv <- file.path(DIR_DIAG, sprintf("%s_residuals_%s.csv",
                                         MODEL_NAME, RUN_ID))
if (!file.exists(resid_csv)) stop("Residuals not found: ", resid_csv)
resid_df <- read.csv(resid_csv, stringsAsFactors = FALSE)
message(sprintf("Loaded residuals: %d rows", nrow(resid_df)))


# =============================================================================
# 1b. Load posterior samples and litter data (for section 5b)
# =============================================================================
# sigma_input is a free parameter (global litter scaling multiplier) estimated
# by the MCMC. Its posterior is a column in the physical-space samples matrix
# saved by save_results(). We use it to show calibration-corrected litter
# trajectories by region.

posterior_phys <- tryCatch({
  f <- file.path("./Calibration_real_data/runs/",
                 sprintf("%s_posterior_%s.rds", MODEL_NAME, RUN_ID))
  if (!file.exists(f)) stop("not found")
  readRDS(f)
}, error = function(e) {
  message(sprintf("WARNING: posterior RDS not found — skipping litter plot (%s)", e$message))
  NULL
})

sigma_input_samples <- if (!is.null(posterior_phys) &&
                           "sigma_input" %in% colnames(posterior_phys)) {
  posterior_phys[, "sigma_input"]
} else {
  message("WARNING: sigma_input not found in posterior — skipping litter plot")
  NULL
}

# Monthly litter inputs (calib_ready plots). Year 2024 excluded: it is the
# carried-forward copy of 2023 used only to place the Komeetta SOC observation.
input_raw <- tryCatch(
  read.csv("./Data/model_inputs/input_raw_monthly.csv"),
  error = function(e) { message("WARNING: input_raw_monthly.csv not found"); NULL }
)


# =============================================================================
# 2.  RF covariate selection
# =============================================================================
# Excluded:
#   OrganicCarbon, OrganicMatter -- collinear with SOC obs by construction
#   soil_depth                   -- replaced by profile_depth_cm (measured)
#   MeanBulkDensity              -- 84% NA in calib_ready (audit)
#   plot_id, year                -- identifiers, not predictors
#   any soc_*, log_*, residual_* -- the targets and derived quantities
#
# ranger handles numeric, integer, and factor columns natively.
# Only character columns need explicit coercion since ranger rejects raw strings.

target <- "residual_log"

candidate_vars <- c(
  # Stratification / categorical
  "KA", "kasvyo_syke", "kasvyo_ahti", "species_code", "soil_code",
  "TexturalClass", "temp_zone", "koppen_class", "ojitustilanne",
  "dev_class_85", "kasvup_tyyppi", "alaryhma",
  # Soil physical
  "clay", "SiltContent", "SandContent", "CoarseFragments",
  "EstimatedBulkDensity", "profile_depth_cm",
  # Organic horizon
  "ofh_lower_cm", "ofh_weight_kgm2",
  # Soil chemistry
  "pH.CaCl2.", "pH.H2O.", "TotalNitrogen", "CN_ratio", "CEC", "base_saturation",
  "ExchangeableAl", "ExchangeableFe", "ExchangeableCa",
  # Climate
  "mean_temp", "mean_precip_annual", "GDD5", "coldest_month_T",
  "warmest_month_T", "T_seasonality", "P_seasonality", "aridity_index",
  "temp_sum_NFI", "elevation_m",
  # Spatial
  "lat_WGS84", "lon_WGS84",
  # Stand & management (1985)
  "stand_age_85", "stand_age_2006_est", "basal_area_85", "mean_height_85_dm",
  "any_cut_85_95", "n_cuts_85_95", "any_trt_85_95", "soil_prep_pre85",
  # Litter quality
  "litter_A_frac", "litter_W_frac", "litter_E_frac", "litter_N_frac",
  "woody_share", "conifer_share",
  # Plot-level summary
  "shallow", "n_soc_obs"
)

available <- intersect(candidate_vars, names(resid_df))
unused    <- setdiff(candidate_vars, names(resid_df))
if (length(unused) > 0)
  message(sprintf("Skipping %d not-found covariates: %s",
                  length(unused), paste(unused, collapse = ", ")))

# Build RF data frame; coerce character columns to factor (ranger requirement)
rf_df <- resid_df[, c(target, available), drop = FALSE]
for (v in available)
  if (is.character(rf_df[[v]])) rf_df[[v]] <- as.factor(rf_df[[v]])

# Drop rows with NA or Inf in the target only.
# Covariates with NA are handled internally by ranger (split on available obs).
complete_rows <- is.finite(rf_df[[target]])
n_dropped     <- sum(!complete_rows)
rf_df         <- rf_df[complete_rows, , drop = FALSE]
message(sprintf("RF data: %d rows x %d covariates (dropped %d non-finite targets)",
                nrow(rf_df), length(available), n_dropped))


# =============================================================================
# 3.  Random Forest with permutation importance
# =============================================================================

set.seed(2025)
rf_fit <- ranger(
  formula      = as.formula(paste(target, "~ .")),
  data         = rf_df,
  num.trees    = 500L,
  importance   = "permutation",
  num.threads  = parallel::detectCores() - 1L
)
oob_r2 <- rf_fit$r.squared   # ranger OOB R² (on residual_log scale)
imp <- sort(rf_fit$variable.importance, decreasing = TRUE)
imp_df <- data.frame(variable = names(imp), importance = unname(imp))
write.csv(imp_df,
          file.path(DIR_DIAG,
                    sprintf("%s_rf_importance_%s.csv", MODEL_NAME, RUN_ID)),
          row.names = FALSE)

rf_summary <- list(
  rf_fit     = rf_fit,
  oob_r2     = oob_r2,
  importance = imp_df,
  n_obs      = nrow(rf_df),
  covariates = available,
  target     = target
)
saveRDS(rf_summary,
        file.path(DIR_DIAG,
                  sprintf("%s_rf_summary_%s.rds", MODEL_NAME, RUN_ID)))


# =============================================================================
# 4.  Plot: RF importance
# =============================================================================

n_show  <- min(20L, nrow(imp_df))
imp_top <- imp_df[seq_len(n_show), ]
imp_top <- imp_top[order(imp_top$importance), ]   # ascending for horizontal bar

png(file.path(DIR_DIAG,
              sprintf("%s_rf_importance_%s.png", MODEL_NAME, RUN_ID)),
    width  = 10L * PX_PER_IN,
    height = max(5L, n_show * 0.3) * PX_PER_IN,
    res    = PX_PER_IN)
par(mar = c(4, 12, 3, 1))
barplot(imp_top$importance, horiz = TRUE,
        names.arg = imp_top$variable, las = 1, cex.names = 0.85,
        col = "steelblue", border = "white",
        xlab = "Permutation importance",
        main = sprintf("%s residual analysis -- top %d covariates  |  OOB R2 = %.3f",
                       MODEL_NAME, n_show, oob_r2))
dev.off()
message(sprintf("RF importance plot: %s",
                file.path(DIR_DIAG,
                          sprintf("%s_rf_importance_%s.png", MODEL_NAME, RUN_ID))))

# =============================================================================
# 5.  Diagnostic plot: 3x3 residual panel
# =============================================================================
# Panel layout:
#   [1] Residuals vs fitted (loess)    [2] By year          [3] By Cajander class (KA)
#   [4] Vs stoniness (loess + rug)     [5] By species       [6] Q-Q plot
#   [7] Vs profile depth (loess + rug) [8] Vs mean temp     [9] Vs CN ratio (loess + rug)
#
# Panels 7-9 are new relative to the original 2x3 layout:
#   [7] profile_depth_cm: tests whether the 1m exponential extrapolation over/
#       underestimates tail SOC for shallow vs deep profiles.
#   [8] mean_temp: direct structural test of the temperature sensitivity function
#       (xi). A trend across Finland's -2 to +6 C gradient is a smoking gun for
#       a mis-specified climate response -- the key structural comparison axis
#       between Yasso variants and RothC.
#   [9] CN_ratio: litter quality proxy. If AWEN fractionation is off for certain
#       site types, or if the model mishandles the quality-decomposition relation,
#       CN ratio should show a residual trend.
# =============================================================================

resid_png <- file.path(DIR_DIAG,
                       sprintf("%s_residuals_diagnostic_%s.png",
                               MODEL_NAME, RUN_ID))
png(resid_png, width = 15L * PX_PER_IN, height = 15L * PX_PER_IN, res = PX_PER_IN)
par(mfrow = c(3, 3), mar = c(4, 4, 3, 1))

resid  <- resid_df$residual_log
fitted <- resid_df$log_hat_mean

# Drop non-finite values (Inf from log(0) if posterior mean <= 0)
finite_mask <- is.finite(fitted) & is.finite(resid)
n_dropped_finite <- sum(!finite_mask)
if (n_dropped_finite > 0)
  message(sprintf("Note: %d non-finite residuals dropped from plots", n_dropped_finite))

fitted_f   <- fitted[finite_mask]
resid_f    <- resid[finite_mask]
resid_df_f <- resid_df[finite_mask, , drop = FALSE]

# --- Helper: scatter with loess + rug ----------------------------------------
scatter_loess <- function(x, y, xlab, main, span = 0.75) {
  mask <- is.finite(x)
  n    <- sum(mask)
  plot(x[mask], y[mask],
       pch = 19, col = adjustcolor("steelblue", 0.45),
       xlab = xlab, ylab = "Residual (log)",
       main = sprintf("%s  (n = %d)", main, n))
  abline(h = 0, lty = 2, col = "grey40")
  if (n > 20) {
    lo  <- loess(y[mask] ~ x[mask], span = span)
    ord <- order(x[mask])
    lines(x[mask][ord], predict(lo)[ord], col = "firebrick", lwd = 2)
  }
  rug(x[mask], col = adjustcolor("steelblue", 0.3), ticksize = 0.02)
}

# Panel 1: residuals vs fitted
plot(fitted_f, resid_f, pch = 19, col = adjustcolor("steelblue", 0.5),
     xlab = "log(SOC predicted)", ylab = "Residual (log)",
     main = "Residuals vs fitted")
abline(h = 0, lty = 2, col = "grey40")
if (length(fitted_f) > 10) {
  lo  <- loess(resid_f ~ fitted_f)
  ord <- order(fitted_f)
  lines(fitted_f[ord], predict(lo)[ord], col = "firebrick", lwd = 2)
}

# Panel 2: by year
boxplot(resid_f ~ resid_df_f$year, col = "lightsteelblue",
        xlab = "Year", ylab = "Residual (log)", main = "Residuals by year")
abline(h = 0, lty = 2, col = "grey40")

# Panel 3: by Cajander fertility class (KA)
# KA 1-4 = mineral fertility classes (OMT down to CT); 11-13 = peatland types
# (excluded from calib_ready, so should not appear here).
if ("KA" %in% names(resid_df_f)) {
  boxplot(resid_f ~ as.factor(resid_df_f$KA), col = "lightsteelblue",
          xlab = "Cajander fertility class (KA)", ylab = "Residual (log)",
          main = "Residuals by Cajander class")
  abline(h = 0, lty = 2, col = "grey40")
} else {
  plot.new(); title("KA not available")
}

# Panel 4: vs stoniness (CoarseFragments)
# CoarseFragments from Biosoil mineral topsoil (M01 + M12 mean).
# Stoniness affects rooting depth, water retention, and bulk density estimates;
# expected to be a moderator of SOC accumulation in Finnish till soils.
if ("CoarseFragments" %in% names(resid_df_f)) {
  scatter_loess(resid_df_f$CoarseFragments, resid_f,
                xlab = "Coarse fragments / stoniness (%)",
                main = "Residuals vs stoniness")
} else {
  plot.new(); title("CoarseFragments not available")
}

# Panel 5: by dominant species
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

# Panel 6: Q-Q plot
qqnorm(resid_f, pch = 19, col = adjustcolor("steelblue", 0.5),
       main = "Q-Q plot of log residuals")
qqline(resid_f, col = "firebrick", lwd = 2)

# Panel 7: vs profile depth (measured Biosoil profile depth, cm)
# A trend here could indicate that the exponential depth model used in
# Data_work.R over/underestimates the tail stock for shallow or deep profiles.
if ("profile_depth_cm" %in% names(resid_df_f)) {
  scatter_loess(resid_df_f$profile_depth_cm, resid_f,
                xlab = "Measured profile depth (cm)",
                main = "Residuals vs profile depth")
} else {
  plot.new(); title("profile_depth_cm not available")
}

# Panel 8: vs mean annual temperature
# Direct structural test of the temperature sensitivity function (xi).
# A systematic trend across the Finnish gradient (~-2 to +6 C) indicates
# a mis-specified climate response -- the key comparison axis across models.
if ("mean_temp" %in% names(resid_df_f)) {
  scatter_loess(resid_df_f$mean_temp, resid_f,
                xlab = "Mean annual temperature (C)",
                main = "Residuals vs mean temperature")
} else {
  plot.new(); title("mean_temp not available")
}

# Panel 9: vs CN ratio
# Litter quality proxy. If AWEN fractionation is off for certain site types,
# or if the model mishandles the quality-decomposition relationship,
# CN ratio should show a residual trend.
if ("CN_ratio" %in% names(resid_df_f)) {
  scatter_loess(resid_df_f$CN_ratio, resid_f,
                xlab = "CN ratio (OrganicCarbon / TotalNitrogen)",
                main = "Residuals vs CN ratio",
                span = 0.85)
} else {
  plot.new(); title("CN_ratio not available")
}

mtext(sprintf("%s residual diagnostics  |  Run: %s", MODEL_NAME, RUN_ID),
      side = 3, outer = TRUE, line = -1.5, cex = 1.0, font = 2)
dev.off()
message(sprintf("Residual diagnostic: %s", resid_png))


# =============================================================================
# 5b. Posterior-corrected litter trajectories by region
# =============================================================================
# sigma_input is a global multiplicative scaling on litter inputs calibrated
# by the MCMC. Its posterior distribution tells us:
#   (a) the direction and magnitude of the input bias the model inferred
#       (sigma_input > 1 = inputs underestimated; < 1 = overestimated)
#   (b) the uncertainty in that correction
#
# This plot shows, for each region:
#   - Raw mean annual litter input trajectory (grey dashed) -- as supplied to
#     the model before any correction
#   - Posterior-corrected trajectory (coloured line) = raw x median(sigma_input)
#   - 95% credible band (coloured ribbon) = raw x quantile(sigma_input, 0.025/0.975)
#
# The corrected trajectory is the model's best estimate of the "true" litter
# forcing after accounting for NFI measurement uncertainty and MUSTIKKA biases.
# Comparing north vs south shows whether the global sigma_input corrects the
# input equally well across the latitudinal gradient.
#
# Year 2024 is excluded: it is the carried-forward copy of 2023 used only to
# anchor the Komeetta SOC observation and does not represent a real litter
# measurement year.
# =============================================================================

if (!is.null(sigma_input_samples) && !is.null(input_raw)) {
  
  si_med  <- median(sigma_input_samples)
  si_lo   <- quantile(sigma_input_samples, 0.025)
  si_hi   <- quantile(sigma_input_samples, 0.975)
  
  message(sprintf("\nsigma_input posterior: median=%.3f  95%% CI=[%.3f, %.3f]",
                  si_med, si_lo, si_hi))
  
  # Annual litter per plot per year (sum all AWEN x size columns x 12 months)
  litter_cols <- grep("^C_(nwl|fwl|cwl)_[AWEN]$", names(input_raw), value = TRUE)
  
  # Filter to calib_ready plots and exclude the 2024 carry-forward year
  calib_plots_ids <- unique(resid_df$plot_id)
  input_calib_lit <- input_raw[
    input_raw$plot_id %in% calib_plots_ids &
      input_raw$year    != 2024L, ]
  
  annual_lit <- input_calib_lit %>%
    group_by(plot_id, year) %>%
    summarise(total_litter = sum(across(all_of(litter_cols)), na.rm = TRUE) * 12,
              .groups = "drop")
  
  # Attach region from resid_df (one value per plot)
  region_lkup <- unique(resid_df[, c("plot_id", "region")])
  annual_lit  <- merge(annual_lit, region_lkup, by = "plot_id", all.x = TRUE)
  annual_lit  <- annual_lit[!is.na(annual_lit$region), ]
  
  regions     <- sort(unique(annual_lit$region))
  n_regions   <- length(regions)
  
  if (n_regions == 0) {
    message("WARNING: no region labels found -- skipping litter plot")
  } else {
    
    # Mean raw litter per year x region
    mean_lit <- annual_lit %>%
      group_by(region, year) %>%
      summarise(mean_raw = mean(total_litter, na.rm = TRUE), .groups = "drop")
    
    # Region colours: use up to 6 distinguishable colours
    reg_cols <- c("#2166ac", "#d6604d", "#1a9641", "#762a83", "#e08214", "#4d4d4d")
    reg_cols <- setNames(reg_cols[seq_len(n_regions)], regions)
    
    litter_png <- file.path(DIR_DIAG,
                            sprintf("%s_litter_posterior_%s.png",
                                    MODEL_NAME, RUN_ID))
    
    png(litter_png,
        width  = max(8L, 5L * n_regions) * PX_PER_IN,
        height = 6L * PX_PER_IN,
        res    = PX_PER_IN)
    par(mfrow = c(1, n_regions), mar = c(4, 4.5, 4, 1))
    
    for (reg in regions) {
      reg_dat <- mean_lit[mean_lit$region == reg, ]
      reg_dat <- reg_dat[order(reg_dat$year), ]
      
      yrs     <- reg_dat$year
      raw_mn  <- reg_dat$mean_raw
      corr_md <- raw_mn * si_med
      corr_lo <- raw_mn * si_lo
      corr_hi <- raw_mn * si_hi
      
      ylim <- range(c(raw_mn, corr_lo, corr_hi), na.rm = TRUE)
      ylim <- ylim + c(-0.05, 0.05) * diff(ylim)
      
      col_reg <- reg_cols[as.character(reg)]
      
      plot(yrs, raw_mn,
           type = "l", lty = 2, lwd = 1.5, col = "grey50",
           ylim = ylim,
           xlab = "Year",
           ylab = "Mean annual litter input (tC/ha/yr)",
           main = sprintf("Region: %s\nsigma_input: %.3f  [%.3f, %.3f]",
                          reg, si_med, si_lo, si_hi))
      
      # 95% CI ribbon for corrected litter
      polygon(c(yrs, rev(yrs)),
              c(corr_lo, rev(corr_hi)),
              col = adjustcolor(col_reg, 0.18), border = NA)
      
      # Corrected median line
      lines(yrs, corr_md, col = col_reg, lwd = 2.5)
      
      # Reference line at sigma_input = 1 (no correction)
      abline(h = mean(raw_mn, na.rm = TRUE), lty = 3, col = "grey70", lwd = 1)
      
      legend("topright",
             legend = c("Raw litter (uncorrected)",
                        sprintf("Posterior median (x %.3f)", si_med),
                        "95% credible band"),
             col    = c("grey50", col_reg, adjustcolor(col_reg, 0.4)),
             lty    = c(2, 1, NA),
             lwd    = c(1.5, 2.5, NA),
             pch    = c(NA, NA, 15),
             pt.cex = c(NA, NA, 2),
             bty    = "n", cex = 0.82)
    }
    
    mtext(sprintf("%s posterior-corrected litter by region  |  Run: %s",
                  MODEL_NAME, RUN_ID),
          side = 3, outer = TRUE, line = -1.5, cex = 1.0, font = 2)
    dev.off()
    message(sprintf("Posterior litter plot: %s", litter_png))
  }
  
} else {
  message("Skipping posterior litter plot (sigma_input or litter data unavailable)")
}


# =============================================================================
# 6.  Append [6] Residual analysis to report
# =============================================================================

report_text <- c(
  "\n",
  "[6] Residual analysis\n",
  sprintf("    Observations used:        %d\n", nrow(rf_df)),
  sprintf("    Covariates considered:    %d\n", length(available)),
  sprintf("    RF Out-Of-Bag R2:         %.3f\n", oob_r2),
  "    Top covariates by permutation importance:\n"
)
for (i in seq_len(min(10, nrow(imp_df)))) {
  report_text <- c(report_text,
                   sprintf("        %-25s  %.4f\n", imp_df$variable[i], imp_df$importance[i]))
}

if (!is.null(sigma_input_samples)) {
  report_text <- c(report_text,
                   sprintf("\n    sigma_input posterior:    median=%.3f  95%% CI=[%.3f, %.3f]\n",
                           median(sigma_input_samples),
                           quantile(sigma_input_samples, 0.025),
                           quantile(sigma_input_samples, 0.975)))
}

report_text <- c(report_text, "\n")
append_to_report(run_config, paste(report_text, collapse = ""))

message("\nResidual analysis complete.")


# =============================================================================
# 7.  Histogram of log-residuals
# =============================================================================

resid_finite <- resid_df$residual_log[is.finite(resid_df$residual_log)]

hist_path <- file.path(DIR_DIAG,
                       sprintf("%s_residuals_histogram_%s.png",
                               MODEL_NAME, RUN_ID))
png(hist_path,
    width = 8L * PX_PER_IN, height = 6L * PX_PER_IN, res = PX_PER_IN)
par(mar = c(4.5, 4.5, 3, 1))

hist(resid_finite,
     breaks = 40,
     col    = "steelblue", border = "white",
     xlab   = "Residual  log(obs) - log(pred)",
     ylab   = "Number of observations",
     main   = sprintf("%s (run %s) -- residual distribution  (n = %d)",
                      MODEL_NAME, RUN_ID, length(resid_finite)))
abline(v = 0,                    lty = 2, lwd = 2, col = "grey30")
abline(v = mean(resid_finite),   lty = 1, lwd = 2, col = "firebrick")
abline(v = median(resid_finite), lty = 1, lwd = 2, col = "darkorange")
legend("topleft",
       legend = c("zero",
                  sprintf("mean   = %.3f", mean(resid_finite)),
                  sprintf("median = %.3f", median(resid_finite))),
       col    = c("grey30", "firebrick", "darkorange"),
       lty    = c(2, 1, 1), lwd = 2, bty = "n")
dev.off()

message(sprintf("Residual histogram: [%s]", hist_path))
message(sprintf("    Residual mean (log):    %+.3f", mean(resid_finite)))
message(sprintf("    Residual median (log):  %+.3f", median(resid_finite)))
message(sprintf("    Residual SD (log):      %.3f",  sd(resid_finite)))