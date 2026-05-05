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
#   ./Calibration/diagnostics/<MODEL>/<MODEL>_residuals_diagnostic_<RUN_ID>.png
#   ./Calibration/diagnostics/<MODEL>/<MODEL>_rf_importance_<RUN_ID>.png
#   ./Calibration/diagnostics/<MODEL>/<MODEL>_rf_summary_<RUN_ID>.rds
#   Appends [6] Residual analysis to the report.
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
  # Pattern matches only the timestamp format (YYYYMMDD_HHMMSS) so we don't
  # accidentally pick up files like "<MODEL>_residuals_diagnostic_<RUN_ID>.png"
  # (which doesn't match the .csv extension anyway, but be defensive).
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
# 1.  Load residuals + site_raw (already merged in predictive script)
# =============================================================================

resid_csv <- file.path(DIR_DIAG, sprintf("%s_residuals_%s.csv",
                                         MODEL_NAME, RUN_ID))
if (!file.exists(resid_csv)) stop("Residuals not found: ", resid_csv)
resid_df <- read.csv(resid_csv, stringsAsFactors = FALSE)
message(sprintf("Loaded residuals: %d rows", nrow(resid_df)))


# =============================================================================
# 2.  RF covariate selection
# =============================================================================
# Excluded:
#   OrganicCarbon, OrganicMatter -- collinear with SOC obs by construction
#   soil_depth                   -- replaced by profile_depth_cm (measured)
#   MeanBulkDensity              -- 84% NA in calib_ready (audit)
#   plot_id, year                -- identifiers, not predictors
#   any soc_*, log_*, residual_* -- the targets and derived quantities

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

available  <- intersect(candidate_vars, names(resid_df))
unused     <- setdiff(candidate_vars, names(resid_df))
if (length(unused) > 0) {
  message(sprintf("Skipping %d not-found covariates: %s",
                  length(unused), paste(unused, collapse = ", ")))
}

# Build RF data: drop rows with NA in target or any selected covariate
rf_df <- resid_df[, c(target, available), drop = FALSE]

# Coerce categorical-looking columns to factor; coerce booleans to factor
for (v in available) {
  x <- rf_df[[v]]
  if (is.character(x) || is.logical(x)) rf_df[[v]] <- as.factor(x)
  if (is.integer(x) && length(unique(na.omit(x))) <= 12)
    rf_df[[v]] <- as.factor(as.character(x))
}

# Drop covariates with zero variance or all-NA after coercion
keep_vars <- vapply(available, function(v) {
  vec <- rf_df[[v]]
  if (all(is.na(vec))) return(FALSE)
  if (is.factor(vec)) return(nlevels(droplevels(vec)) > 1)
  length(unique(vec[!is.na(vec)])) > 1
}, logical(1))
available <- available[keep_vars]
rf_df     <- rf_df[, c(target, available), drop = FALSE]

# Drop rows with NA or Inf in the target (residual_log can be Inf if a
# posterior-mean prediction was 0 or negative). complete.cases handles NA
# but not Inf, so check both.
complete_rows <- complete.cases(rf_df) & is.finite(rf_df[[target]])
n_dropped <- sum(!complete_rows)
rf_df     <- rf_df[complete_rows, , drop = FALSE]
message(sprintf("RF data: %d rows × %d covariates (dropped %d incomplete/non-finite)",
                nrow(rf_df), length(available), n_dropped))


# =============================================================================
# 3.  Random Forest with permutation importance
# =============================================================================

set.seed(2025)   # standardised seed
rf_fit <- ranger(
  formula     = as.formula(paste(target, "~ .")),
  data        = rf_df,
  num.trees   = 500L,
  importance  = "permutation",
  num.threads = parallel::detectCores() - 1L
)

oob_r2 <- 1 - rf_fit$prediction.error / var(rf_df[[target]])
message(sprintf("\nRandom Forest OOB R²: %.3f", oob_r2))
message("(fraction of residual log-SOC variance explained by site covariates)")

imp <- sort(rf_fit$variable.importance, decreasing = TRUE)
imp_df <- data.frame(variable = names(imp), importance = unname(imp))
write.csv(imp_df,
          file.path(DIR_DIAG,
                    sprintf("%s_rf_importance_%s.csv", MODEL_NAME, RUN_ID)),
          row.names = FALSE)

# Save full RF object + summary
rf_summary <- list(
  rf_fit       = rf_fit,
  oob_r2       = oob_r2,
  importance   = imp_df,
  n_obs        = nrow(rf_df),
  covariates   = available,
  target       = target
)
saveRDS(rf_summary,
        file.path(DIR_DIAG,
                  sprintf("%s_rf_summary_%s.rds", MODEL_NAME, RUN_ID)))


# =============================================================================
# 4.  Plot: RF importance
# =============================================================================

n_show <- min(20L, nrow(imp_df))
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
        main = sprintf("%s residual analysis -- top %d covariates  |  OOB R² = %.3f",
                       MODEL_NAME, n_show, oob_r2))
dev.off()
message(sprintf("RF importance plot: %s",
                file.path(DIR_DIAG,
                          sprintf("%s_rf_importance_%s.png", MODEL_NAME, RUN_ID))))


# =============================================================================
# 5.  Diagnostic plot: residuals by year, by KA, by stoniness, by species; QQ
# =============================================================================
# Six panels in one PNG.

resid_png <- file.path(DIR_DIAG,
                       sprintf("%s_residuals_diagnostic_%s.png",
                               MODEL_NAME, RUN_ID))
png(resid_png, width = 12L * PX_PER_IN, height = 8L * PX_PER_IN, res = PX_PER_IN)
par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))

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

# Panel 1: residuals vs fitted with loess
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

# Panel 4: residuals vs CoarseFragments (stoniness, %)
# CoarseFragments from Biosoil mineral topsoil (M01 + M12 mean).
# Stoniness affects rooting depth, water retention, and bulk density estimates;
# expected to be a moderator of SOC accumulation in Finnish till soils.
# Loess added when enough non-NA obs are available.
if ("CoarseFragments" %in% names(resid_df_f)) {
  stone <- resid_df_f$CoarseFragments
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
  
  # Marginal rug on x-axis to show stoniness density
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

# Panel 6: QQ plot
qqnorm(resid_f, pch = 19, col = adjustcolor("steelblue", 0.5),
       main = "Q-Q plot of log residuals")
qqline(resid_f, col = "firebrick", lwd = 2)

mtext(sprintf("%s residual diagnostics  |  Run: %s", MODEL_NAME, RUN_ID),
      side = 3, outer = TRUE, line = -1.5, cex = 1.0, font = 2)
dev.off()
message(sprintf("Residual diagnostic: %s", resid_png))





# =============================================================================
# 6.  Append [6] Residual analysis to report
# =============================================================================

report_text <- c(
  "\n",
  "[6] Residual analysis\n",
  sprintf("    Observations used:        %d\n", nrow(rf_df)),
  sprintf("    Covariates considered:    %d\n", length(available)),
  sprintf("    RF Out-Of-Bag R²:         %.3f\n", oob_r2),
  "    Top covariates by permutation importance:\n"
)
for (i in seq_len(min(10, nrow(imp_df)))) {
  report_text <- c(report_text,
                   sprintf("        %-25s  %.4f\n", imp_df$variable[i], imp_df$importance[i]))
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
     main   = sprintf("%s (run %s) — residual distribution  (n = %d)",
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
message(sprintf("    Residual mean (log):    %+.3f\n", mean(resid_finite)),
        sprintf("    Residual median (log):  %+.3f\n", median(resid_finite)),
        sprintf("    Residual SD (log):      %.3f\n",  sd(resid_finite)))

