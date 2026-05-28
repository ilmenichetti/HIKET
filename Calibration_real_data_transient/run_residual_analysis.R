# =============================================================================
# run_residual_analysis.R
#
# Model-agnostic residual analysis stage of the HIKET pipeline.
# Reads the residuals CSV produced by run_<MODEL>_predictive.R, runs Random
# Forest variable importance and basic by-stratum residual visualisations,
# and appends a [6] Residual analysis section to the calibration report.
#
# The analysis is run twice in parallel structure:
#   - Calibration set  (is_holdout == FALSE)
#   - Holdout set      (is_holdout == TRUE)  -- independent validation
# Output files are identical in format; holdout versions carry a _holdout
# suffix in the filename and "Independent validation" in the title.
#
# USAGE:
#   Rscript run_residual_analysis.R <MODEL_NAME> <RUN_ID>
#   e.g.  Rscript run_residual_analysis.R Yasso07 20260427_213045
#
# OUTPUTS:
#   ./Calibration_real_data_transient/diagnostics/<MODEL>/
#     <MODEL>_residuals_diagnostic_<RUN_ID>.png         -- 3x3 panel, calibration
#     <MODEL>_residuals_diagnostic_holdout_<RUN_ID>.png -- 3x3 panel, holdout
#     <MODEL>_litter_posterior_<RUN_ID>.png             -- posterior-corrected litter
#     <MODEL>_residuals_histogram_<RUN_ID>.png          -- calibration
#     <MODEL>_residuals_histogram_holdout_<RUN_ID>.png  -- holdout
#     <MODEL>_rf_importance_<RUN_ID>.png/.csv           -- calibration
#     <MODEL>_rf_importance_holdout_<RUN_ID>.png/.csv   -- holdout
#     <MODEL>_rf_summary_<RUN_ID>.rds                   -- calibration
#     <MODEL>_rf_summary_holdout_<RUN_ID>.rds           -- holdout
#   Appends [6] Residual analysis to the calibration report.
# =============================================================================

library(ranger)
library(dplyr)

source("./Calibration_real_data_transient/calibration_engine_transient.R")


# =============================================================================
# 0.  Configuration / argument parsing
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  diag_root <- "./Calibration_real_data_transient/diagnostics"
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

DIR_DIAG <- file.path("./Calibration_real_data_transient/diagnostics", MODEL_NAME)
run_config <- list(MODEL_NAME = MODEL_NAME, RUN_ID = RUN_ID, DIR_DIAG = DIR_DIAG)

PX_PER_IN <- 150L

message("=============================================================")
message(sprintf("  %s Residual Analysis  |  Run: %s", MODEL_NAME, RUN_ID))
message("=============================================================\n")


# =============================================================================
# 1.  Load residuals
# =============================================================================

resid_csv <- file.path(DIR_DIAG, sprintf("%s_residuals_%s.csv", MODEL_NAME, RUN_ID))
if (!file.exists(resid_csv)) stop("Residuals not found: ", resid_csv)
resid_df <- read.csv(resid_csv, stringsAsFactors = FALSE)
message(sprintf("Loaded residuals: %d rows", nrow(resid_df)))

# Ensure is_holdout is logical; default FALSE if absent (old bundles)
if (!"is_holdout" %in% names(resid_df)) {
  message("WARNING: is_holdout column not found -- treating all rows as calibration")
  resid_df$is_holdout <- FALSE
}
resid_df$is_holdout <- as.logical(resid_df$is_holdout)

res_calib   <- resid_df[!resid_df$is_holdout, ]
res_holdout <- resid_df[ resid_df$is_holdout, ]

message(sprintf("  Calibration rows: %d  |  Holdout rows: %d",
                nrow(res_calib), nrow(res_holdout)))


# =============================================================================
# 1b. Load posterior samples and litter data (for litter plot, Section 5b)
# =============================================================================

posterior_phys <- tryCatch({
  f <- file.path("./Calibration_real_data_transient/runs/",
                 sprintf("%s_posterior_%s.rds", MODEL_NAME, RUN_ID))
  if (!file.exists(f)) stop("not found")
  readRDS(f)
}, error = function(e) {
  message(sprintf("WARNING: posterior RDS not found -- skipping litter plot (%s)",
                  e$message))
  NULL
})

sigma_input_samples <- if (!is.null(posterior_phys) &&
                           "sigma_input" %in% colnames(posterior_phys)) {
  posterior_phys[, "sigma_input"]
} else {
  message("WARNING: sigma_input not found in posterior -- skipping litter plot")
  NULL
}

input_raw <- tryCatch(
  read.csv("./Data/model_inputs/input_raw_monthly.csv"),
  error = function(e) { message("WARNING: input_raw_monthly.csv not found"); NULL }
)


# =============================================================================
# 2.  Covariate list (shared between calibration and holdout RF runs)
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


# =============================================================================
# 3.  Helper: run RF + produce importance plot + 3x3 panel + histogram
# =============================================================================
# suffix : ""          for calibration outputs
#          "_holdout"  for holdout outputs
# label  : "Calibration" or "Independent validation" (used in plot titles)

run_rf_and_plots <- function(df, suffix, label) {

  tag <- if (nchar(suffix) > 0) sprintf(" [%s]", label) else ""
  message(sprintf("\n--- RF analysis%s ---", tag))

  # --- Covariate selection ---------------------------------------------------
  available <- intersect(candidate_vars, names(df))
  unused    <- setdiff(candidate_vars, names(df))
  if (length(unused) > 0)
    message(sprintf("Skipping %d not-found covariates: %s",
                    length(unused), paste(unused, collapse = ", ")))

  rf_df <- df[, c(target, available), drop = FALSE]
  for (v in available)
    if (is.character(rf_df[[v]])) rf_df[[v]] <- as.factor(rf_df[[v]])

  complete_rows <- is.finite(rf_df[[target]])
  n_dropped     <- sum(!complete_rows)
  rf_df         <- rf_df[complete_rows, , drop = FALSE]
  message(sprintf("RF data: %d rows x %d covariates (dropped %d non-finite targets)",
                  nrow(rf_df), length(available), n_dropped))

  # Drop columns with >10% missing
  na_frac <- colMeans(is.na(rf_df))
  rf_df   <- rf_df[, na_frac < 0.10]
  rf_df   <- rf_df[complete.cases(rf_df), ]
  message(sprintf("RF data after NA removal: %d rows x %d cols",
                  nrow(rf_df), ncol(rf_df)))

  if (nrow(rf_df) < 10) {
    message(sprintf("WARNING: too few rows (%d) for RF%s -- skipping",
                    nrow(rf_df), tag))
    return(NULL)
  }

  # --- RF fit ----------------------------------------------------------------
  set.seed(2025)
  rf_fit <- ranger(
    formula     = as.formula(paste(target, "~ .")),
    data        = rf_df,
    num.trees   = 500L,
    importance  = "permutation",
    num.threads = parallel::detectCores() - 1L
  )
  oob_r2 <- rf_fit$r.squared

  imp    <- sort(rf_fit$variable.importance, decreasing = TRUE)
  imp_df <- data.frame(variable = names(imp), importance = unname(imp))

  write.csv(imp_df,
            file.path(DIR_DIAG,
                      sprintf("%s_rf_importance%s_%s.csv",
                              MODEL_NAME, suffix, RUN_ID)),
            row.names = FALSE)

  rf_summary <- list(rf_fit     = rf_fit,
                     oob_r2     = oob_r2,
                     importance = imp_df,
                     n_obs      = nrow(rf_df),
                     covariates = available,
                     target     = target,
                     subset     = label)
  saveRDS(rf_summary,
          file.path(DIR_DIAG,
                    sprintf("%s_rf_summary%s_%s.rds",
                            MODEL_NAME, suffix, RUN_ID)))
  message(sprintf("RF OOB R²%s: %.3f", tag, oob_r2))

  # --- RF importance plot ----------------------------------------------------
  n_show  <- min(20L, nrow(imp_df))
  imp_top <- imp_df[seq_len(n_show), ]
  imp_top <- imp_top[order(imp_top$importance), ]

  png(file.path(DIR_DIAG,
                sprintf("%s_rf_importance%s_%s.png",
                        MODEL_NAME, suffix, RUN_ID)),
      width  = 10L * PX_PER_IN,
      height = max(5L, n_show * 0.3) * PX_PER_IN,
      res    = PX_PER_IN)
  par(mar = c(4, 12, 3, 1))
  barplot(imp_top$importance, horiz = TRUE,
          names.arg = imp_top$variable, las = 1, cex.names = 0.85,
          col = "steelblue", border = "white",
          xlab = "Permutation importance",
          main = sprintf("%s residual RF — %s — top %d  |  OOB R² = %.3f",
                         MODEL_NAME, label, n_show, oob_r2))
  dev.off()

  # --- 3x3 diagnostic panel --------------------------------------------------
  resid_png <- file.path(DIR_DIAG,
                         sprintf("%s_residuals_diagnostic%s_%s.png",
                                 MODEL_NAME, suffix, RUN_ID))
  png(resid_png, width = 15L * PX_PER_IN, height = 15L * PX_PER_IN, res = PX_PER_IN)
  par(mfrow = c(3, 3), mar = c(4, 4, 3, 1))

  resid  <- df$residual_log
  fitted <- df$log_hat_mean

  finite_mask <- is.finite(fitted) & is.finite(resid)
  fitted_f    <- fitted[finite_mask]
  resid_f     <- resid[finite_mask]
  df_f        <- df[finite_mask, , drop = FALSE]

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
  boxplot(resid_f ~ df_f$year, col = "lightsteelblue",
          xlab = "Year", ylab = "Residual (log)", main = "Residuals by year")
  abline(h = 0, lty = 2, col = "grey40")

  # Panel 3: by Cajander fertility class
  if ("KA" %in% names(df_f)) {
    boxplot(resid_f ~ as.factor(df_f$KA), col = "lightsteelblue",
            xlab = "Cajander fertility class (KA)", ylab = "Residual (log)",
            main = "Residuals by Cajander class")
    abline(h = 0, lty = 2, col = "grey40")
  } else { plot.new(); title("KA not available") }

  # Panel 4: vs stoniness
  if ("CoarseFragments" %in% names(df_f)) {
    scatter_loess(df_f$CoarseFragments, resid_f,
                  xlab = "Coarse fragments / stoniness (%)",
                  main = "Residuals vs stoniness")
  } else { plot.new(); title("CoarseFragments not available") }

  # Panel 5: by dominant species
  if ("species_code" %in% names(df_f)) {
    sp_labels <- c("1" = "Scots pine", "2" = "Nor. spruce", "3" = "Birch",
                   "4" = "Aspen",      "5" = "Grey alder",  "6" = "Blk. alder",
                   "7" = "Other brdlv.")
    sp_factor <- as.factor(df_f$species_code)
    levels(sp_factor) <- sp_labels[levels(sp_factor)]
    boxplot(resid_f ~ sp_factor, col = "lightsteelblue", las = 2,
            xlab = "", ylab = "Residual (log)",
            main = "Residuals by dominant species", cex.axis = 0.78)
    abline(h = 0, lty = 2, col = "grey40")
  } else { plot.new(); title("species_code not available") }

  # Panel 6: Q-Q plot
  qqnorm(resid_f, pch = 19, col = adjustcolor("steelblue", 0.5),
         main = "Q-Q plot of log residuals")
  qqline(resid_f, col = "firebrick", lwd = 2)

  # Panel 7: vs profile depth
  if ("profile_depth_cm" %in% names(df_f)) {
    scatter_loess(df_f$profile_depth_cm, resid_f,
                  xlab = "Measured profile depth (cm)",
                  main = "Residuals vs profile depth")
  } else { plot.new(); title("profile_depth_cm not available") }

  # Panel 8: vs mean annual temperature
  if ("mean_temp" %in% names(df_f)) {
    scatter_loess(df_f$mean_temp, resid_f,
                  xlab = "Mean annual temperature (C)",
                  main = "Residuals vs mean temperature")
  } else { plot.new(); title("mean_temp not available") }

  # Panel 9: vs CN ratio
  if ("CN_ratio" %in% names(df_f)) {
    scatter_loess(df_f$CN_ratio, resid_f,
                  xlab = "CN ratio (OrganicCarbon / TotalNitrogen)",
                  main = "Residuals vs CN ratio", span = 0.85)
  } else { plot.new(); title("CN_ratio not available") }

  mtext(sprintf("%s residual diagnostics — %s  |  Run: %s",
                MODEL_NAME, label, RUN_ID),
        side = 3, outer = TRUE, line = -1.5, cex = 1.0, font = 2)
  dev.off()
  message(sprintf("Diagnostic panel%s: %s", tag, resid_png))

  # --- Histogram -------------------------------------------------------------
  resid_finite <- resid[is.finite(resid)]

  hist_path <- file.path(DIR_DIAG,
                         sprintf("%s_residuals_histogram%s_%s.png",
                                 MODEL_NAME, suffix, RUN_ID))
  png(hist_path, width = 8L * PX_PER_IN, height = 6L * PX_PER_IN, res = PX_PER_IN)
  par(mar = c(4.5, 4.5, 3, 1))
  hist(resid_finite,
       breaks = 40,
       col    = "steelblue", border = "white",
       xlab   = "Residual  log(obs) - log(pred)",
       ylab   = "Number of observations",
       main   = sprintf("%s — %s — residual distribution  (n = %d, run %s)",
                        MODEL_NAME, label, length(resid_finite), RUN_ID))
  abline(v = 0,                    lty = 2, lwd = 2, col = "grey30")
  abline(v = mean(resid_finite),   lty = 1, lwd = 2, col = "firebrick")
  abline(v = median(resid_finite), lty = 1, lwd = 2, col = "darkorange")
  legend("topleft",
         legend = c("zero",
                    sprintf("mean   = %.3f", mean(resid_finite)),
                    sprintf("median = %.3f", median(resid_finite))),
         col = c("grey30", "firebrick", "darkorange"),
         lty = c(2, 1, 1), lwd = 2, bty = "n")
  dev.off()
  message(sprintf("Histogram%s: mean=%.3f  median=%.3f  sd=%.3f",
                  tag, mean(resid_finite), median(resid_finite),
                  sd(resid_finite)))

  list(oob_r2 = oob_r2, imp_df = imp_df, n_obs = nrow(rf_df))
}


# =============================================================================
# 4.  Run RF + plots for calibration and holdout subsets
# =============================================================================

rf_calib   <- run_rf_and_plots(res_calib,   suffix = "",         label = "Calibration")
rf_holdout <- run_rf_and_plots(res_holdout, suffix = "_holdout", label = "Independent validation")


# =============================================================================
# 5.  Posterior-corrected litter trajectories by region (calibration only)
#     No holdout mirror: sigma_input is a global model parameter, not
#     specific to either subset.
# =============================================================================

if (!is.null(sigma_input_samples) && !is.null(input_raw)) {

  si_med <- median(sigma_input_samples)
  si_lo  <- quantile(sigma_input_samples, 0.025)
  si_hi  <- quantile(sigma_input_samples, 0.975)

  message(sprintf("\nsigma_input posterior: median=%.3f  95%% CI=[%.3f, %.3f]",
                  si_med, si_lo, si_hi))

  litter_cols <- grep("^C_(nwl|fwl|cwl)_[AWEN]$", names(input_raw), value = TRUE)

  calib_plot_ids  <- unique(res_calib$plot_id)
  input_calib_lit <- input_raw[
    input_raw$plot_id %in% calib_plot_ids &
    input_raw$year    != 2024L, ]

  annual_lit <- input_calib_lit %>%
    group_by(plot_id, year) %>%
    summarise(total_litter = sum(across(all_of(litter_cols)), na.rm = TRUE) * 12,
              .groups = "drop")

  region_lkup <- unique(res_calib[, c("plot_id", "region")])
  annual_lit  <- merge(annual_lit, region_lkup, by = "plot_id", all.x = TRUE)
  annual_lit  <- annual_lit[!is.na(annual_lit$region), ]

  regions   <- sort(unique(annual_lit$region))
  n_regions <- length(regions)

  if (n_regions == 0) {
    message("WARNING: no region labels found -- skipping litter plot")
  } else {
    mean_lit <- annual_lit %>%
      group_by(region, year) %>%
      summarise(mean_raw = mean(total_litter, na.rm = TRUE), .groups = "drop")

    reg_cols <- c("#2166ac", "#d6604d", "#1a9641", "#762a83", "#e08214", "#4d4d4d")
    reg_cols <- setNames(reg_cols[seq_len(n_regions)], regions)

    litter_png <- file.path(DIR_DIAG,
                            sprintf("%s_litter_posterior_%s.png", MODEL_NAME, RUN_ID))
    png(litter_png,
        width  = max(8L, 5L * n_regions) * PX_PER_IN,
        height = 6L * PX_PER_IN, res = PX_PER_IN)
    par(mfrow = c(1, n_regions), mar = c(4, 4.5, 4, 1))

    for (reg in regions) {
      reg_dat  <- mean_lit[mean_lit$region == reg, ]
      reg_dat  <- reg_dat[order(reg_dat$year), ]
      yrs      <- reg_dat$year
      raw_mn   <- reg_dat$mean_raw
      corr_md  <- raw_mn * si_med
      corr_lo  <- raw_mn * si_lo
      corr_hi  <- raw_mn * si_hi
      ylim     <- range(c(raw_mn, corr_lo, corr_hi), na.rm = TRUE)
      ylim     <- ylim + c(-0.05, 0.05) * diff(ylim)
      col_reg  <- reg_cols[as.character(reg)]

      plot(yrs, raw_mn, type = "l", lty = 2, lwd = 1.5, col = "grey50",
           ylim = ylim, xlab = "Year",
           ylab = "Mean annual litter input (tC/ha/yr)",
           main = sprintf("Region: %s\nsigma_input: %.3f  [%.3f, %.3f]",
                          reg, si_med, si_lo, si_hi))
      polygon(c(yrs, rev(yrs)), c(corr_lo, rev(corr_hi)),
              col = adjustcolor(col_reg, 0.18), border = NA)
      lines(yrs, corr_md, col = col_reg, lwd = 2.5)
      abline(h = mean(raw_mn, na.rm = TRUE), lty = 3, col = "grey70", lwd = 1)
      legend("topright",
             legend = c("Raw litter (uncorrected)",
                        sprintf("Posterior median (x %.3f)", si_med),
                        "95% credible band"),
             col  = c("grey50", col_reg, adjustcolor(col_reg, 0.4)),
             lty  = c(2, 1, NA), lwd = c(1.5, 2.5, NA),
             pch  = c(NA, NA, 15), pt.cex = c(NA, NA, 2),
             bty  = "n", cex = 0.82)
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
# 6.  Append to calibration report
# =============================================================================

oob_calib   <- if (!is.null(rf_calib))   rf_calib$oob_r2   else NA
oob_holdout <- if (!is.null(rf_holdout)) rf_holdout$oob_r2 else NA

report_text <- c(
  "\n",
  "[6] Residual analysis\n",
  sprintf("    Calibration obs:          %d\n", nrow(res_calib)),
  sprintf("    Holdout obs:              %d\n", nrow(res_holdout)),
  sprintf("    RF OOB R² (calibration):  %.3f\n", oob_calib),
  sprintf("    RF OOB R² (holdout):      %.3f\n", oob_holdout),
  "    Top calibration covariates (permutation importance):\n"
)

if (!is.null(rf_calib)) {
  for (i in seq_len(min(10, nrow(rf_calib$imp_df))))
    report_text <- c(report_text,
                     sprintf("        %-25s  %.4f\n",
                             rf_calib$imp_df$variable[i],
                             rf_calib$imp_df$importance[i]))
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
