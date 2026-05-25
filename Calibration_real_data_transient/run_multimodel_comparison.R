# =============================================================================
# run_multimodel_comparison.R
#
# Multi-model comparison stage of the HIKET pipeline (transient version).
# Loads the posterior predictive bundles from all five models (SP1, TP2,
# Yasso07, Yasso15, Yasso20) and produces a single comparison report and
# figure set.
#
# USAGE:
#   Rscript run_multimodel_comparison.R \
#     SP1:<RUN_ID_SP1> TP2:<RUN_ID_TP2> \
#     Yasso07:<RUN_ID_07> Yasso15:<RUN_ID_15> Yasso20:<RUN_ID_20>
#
#   or with no arguments to auto-detect the most recent run for each model.
#
# INPUTS (one per model):
#   ./Calibration_real_data_transient/runs/<MODEL>_posterior_predictive_<RUN_ID>.rds
#
# OUTPUTS:
#   ./Calibration_real_data_transient/diagnostics/multimodel/
#     multimodel_metrics_<COMP_ID>.txt          -- formatted metrics table
#     multimodel_obs_vs_pred_<COMP_ID>.png      -- 2x3 scatter + bias panel
#     multimodel_delta_soc_<COMP_ID>.png        -- annual ΔSOC rate trajectories
#     multimodel_residuals_<COMP_ID>.png        -- residual distributions + by KA
#     multimodel_summary_<COMP_ID>.rds          -- metrics list for downstream use
#
#   where COMP_ID = <RUN_ID_SP1>_<RUN_ID_TP2>_<RUN_ID_07>_<RUN_ID_15>_<RUN_ID_20>
# =============================================================================

library(dplyr)

MODELS <- c("SP1", "TP2", "Yasso07", "Yasso15", "Yasso20")
DIR_RUNS <- "./Calibration_real_data_transient/runs/"
DIR_OUT  <- "./Calibration_real_data_transient/diagnostics/multimodel"
dir.create(DIR_OUT, showWarnings = FALSE, recursive = TRUE)
DIR_DIAG <- "./Calibration_real_data_transient/diagnostics"

PX_PER_IN <- 150L

# Model colours used consistently across all comparison plots.
# SP1/TP2 use brown tones to visually distinguish the simpler benchmark
# models from the three YASSO variants (blue/red/green).
MODEL_COLS <- c(
  SP1     = "#8c510a",   # brown  -- 1-pool baseline
  TP2     = "#bf812d",   # tan    -- 2-pool intermediate
  Yasso07 = "#2166ac",   # blue
  Yasso15 = "#d6604d",   # red
  Yasso20 = "#1a9641"    # green
)


# =============================================================================
# 0.  Parse arguments / auto-detect RUN_IDs
# =============================================================================

detect_latest <- function(model) {
  pat <- sprintf("^%s_posterior_predictive_[0-9]{8}_[0-9]{6}\\.rds$", model)
  fns <- list.files(DIR_RUNS, pattern = pat)
  if (length(fns) == 0)
    stop(sprintf("No posterior predictive found for %s in %s", model, DIR_RUNS))
  fn <- sort(fns, decreasing = TRUE)[1]
  sub(sprintf("^%s_posterior_predictive_(.+)\\.rds$", model), "\\1", fn)
}

args <- commandArgs(trailingOnly = TRUE)
run_ids <- setNames(vector("list", length(MODELS)), MODELS)

if (length(args) >= 1) {
  # Parse "Model:RUN_ID" pairs from command line
  for (a in args) {
    parts <- strsplit(a, ":", fixed = TRUE)[[1]]
    if (length(parts) == 2 && parts[1] %in% MODELS)
      run_ids[[parts[1]]] <- parts[2]
  }
}

# Fall back to auto-detect for any model not supplied
for (m in MODELS) {
  if (is.null(run_ids[[m]])) {
    run_ids[[m]] <- detect_latest(m)
    message(sprintf("Auto-detected %s RUN_ID: %s", m, run_ids[[m]]))
  }
}

COMP_ID <- paste(run_ids, collapse = "_")

message("=============================================================")
message("  HIKET Multi-Model Comparison")
for (m in MODELS) message(sprintf("    %-8s  %s", m, run_ids[[m]]))
message(sprintf("  COMP_ID: %s", COMP_ID))
message("=============================================================\n")


# =============================================================================
# 1.  Load posterior predictive bundles
# =============================================================================

pp <- lapply(MODELS, function(m) {
  f <- file.path(DIR_RUNS,
                 sprintf("%s_posterior_predictive_%s.rds", m, run_ids[[m]]))
  if (!file.exists(f)) stop(sprintf("File not found: %s", f))
  message(sprintf("Loading %s ...", f))
  readRDS(f)
})
names(pp) <- MODELS

# Quick check: each bundle has the expected structure
for (m in MODELS) {
  stopifnot(!is.null(pp[[m]]$metrics),
            !is.null(pp[[m]]$residuals_df),
            !is.null(pp[[m]]$posterior_predictions))
}
message("All bundles loaded and validated.\n")


# =============================================================================
# 2.  Metrics table
# =============================================================================

metrics_df <- do.call(rbind, lapply(MODELS, function(m) {
  mt <- pp[[m]]$metrics
  data.frame(
    Model       = m,
    R2          = mt$R2,
    RMSE_median = mt$RMSE_median,
    Bias_median = mt$bias_median,
    RMSE_mean   = mt$RMSE_mean,
    Bias_mean   = mt$bias_mean,
    Coverage_95 = mt$coverage_95,
    N_obs       = nrow(pp[[m]]$residuals_df),
    stringsAsFactors = FALSE
  )
}))

metrics_txt <- file.path(DIR_OUT, sprintf("multimodel_metrics_%s.txt", COMP_ID))
sink(metrics_txt)
cat("HIKET Multi-Model Comparison — Predictive Metrics\n")
cat(sprintf("Generated: %s\n\n", format(Sys.time())))
for (m in MODELS)
  cat(sprintf("  %-8s  RUN_ID: %s\n", m, run_ids[[m]]))
cat("\n")
cat(sprintf("%-10s  %6s  %10s  %10s  %10s\n",
            "Model", "R²", "RMSE (med)", "Bias (med)", "95% cov"))
cat(strrep("-", 54), "\n")
for (i in seq_len(nrow(metrics_df))) {
  cat(sprintf("%-10s  %6.3f  %10.2f  %+10.2f  %10.3f\n",
              metrics_df$Model[i],
              metrics_df$R2[i],
              metrics_df$RMSE_median[i],
              metrics_df$Bias_median[i],
              metrics_df$Coverage_95[i]))
}
cat("\n")
sink()
message(sprintf("Metrics table: %s", metrics_txt))
print(metrics_df[, c("Model","R2","RMSE_median","Bias_median","Coverage_95")])


# =============================================================================
# 3.  Plot A: Obs vs predicted — 2x3 grid (5 model panels + bias barplot)
# =============================================================================

obs_pred_png <- file.path(DIR_OUT,
                          sprintf("multimodel_obs_vs_pred_%s.png", COMP_ID))

# Shared axis limits across all five models for comparability
ax_max <- max(sapply(MODELS, function(m) {
  max(c(pp[[m]]$residuals_df$soc_obs_tCha,
        pp[[m]]$residuals_df$soc_q975), na.rm = TRUE)
})) * 1.05
ax_lim <- c(0, ax_max)

png(obs_pred_png,
    width  = 18L * PX_PER_IN,
    height = 14L * PX_PER_IN,
    res    = PX_PER_IN)
par(mfrow = c(2, 3), mar = c(4.5, 4.5, 3.5, 1))

for (m in MODELS) {
  rd <- pp[[m]]$residuals_df
  mt <- pp[[m]]$metrics
  
  ka_col <- c("1" = "#2166ac","2" = "#74add1",
              "3" = "#f4a582","4" = "#d6604d")[as.character(rd$KA)]
  ka_col[is.na(ka_col)] <- "grey60"
  
  plot(NA, xlim = ax_lim, ylim = ax_lim,
       xlab = "Observed SOC (tC/ha)",
       ylab = "Predicted SOC — posterior median (tC/ha)",
       main = m)
  abline(0, 1, lty = 2, col = "grey40", lwd = 1.5)
  segments(rd$soc_obs_tCha, rd$soc_q025, rd$soc_obs_tCha, rd$soc_q975,
           col = adjustcolor(ka_col, 0.25), lwd = 0.7)
  points(rd$soc_obs_tCha, rd$soc_median,
         pch = 16, cex = 0.8, col = adjustcolor(ka_col, 0.65))
  legend("topleft",
         legend = c(sprintf("R² = %.3f",      mt$R2),
                    sprintf("RMSE = %.1f",     mt$RMSE_median),
                    sprintf("Bias = %+.1f",    mt$bias_median),
                    sprintf("Cov = %.2f",      mt$coverage_95)),
         bty = "n", cex = 0.8)
}

# Panel 6: bias barplot across all five models.
# Bias_median = median(log(obs) - log(pred)) across all observations for that
# model. Positive = systematic under-prediction; negative = over-prediction.
# Reference line at 0 = unbiased. Models ordered by complexity (SP1 -> Yasso20)
# so complexity-vs-bias trends are immediately visible.
bias_vals <- sapply(MODELS, function(m) pp[[m]]$metrics$bias_median)
ylim_bias <- c(min(bias_vals, 0), max(bias_vals, 0)) *
  c(ifelse(min(bias_vals) < 0, 1.25, 1), ifelse(max(bias_vals) > 0, 1.25, 1))

bp <- barplot(bias_vals,
              col       = MODEL_COLS[MODELS],
              names.arg = MODELS,
              ylab      = "Bias — posterior median  log(obs/pred)",
              main      = "Systematic bias by model",
              border    = "white",
              ylim      = ylim_bias,
              las       = 1)
abline(h = 0, lty = 2, col = "grey40", lwd = 1.5)
text(x      = bp,
     y      = bias_vals + sign(bias_vals) * diff(ylim_bias) * 0.04,
     labels = sprintf("%+.3f", bias_vals),
     cex    = 0.9, font = 2)

mtext("HIKET — Observed vs Predicted SOC (posterior median)",
      side = 3, outer = TRUE, line = -1.5, cex = 1.1, font = 2)
dev.off()
message(sprintf("Obs vs pred: %s", obs_pred_png))


# =============================================================================
# 4.  Plot B: ΔSOC trajectories overlaid
# =============================================================================
# Both the model trajectories and the observed reference are expressed as an
# annual rate of change relative to each plot's first simulation year:
#
#   annual_rate(t) = (SOC(t) - SOC(t0)) / (t - t0)    [tC/ha/yr]
#
# This is identical to the formula used in the individual predictive scripts
# (run_Yasso07_predictive.R, etc.) and makes the model trajectories directly
# comparable to the observed metric, which is also an annual rate computed
# over each plot's monitoring window.  The first year is dropped (0/0
# undefined); the trajectory therefore starts one year after t0 and is
# expected to converge toward the long-run observed rate as time increases.

# Build trajectory summary per model
traj_list <- lapply(MODELS, function(m) {
  pp[[m]]$posterior_predictions %>%
    group_by(plot_id, draw) %>%
    arrange(year, .by_group = TRUE) %>%
    mutate(
      first_soc  = first(total_soc),
      first_year = first(year),
      annual_rate = (total_soc - first_soc) / (year - first_year)
    ) %>%
    filter(year != first_year) %>%   # drop t=0 (rate is 0/0 = undefined)
    ungroup() %>%
    group_by(draw, year) %>%
    summarise(mean_annual_rate = mean(annual_rate, na.rm = TRUE),
              .groups = "drop") %>%
    group_by(year) %>%
    summarise(
      median_delta = median(mean_annual_rate, na.rm = TRUE),
      q025         = quantile(mean_annual_rate, 0.025, na.rm = TRUE),
      q975         = quantile(mean_annual_rate, 0.975, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(year) %>%
    mutate(model = m)
})
names(traj_list) <- MODELS

# Observed annual rate per plot: (SOC_last - SOC_first) / (year_last - year_first).
# Uses the first model's bundle since observations are the same across models.
obs_ref <- pp[[MODELS[1]]]$posterior_summary
obs_delta <- obs_ref %>%
  filter(!is.na(soc_obs_tCha)) %>%
  group_by(plot_id) %>%
  arrange(year, .by_group = TRUE) %>%
  summarise(
    delta_obs  = (last(soc_obs_tCha) - first(soc_obs_tCha)) /
      (last(year) - first(year)),
    year_first = first(year), year_last = last(year),
    .groups = "drop"
  ) %>% filter(year_first != year_last)

obs_mean_delta <- mean(obs_delta$delta_obs, na.rm = TRUE)
obs_se_delta   <- sd(obs_delta$delta_obs, na.rm = TRUE) / sqrt(nrow(obs_delta))

# Shared y limits
all_vals <- c(unlist(lapply(traj_list, function(t) c(t$q025, t$q975))),
              obs_mean_delta + c(-1.96, 1.96) * obs_se_delta, 0)
ylim_r <- range(all_vals, na.rm = TRUE)
ylim_r <- ylim_r + c(-0.05, 0.05) * diff(ylim_r)

all_years <- sort(unique(unlist(lapply(traj_list, `[[`, "year"))))

delta_png <- file.path(DIR_OUT, sprintf("multimodel_delta_soc_%s.png", COMP_ID))
png(delta_png, width = 13L * PX_PER_IN, height = 7L * PX_PER_IN, res = PX_PER_IN)
par(mar = c(4, 5, 3.5, 1))

plot(NA, xlim = range(all_years), ylim = ylim_r,
     xlab = "Year",
     ylab = expression(paste("Mean annual ",
                             Delta, "SOC / years since t"[0],
                             "  (tC/ha/yr)")),
     main = "HIKET — Mean annual ΔSOC rate across models")
abline(h = 0, lty = 3, col = "grey50")

# 95% CI ribbons then median lines per model
for (m in MODELS) {
  tr  <- traj_list[[m]]
  col <- MODEL_COLS[m]
  polygon(c(tr$year, rev(tr$year)),
          c(tr$q025, rev(tr$q975)),
          col = adjustcolor(col, 0.12), border = NA)
}
for (m in MODELS) {
  tr  <- traj_list[[m]]
  col <- MODEL_COLS[m]
  lines(tr$year, tr$median_delta, col = col, lwd = 2.5)
}

# Observed mean annual rate as horizontal reference line spanning the obs window
abline(h = obs_mean_delta, col = "grey20", lwd = 2, lty = 2)
rect(xleft   = min(obs_delta$year_first),
     xright  = max(obs_delta$year_last),
     ybottom = obs_mean_delta - 1.96 * obs_se_delta,
     ytop    = obs_mean_delta + 1.96 * obs_se_delta,
     col = adjustcolor("grey20", 0.10), border = NA)

legend("topleft",
       legend = c(MODELS,
                  sprintf("Observed mean annual ΔSOC ± 95%% CI  (n=%d)", nrow(obs_delta))),
       col    = c(MODEL_COLS[MODELS], "grey20"),
       lwd    = c(rep(2.5, length(MODELS)), 2),
       lty    = c(rep(1,   length(MODELS)), 2),
       bty = "n", cex = 0.9)
dev.off()
message(sprintf("ΔSOC trajectory: %s", delta_png))


# =============================================================================
# 5.  Plot C: Residual distributions + by KA  (2 rows, 2 cols)
# =============================================================================

resid_png <- file.path(DIR_OUT, sprintf("multimodel_residuals_%s.png", COMP_ID))
png(resid_png, width = 14L * PX_PER_IN, height = 10L * PX_PER_IN, res = PX_PER_IN)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

# Panel 1: overlaid log-residual density curves
finite_resids <- lapply(MODELS, function(m) {
  r <- pp[[m]]$residuals_df$residual_log
  r[is.finite(r)]
})
names(finite_resids) <- MODELS

x_range <- range(unlist(finite_resids), na.rm = TRUE)
x_range <- c(max(x_range[1], -3), min(x_range[2], 3))  # clip extreme tails

y_max <- 0
dens_list <- lapply(MODELS, function(m) {
  d <- density(finite_resids[[m]], from = x_range[1], to = x_range[2])
  y_max <<- max(y_max, max(d$y))
  d
})
# Rescale y axis then draw
plot(NA, xlim = x_range, ylim = c(0, y_max * 1.1),
     xlab = "Residual  log(obs) - log(pred)",
     ylab = "Density", main = "Residual distributions")
abline(v = 0, lty = 2, col = "grey40")
for (m in MODELS)
  lines(dens_list[[m]], col = MODEL_COLS[m], lwd = 2)
legend("topright",
       legend = sprintf("%s  (mean=%+.3f)", MODELS,
                        sapply(finite_resids, mean)),
       col = MODEL_COLS[MODELS], lwd = 2, bty = "n", cex = 0.85)

# Panel 2: residual bias by KA class (grouped boxplot)
# Build combined data frame with model column
resid_combined <- do.call(rbind, lapply(MODELS, function(m) {
  rd <- pp[[m]]$residuals_df
  data.frame(model = m, KA = rd$KA,
             residual_log = rd$residual_log,
             stringsAsFactors = FALSE)
}))
resid_combined <- resid_combined[is.finite(resid_combined$residual_log), ]

ka_levels <- sort(unique(na.omit(resid_combined$KA)))
if (length(ka_levels) > 0) {
  n_ka  <- length(ka_levels)
  n_mod <- length(MODELS)
  width <- 0.14   # narrower to keep five models legible within each KA cluster
  
  x_centres <- seq_len(n_ka)
  offsets    <- seq(-(n_mod - 1)/2, (n_mod - 1)/2, length.out = n_mod) * width * 1.5
  
  plot(NA,
       xlim = c(0.5, n_ka + 0.5),
       ylim = range(resid_combined$residual_log, na.rm = TRUE),
       xlab = "Cajander fertility class (KA)",
       ylab = "Residual (log)",
       main = "Residuals by KA class",
       xaxt = "n")
  axis(1, at = x_centres, labels = paste("KA", ka_levels))
  abline(h = 0, lty = 2, col = "grey40")
  
  for (j in seq_along(MODELS)) {
    m <- MODELS[j]
    for (i in seq_along(ka_levels)) {
      ka  <- ka_levels[i]
      vals <- resid_combined$residual_log[
        resid_combined$model == m & resid_combined$KA == ka &
          is.finite(resid_combined$residual_log)]
      if (length(vals) < 3) next
      xc <- x_centres[i] + offsets[j]
      boxplot(vals, at = xc, add = TRUE, boxwex = width,
              col = adjustcolor(MODEL_COLS[m], 0.6),
              border = MODEL_COLS[m], outline = FALSE, axes = FALSE)
    }
  }
  legend("topright", legend = MODELS, fill = adjustcolor(MODEL_COLS[MODELS], 0.6),
         border = MODEL_COLS[MODELS], bty = "n", cex = 0.85)
} else {
  plot.new(); title("KA not available")
}

# Panel 3: residuals vs fitted (all five models, same axes)
plot(NA,
     xlim = range(unlist(lapply(MODELS, function(m)
       pp[[m]]$residuals_df$log_hat_mean[
         is.finite(pp[[m]]$residuals_df$log_hat_mean)])), na.rm = TRUE),
     ylim = x_range,
     xlab = "log(SOC predicted)",
     ylab = "Residual (log)",
     main = "Residuals vs fitted")
abline(h = 0, lty = 2, col = "grey40")
for (m in MODELS) {
  rd <- pp[[m]]$residuals_df
  mask <- is.finite(rd$log_hat_mean) & is.finite(rd$residual_log)
  points(rd$log_hat_mean[mask], rd$residual_log[mask],
         pch = 16, cex = 0.5, col = adjustcolor(MODEL_COLS[m], 0.3))
  if (sum(mask) > 10) {
    lo  <- loess(rd$residual_log[mask] ~ rd$log_hat_mean[mask])
    ord <- order(rd$log_hat_mean[mask])
    lines(rd$log_hat_mean[mask][ord], predict(lo)[ord],
          col = MODEL_COLS[m], lwd = 2)
  }
}
legend("topright", legend = MODELS, col = MODEL_COLS[MODELS],
       lwd = 2, bty = "n", cex = 0.85)

# Panel 4: RMSE barplot across all five models
rmse_vals <- sapply(MODELS, function(m) pp[[m]]$metrics$RMSE_median)
barplot(rmse_vals,
        col    = MODEL_COLS[MODELS],
        names.arg = MODELS,
        ylab   = "RMSE — posterior median (tC/ha)",
        main   = "RMSE comparison",
        border = "white",
        ylim   = c(0, max(rmse_vals) * 1.2))
text(x   = seq(0.7, by = 1.2, length.out = length(MODELS)),
     y   = rmse_vals + max(rmse_vals) * 0.03,
     labels = sprintf("%.2f", rmse_vals), cex = 0.9)

mtext(sprintf("HIKET multi-model residuals  |  comp %s", COMP_ID),
      side = 3, outer = TRUE, line = -1.5, cex = 1.0, font = 2)
dev.off()
message(sprintf("Residual comparison: %s", resid_png))



# =============================================================================
# 6.  Plot D: RF importance heatmap (cross-model residual structure)
#
# Reads the per-model RF importance CSVs written by run_residual_analysis.R
# and produces a single heatmap where:
#   - rows  = union of the top-15 predictors across all models
#   - cols  = models
#   - colour = relative importance (within-model max = 1)
#
# WHY RELATIVE: raw permutation importance scales with each model's residual
# variance and fit quality, so raw values are not cross-model comparable.
# Dividing by the within-model maximum maps every model to [0, 1]: the colour
# answers "how consistently does this predictor structure the residuals across
# models?" rather than "how large is the number for this model?"
#
# Rows sorted by mean relative importance (ascending so the most important
# predictor appears at the top — heatmap() renders row 1 at the bottom).
#
# Requires (per model):
#   <DIR_DIAG>/<MODEL>/<MODEL>_rf_importance_<RUN_ID>.csv
#   <DIR_DIAG>/<MODEL>/<MODEL>_rf_summary_<RUN_ID>.rds   (for OOB R²)
#
# If any model's CSV is missing a warning is issued and that model's column
# will be absent from the heatmap.
# =============================================================================

TOP_N    <- 15L
imp_list <- setNames(vector("list", length(MODELS)), MODELS)
oob_r2   <- setNames(rep(NA_real_,  length(MODELS)), MODELS)

for (m in MODELS) {
  csv_path <- file.path(DIR_DIAG, m,
                        sprintf("%s_rf_importance_%s.csv", m, run_ids[[m]]))
  rds_path <- file.path(DIR_DIAG, m,
                        sprintf("%s_rf_summary_%s.rds",   m, run_ids[[m]]))
  if (!file.exists(csv_path)) {
    warning(sprintf("[heatmap] RF importance CSV not found for %s: %s", m, csv_path))
    next
  }
  imp_list[[m]] <- read.csv(csv_path, stringsAsFactors = FALSE)
  if (file.exists(rds_path))
    oob_r2[m] <- readRDS(rds_path)$oob_r2
}

ok_models <- MODELS[!sapply(imp_list, is.null)]

if (length(ok_models) == 0) {
  warning("[heatmap] No RF importance data found; skipping heatmap. ",
          "Run run_residual_analysis.R for each model first.")
} else {
  
  if (length(ok_models) < length(MODELS))
    message(sprintf("[heatmap] Missing RF data for: %s",
                    paste(setdiff(MODELS, ok_models), collapse = ", ")))
  message(sprintf("[heatmap] Building importance matrix for: %s",
                  paste(ok_models, collapse = ", ")))
  
  # Relative importance: divide each model's vector by its own maximum
  rel_list <- lapply(ok_models, function(m) {
    df     <- imp_list[[m]]
    df$rel <- df$importance / max(df$importance, na.rm = TRUE)
    df
  })
  names(rel_list) <- ok_models
  
  # Union of top-N predictors across all models
  top_preds <- unique(unlist(lapply(ok_models, function(m) {
    df <- rel_list[[m]][order(rel_list[[m]]$importance, decreasing = TRUE), ]
    head(df$variable, TOP_N)
  })))
  
  # Matrix: rows = predictors, cols = models (absent predictors get 0)
  imp_mat <- matrix(0, nrow = length(top_preds), ncol = length(ok_models),
                    dimnames = list(top_preds, ok_models))
  for (m in ok_models) {
    df  <- rel_list[[m]]
    idx <- match(top_preds, df$variable)
    has <- !is.na(idx)
    imp_mat[has, m] <- df$rel[idx[has]]
  }
  
  # Sort ascending: heatmap() renders row 1 at the bottom, so the most
  # important predictor (last row) ends up at the top of the image.
  imp_mat <- imp_mat[order(rowMeans(imp_mat), decreasing = FALSE), , drop = FALSE]
  
  # Column labels with OOB R² where available
  col_labels <- sapply(ok_models, function(m) {
    if (!is.na(oob_r2[m])) sprintf("%s (OOB R2=%.2f)", m, oob_r2[m]) else m
  })
  
  heat_pal <- colorRampPalette(
    c("white", "#deebf7", "#9ecae1", "#3182bd", "#08306b")
  )(100)
  
  n_pred   <- nrow(imp_mat)
  plot_h   <- max(7L, round(n_pred * 0.45 + 4L))
  heat_png <- file.path(DIR_OUT,
                        sprintf("multimodel_rf_importance_heatmap_%s.png", COMP_ID))
  
  png(heat_png,
      width  = 10L * 300L,
      height = plot_h * 300L,
      res    = 300L)
  
  heatmap(imp_mat,
          Rowv    = NA,
          Colv    = NA,
          scale   = "none",
          col     = heat_pal,
          labCol  = col_labels,
          margins = c(8, 14),
          main    = "RF residual predictor importance (relative, within-model)",
          cexRow  = 0.85,
          cexCol  = 0.80)
  
  dev.off()
  message(sprintf("[heatmap] Written: %s", heat_png))
}



# =============================================================================
# 6b.  Plot E: 60-year SOC projection by NFI region
# =============================================================================
# For each model: aggregate projection_predictions by NFI region (North/South),
# compute the cross-draw median and 95% CI of mean regional SOC, and plot as
# a continuous trajectory spanning both the historical calibration period and
# the 60-year forward projection.
#
# nfi_region is read from site_raw.csv (written by assign_nfi_regions.R).
# Northern Finland = Lappi, Pohjois-Pohjanmaa, Kainuu (NFI convention).
#
# If projection_predictions is absent from a bundle (old run without Section 4b),
# that model's projection panel is left blank with a warning.
#
# Y-axis is capped at Y_MAX_PROJ (tC/ha). Draws that exceed this produce CI
# polygons clipped to Y_MAX_PROJ -- this is intentional: the clipping itself
# visualises the structural instability of models whose posterior spans
# physically impossible SOC values in the projection window.
# =============================================================================

Y_MAX_PROJ <- 300L   # tC/ha display ceiling; Finnish forest SOC < 250 tC/ha

# --- Load region lookup -------------------------------------------------------
site_raw_proj <- tryCatch(
  read.csv("./Data/model_inputs/site_raw.csv", stringsAsFactors = FALSE),
  error = function(e) {
    warning("[projection] site_raw.csv not found -- skipping projection plot")
    NULL
  })

if (!is.null(site_raw_proj) && "nfi_region" %in% names(site_raw_proj)) {
  
  site_raw_proj$plot_id <- as.character(site_raw_proj$plot_id)
  region_lookup <- site_raw_proj[!is.na(site_raw_proj$nfi_region),
                                 c("plot_id", "nfi_region")]
  
  # --- Aggregate historical and projection per model -------------------------
  aggregate_regional <- function(pred_df) {
    if (is.null(pred_df) || nrow(pred_df) == 0) return(NULL)
    pred_df$plot_id <- as.character(pred_df$plot_id)
    pred_df <- merge(pred_df, region_lookup, by = "plot_id", all.x = FALSE)
    if (nrow(pred_df) == 0) return(NULL)
    
    pred_df %>%
      group_by(draw, nfi_region, year) %>%
      summarise(mean_soc = mean(total_soc, na.rm = TRUE), .groups = "drop") %>%
      group_by(nfi_region, year) %>%
      summarise(
        soc_median = median(mean_soc, na.rm = TRUE),
        soc_q025   = quantile(mean_soc, 0.025, na.rm = TRUE),
        soc_q975   = quantile(mean_soc, 0.975, na.rm = TRUE),
        .groups = "drop"
      )
  }
  
  hist_reg <- lapply(MODELS, function(m)
    aggregate_regional(pp[[m]]$posterior_predictions))
  names(hist_reg) <- MODELS
  
  proj_reg <- lapply(MODELS, function(m) {
    if (is.null(pp[[m]]$projection_predictions)) {
      message(sprintf("[projection] No projection_predictions in %s bundle", m))
      return(NULL)
    }
    aggregate_regional(pp[[m]]$projection_predictions)
  })
  names(proj_reg) <- MODELS
  
  # --- Shared x range (historical + projection) ------------------------------
  all_years <- sort(unique(c(
    unlist(lapply(hist_reg, function(d) if (!is.null(d)) d$year)),
    unlist(lapply(proj_reg, function(d) if (!is.null(d)) d$year))
  )))
  x_range <- range(all_years, na.rm = TRUE)
  
  # --- Plot ------------------------------------------------------------------
  REGION_COLS <- c(North = "#2166ac", South = "#d6604d")
  REGION_LTY  <- c(North = 1,        South  = 2)
  
  proj_png <- file.path(DIR_OUT,
                        sprintf("multimodel_projection_%s.png", COMP_ID))
  png(proj_png,
      width  = 20L * PX_PER_IN,
      height =  7L * PX_PER_IN,
      res    = PX_PER_IN)
  par(mfrow = c(1, 5),
      mar   = c(4, 4, 3.5, 0.5),
      oma   = c(0, 0, 2.5, 0))
  
  for (m in MODELS) {
    hd <- hist_reg[[m]]
    pd <- proj_reg[[m]]
    
    plot(NA,
         xlim = x_range,
         ylim = c(0, Y_MAX_PROJ),
         xlab = "Year",
         ylab = if (m == MODELS[1]) "Mean SOC (tC/ha)" else "",
         main = m,
         yaxt = if (m == MODELS[1]) "s" else "n")
    
    # Vertical reference at last observation year
    abline(v = 2006, lty = 3, col = "grey55", lwd = 1.2)
    
    for (reg in c("North", "South")) {
      col <- REGION_COLS[reg]
      lty <- REGION_LTY[reg]
      
      # Historical ribbon + median
      if (!is.null(hd)) {
        hd_r <- hd[hd$nfi_region == reg, ]
        hd_r <- hd_r[order(hd_r$year), ]
        if (nrow(hd_r) > 0) {
          polygon(c(hd_r$year, rev(hd_r$year)),
                  pmin(c(hd_r$soc_q025, rev(hd_r$soc_q975)), Y_MAX_PROJ),
                  col = adjustcolor(col, 0.12), border = NA)
          lines(hd_r$year, pmin(hd_r$soc_median, Y_MAX_PROJ),
                col = col, lwd = 1.5, lty = lty)
        }
      }
      
      # Projection ribbon + median (slightly darker alpha to distinguish period)
      if (!is.null(pd)) {
        pd_r <- pd[pd$nfi_region == reg, ]
        pd_r <- pd_r[order(pd_r$year), ]
        if (nrow(pd_r) > 0) {
          polygon(c(pd_r$year, rev(pd_r$year)),
                  pmin(c(pd_r$soc_q025, rev(pd_r$soc_q975)), Y_MAX_PROJ),
                  col = adjustcolor(col, 0.20), border = NA)
          lines(pd_r$year, pmin(pd_r$soc_median, Y_MAX_PROJ),
                col = col, lwd = 2.2, lty = lty)
        }
      }
    }
    
    # Note if 95% CI exceeds display ceiling (scientifically meaningful)
    if (!is.null(pd)) {
      max_ci <- max(pd$soc_q975, na.rm = TRUE)
      if (max_ci > Y_MAX_PROJ)
        mtext(sprintf("95%% CI max: %.0f tC/ha (clipped)", max_ci),
              side = 1, line = -1.5, cex = 0.62, col = "grey35", adj = 0.05)
    }
    
    # Legend on first panel only
    if (m == MODELS[1])
      legend("topleft",
             legend  = c("North — median (95% CI)",
                         "South — median (95% CI)",
                         "Last obs. year (2006)"),
             col     = c(REGION_COLS, "grey55"),
             lwd     = c(2.2, 2.2, 1.2),
             lty     = c(REGION_LTY, 3),
             bty     = "n", cex = 0.78)
  }
  
  mtext(
    sprintf("HIKET — 60-year SOC projection by NFI region  |  y-axis capped at %d tC/ha  |  comp %s",
            Y_MAX_PROJ, COMP_ID),
    side = 3, outer = TRUE, line = 0.8, cex = 0.92, font = 2)
  dev.off()
  message(sprintf("Projection plot: %s", proj_png))
  
} else {
  message("[projection] nfi_region not found in site_raw.csv -- ",
          "run assign_nfi_regions.R first; skipping projection plot.")
}

# =============================================================================
# 7.  Save summary bundle + print metrics
# =============================================================================

summary_out <- list(
  run_ids    = run_ids,
  comp_id    = COMP_ID,
  metrics    = metrics_df,
  timestamp  = Sys.time()
)
summary_rds <- file.path(DIR_OUT,
                         sprintf("multimodel_summary_%s.rds", COMP_ID))
saveRDS(summary_out, summary_rds)
message(sprintf("\nSummary saved: %s", summary_rds))

message("\n=============================================================")
message("  Multi-Model Comparison Complete")
message(sprintf("  COMP_ID: %s", COMP_ID))
message("=============================================================")
message("\nMetrics summary:")
print(metrics_df[, c("Model","R2","RMSE_median","Bias_median","Coverage_95","N_obs")],
      row.names = FALSE)