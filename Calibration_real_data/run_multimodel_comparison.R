# =============================================================================
# run_multimodel_comparison.R
#
# Multi-model comparison stage of the HIKET pipeline.
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
#   ./Calibration_real_data/runs/<MODEL>_posterior_predictive_<RUN_ID>.rds
#
# OUTPUTS:
#   ./Calibration_real_data/diagnostics/multimodel/
#     multimodel_metrics_<COMP_ID>.txt               -- formatted metrics table
#     multimodel_obs_vs_pred_<COMP_ID>.png           -- 2x3 scatter + bias panel
#     multimodel_delta_soc_<COMP_ID>.png             -- ΔSOC trajectories overlaid
#     multimodel_residuals_<COMP_ID>.png             -- residual distributions + by KA
#     multimodel_rf_importance_heatmap_<COMP_ID>.png -- cross-model RF importance
#     multimodel_summary_<COMP_ID>.rds               -- metrics list for downstream use
#
#   where COMP_ID = <RUN_ID_SP1>_<RUN_ID_TP2>_<RUN_ID_07>_<RUN_ID_15>_<RUN_ID_20>
# =============================================================================

library(dplyr)

MODELS <- c("SP1", "TP2", "Yasso07", "Yasso15", "Yasso20")
DIR_RUNS     <- "./Calibration_real_data/runs/"
DIR_OUT      <- "./Calibration_real_data/diagnostics/multimodel"
DIR_DIAG     <- "./Calibration_real_data/diagnostics"
dir.create(DIR_OUT, showWarnings = FALSE, recursive = TRUE)

PX_PER_IN <- 150L

# Model colours used consistently across all comparison plots.
# Ordered from simplest (SP1, 1 pool) to most complex (Yasso20, 5 pools)
# so that hue progression mirrors structural complexity throughout.
MODEL_COLS <- c(
  SP1     = "#762a83",   # purple  -- simplest (1 pool)
  TP2     = "#e08214",   # amber   -- intermediate (2 pools)
  Yasso07 = "#2166ac",   # blue
  Yasso15 = "#d6604d",   # red-orange
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

args    <- commandArgs(trailingOnly = TRUE)
run_ids <- setNames(vector("list", length(MODELS)), MODELS)

if (length(args) >= 1) {
  for (a in args) {
    parts <- strsplit(a, ":", fixed = TRUE)[[1]]
    if (length(parts) == 2 && parts[1] %in% MODELS)
      run_ids[[parts[1]]] <- parts[2]
  }
}

# Fall back to auto-detect for any model not supplied on the command line
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
            "Model", "R\u00b2", "RMSE (med)", "Bias (med)", "95% cov"))
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
# 3.  Plot A: Obs vs predicted — 2x3 grid (5 scatter + bias barplot)
# =============================================================================

obs_pred_png <- file.path(DIR_OUT,
                          sprintf("multimodel_obs_vs_pred_%s.png", COMP_ID))

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
  
  ka_col <- c("1" = "#2166ac", "2" = "#74add1",
              "3" = "#f4a582", "4" = "#d6604d")[as.character(rd$KA)]
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
         legend = c(sprintf("R\u00b2 = %.3f",   mt$R2),
                    sprintf("RMSE = %.1f",        mt$RMSE_median),
                    sprintf("Bias = %+.1f",       mt$bias_median),
                    sprintf("Cov = %.2f",         mt$coverage_95)),
         bty = "n", cex = 0.8)
}

# Panel 6: bias barplot across all five models.
# Bias_median = median(log(obs) - log(pred)) across all observations for that
# model. Positive = systematic under-prediction; negative = over-prediction.
# Sorted by model complexity (SP1 -> Yasso20) to reveal complexity-vs-bias trends.
bias_vals <- sapply(MODELS, function(m) pp[[m]]$metrics$bias_median)
bar_cols  <- MODEL_COLS[MODELS]
ylim_bias <- c(min(bias_vals, 0), max(bias_vals, 0)) *
  c(ifelse(min(bias_vals) < 0, 1.25, 1), ifelse(max(bias_vals) > 0, 1.25, 1))

bp <- barplot(bias_vals,
              col       = bar_cols,
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

# Build trajectory summary per model.
# Annual rate = year-on-year difference in mean SOC across plots (tC/ha/yr),
# so model and observed are on the same axis. The first year of each plot
# produces NA (no lag) and is dropped by na.rm = TRUE in summarise.
traj_list <- lapply(MODELS, function(m) {
  preds <- pp[[m]]$posterior_predictions
  preds %>%
    group_by(plot_id, draw) %>%
    arrange(year, .by_group = TRUE) %>%
    mutate(delta_soc = total_soc - lag(total_soc)) %>%
    ungroup() %>%
    group_by(draw, year) %>%
    summarise(mean_delta = mean(delta_soc, na.rm = TRUE), .groups = "drop") %>%
    group_by(year) %>%
    summarise(
      median_delta = median(mean_delta, na.rm = TRUE),
      q025         = quantile(mean_delta, 0.025, na.rm = TRUE),
      q975         = quantile(mean_delta, 0.975, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(year) %>%
    mutate(model = m)
})
names(traj_list) <- MODELS

# Observed mean annual rate: (SOC_last - SOC_first) / (year_last - year_first)
# per plot, then averaged across plots. Same units as model trajectories.
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
  ) %>%
  filter(year_first != year_last)

obs_mean_delta <- mean(obs_delta$delta_obs, na.rm = TRUE)
obs_se_delta   <- sd(obs_delta$delta_obs, na.rm = TRUE) / sqrt(nrow(obs_delta))

all_vals  <- c(unlist(lapply(traj_list, function(t) c(t$q025, t$q975))),
               obs_mean_delta + c(-1.96, 1.96) * obs_se_delta, 0)
ylim_r    <- range(all_vals, na.rm = TRUE)
ylim_r    <- ylim_r + c(-0.05, 0.05) * diff(ylim_r)
all_years <- sort(unique(unlist(lapply(traj_list, `[[`, "year"))))

delta_png <- file.path(DIR_OUT, sprintf("multimodel_delta_soc_%s.png", COMP_ID))
png(delta_png, width = 13L * PX_PER_IN, height = 7L * PX_PER_IN, res = PX_PER_IN)
par(mar = c(4, 4.5, 3.5, 1))

plot(NA, xlim = range(all_years), ylim = ylim_r,
     xlab = "Year",
     ylab = expression(paste(Delta, "SOC annual rate (tC/ha/yr)")),
     main = "HIKET — Mean ΔSOC across models")
abline(h = 0, lty = 3, col = "grey50")

for (m in MODELS) {
  tr  <- traj_list[[m]]
  col <- MODEL_COLS[m]
  polygon(c(tr$year, rev(tr$year)),
          c(tr$q025, rev(tr$q975)),
          col = adjustcolor(col, 0.12), border = NA)
}
for (m in MODELS) {
  tr  <- traj_list[[m]]
  lines(tr$year, tr$median_delta, col = MODEL_COLS[m], lwd = 2.5)
}

abline(h = obs_mean_delta, col = "grey20", lwd = 2, lty = 2)
rect(xleft   = min(obs_delta$year_first),
     xright  = max(obs_delta$year_last),
     ybottom = obs_mean_delta - 1.96 * obs_se_delta,
     ytop    = obs_mean_delta + 1.96 * obs_se_delta,
     col = adjustcolor("grey20", 0.10), border = NA)

legend("topleft",
       legend = c(MODELS,
                  sprintf("Observed mean \u0394SOC \u00b1 95%% CI  (n=%d)", nrow(obs_delta))),
       col    = c(MODEL_COLS[MODELS], "grey20"),
       lwd    = c(rep(2.5, length(MODELS)), 2),
       lty    = c(rep(1,   length(MODELS)), 2),
       bty    = "n", cex = 0.9)
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
x_range <- c(max(x_range[1], -3), min(x_range[2], 3))

y_max     <- 0
dens_list <- lapply(MODELS, function(m) {
  d     <- density(finite_resids[[m]], from = x_range[1], to = x_range[2])
  y_max <<- max(y_max, max(d$y))
  d
})

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
  width <- 0.2
  
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
      vals <- resid_combined$residual_log[
        resid_combined$model == m & resid_combined$KA == ka_levels[i] &
          is.finite(resid_combined$residual_log)]
      if (length(vals) < 3) next
      boxplot(vals, at = x_centres[i] + offsets[j], add = TRUE,
              boxwex = width,
              col    = adjustcolor(MODEL_COLS[m], 0.6),
              border = MODEL_COLS[m], outline = FALSE, axes = FALSE)
    }
  }
  legend("topright", legend = MODELS,
         fill   = adjustcolor(MODEL_COLS[MODELS], 0.6),
         border = MODEL_COLS[MODELS], bty = "n", cex = 0.85)
} else {
  plot.new(); title("KA not available")
}

# Panel 3: residuals vs fitted (all models, same axes)
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
  rd   <- pp[[m]]$residuals_df
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

# Panel 4: RMSE bar chart
rmse_vals <- sapply(MODELS, function(m) pp[[m]]$metrics$RMSE_median)
barplot(rmse_vals,
        col       = MODEL_COLS[MODELS],
        names.arg = MODELS,
        ylab      = "RMSE — posterior median (tC/ha)",
        main      = "RMSE comparison",
        border    = "white",
        ylim      = c(0, max(rmse_vals) * 1.2))
text(x      = seq(0.7, by = 1.2, length.out = length(MODELS)),
     y      = rmse_vals + max(rmse_vals) * 0.03,
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