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
#     multimodel_metrics_<COMP_ID>.txt          -- formatted metrics table
#     multimodel_obs_vs_pred_<COMP_ID>.png      -- 5-panel scatter
#     multimodel_delta_soc_<COMP_ID>.png        -- annual ΔSOC rate trajectories
#     multimodel_residuals_<COMP_ID>.png        -- residual distributions + by KA
#     multimodel_summary_<COMP_ID>.rds          -- metrics list for downstream use
#
#   where COMP_ID = <RUN_ID_SP1>_<RUN_ID_TP2>_<RUN_ID_07>_<RUN_ID_15>_<RUN_ID_20>
# =============================================================================

library(dplyr)

MODELS <- c("SP1", "TP2", "Yasso07", "Yasso15", "Yasso20")
DIR_RUNS <- "./Calibration_real_data/runs/"
DIR_OUT  <- "./Calibration_real_data/diagnostics/multimodel"
dir.create(DIR_OUT, showWarnings = FALSE, recursive = TRUE)

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

# Quick structural check: each bundle must carry the three expected components
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
# 3.  Plot A: Obs vs predicted — 5 panels side by side
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
    width  = 28L * PX_PER_IN,   # 5 panels at ~5.6 in each
    height =  7L * PX_PER_IN,
    res    = PX_PER_IN)
par(mfrow = c(1, 5), mar = c(4.5, 4.5, 3.5, 1))

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
         legend = c(sprintf("R² = %.3f",   mt$R2),
                    sprintf("RMSE = %.1f",  mt$RMSE_median),
                    sprintf("Bias = %+.1f", mt$bias_median),
                    sprintf("Cov = %.2f",   mt$coverage_95)),
         bty = "n", cex = 0.8)
}

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
# Observations are the same across all model bundles; use the first model's
# posterior_summary to avoid any accidental differences from bundle construction.
obs_ref   <- pp[[MODELS[1]]]$posterior_summary
obs_delta <- obs_ref %>%
  filter(!is.na(soc_obs_tCha)) %>%
  group_by(plot_id) %>%
  arrange(year, .by_group = TRUE) %>%
  summarise(
    delta_obs  = (last(soc_obs_tCha) - first(soc_obs_tCha)) /
      (last(year) - first(year)),
    year_first = first(year),
    year_last  = last(year),
    .groups = "drop"
  ) %>%
  filter(year_first != year_last)

obs_mean_delta <- mean(obs_delta$delta_obs, na.rm = TRUE)
obs_se_delta   <- sd(obs_delta$delta_obs,   na.rm = TRUE) / sqrt(nrow(obs_delta))

# Shared y limits: union of all model CI bands + observed CI + zero
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

# 95% CI ribbons (drawn first so lines sit on top)
for (m in MODELS) {
  tr  <- traj_list[[m]]
  polygon(c(tr$year, rev(tr$year)),
          c(tr$q025, rev(tr$q975)),
          col = adjustcolor(MODEL_COLS[m], 0.12), border = NA)
}
# Median lines
for (m in MODELS) {
  tr <- traj_list[[m]]
  lines(tr$year, tr$median_delta, col = MODEL_COLS[m], lwd = 2.5)
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
                  sprintf("Observed mean annual ΔSOC ± 95%% CI  (n=%d)",
                          nrow(obs_delta))),
       col = c(MODEL_COLS[MODELS], "grey20"),
       lwd = c(rep(2.5, length(MODELS)), 2),
       lty = c(rep(1,   length(MODELS)), 2),
       bty = "n", cex = 0.9)
dev.off()
message(sprintf("ΔSOC trajectory: %s", delta_png))


# =============================================================================
# 5.  Plot C: Residual distributions + by KA  (2 rows × 2 cols)
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
x_range <- c(max(x_range[1], -3), min(x_range[2], 3))   # clip extreme tails

y_max     <- 0
dens_list <- lapply(MODELS, function(m) {
  d <- density(finite_resids[[m]], from = x_range[1], to = x_range[2])
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

# Panel 2: residual bias by KA class (grouped boxplot, one cluster per class)
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
  # Slightly narrower boxes to keep five models legible within each KA cluster
  width   <- 0.14
  offsets <- seq(-(n_mod - 1) / 2, (n_mod - 1) / 2,
                 length.out = n_mod) * width * 1.6
  
  x_centres <- seq_len(n_ka)
  
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
      ka   <- ka_levels[i]
      vals <- resid_combined$residual_log[
        resid_combined$model == m & resid_combined$KA == ka]
      if (length(vals) < 3) next
      xc <- x_centres[i] + offsets[j]
      boxplot(vals, at = xc, add = TRUE, boxwex = width,
              col    = adjustcolor(MODEL_COLS[m], 0.6),
              border = MODEL_COLS[m],
              outline = FALSE, axes = FALSE)
    }
  }
  legend("topright", legend = MODELS,
         fill   = adjustcolor(MODEL_COLS[MODELS], 0.6),
         border = MODEL_COLS[MODELS],
         bty = "n", cex = 0.85)
} else {
  plot.new(); title("KA not available")
}

# Panel 3: residuals vs fitted (all five models overlaid, same axes)
log_hat_range <- range(unlist(lapply(MODELS, function(m)
  pp[[m]]$residuals_df$log_hat_mean[
    is.finite(pp[[m]]$residuals_df$log_hat_mean)])), na.rm = TRUE)

plot(NA,
     xlim = log_hat_range,
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

# Panel 4: RMSE bar chart (one bar per model, sorted by complexity)
rmse_vals <- sapply(MODELS, function(m) pp[[m]]$metrics$RMSE_median)
bp <- barplot(rmse_vals,
              col       = MODEL_COLS[MODELS],
              names.arg = MODELS,
              ylab      = "RMSE — posterior median (tC/ha)",
              main      = "RMSE comparison",
              border    = "white",
              ylim      = c(0, max(rmse_vals) * 1.2))
text(x      = bp,
     y      = rmse_vals + max(rmse_vals) * 0.03,
     labels = sprintf("%.2f", rmse_vals),
     cex    = 0.9)

mtext(sprintf("HIKET multi-model residuals  |  comp %s", COMP_ID),
      side = 3, outer = TRUE, line = -1.5, cex = 1.0, font = 2)
dev.off()
message(sprintf("Residual comparison: %s", resid_png))


# =============================================================================
# 6.  Save summary bundle + print metrics
# =============================================================================

summary_out <- list(
  run_ids   = run_ids,
  comp_id   = COMP_ID,
  metrics   = metrics_df,
  timestamp = Sys.time()
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
print(metrics_df[, c("Model", "R2", "RMSE_median", "Bias_median",
                     "Coverage_95", "N_obs")],
      row.names = FALSE)