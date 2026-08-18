source("manuscript/figures/run_ids.R")   # auto-selects current RUN_IDs
setwd("/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling")
# =============================================================================
# S14 -- POSTERIOR DISTRIBUTION OF RMSE, one panel per model.
# (Supplementary numbering; the surrounding appendix module is
#  manuscript/appendices/appendix_rmse_posterior.tex.)
#
# THE POINT: a cross-run comparison instrument. The multimodel metrics table
# reports ONE RMSE per model (computed from the posterior-mean prediction), which
# says nothing about how tightly the fit is pinned. Here RMSE is recomputed FOR
# EVERY POSTERIOR DRAW, so each panel is a distribution: its LOCATION is the fit,
# its WIDTH is how much the posterior disagrees with itself about the fit.
#
# WHY IT IS WORTH ITS PAGE. When priors are tightened (2026-08-17: sigma_init /
# sigma_input), the expectation is that fit DEGRADES -- that is the intended cost
# of a physically-bounded prior. This figure is how that cost is read off, and the
# distribution matters as much as the point: a prior that moves the location while
# leaving the width alone is doing something different from one that does both.
#
# CALIBRATION vs HOLDOUT are drawn in the same panel, because the gap between them
# is the part that is not just a fit statistic.
#
# ⚠ COMPARABILITY ACROSS RUNS. The x-axis is FIXED (RMSE_XLIM below), not fitted to
# the data, so two runs' figures can be laid side by side. If a future run falls
# outside it, widen it ONCE and rebuild BOTH, rather than letting each run pick its
# own range. The per-draw quantiles are also written to a CSV next to the PNG so
# the comparison can be numerical, not visual.
#
# Run from repo root:  Rscript manuscript/figures/build_S14_rmse_posterior.R
# =============================================================================

source("manuscript/figures/model_palette.R")

RMSE_XLIM <- c(40, 60)          # FIXED -- see comparability note above
OUT_PNG   <- "manuscript/figures/S14_rmse_posterior.png"
OUT_CSV   <- "manuscript/figures/S14_rmse_posterior.csv"

rid <- as.list(RID)
res <- list(); tab <- list()

for (m in FIG_MODELS) {
  b  <- readRDS(sprintf("Calibration_real_data_transient/runs/%s_posterior_predictive_%s.rds",
                        m, rid[[m]]))
  obs <- as.data.frame(b$residuals_df)[, c("plot_id", "year", "soc_obs_tCha")]
  pp  <- b$posterior_predictions[, c("plot_id", "year", "draw", "total_soc")]
  # keep only prediction rows that land on an observation
  j <- merge(pp, obs, by = c("plot_id", "year"))
  ho <- b$holdout_plots
  j$grp <- ifelse(j$plot_id %in% ho, "holdout", "calibration")
  j$sq  <- (j$total_soc - j$soc_obs_tCha)^2
  # RMSE per draw, per group
  agg <- aggregate(sq ~ draw + grp, j, function(z) sqrt(mean(z)))
  res[[m]] <- agg
  for (g in unique(agg$grp)) {
    v <- agg$sq[agg$grp == g]
    tab[[length(tab) + 1]] <- data.frame(
      model = m, run_id = rid[[m]], set = g, n_draws = length(v),
      q05 = quantile(v, .05), median = median(v), q95 = quantile(v, .95),
      n_obs = sum(j$grp == g) / length(unique(j$draw)))
  }
  message(sprintf("  %-8s calib %.2f | holdout %.2f",
                  m, median(agg$sq[agg$grp == "calibration"]),
                  median(agg$sq[agg$grp == "holdout"])))
}

tabd <- do.call(rbind, tab); rownames(tabd) <- NULL
write.csv(tabd, OUT_CSV, row.names = FALSE)

png(OUT_PNG, width = 10.5, height = 6.4, units = "in", res = 200)
par(mfrow = c(2, 3), mar = c(4.0, 4.2, 2.6, 1.0), mgp = c(2.4, 0.7, 0), las = 1)

GCOL <- c(calibration = "#2F6F9F", holdout = "#C26B51")
for (m in FIG_MODELS) {
  agg <- res[[m]]
  plot(NA, xlim = RMSE_XLIM, ylim = c(0, 1.08), xlab = "RMSE (tC/ha)",
       ylab = "posterior density (scaled)", main = m, font.main = 1)
  for (g in names(GCOL)) {
    v <- agg$sq[agg$grp == g]; if (!length(v)) next
    d <- density(v); d$y <- d$y / max(d$y)
    polygon(c(d$x, rev(d$x)), c(d$y, rep(0, length(d$y))),
            col = adjustcolor(GCOL[g], 0.30), border = NA)
    lines(d$x, d$y, col = GCOL[g], lwd = 2)
    segments(median(v), 0, median(v), 1.0, col = GCOL[g], lwd = 1.6, lty = 2)
  }
  if (m == FIG_MODELS[1])
    legend("topright", bty = "n", cex = 0.85, fill = adjustcolor(GCOL, 0.30),
           border = NA, legend = names(GCOL))
  mtext(sprintf("run %s", rid[[m]]), side = 3, line = -0.1, cex = 0.55, col = "grey45")
}
dev.off()
cat("wrote", OUT_PNG, "and", OUT_CSV, "\n")
print(tabd, row.names = FALSE, digits = 4)
