source("manuscript/figures/run_ids.R")   # auto-selects current RUN_IDs
setwd("/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling")
# F6 (redesign) -- observed vs predicted SOC, HOLDOUT plots only (calibration -> supp),
# points classed by stand BASAL AREA (low/med/high). Puts the #1 residual driver
# (F11) directly in the obs-vs-pred cloud: if skill tracks productivity, the classes
# separate. 2x3, one panel per model; holdout R^2 from each bundle's metrics_holdout.
# Data: <MODEL> predictive bundle -> residuals_df (soc_obs_tCha, soc_median,
# basal_area_85, is_holdout) + metrics_holdout$R2.

runs <- "Calibration_real_data_transient/runs"
bundles <- c(
  SP1     = sprintf("SP1_posterior_predictive_%s.rds", RID[["SP1"]]),
  TP2     = sprintf("TP2_posterior_predictive_%s.rds", RID[["TP2"]]),
  TP3     = sprintf("TP3_posterior_predictive_%s.rds", RID[["TP3"]]),
  Yasso07 = sprintf("Yasso07_posterior_predictive_%s.rds", RID[["Yasso07"]]),
  Yasso15 = sprintf("Yasso15_posterior_predictive_%s.rds", RID[["Yasso15"]]),
  Yasso20 = sprintf("Yasso20_posterior_predictive_%s.rds", RID[["Yasso20"]]))

source("manuscript/figures/model_palette.R")   # shared basal-area quintile palette
cls_col <- BASAL_COL

dat <- list(); r2h <- c()
for (m in names(bundles)) {
  b  <- readRDS(file.path(runs, bundles[[m]]))
  rd <- b$residuals_df
  rd <- rd[rd$is_holdout %in% TRUE, ]
  rd <- rd[is.finite(rd$soc_obs_tCha) & is.finite(rd$soc_median) &
             is.finite(rd$basal_area_85), ]
  dat[[m]] <- rd[, c("soc_obs_tCha", "soc_median", "basal_area_85")]
  r2h[m]   <- b$metrics_holdout$R2
  rm(b); gc()
}

# common basal-area QUINTILE classes across all models (shared helper -> matches S2)
allba <- unlist(lapply(dat, function(d) d$basal_area_85))
bc    <- basal_classes(allba); classify_ba <- bc$classify

lim <- range(unlist(lapply(dat, function(d) c(d$soc_obs_tCha, d$soc_median))), na.rm = TRUE)

# Print size (0.95 x 6.53 in text width, 9 pt), so nothing prints below ~7 pt.
png("manuscript/figures/F6_obs_vs_pred_holdout.png", width = 6.2, height = 4.9, units = "in", res = 300, pointsize = 9)
# six panels plus a legend strip underneath (inside a panel it covered the points)
layout(rbind(matrix(1:6, 2, byrow = TRUE), 7), heights = c(1, 1, 0.3))
par(cex = 1)   # layout/mfrow shrink text to 0.66 (or 0.83 in 2x2); keep the 9 pt
par(mar = c(3.0, 3.4, 1.6, 0.4), mgp = c(2.0, 0.5, 0), las = 1, tcl = -0.3)
for (m in names(bundles)) {
  k  <- match(m, names(bundles))
  d  <- dat[[m]]; cl <- classify_ba(d$basal_area_85)
  plot(d$soc_obs_tCha, d$soc_median, xlim = lim, ylim = lim, pch = 19, cex = 0.6,
       col = adjustcolor(cls_col[as.character(cl)], 0.75),
       xlab = "Observed (tC/ha)", ylab = if (k %in% c(1, 4)) "Predicted (tC/ha)" else "")
  mtext(sprintf("(%s) %s, R\u00b2 %.3f", letters[k], m, r2h[m]), 3, 0.3, adj = 0, font = 2)
  # See build_S2: dashed = 1:1, dotted = predicting the mean (no skill), red = OLS fit.
  abline(0, 1, col = "grey40", lwd = 1.2, lty = 2)
  abline(h = mean(d$soc_median), col = "steelblue4", lwd = 1.2, lty = 3)
  abline(lm(soc_median ~ soc_obs_tCha, d), col = "firebrick", lwd = 1.4)
}
par(mar = c(0, 0, 0, 0)); plot.new()
legend("left", bty = "n", cex = 0.85, pt.cex = 1.1, ncol = 3, pch = 19, col = cls_col[seq_len(bc$nb - 1)],
       title = "Basal area 1985 (m2/ha)", legend = bc$labels)
legend("right", bty = "n", cex = 0.85, lwd = c(1.2, 1.2, 1.4), lty = c(2, 3, 1),
       col = c("grey40", "steelblue4", "firebrick"),
       legend = c("1:1", "mean prediction", "regression"))
dev.off()
cat("wrote manuscript/figures/F6_obs_vs_pred_holdout.png\n")
cat("holdout R2:", paste(sprintf("%s %.3f", names(r2h), r2h), collapse = "  "), "\n")
