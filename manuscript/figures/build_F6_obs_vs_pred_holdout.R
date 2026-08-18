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

png("manuscript/figures/F6_obs_vs_pred_holdout.png", width = 10.5, height = 7, units = "in", res = 200)
par(mfrow = c(2, 3), mar = c(4.2, 4.2, 2.8, 1), mgp = c(2.4, 0.7, 0), las = 1,
    cex.axis = 1.05, cex.lab = 1.25)
for (m in names(bundles)) {
  d  <- dat[[m]]; cl <- classify_ba(d$basal_area_85)
  plot(d$soc_obs_tCha, d$soc_median, xlim = lim, ylim = lim, pch = 19, cex = 0.8,
       col = adjustcolor(cls_col[as.character(cl)], 0.75),
       xlab = "Observed SOC (tC/ha)", ylab = "Predicted SOC (tC/ha)",
       main = sprintf("%s   (holdout R2 = %.3f)", m, r2h[m]), font.main = 1, cex.main = 1.2)
  # See build_S2: dashed = 1:1, dotted = predicting the mean (no skill), red = OLS fit.
  abline(0, 1, col = "grey40", lwd = 1.4, lty = 2)
  abline(h = mean(d$soc_median), col = "steelblue4", lwd = 1.4, lty = 3)
  abline(lm(soc_median ~ soc_obs_tCha, d), col = "firebrick", lwd = 1.6)
}
# shared legend in the last cell margin
legend("bottomright", inset = c(0.02, 0.02), bty = "n", cex = 0.92, pt.cex = 1.3,
       pch = 19, col = cls_col[seq_len(bc$nb - 1)],
       title = "Stand basal area (quintiles, m2/ha)", legend = bc$labels)
legend("topleft", inset = c(0.02, 0.02), bty = "n", cex = 0.88,
       lwd = c(1.4, 1.4, 1.6), lty = c(2, 3, 1),
       col = c("grey40", "steelblue4", "firebrick"),
       legend = c("1:1 (perfect)", "mean prediction (no skill)", "OLS fit"))
dev.off()
cat("wrote manuscript/figures/F6_obs_vs_pred_holdout.png\n")
cat("holdout R2:", paste(sprintf("%s %.3f", names(r2h), r2h), collapse = "  "), "\n")
