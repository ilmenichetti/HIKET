source("manuscript/figures/run_ids.R")   # auto-selects current RUN_IDs
setwd("/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling")
# S2 (supplement) -- calibration obs-vs-pred, the in-sample companion to the main-text
# holdout F6. COLOUR CODING (per user, 2026-07-15):
#   * the six scatterplots -> coloured by stand basal-area quintile, SAME as F6 (Green-Gold);
#   * the systematic-bias barplot -> the shared per-model palette (Temperature Diverging).
# Data: <MODEL> predictive bundle -> residuals_df (calibration rows) + metrics_calib.
source("manuscript/figures/model_palette.R")

runs <- "Calibration_real_data_transient/runs"
bundles <- c(
  SP1     = sprintf("SP1_posterior_predictive_%s.rds", RID[["SP1"]]),
  TP2     = sprintf("TP2_posterior_predictive_%s.rds", RID[["TP2"]]),
  TP3     = sprintf("TP3_posterior_predictive_%s.rds", RID[["TP3"]]),
  Yasso07 = sprintf("Yasso07_posterior_predictive_%s.rds", RID[["Yasso07"]]),
  Yasso15 = sprintf("Yasso15_posterior_predictive_%s.rds", RID[["Yasso15"]]),
  Yasso20 = sprintf("Yasso20_posterior_predictive_%s.rds", RID[["Yasso20"]]))

dat <- list(); r2c <- c(); bias <- c()
for (m in names(bundles)) {
  b  <- readRDS(file.path(runs, bundles[[m]]))
  rd <- b$residuals_df
  rd <- rd[rd$is_holdout %in% FALSE, ]                       # CALIBRATION plots
  rd <- rd[is.finite(rd$soc_obs_tCha) & is.finite(rd$soc_median) &
             is.finite(rd$basal_area_85), ]
  dat[[m]] <- rd[, c("soc_obs_tCha", "soc_median", "basal_area_85")]
  r2c[m]   <- b$metrics_calib$R2
  bias[m]  <- b$metrics_calib$bias_mean
  rm(b); gc()
}

allba <- unlist(lapply(dat, function(d) d$basal_area_85))
bc    <- basal_classes(allba)                                # shared quintile scheme
lim   <- range(unlist(lapply(dat, function(d) c(d$soc_obs_tCha, d$soc_median))), na.rm = TRUE)

png("manuscript/figures/S2_obs_vs_pred_calib.png", width = 12.5, height = 6.6, units = "in", res = 190)
layout(matrix(1:8, nrow = 2, byrow = TRUE))
par(mar = c(4, 4, 2.4, 1), mgp = c(2.4, 0.7, 0), las = 1)

# --- six calibration scatterplots (basal-area quintile colours, like F6) ---
for (m in names(bundles)) {
  d  <- dat[[m]]; cl <- bc$classify(d$basal_area_85)
  plot(d$soc_obs_tCha, d$soc_median, xlim = lim, ylim = lim, pch = 19, cex = 0.7,
       col = adjustcolor(BASAL_COL[as.character(cl)], 0.75),
       xlab = "Observed SOC (tC/ha)", ylab = "Predicted SOC (tC/ha)",
       main = sprintf("%s   (calib R2 = %.3f)", m, r2c[m]), font.main = 1, cex.main = 1.0)
  # Three references. Dashed grey = 1:1 (perfect). Dotted blue = predicting the mean
  # regardless of the plot (NO skill). Red = the actual OLS fit -- it lies ~80-85% of
  # the way from 1:1 down to flat, which is what R2 ~ 0.02 looks like. NB the cloud's
  # major axis IS near 1:1 (SMA slope ~1.15); OLS is flatter by the factor r ~ 0.14,
  # so "centred on 1:1" and "regresses flat" are both true and not in conflict.
  abline(0, 1, col = "grey40", lwd = 1.3, lty = 2)
  abline(h = mean(d$soc_median), col = "steelblue4", lwd = 1.3, lty = 3)
  abline(lm(soc_median ~ soc_obs_tCha, d), col = "firebrick", lwd = 1.5)
}

# --- systematic-bias barplot (per-model palette) ---
par(mar = c(4.2, 4.4, 2.4, 1))
# Bias may be of EITHER sign (it flipped negative on the corrected target), so the
# axis must span 0 and the bars both ways -- ylim = c(0, max(bias)*1.15) silently
# inverted the panel and clipped five of six bars once every bias went negative.
.bl <- range(c(0, bias * 1.15))
bp <- barplot(bias[MODEL_ORDER], col = MODEL_COL[MODEL_ORDER], border = NA, las = 2,
              cex.names = 0.8, ylab = "Mean bias (tC/ha)",
              main = "Systematic bias by model", font.main = 1, ylim = .bl)
abline(h = 0, col = "grey40", lwd = 1)
text(bp, bias[MODEL_ORDER], sprintf("%+.1f", bias[MODEL_ORDER]),
     pos = ifelse(bias[MODEL_ORDER] >= 0, 3, 1), cex = 0.75, xpd = NA)

# --- basal-area legend in the 8th cell ---
plot.new()
legend("center", bty = "n", cex = 1.05, pt.cex = 1.6, pch = 19,
       col = BASAL_COL[seq_len(bc$nb - 1)], title = "Stand basal area (quintiles, m2/ha)\nscatterplot point colour",
       legend = bc$labels)
legend("bottom", bty = "n", cex = 0.9, lwd = c(1.3, 1.3, 1.5), lty = c(2, 3, 1),
       col = c("grey40", "steelblue4", "firebrick"),
       legend = c("1:1 (perfect)", "mean prediction (no skill)", "OLS fit"))
dev.off()
cat("wrote manuscript/figures/S2_obs_vs_pred_calib.png\n")
cat("calib R2:", paste(sprintf("%s %.3f", names(r2c), r2c), collapse = "  "), "\n")
cat("bias:", paste(sprintf("%s %+.1f", names(bias), bias), collapse = "  "), "\n")
