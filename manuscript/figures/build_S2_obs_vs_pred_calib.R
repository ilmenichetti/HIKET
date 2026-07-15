setwd("/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling")
# S2 (supplement) -- calibration obs-vs-pred, the in-sample companion to the main-text
# holdout F6. COLOUR CODING (per user, 2026-07-15):
#   * the six scatterplots -> coloured by stand basal-area quintile, SAME as F6 (Green-Gold);
#   * the systematic-bias barplot -> the shared per-model palette (Temperature Diverging).
# Data: <MODEL> predictive bundle -> residuals_df (calibration rows) + metrics_calib.
source("manuscript/figures/model_palette.R")

runs <- "Calibration_real_data_transient/runs"
bundles <- c(
  SP1     = "SP1_posterior_predictive_20260710_104903.rds",
  TP2     = "TP2_posterior_predictive_20260710_104904.rds",
  TP3     = "TP3_posterior_predictive_20260710_104904.rds",
  Yasso07 = "Yasso07_posterior_predictive_20260710_104902.rds",
  Yasso15 = "Yasso15_posterior_predictive_20260710_104902.rds",
  Yasso20 = "Yasso20_posterior_predictive_20260710_102431.rds")

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
  abline(0, 1, col = "grey40", lwd = 1.3, lty = 2)
  abline(lm(soc_median ~ soc_obs_tCha, d), col = "firebrick", lwd = 1.5)
}

# --- systematic-bias barplot (per-model palette) ---
par(mar = c(4.2, 4.4, 2.4, 1))
bp <- barplot(bias[MODEL_ORDER], col = MODEL_COL[MODEL_ORDER], border = NA, las = 2,
              cex.names = 0.8, ylab = "Mean bias (tC/ha)",
              main = "Systematic bias by model", font.main = 1,
              ylim = c(0, max(bias) * 1.15))
text(bp, bias[MODEL_ORDER], sprintf("%+.1f", bias[MODEL_ORDER]), pos = 3, cex = 0.75, xpd = NA)

# --- basal-area legend in the 8th cell ---
plot.new()
legend("center", bty = "n", cex = 1.05, pt.cex = 1.6, pch = 19,
       col = BASAL_COL[seq_len(bc$nb - 1)], title = "Stand basal area (quintiles, m2/ha)\nscatterplot point colour",
       legend = bc$labels)
dev.off()
cat("wrote manuscript/figures/S2_obs_vs_pred_calib.png\n")
cat("calib R2:", paste(sprintf("%s %.3f", names(r2c), r2c), collapse = "  "), "\n")
cat("bias:", paste(sprintf("%s %+.1f", names(bias), bias), collapse = "  "), "\n")
