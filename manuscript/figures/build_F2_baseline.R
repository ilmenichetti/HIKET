setwd("/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling")
# F2 (redesign) -- the uncalibrated baseline, GENERAL (cross-plot mean, all 447 plots),
# built like F3 rather than 4 example plots. Yasso20 run forward at published defaults
# with the conventional equilibrium (steady-state) start. Two failures in one glance:
#   (1) FLAT  -- an equilibrium start is already at its steady state -> no accumulation;
#   (2) LOW   -- boreal-blind published defaults sit far below the observed stocks.
# Sources: baseline predictions (defaults, steady-state) + observed campaign means.

pred <- readRDS("Calibration_real_data_transient/diagnostics/Yasso20_baseline/Yasso20_baseline_predictions.rds")
traj <- aggregate(total_soc ~ year, pred, function(x)
          c(m = mean(x), se = sd(x) / sqrt(length(x))))
tj <- data.frame(year = traj$year, m = traj$total_soc[, "m"], se = traj$total_soc[, "se"])

om <- readRDS("Data/model_inputs/Yasso20_inputs_20260710_102431.rds")$obs_meta
obs <- do.call(rbind, lapply(names(om), function(pid) {
  z <- om[[pid]]; if (length(z$soc_obs) == 0) return(NULL)
  data.frame(year = 1984L + z$idx, soc = z$soc_obs)
}))
cm <- aggregate(soc ~ year, obs, function(x)
        c(m = mean(x), lo = mean(x) - 1.96*sd(x)/sqrt(length(x)),
          hi = mean(x) + 1.96*sd(x)/sqrt(length(x))))
cm <- data.frame(year = cm$year, m = cm$soc[,"m"], lo = cm$soc[,"lo"], hi = cm$soc[,"hi"])

png("manuscript/figures/F2_baseline.png", width = 8.4, height = 5.2, units = "in", res = 200)
par(mar = c(4.0, 4.6, 3.0, 1.2), mgp = c(2.6, 0.7, 0), las = 1)
plot(NA, xlim = c(min(tj$year), max(tj$year) + 2), ylim = c(40, 118), xlab = "Year",
     ylab = "Mean SOC across plots (tC/ha)",
     main = "Uncalibrated baseline: published defaults miss the accumulation")

# under-prediction gap at each campaign (grey droplines)
mb <- approx(tj$year, tj$m, cm$year)$y
lab_side <- ifelse(cm$year == max(cm$year), 2, 4)   # last campaign labels to the left
segments(cm$year, mb, cm$year, cm$m, col = "grey65", lwd = 8, lend = 1)
text(cm$year, (mb + cm$m)/2, sprintf("-%.0f", cm$m - mb), pos = lab_side, offset = 0.5,
     cex = 0.8, col = "grey35", font = 2)

# baseline mean trajectory + cross-plot SE band
polygon(c(tj$year, rev(tj$year)), c(tj$m - tj$se, rev(tj$m + tj$se)),
        col = adjustcolor("steelblue", 0.30), border = NA)
lines(tj$year, tj$m, col = "steelblue", lwd = 2.6)

# observed campaign means + 95% CI
arrows(cm$year, cm$lo, cm$year, cm$hi, angle = 90, code = 3, length = 0.05, col = "firebrick", lwd = 1.8)
points(cm$year, cm$m, pch = 19, col = "firebrick", cex = 1.5)
text(cm$year, cm$hi, c("VMI8", "Biosoil", "Komeetta"),
     pos = ifelse(cm$year == max(cm$year), 2, 4), offset = 0.6, cex = 0.78, col = "firebrick")

# annotations for the two failures
text(2004, 54, "equilibrium start -> flat: no accumulation captured", col = "steelblue", font = 3, cex = 0.82, pos = 3)
legend("topleft", inset = c(0.01, 0.02), bty = "n", cex = 0.85,
       legend = c("Yasso20 at published defaults (mean +/- SE, 447 plots)",
                  "Observed campaign mean +/- 95% CI"),
       pch = c(NA, 19), lwd = c(2.6, NA), col = c("steelblue", "firebrick"))
dev.off()
cat("Wrote manuscript/figures/F2_baseline.png\n")
cat(sprintf("Baseline mean SOC: %.1f (1985) .. %.1f (2024);  obs %.0f -> %.0f -> %.0f\n",
            tj$m[1], tj$m[nrow(tj)], cm$m[1], cm$m[2], cm$m[3]))
