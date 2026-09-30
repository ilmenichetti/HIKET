source("manuscript/figures/run_ids.R")   # auto-selects current RUN_IDs
setwd("/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling")
# F3_mean_SOCalt -- the OPENING figure: the observed campaign means ALONE, no models.
#
# WHY A SECOND VERSION OF F3. F3_mean_soc answers "do the six models track the
# campaigns?", which is a Results question; opening the Introduction with it puts the
# models in the reader's head before the problem has been stated. This one carries only
# the datum the Introduction needs: the inventory measured the stock three times, it
# rose, and the second rise is much smaller than the first. Nothing here is model-derived.
#
# SAME BASIS AS F3 (obs_basis.R): balanced plot set (observed in all three campaigns),
# unweighted, whole profile, true observation years -- so the two figures are directly
# comparable and F3 can be shown later as "and now the models".
#
# THE TWO SEGMENTS ARE THE STORY. Their slopes are printed on them, with 95% intervals
# from the PAIRED plot-level difference (the campaigns share plots, so the paired SE is
# the right one -- it is roughly half the unpaired one and is what the calibration sees).
#
# AREA-WEIGHTED SERIES (added 2026-09-30, Aleksi's question). The unweighted series is
# the average SAMPLED PLOT, which is what the unweighted model means are compared with.
# The network samples the North at 1/3 the southern density (design weight 3 vs 1), so
# the average HECTARE needs weights; it is drawn alongside, open markers and dashed,
# on the SAME plots. The calibration itself stays unweighted: each plot is one
# observation, and a weight of 3 would count the northern plots as three measurements.
#
# ⚠ Deliberately NO sub-year means inside VMI8: those groups are different plot SUBSETS
# and differ geographically, not temporally (see the warning in obs_basis.R). The
# decade-wide sampling span is shown as a bar instead, which is the honest form.
source("manuscript/figures/model_palette.R")
source("manuscript/figures/obs_basis.R")

om  <- readRDS(sprintf("Data/model_inputs/Yasso20_inputs_%s.rds", RID[["Yasso20"]]))$obs_meta
BAL <- balanced_plots(om)
obs <- do.call(rbind, lapply(as.character(BAL), function(p) {
  z <- om[[p]]
  data.frame(plot = as.integer(p), year = 1984L + z$idx, soc = z$soc_obs)
}))
obs$camp <- campaign_of(obs$year)
cm <- obs_campaigns(om, BAL)          # one marker per campaign, at its true mean year
message("F3_mean_SOCalt ", basis_note(BAL))

# design weights (3 North, 1 South) from the homogenised baseline, one per plot
.P  <- read.csv("Data/SOC_homogeneized/soc_homogenized_plot.csv")
WT  <- tapply(.P$weight, .P$plot_id, function(w) w[1])
obs$w <- as.numeric(WT[as.character(obs$plot)])
if (anyNA(obs$w)) stop("F3_mean_SOCalt: design weight missing for ", sum(is.na(obs$w)), " rows")
wmean_ci <- function(x, w) {                           # weighted mean, linearised SE
  m <- sum(w * x) / sum(w); n <- length(x)
  se <- sqrt(n / (n - 1) * sum(w^2 * (x - m)^2)) / sum(w)
  c(m = m, lo = m - 1.96 * se, hi = m + 1.96 * se)
}
cw <- do.call(rbind, lapply(1:3, function(k) wmean_ci(obs$soc[obs$camp == k], obs$w[obs$camp == k])))

CCOL <- unname(CAMPAIGN_COL[CAMPAIGN_ORDER])           # light -> dark with time
rate <- diff(cm$m) / diff(cm$year)                     # tC/ha/yr per interval
wide <- reshape(obs[, c("plot", "camp", "soc")], idvar = "plot",
                timevar = "camp", direction = "wide")
rate_ci <- t(sapply(1:2, function(k) {                 # PAIRED difference CI
  d  <- wide[[paste0("soc.", k + 1)]] - wide[[paste0("soc.", k)]]
  se <- sd(d) / sqrt(length(d))
  c(mean(d) - 1.96 * se, mean(d) + 1.96 * se) / (cm$year[k + 1] - cm$year[k])
}))
wwide <- merge(wide, unique(obs[, c("plot", "w")]), by = "plot")
rate_w <- t(sapply(1:2, function(k) {                  # weighted PAIRED rate + CI
  d <- wwide[[paste0("soc.", k + 1)]] - wwide[[paste0("soc.", k)]]
  wmean_ci(d, wwide$w) / (cm$year[k + 1] - cm$year[k])
}))

png("manuscript/figures/F3_mean_SOCalt.png", width = 8.6, height = 5.6, units = "in", res = 200)
par(mar = c(5.8, 5.2, 3.2, 1.6), mgp = c(3.0, 0.8, 0), las = 1,
    cex.axis = 1.15, cex.lab = 1.3)
yl <- range(cm$lo, cm$hi, cw[, "lo"], cw[, "hi"]) + c(-3.0, 3.5)
plot(NA, xlim = c(1985, 2027), ylim = yl, xlab = "", ylab = "",
     main = "Observed SOC: the mean stock across the three campaigns")
title(xlab = "Year", line = 2.6, cex.lab = 1.3)
title(ylab = "Mean SOC across plots (tC/ha)", line = 3.3, cex.lab = 1.3)
abline(h = pretty(yl), col = "grey93")

segments(cm$year[-3], cm$m[-3], cm$year[-1], cm$m[-1], lwd = 3.4, col = "grey25")
# VMI8 sampling span, drawn low so it cannot be read as part of the trajectory
.sp <- range(obs$year[obs$camp == 1]); .yb <- yl[1] + 0.05 * diff(yl)
segments(.sp[1], .yb, .sp[2], .yb, lwd = 2.2, col = CCOL[1], lend = 1)
segments(.sp, .yb - 0.013 * diff(yl), .sp, .yb + 0.013 * diff(yl), lwd = 2.2, col = CCOL[1])
text(mean(.sp), .yb, "VMI8 sampled 1986-1995", pos = 1, offset = 0.5, cex = 0.85, col = CCOL[2])

XW <- cm$year + 0.9                                   # weighted series, nudged right
segments(XW[-3], cw[-3, "m"], XW[-1], cw[-1, "m"], lwd = 2.4, lty = 2, col = "grey45")
arrows(XW, cw[, "lo"], XW, cw[, "hi"], angle = 90, code = 3, length = 0.045,
       col = "grey45", lwd = 1.6)
points(XW, cw[, "m"], pch = 21, bg = "white", col = CCOL, cex = 1.8, lwd = 2.2)
for (k in 1:2)
  text(mean(XW[k:(k + 1)]), mean(cw[k:(k + 1), "m"]),
       sprintf("weighted %+.2f [%+.2f, %+.2f]", rate_w[k, "m"], rate_w[k, "lo"], rate_w[k, "hi"]),
       pos = 1, offset = c(2.2, 1.4)[k], cex = 0.8, font = 3, col = "grey40")
legend("topleft", bty = "n", cex = 0.85, inset = c(0.01, 0.0),
       pch = 21, pt.bg = c("grey60", "white"), col = "grey15", lty = c(1, 2), lwd = c(3, 2.4),
       legend = c("unweighted: the average sampled plot",
                  "area-weighted (North x3): the average hectare"))

arrows(cm$year, cm$lo, cm$year, cm$hi, angle = 90, code = 3, length = 0.055,
       col = "grey25", lwd = 2.2)
points(cm$year, cm$m, pch = 21, bg = CCOL, col = "grey15", cex = 2.3, lwd = 1.5)
text(cm$year, cm$hi, sprintf("%s\n%.1f", CAMPAIGN_ORDER, cm$m),
     pos = c(4, 3, 2), offset = 0.85, cex = 0.95, font = 2)
for (k in 1:2) {
  xm <- mean(cm$year[k:(k + 1)]); ym <- mean(cm$m[k:(k + 1)])
  text(xm, ym, sprintf("%+.2f tC/ha/yr\n[%+.2f, %+.2f]", rate[k], rate_ci[k, 1], rate_ci[k, 2]),
       pos = c(3, 1)[k], offset = c(1.5, 1.0)[k], cex = 0.92, col = "grey25")
}
mtext(sprintf("plots measured in all three campaigns (n = %d); bars = 95%% CI of the mean",
              length(BAL)), side = 1, line = 3.5, cex = 0.8, col = "grey35")
mtext("organic + mineral soil to 1 m or auger refusal: measured to 40 cm, extrapolated below",
      side = 1, line = 4.3, cex = 0.8, col = "grey35")
dev.off()
cat(sprintf("wrote F3_mean_SOCalt.png  means %.1f -> %.1f -> %.1f; rates %+.3f, %+.3f\n",
            cm$m[1], cm$m[2], cm$m[3], rate[1], rate[2]))
cat(sprintf("  weighted: means %.1f -> %.1f -> %.1f; rates %+.3f, %+.3f; full %+.3f\n",
            cw[1, "m"], cw[2, "m"], cw[3, "m"], rate_w[1, "m"], rate_w[2, "m"],
            (cw[3, "m"] - cw[1, "m"]) / (cm$year[3] - cm$year[1])))
cat(sprintf("  unweighted full %+.3f; North share of plots %.3f\n",
            (cm$m[3] - cm$m[1]) / (cm$year[3] - cm$year[1]), mean(unique(obs[, c("plot","w")])$w == 3)))
