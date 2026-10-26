# =============================================================================
# F3_mean_SOCalt_by_region.R   (2026-10-26)
#
# Fig 2 of draft v2 (manuscript/figures/build_F3_mean_SOCalt.R) split into South and
# North. Prompted by doublechecks/ghgi_steady_state_vs_obs.R: the far north LOSES
# carbon on the balanced plots while the South gains, so the national curve averages
# two opposite trajectories.
#
# Same basis as Fig 2: balanced plot set (all three campaigns), whole profile, true
# observation years, campaign markers at their mean sampling year, rates from the
# PAIRED plot-level differences. "North" = design weight 3 (as Fig 2's weighting);
# it matches site_attributes$region == 2 except for one plot. Within a region every
# plot has the same weight, so each panel is both the plot mean and the area mean.
# Shared y-axis span so the two slopes compare by eye.
#
# Run from repo root:  Rscript doublechecks/F3_mean_SOCalt_by_region.R
# Output: doublechecks/figures/F3_mean_SOCalt_by_region.png
# =============================================================================
suppressMessages({
  source("manuscript/figures/run_ids.R"); source("manuscript/figures/obs_basis.R")
  source("manuscript/figures/model_palette.R")
})
om  <- readRDS(sprintf("Data/model_inputs/Yasso20_inputs_%s.rds", RID[["Yasso20"]]))$obs_meta
BAL <- balanced_plots(om)
P   <- read.csv("Data/SOC_homogeneized/soc_homogenized_plot.csv")
WT  <- tapply(P$weight, P$plot_id, function(w) w[1])
REG <- list(South = BAL[WT[as.character(BAL)] == 1], North = BAL[WT[as.character(BAL)] == 3])
stopifnot(!anyNA(unlist(REG)), sum(lengths(REG)) == length(BAL))

summ_region <- function(pl) {
  cm  <- obs_campaigns(om, pl)
  obs <- do.call(rbind, lapply(as.character(pl), function(p) {
    z <- om[[p]]; data.frame(plot = p, year = 1984L + z$idx, soc = z$soc_obs) }))
  obs$camp <- campaign_of(obs$year)
  wide <- reshape(obs[, c("plot", "camp", "soc")], idvar = "plot", timevar = "camp", direction = "wide")
  rate <- t(sapply(list(c(1, 2), c(2, 3), c(1, 3)), function(ab) {
    d <- wide[[paste0("soc.", ab[2])]] - wide[[paste0("soc.", ab[1])]]; dt <- cm$year[ab[2]] - cm$year[ab[1]]
    se <- sd(d) / sqrt(length(d)); c(r = mean(d), lo = mean(d) - 1.96 * se, hi = mean(d) + 1.96 * se) / dt }))
  rownames(rate) <- c("1985-2006", "2006-2024", "1985-2024")
  list(cm = cm, rate = rate, n = length(pl))
}
S <- lapply(REG, summ_region)
for (k in names(S)) {
  cat(sprintf("%s (n = %d): means %s\n", k, S[[k]]$n, paste(round(S[[k]]$cm$m, 1), collapse = " -> ")))
  print(round(S[[k]]$rate, 3))
}

# === Figure (print size, 9 pt) ===
png("doublechecks/figures/F3_mean_SOCalt_by_region.png", width = 6.5, height = 3.6, units = "in",
    res = 300, pointsize = 9)
par(mfrow = c(1, 2)); par(cex = 1)
par(mgp = c(2.2, 0.5, 0), las = 1, tcl = -0.3)
CCOL <- unname(CAMPAIGN_COL[CAMPAIGN_ORDER])
span <- max(sapply(S, function(s) diff(range(s$cm$lo, s$cm$hi)))) + 8
for (k in seq_along(S)) {
  s <- S[[k]]; cm <- s$cm
  mid <- mean(range(cm$lo, cm$hi)); yl <- mid + c(-0.5, 0.5) * span
  par(mar = c(3.0, if (k == 1) 3.8 else 2.6, 1.8, 0.6))
  plot(NA, xlim = c(1983, 2030), ylim = yl, xlab = "",
       ylab = if (k == 1) expression("Mean stock  (tC ha"^-1*")") else "")
  mtext(sprintf("(%s) %s, %d plots", letters[k], names(S)[k], s$n), 3, 0.4, adj = 0, font = 2)
  abline(h = pretty(yl), col = "grey93")
  segments(cm$year[-3], cm$m[-3], cm$year[-1], cm$m[-1], lwd = 3, col = "grey25")
  arrows(cm$year, cm$lo, cm$year, cm$hi, angle = 90, code = 3, length = 0.04, col = "grey25", lwd = 1.8)
  points(cm$year, cm$m, pch = 21, bg = CCOL, col = "grey15", cex = 1.8, lwd = 1.3)
  text(cm$year, cm$hi, sprintf("%.1f", cm$m), pos = 3, offset = 0.4, cex = 0.85, font = 2)
  for (j in 1:2) {
    xm <- mean(cm$year[j:(j + 1)]); ym <- mean(cm$m[j:(j + 1)]); r <- s$rate[j, ]
    text(xm, ym, sprintf("%+.2f\n[%+.2f, %+.2f]", r["r"], r["lo"], r["hi"]),
         pos = 1, offset = 1.6, cex = 0.8, col = "grey25")
  }
  r <- s$rate["1985-2024", ]
  mtext(sprintf("whole span %+.2f [%+.2f, %+.2f] tC/ha/yr", r["r"], r["lo"], r["hi"]),
        1, 1.9, cex = 0.85, col = "grey25")
}
dev.off()
cat("Wrote doublechecks/figures/F3_mean_SOCalt_by_region.png\n")
