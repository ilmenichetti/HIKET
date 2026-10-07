# =============================================================================
# residuals_basal_area_fertility.R — the 1985 stand AND site fertility together
#
# WHY. Basal area 1985 is the main residual predictor; forest site type (the
# NFI field Cajander class, i.e. fertility) predicts nothing alone (R² 0.04),
# but at a GIVEN basal area a monotone gradient appears (fertile under-predicted,
# xeric over-predicted; ΔR² 0.03, p < 0.001). Poor sites carry small stands, so
# the two effects mask each other unless shown together.
#
# WHAT. Residual = mean of Yasso07 and Yasso15 plot-mean log residuals (one point
# per plot; + = model under-predicts). Fertility groups: fertile (classes 1-2;
# class 1 has 5 plots), mesic (3), sub-xeric (4), xeric (5).
#   Left : residual vs basal area 1985, points by fertility; curves from ONE
#          model, r ~ ns(basal area, 3) + fertility, i.e. a common shape shifted
#          per group (parallel curves)
#   Right: fertility effect at fixed basal area = mean residual after removing
#          the basal-area curve, per group, with 95% CI
#
# Run from repo root:  Rscript doublechecks/residuals_basal_area_fertility.R
# Writes doublechecks/figures/residuals_basal_area_fertility.png
# =============================================================================

suppressPackageStartupMessages(library(splines))
OUT <- "doublechecks/figures"
B   <- readRDS(file.path(OUT, "residuals_rf_everything_Yasso07_Yasso15.rds"))
d   <- data.frame(r = B$X$resid, ba = B$X$wb.ppa_kaikki_85,
                  ft = suppressWarnings(as.integer(as.character(B$X$cov.kasvup_tyyppi))))
d   <- d[is.finite(d$r) & is.finite(d$ba) & d$ft %in% 1:5, ]
GRP <- c("Fertile (1-2)", "Mesic (3)", "Sub-xeric (4)", "Xeric (5)")
d$g <- factor(GRP[c(1, 1, 2, 3, 4)][d$ft], levels = GRP)
COL <- setNames(c("#6da7ec", "#2a78d6", "#1c5cab", "#0d366b")[4:1], GRP)   # dark = fertile
INK <- "#3d3d3a"; GRID <- "#e6e5df"

m_ba  <- lm(r ~ ns(ba, 3), d)
m_all <- lm(r ~ ns(ba, 3) + g, d)
cat(sprintf("n %d | R² basal area %.3f | + fertility %.3f (ΔR² %.3f, p %.2g)\n", nrow(d),
            summary(m_ba)$r.squared, summary(m_all)$r.squared,
            summary(m_all)$r.squared - summary(m_ba)$r.squared, anova(m_ba, m_all)[2, 6]))

png(file.path(OUT, "residuals_basal_area_fertility.png"), 13, 5.6, units = "in", res = 150)
layout(matrix(1:2, 1), widths = c(1.55, 1))
par(mar = c(4.2, 4.2, 3.2, 1), mgp = c(2.4, 0.6, 0), tcl = -0.3,
    col.axis = INK, col.lab = INK, fg = INK)
YL <- quantile(d$r, c(.005, .995))

# --- left: scatter + parallel curves -----------------------------------------
plot(NA, xlim = c(0, max(d$ba)), ylim = YL, bty = "l",
     xlab = "Basal area 1985 (m2/ha)", ylab = "Log residual (obs/pred), plot mean")
abline(h = 0, col = GRID, lwd = 2)
points(d$ba, d$r, pch = 16, cex = 0.75, col = adjustcolor(COL[d$g], 0.6))
for (k in GRP) {
  rg <- range(d$ba[d$g == k]); xs <- seq(rg[1], rg[2], length.out = 80)
  lines(xs, predict(m_all, data.frame(ba = xs, g = factor(k, levels = GRP))), col = "white", lwd = 5)
  lines(xs, predict(m_all, data.frame(ba = xs, g = factor(k, levels = GRP))), col = COL[[k]], lwd = 2.5)
}
legend("topright", sprintf("%s  n %d", GRP, as.vector(table(d$g))), col = COL, lwd = 2.5, pch = 16,
       bty = "n", cex = 0.85, title = "Site fertility (NFI site type)")
mtext("Residual vs the 1985 stand, by site fertility", 3, 1.4, adj = 0, font = 2)
mtext(sprintf("Curves: one fit, common basal-area shape + a shift per fertility group (R² %.2f; basal area alone %.2f)",
              summary(m_all)$r.squared, summary(m_ba)$r.squared), 3, 0.3, adj = 0, cex = 0.75)

# --- right: fertility effect at fixed basal area ------------------------------
pr <- resid(m_ba)
st <- t(sapply(GRP, function(k) { x <- pr[d$g == k]; c(m = mean(x), se = sd(x) / sqrt(length(x)), n = length(x)) }))
plot(NA, xlim = c(0.5, 4.5), ylim = range(c(st[, "m"] - 2.2 * st[, "se"], st[, "m"] + 2.2 * st[, "se"], 0)),
     bty = "l", xaxt = "n", xlab = "", ylab = "Mean residual after removing basal area")
axis(1, 1:4, sub(" \\(", "\n(", GRP), padj = 0.5, cex.axis = 0.85)
abline(h = 0, col = GRID, lwd = 2)
segments(1:4, st[, "m"] - 1.96 * st[, "se"], 1:4, st[, "m"] + 1.96 * st[, "se"], col = COL, lwd = 2.5)
points(1:4, st[, "m"], pch = 21, bg = COL, col = "white", cex = 2, lwd = 2)
text(1:4, st[, "m"], sprintf("%+.2f", st[, "m"]), pos = 4, offset = 0.8, cex = 0.8, col = INK)
mtext("Fertility effect at a given basal area", 3, 1.4, adj = 0, font = 2)
mtext(sprintf("Means ± 95%% CI; + = more carbon than predicted (ΔR² %.3f, p %.1g)",
              summary(m_all)$r.squared - summary(m_ba)$r.squared, anova(m_ba, m_all)[2, 6]),
      3, 0.3, adj = 0, cex = 0.75)
dev.off()
print(round(st, 3))
cat("Wrote", file.path(OUT, "residuals_basal_area_fertility.png"), "\n")
