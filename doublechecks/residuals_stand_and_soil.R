# =============================================================================
# residuals_stand_and_soil.R — the residual in two parts: the 1985 stand and
# the soil
#
# WHY. Basal area 1985 explains most of the local residual; fertility, texture
# and stoniness each add a little at a given basal area, and they overlap
# (fine soils are the fertile ones). One JOINT model estimates each soil effect
# with the others and basal area held fixed, so nothing is double-counted.
#
# WHAT. Residual = mean of Yasso07 and Yasso15 plot-mean log residuals.
#   model: r ~ ns(basal area 1985, 3) + fertility + texture + stoniness
#   fertility  NFI site type: fertile (1-2), mesic (3), sub-xeric (4), xeric (5)
#   texture    1985 field soil type (maap_laatu): coarse = coarse moraine,
#              gravel, sand (3,5,6); fine moraine (4); fine = fine sand, silt,
#              clay (7-9); bedrock/rocky/peat (0-2, 12 plots) left out
#   stoniness  BioSoil coarse fragments: < 20, 20-50, > 50 % vol.
# The regression is additive, so each plot's predicted residual decomposes as
#   average residual + basal-area part + fertility part + texture part + stoniness part
# with every part measured FROM THE AVERAGE PLOT (class k: every plot set to k,
# others as observed, minus the average plot). Left = basal-area part over a
# grid, right = the soil parts. 95% CI by the delta method (linear in the
# coefficients). +0.2 ≈ +22% carbon.
#   Left : basal-area effect (same construction over a grid of basal area)
#   Right: the soil effects, dot-and-whisker, rows grouped by variable
#
# Run from repo root:  Rscript doublechecks/residuals_stand_and_soil.R
# Writes doublechecks/figures/residuals_stand_and_soil.png
# =============================================================================

suppressPackageStartupMessages(library(splines))
OUT <- "doublechecks/figures"
B   <- readRDS(file.path(OUT, "residuals_rf_everything_Yasso07_Yasso15.rds")); X <- B$X
num <- function(x) suppressWarnings(as.integer(as.character(x)))
ft  <- num(X$cov.kasvup_tyyppi); ml <- num(X$k85.maap_laatu); st <- X$cov.CoarseFragments
d <- data.frame(
  r  = X$resid, ba = X$wb.ppa_kaikki_85,
  fert  = factor(c("Fertile (1-2)", "Fertile (1-2)", "Mesic (3)", "Sub-xeric (4)", "Xeric (5)")[ft],
                 levels = c("Fertile (1-2)", "Mesic (3)", "Sub-xeric (4)", "Xeric (5)")),
  tex   = factor(ifelse(ml %in% c(3, 5, 6), "Coarse (sand, gravel, coarse moraine)",
                 ifelse(ml == 4, "Fine moraine", ifelse(ml %in% 7:9, "Fine (fine sand, silt, clay)", NA))),
                 levels = c("Coarse (sand, gravel, coarse moraine)", "Fine moraine", "Fine (fine sand, silt, clay)")),
  stone = cut(st, c(-Inf, 20, 50, Inf), labels = c("< 20 %", "20-50 %", "> 50 %")))
d <- d[complete.cases(d), ]
m <- lm(r ~ ns(ba, 3) + fert + tex + stone, d)
V <- vcov(m); bh <- coef(m)
r2 <- function(f) summary(lm(f, d))$r.squared
R2 <- c(ba = r2(r ~ ns(ba, 3)), full = summary(m)$r.squared,
        drop_fert = r2(r ~ ns(ba, 3) + tex + stone), drop_tex = r2(r ~ ns(ba, 3) + fert + stone),
        drop_stone = r2(r ~ ns(ba, 3) + fert + tex))
cat(sprintf("n %d | R² basal area %.3f | full %.3f | unique: fertility %.3f, texture %.3f, stoniness %.3f\n",
            nrow(d), R2["ba"], R2["full"], R2["full"] - R2["drop_fert"], R2["full"] - R2["drop_tex"],
            R2["full"] - R2["drop_stone"]))

# deviation of class k of variable v from the average plot, with delta-method SE
Xobs <- colMeans(model.matrix(m))                  # the average plot
dev <- function(v, k) {   # every plot set to class k, others as observed, minus the average plot
  nd <- d; nd[[v]] <- factor(k, levels = levels(d[[v]]))
  a <- colMeans(model.matrix(delete.response(terms(m)), nd)) - Xobs
  c(est = sum(a * bh), se = sqrt(drop(t(a) %*% V %*% a)))
}
rows <- do.call(rbind, lapply(c("fert", "tex", "stone"), function(v)
  data.frame(var = v, cls = levels(d[[v]]), n = as.vector(table(d[[v]])),
             t(sapply(levels(d[[v]]), function(k) dev(v, k))))))
VAR_LAB <- c(fert = "Site fertility (NFI site type)", tex = "Soil texture (1985 field record)",
             stone = "Stoniness (coarse fragments, % vol.)")
print(transform(rows, est = round(est, 3), se = round(se, 3)), row.names = FALSE)

# basal-area curve as deviation from the average plot, over a grid
grid <- seq(quantile(d$ba, .02), quantile(d$ba, .98), length.out = 60)
bacurve <- t(sapply(grid, function(b) {
  nd <- d; nd$ba <- b
  a <- colMeans(model.matrix(delete.response(terms(m)), nd)) - Xobs   # basal-area part, from the average plot
  c(est = sum(a * bh), se = sqrt(drop(t(a) %*% V %*% a)))
}))

# --- draw ------------------------------------------------------------------------
INK <- "#3d3d3a"; GRID <- "#e6e5df"; MUTED <- "#8a8983"
# One hue per soil variable (categorical slots 1-3: blue, orange, aqua), classes
# shaded light -> dark within a block, DARK = more carbon than predicted
# (fertile, fine, stone-poor). Lightest shade kept >= ~2:1 against white.
HUE <- list(fert = c("#86b6ef", "#3987e5", "#1c5cab", "#0d366b"),       # blue ramp 250/400/550/700
            tex  = colorRampPalette(c("#f0a07e", "#eb6834", "#9c3d14"))(3),
            stone = colorRampPalette(c("#6fcfa8", "#1baf7a", "#0b6a48"))(3))
DARK_FIRST <- c(fert = TRUE, tex = FALSE, stone = TRUE)   # row order: is the first row the "more carbon" end?
rowcol <- unlist(lapply(c("fert", "tex", "stone"), function(v) {
  k <- levels(d[[v]]); h <- HUE[[v]][seq_along(k)]
  if (DARK_FIRST[[v]]) rev(h) else h }))
blockcol <- c(fert = HUE$fert[3], tex = HUE$tex[2], stone = HUE$stone[3])
pct <- function(x) sprintf("%+d%%", round(100 * (exp(x) - 1)))
png(file.path(OUT, "residuals_stand_and_soil.png"), 14.5, 7.4, units = "in", res = 150)
layout(matrix(1:2, 1), widths = c(1, 1.15))
par(oma = c(3.2, 0, 0, 0))
par(mgp = c(2.4, 0.6, 0), tcl = -0.3, col.axis = INK, col.lab = INK, fg = INK)

# left: stand
par(mar = c(4.6, 5.4, 4.6, 1.2))
# partial residuals: model residual + the plot's own basal-area term (centred on
# the data mean, i.e. the same deviation-from-average construction as the curve)
# partial residuals: regression residual + the plot's own basal-area part
part <- resid(m) + predict(m, type = "terms")[, "ns(ba, 3)"]
YL <- range(c(quantile(part, c(.01, .99)), bacurve[, "est"] + 2 * bacurve[, "se"]))
plot(NA, xlim = range(grid), ylim = YL, bty = "l", xlab = "Basal area 1985 (m2/ha)",
     ylab = "Basal-area part of the residual")
abline(h = 0, col = GRID, lwd = 2)
points(d$ba, part, pch = 16, cex = 0.55, col = adjustcolor(MUTED, 0.45))
polygon(c(grid, rev(grid)), c(bacurve[, "est"] - 1.96 * bacurve[, "se"], rev(bacurve[, "est"] + 1.96 * bacurve[, "se"])),
        col = adjustcolor(INK, 0.15), border = NA)
lines(grid, bacurve[, "est"], col = INK, lwd = 3)
# direction notes on the left, outside the axis title (text runs upward,
# so the arrows point up/down), split at zero
mtext("\u2190 overpredicting", 2, 3.6, at = -0.1, adj = 1, cex = 0.75, col = MUTED)
mtext("underpredicting \u2192", 2, 3.6, at = 0.1, adj = 0, cex = 0.75, col = MUTED)
mtext("Basal-area part of the residual (stand in 1985)", 3, 2.2, adj = 0, font = 2)
mtext(sprintf("Curve ± 95%% CI; points: plots (partial residuals). Basal area alone R² %.2f",
              R2["ba"]), 3, 1.0, adj = 0, cex = 0.72)

# right: soil dot-and-whisker
par(mar = c(5.2, 17, 4.6, 3.5))
gap <- 0.8; y <- numeric(nrow(rows)); pos <- 0; yv <- c()
for (v in c("fert", "tex", "stone")) { k <- which(rows$var == v)
  y[k] <- pos - seq_along(k) + 1; yv[v] <- max(y[k]) + 0.75; pos <- min(y[k]) - 1 - gap }
lo <- rows$est - 1.96 * rows$se; hi <- rows$est + 1.96 * rows$se
XL <- range(c(lo, hi, 0)) * 1.1
plot(NA, xlim = XL, ylim = range(y) + c(-0.6, 1.0), bty = "l", yaxt = "n",
     xlab = "Soil part of the residual", ylab = "")
abline(v = 0, col = GRID, lwd = 2)
for (v in names(yv)[-1]) abline(h = max(y[rows$var == v]) + (1 + gap) / 2 + 0.1, col = GRID, lwd = 1)
segments(lo, y, hi, y, col = rowcol, lwd = 2.5)
points(rows$est, y, pch = 21, bg = rowcol, col = "white", cex = 1.9, lwd = 1.5)
axis(2, y, sprintf("%s  (n %d)", rows$cls, rows$n), las = 1, tick = FALSE, cex.axis = 0.8)
for (v in names(yv)) mtext(VAR_LAB[[v]], 2, 16, at = yv[[v]], las = 1, adj = 0, font = 2, cex = 0.9, col = blockcol[[v]])
text(hi, y, pct(rows$est), pos = 4, cex = 0.72, col = INK, xpd = NA)
mtext("Soil parts of the residual, at the same basal area", 3, 2.2, adj = 0, font = 2)
mtext(sprintf("Each class ± 95%% CI; a plot takes one row from each block.\nSoil raises R² %.2f -> %.2f. Labels: SOC relative to the average plot.",
              R2["ba"], R2["full"]), 3, 0.3, adj = 0, cex = 0.72)
mtext("\u2190 overpredicting", 1, 3.8, at = XL[1], adj = 0, cex = 0.7, col = MUTED)
mtext("underpredicting \u2192", 1, 3.8, at = XL[2], adj = 1, cex = 0.7, col = MUTED)
mtext(sprintf("Residual = log(observed / predicted SOC), plot mean of Yasso07 and Yasso15 (n = %d). Linear regression on 1985 basal area (spline) + fertility + texture + stoniness, so that\npredicted residual = average (%.2f) + basal-area part (left) + fertility, texture and stoniness parts (right, one row each). Positive = models predict too little SOC.",
              nrow(d), mean(d$r)), 1, 0.2, outer = TRUE, cex = 0.72, col = INK, adj = 0.03)
dev.off()
cat("Wrote", file.path(OUT, "residuals_stand_and_soil.png"), "\n")
