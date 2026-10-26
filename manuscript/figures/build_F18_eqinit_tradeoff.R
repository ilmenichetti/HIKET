# =============================================================================
# build_F18_eqinit_tradeoff.R   (2026-10-08)   -- SUPPLEMENT
#
# F18 -- WHERE EACH MODEL SITS ON THE INPUT-TURNOVER TRADE-OFF, per start.
# One panel per model. x = effective litter input (sigma_input x J_bar, tC/ha/yr),
# y = intrinsic mean transit time (yr), both log scale. Posterior draws of the
# transient start in model colour, of the equilibrium start in grey; the arrow
# joins the two posterior medians. Vertical dotted line: the tree-litter product
# alone (sigma_input = 1).
#
# Pre-registered reading: without the inherited deficit, the equilibrium start
# moves towards MORE input and/or SHORTER transit time; SP1 moves least.
# Transit time: as F14/S13 (unit input, dataset-mean reference) -- eqinit_draws.R.
#
# Usage:  Rscript manuscript/figures/eqinit_draws.R      (once per run pair)
#         Rscript manuscript/figures/build_F18_eqinit_tradeoff.R
# =============================================================================

setwd("/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling")
source("manuscript/figures/run_ids_eqinit.R")
source("manuscript/figures/model_palette.R")
D <- readRDS("manuscript/figures/eqinit_draws.rds")
if (!identical(D$stamp, list(prod = RID[EQ_MODELS], eq = RID_EQ[EQ_MODELS])))
  stop("eqinit_draws.rds was built from other RUN_IDs -- rerun manuscript/figures/eqinit_draws.R")

# Print size (0.95 x 6.53 in, 9 pt); title and run note live in the caption (fig:F18).
png("manuscript/figures/F18_eqinit_tradeoff.png", width = 6.2, height = 4.7, units = "in", res = 300, pointsize = 9)
par(mfrow = c(2, 3), mar = c(3.2, 3.6, 1.6, 0.4), oma = c(0, 0, 0, 0), mgp = c(2.0, 0.5, 0), las = 1, tcl = -0.3)
par(cex = 1)   # layout/mfrow shrink text to 0.66 (or 0.83 in 2x2); keep the 9 pt
med <- list()
for (M in FIG_MODELS) {
  if (!M %in% D$models) { plot.new(); title(main = M); text(0.5, 0.5, "balanced-start arm\nnot available", col = "grey40"); next }
  d  <- D$d[[M]]; Jb <- d$prod$J_bar
  pr <- d$prod$thin; eq <- d$eq$thin
  pr <- pr[is.finite(pr$mtt), ]; eq <- eq[is.finite(eq$mtt), ]
  xr <- range(c(pr$si, eq$si) * Jb, Jb); yr <- range(c(pr$mtt, eq$mtt))
  k <- match(M, FIG_MODELS)
  plot(NA, xlim = xr, ylim = yr, log = "xy", main = "",
       xlab = expression("Input  (tC ha"^-1*" yr"^-1*")"), ylab = if (k %in% c(1, 4)) "Mean transit time (yr)" else "")
  mtext(sprintf("(%s) %s", letters[k], M), 3, 0.3, adj = 0, font = 2)
  abline(v = Jb, lty = 3, col = "grey40")
  points(pr$si * Jb, pr$mtt, pch = 16, cex = 0.5, col = adjustcolor(MODEL_COL[[M]], 0.45))
  points(eq$si * Jb, eq$mtt, pch = 16, cex = 0.5, col = adjustcolor("grey30", 0.35))
  m0 <- c(median(pr$si) * Jb, median(pr$mtt)); m1 <- c(median(eq$si) * Jb, median(eq$mtt))
  if (max(abs(log(m1 / m0))) > 1e-6)                  # identical arms (tests): no arrow
    arrows(m0[1], m0[2], m1[1], m1[2], length = 0.09, lwd = 2.2, col = "grey10")
  points(m0[1], m0[2], pch = 21, bg = MODEL_COL[[M]], col = "grey10", cex = 1.4, lwd = 1.3)
  points(m1[1], m1[2], pch = 21, bg = "white", col = "grey10", cex = 1.4, lwd = 1.3)
  med[[M]] <- data.frame(model = M, flux_prod = m0[1], flux_eq = m1[1], mtt_prod = m0[2], mtt_eq = m1[2])
  if (M == D$models[1])
    legend("topright", bty = "n", cex = 0.85, legend = c("transient start", "balanced start", "tree litter alone"),
           pch = c(21, 21, NA), pt.bg = c(MODEL_COL[[M]], "white", NA), col = c("grey10", "grey10", "grey40"),
           lty = c(NA, NA, 3), pt.cex = 1.2)
}
dev.off()
cat("\nWrote manuscript/figures/F18_eqinit_tradeoff.png\n\n")
print(do.call(rbind, med), row.names = FALSE, digits = 3)
