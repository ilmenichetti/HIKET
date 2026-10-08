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

png("manuscript/figures/F18_eqinit_tradeoff.png", width = 11, height = 7.2, units = "in", res = 200)
par(mfrow = c(2, 3), mar = c(4, 4.4, 2.4, 0.8), oma = c(1.4, 0, 2.2, 0), mgp = c(2.5, 0.6, 0), las = 1)
med <- list()
for (M in FIG_MODELS) {
  if (!M %in% D$models) { plot.new(); title(main = M); text(0.5, 0.5, "equilibrium arm\nnot available", col = "grey40"); next }
  d  <- D$d[[M]]; Jb <- d$prod$J_bar
  pr <- d$prod$thin; eq <- d$eq$thin
  pr <- pr[is.finite(pr$mtt), ]; eq <- eq[is.finite(eq$mtt), ]
  xr <- range(c(pr$si, eq$si) * Jb, Jb); yr <- range(c(pr$mtt, eq$mtt))
  plot(NA, xlim = xr, ylim = yr, log = "xy", main = M,
       xlab = expression("Effective litter input  (tC ha"^-1*" yr"^-1*")"), ylab = "Mean transit time (yr)")
  abline(v = Jb, lty = 3, col = "grey40")
  points(pr$si * Jb, pr$mtt, pch = 16, cex = 0.5, col = adjustcolor(MODEL_COL[[M]], 0.45))
  points(eq$si * Jb, eq$mtt, pch = 16, cex = 0.5, col = adjustcolor("grey30", 0.35))
  m0 <- c(median(pr$si) * Jb, median(pr$mtt)); m1 <- c(median(eq$si) * Jb, median(eq$mtt))
  if (max(abs(log(m1 / m0))) > 1e-6)                  # identical arms (tests): no arrow
    arrows(m0[1], m0[2], m1[1], m1[2], length = 0.09, lwd = 2.2, col = "grey10")
  points(m0[1], m0[2], pch = 21, bg = MODEL_COL[[M]], col = "grey10", cex = 1.7, lwd = 1.5)
  points(m1[1], m1[2], pch = 21, bg = "white", col = "grey10", cex = 1.7, lwd = 1.5)
  med[[M]] <- data.frame(model = M, flux_prod = m0[1], flux_eq = m1[1], mtt_prod = m0[2], mtt_eq = m1[2])
  if (M == D$models[1])
    legend("topright", bty = "n", cex = 0.8, legend = c("transient start", "equilibrium start", "tree litter alone"),
           pch = c(21, 21, NA), pt.bg = c(MODEL_COL[[M]], "white", NA), col = c("grey10", "grey10", "grey40"),
           lty = c(NA, NA, 3), pt.cex = 1.4)
}
mtext("Litter input against transit time: transient vs equilibrium start", outer = TRUE,
      line = 0.4, cex = 1.05, font = 2)
mtext(sprintf("posterior draws (thinned) and medians; transit time at unit input and dataset-mean climate; runs %s / %s",
              paste(unique(substr(D$stamp$prod, 1, 8)), collapse = "/"), paste(unique(substr(D$stamp$eq, 1, 8)), collapse = "/")),
      side = 1, outer = TRUE, line = 0.3, cex = 0.62, col = "grey35", adj = 0.01)
dev.off()
cat("\nWrote manuscript/figures/F18_eqinit_tradeoff.png\n\n")
print(do.call(rbind, med), row.names = FALSE, digits = 3)
