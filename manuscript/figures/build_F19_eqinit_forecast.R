# =============================================================================
# build_F19_eqinit_forecast.R   (2026-10-08)
#
# F19 -- WHAT THE START ASSUMPTION DOES TO THE FORECAST.
# The projection holds litter at 2024 and recycles the 2005-2024 climate (no
# climate change), so its slope is the disequilibrium each arm carries out of
# 2024. That is exactly what the initialisation decides, which makes this the
# inventory-relevant comparison.
#   (a) projected change since 2024, per model: solid = transient, dashed =
#       equilibrium start; bold = ensemble (mean of model medians)
#   (b) projected sink over the first 20-yr climate cycle (2025-2044)
#   (c) headroom: equilibrium stock at 2024 litter / 2024 stock - 1
# (b) and (c) use the definitions of build_forward_scenarios.R; rates over WHOLE
# 20-yr cycles because the recycled climate alternates by decade.
# Intervals: 90% over posterior draws of the national mean.
#
# Usage:  Rscript manuscript/figures/build_F19_eqinit_forecast.R
# =============================================================================

setwd("/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling")
source("manuscript/figures/eqinit_common.R")
z  <- load_eqinit_comparison()
MS <- z$models
yr <- as.integer(colnames(z$arms[[MS[1]]]$prod$traj)); py <- yr[yr >= 2024]
OFF <- c(prod = 0.17, eq = -0.17)

dC <- function(m, a) {                                    # median change since 2024
  t <- z$arms[[m]][[a]]$traj
  apply(t[, as.character(py), drop = FALSE] - t[, "2024"], 2, median)
}
fs <- function(m, a, v) q90(z$arms[[m]][[a]]$fore[[v]])

dot_panel <- function(v, xlab, main, scale = 1) {
  s <- do.call(rbind, lapply(MS, function(m) do.call(rbind, lapply(c("prod", "eq"), function(a)
         data.frame(model = m, arm = a, t(fs(m, a, v) * scale))))))
  yy <- setNames(rev(seq_along(FIG_MODELS)), FIG_MODELS)
  xr <- range(c(s$lo, s$hi, 0)); xr <- xr + c(-0.06, 0.06) * diff(xr)
  plot(NA, xlim = xr, ylim = c(0.4, length(FIG_MODELS) + 0.6), yaxt = "n", xlab = xlab, ylab = "", main = main)
  rect(xr[1] - 1e3, 0, 0, length(FIG_MODELS) + 1, col = adjustcolor("#AA3333", 0.06), border = NA)
  abline(v = 0, col = "grey35", lty = 2)
  for (k in seq_len(nrow(s))) {
    y <- yy[[s$model[k]]] + OFF[[s$arm[k]]]; col <- MODEL_COL[[s$model[k]]]
    segments(s$lo[k], y, s$hi[k], y, col = col, lwd = 2.4, lty = ARM_LTY[[s$arm[k]]])
    points(s$med[k], y, pch = 21, cex = 1.5, lwd = 1.5, col = "grey20",
           bg = if (s$arm[k] == "prod") col else "white")
  }
  axis(2, at = yy, labels = ifelse(FIG_MODELS %in% MS, FIG_MODELS, paste(FIG_MODELS, "(n/a)")),
       tick = FALSE)
  s
}

png("manuscript/figures/F19_eqinit_forecast.png", width = 12.5, height = 5, units = "in", res = 200)
layout(matrix(1:3, 1), widths = c(1.35, 1, 1))
par(oma = c(1.4, 0, 2.4, 0), mgp = c(2.5, 0.6, 0), las = 1)

par(mar = c(4.2, 4.6, 2.4, 1))
yl <- range(unlist(lapply(MS, function(m) c(dC(m, "prod"), dC(m, "eq")))), 0)
plot(NA, xlim = range(py), ylim = yl, xlab = "Year",
     ylab = expression("Change since 2024  (tC ha"^-1*")"), main = "(a) Projected stock change")
abline(h = 0, col = "grey35", lty = 2)
for (m in MS) for (a in c("prod", "eq")) lines(py, dC(m, a), col = MODEL_COL[[m]], lwd = 1.4, lty = ARM_LTY[[a]])
for (a in c("prod", "eq")) lines(py, rowMeans(sapply(MS, dC, a = a)), col = "grey10", lwd = 3, lty = ARM_LTY[[a]])
legend("topleft", bty = "n", cex = 0.78, legend = c(MS, ARM_LAB),
       col = c(MODEL_COL[MS], "grey10", "grey10"), lty = c(rep(1, length(MS)), 1, 2),
       lwd = c(rep(1.6, length(MS)), 2.6, 2.6))

par(mar = c(4.2, 6.2, 2.4, 1))
s_sink <- dot_panel("sink_first20", expression("tC ha"^-1*" yr"^-1), "(b) Projected sink, 2025-2044")
s_head <- dot_panel("headroom", "% of the 2024 stock", "(c) Headroom to equilibrium", scale = 100)
legend("bottomright", bty = "n", cex = 0.8, legend = ARM_LAB, pch = 21, pt.bg = c("grey50", "white"),
       col = "grey20", lty = c(1, 2), lwd = 2, pt.cex = 1.3)

mtext("Forecast under a transient vs an equilibrium start (litter held at 2024, recent climate recycled)",
      outer = TRUE, line = 0.5, cex = 1.05, font = 2)
mtext(sprintf("%s; intervals: 90%% over posterior draws of the national mean", basis_note(seq_len(z$obs$n))),
      side = 1, outer = TRUE, line = 0.3, cex = 0.62, col = "grey35", adj = 0.01)
dev.off()

cat("\nWrote manuscript/figures/F19_eqinit_forecast.png\n\n")
tab <- merge(s_sink[, c("model", "arm", "med")], s_head[, c("model", "arm", "med")],
             by = c("model", "arm"), suffixes = c(".sink_2025_44", ".headroom_pct"))
print(transform(tab, med.sink_2025_44 = round(med.sink_2025_44, 3), med.headroom_pct = round(med.headroom_pct, 1)),
      row.names = FALSE)
