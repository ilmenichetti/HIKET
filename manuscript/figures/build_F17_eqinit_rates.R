# =============================================================================
# build_F17_eqinit_rates.R   (2026-10-08)
#
# F17 -- RATES OF CHANGE, transient vs equilibrium start (built; use undecided).
# Three columns: 1985->2006, 2006->2024, 1985->2024. One row per model, two dots
# per row (filled = transient start, open = equilibrium start), median and 90%
# posterior interval of the national-mean rate. Observed rate as a vertical band
# (95% CI over plots).
#
# Same basis as F5: balanced plots, each plot on its OWN campaign years, so the
# denominators are the real spans (~35 yr for 1985->2024, not 39). Uncertainty on
# the model side is POSTERIOR (over draws), not plot scatter -- so the intervals
# are narrower than F5's and comparable between the two arms.
#
# Pre-registered reading (NEXT_RUN_equilibrium_init.md): the equilibrium arm
# under-produces the rates, most visibly 2006->2024; SP1 changes least.
#
# Usage:  Rscript manuscript/figures/build_F17_eqinit_rates.R
# =============================================================================

setwd("/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling")
source("manuscript/figures/eqinit_common.R")
z  <- load_eqinit_comparison()
MS <- z$models
OFF <- c(prod = 0.17, eq = -0.17)

stats <- do.call(rbind, lapply(names(INTERVALS), function(iv)
  do.call(rbind, lapply(MS, function(m) do.call(rbind, lapply(c("prod", "eq"), function(a) {
    r <- z$arms[[m]][[a]]$rates[, iv]
    data.frame(interval = iv, model = m, arm = a, med = median(r),
               lo = quantile(r, .05), hi = quantile(r, .95))
  }))))))

# Print size (0.95 x 6.53 in, 9 pt); model names on the first panel only; title and basis
# note live in the caption (fig:S_eqinit_rates).
png("manuscript/figures/F17_eqinit_rates.png", width = 6.2, height = 3.4, units = "in", res = 300, pointsize = 9)
layout(matrix(1:3, 1), widths = c(1.42, 1, 1))
par(cex = 1)   # layout/mfrow shrink text to 0.66 (or 0.83 in 2x2); keep the 9 pt
par(oma = c(0, 0, 0, 0), mgp = c(2.0, 0.5, 0), las = 1, tcl = -0.3)
yy <- setNames(rev(seq_along(FIG_MODELS)), FIG_MODELS)
for (iv in names(INTERVALS)) {
  s <- stats[stats$interval == iv, ]; o <- z$obs$rates[z$obs$rates$interval == iv, ]
  xr <- range(c(s$lo, s$hi, o$lo, o$hi, 0)); xr <- xr + c(-0.06, 0.06) * diff(xr)
  kp <- match(iv, names(INTERVALS))
  par(mar = c(3.3, if (kp == 1) 5.0 else 0.6, 1.8, 0.5))
  plot(NA, xlim = xr, ylim = c(0.4, length(FIG_MODELS) + 1.75), yaxt = "n",
       xlab = expression("Change  (tC ha"^-1*" yr"^-1*")"), ylab = "", main = "")
  mtext(sprintf("(%s) %s", letters[kp], sub("-", " \u2192 ", iv)), 3, 0.4, adj = 0, font = 2)
  rect(xr[1] - 1, 0, 0, length(FIG_MODELS) + 2, col = adjustcolor("#AA3333", 0.06), border = NA)
  rect(o$lo, 0, o$hi, length(FIG_MODELS) + 2, col = adjustcolor("#4477AA", 0.13), border = NA)
  abline(v = o$obs, col = "#4477AA", lwd = 2.2); abline(v = 0, col = "grey35", lty = 2)
  for (k in seq_len(nrow(s))) {
    y <- yy[[s$model[k]]] + OFF[[s$arm[k]]]; col <- MODEL_COL[[s$model[k]]]
    segments(s$lo[k], y, s$hi[k], y, col = col, lwd = 2.4, lty = ARM_LTY[[s$arm[k]]])
    points(s$med[k], y, pch = 21, cex = 1.2, lwd = 1.3, col = "grey20",
           bg = if (s$arm[k] == "prod") col else "white")
  }
  if (kp == 1) axis(2, at = yy, labels = ifelse(FIG_MODELS %in% MS, FIG_MODELS, paste(FIG_MODELS, "(n/a)")),
       tick = FALSE, cex.axis = 1.0)
  if (iv == names(INTERVALS)[1])
    legend("topleft", bty = "n", cex = 0.8, bg = "white",
           legend = c(ARM_LAB, "observed (95% CI)"),
           pch = c(21, 21, NA), pt.bg = c("grey50", "white", NA), col = c("grey20", "grey20", "#4477AA"),
           lty = c(1, 2, 1), lwd = c(2, 2, 2.2), pt.cex = 1.1)
}
dev.off()

cat("\nWrote manuscript/figures/F17_eqinit_rates.png\n\n")
w <- reshape(stats[, c("interval", "model", "arm", "med")], idvar = c("interval", "model"),
             timevar = "arm", direction = "wide")
w$obs <- z$obs$rates$obs[match(w$interval, z$obs$rates$interval)]
print(transform(w, med.prod = round(med.prod, 3), med.eq = round(med.eq, 3), obs = round(obs, 3)),
      row.names = FALSE)
