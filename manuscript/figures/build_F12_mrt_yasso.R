suppressWarnings(suppressMessages(library(BayesianTools)))
set.seed(2025)
source("manuscript/figures/model_palette.R")

# Intrinsic model MRT: computed by doublechecks/intrinsic_mrt.R, which feeds each
# model a UNIT litter input at a fixed reference condition (dataset-mean climate
# and dataset-mean AWEN x size composition) and takes sum(C_ss) from the PURE
# steady-state routine. No sigma_input, no sigma_init, no SOC observations, no
# pre-run -- so this is a property of the model's generator, identical method for
# our posterior and for the published MCMC samples.
#
# BUT: MRT as an INFERRED quantity IS conditional on the sigma_input prior. The
# likelihood mainly pins the PRODUCT MRT x sigma_input x J_raw (= observed stock),
# so the sigma_input prior decides how that product splits. Yasso15: published
# kinetics need sigma_input 1.305, ours 2.432 -- ratios inverse to ~8%. So
# 'MRT too short' and 'sigma_input too high' are ONE finding from two ends.
# Read this figure against the sigma_input posterior, and say so in the caption.
RES <- readRDS("doublechecks/intrinsic_mrt.rds")
MODELS <- c("Yasso07","Yasso15","Yasso20")
R <- RES$out; ref <- RES$ref

# --- sanity: MRT must be invariant to sigma_input (linear system) -------------
cat("\n--- intrinsic MRT (yr) ---\n")
for (m in MODELS) cat(sprintf("%-8s published %6.2f%s | ours %6.2f [%.2f, %.2f]\n",
  m, if (length(R[[m]]$pub) > 1) median(R[[m]]$pub) else R[[m]]$pub_pt,
  if (length(R[[m]]$pub) > 1) " (posterior)" else " (point)    ",
  median(R[[m]]$ours), quantile(R[[m]]$ours,.025), quantile(R[[m]]$ours,.975)))

# =============================== figure ======================================
png("manuscript/figures/F12_mrt_yasso.png", width = 7.2, height = 5.4, units = "in", res = 300)
par(mar = c(4.6, 4.6, 3.2, 1.2))

pub_of <- function(m) if (length(R[[m]]$pub) > 1) R[[m]]$pub else R[[m]]$pub_pt
allv <- unlist(c(lapply(MODELS, function(m) R[[m]]$ours), lapply(MODELS, pub_of)))
yl <- range(pretty(c(min(allv) * 0.90, max(allv) * 1.06)))

plot(NA, xlim = c(0.5, length(MODELS) + 0.5), ylim = yl, axes = FALSE, xlab = "", ylab = "")
abline(h = axTicks(2), col = "grey92", lwd = 0.8)

# GREY, wide, behind: the published parameterisation. A box where a published
# MCMC sample exists (Yasso15/20); a line where only a point value does (Yasso07).
for (i in seq_along(MODELS)) {
  v <- pub_of(MODELS[i])
  if (length(v) > 1) {
    boxplot(v, at = i, add = TRUE, axes = FALSE, outline = FALSE, boxwex = 0.74,
            whisklty = 1, staplewex = 0.35, col = "grey86", border = "grey45", medlwd = 2.6)
  } else {
    segments(i - 0.37, v, i + 0.37, v, col = "grey35", lwd = 3)
  }
}

# COLOURED, narrow, in front: this study's posterior
boxplot(lapply(MODELS, function(m) R[[m]]$ours), at = seq_along(MODELS), add = TRUE,
        axes = FALSE, outline = FALSE, boxwex = 0.34, whisklty = 1, staplewex = 0.28,
        col = adjustcolor(MODEL_COL[MODELS], alpha.f = 0.9),
        border = MODEL_COL[MODELS], medlwd = 2.4)

axis(1, at = seq_along(MODELS), labels = MODELS, tick = FALSE, line = -0.4)
for (i in seq_along(MODELS)) {
  v <- pub_of(MODELS[i]); ref_med <- if (length(v) > 1) median(v) else v
  mtext(sprintf("%.2f×", median(R[[MODELS[i]]]$ours) / ref_med),
        side = 1, at = i, line = 1.15, cex = 0.86, col = "grey25", font = 2)
}
axis(2, las = 1)
mtext("Intrinsic mean residence time (yr)", side = 2, line = 3.1, cex = 1.02)
box(col = "grey55")

title(main = "Intrinsic mean residence time, Yasso family", cex.main = 1.04, line = 1.9)
mtext(sprintf("unit litter input at a fixed reference (T = %.1f\u00b0C, P = %.0f mm); independent of \u03c3_input and \u03c3_init",
              ref$clim$temp_mean, ref$clim$precip),
      side = 3, line = 0.55, cex = 0.76, col = "grey30")

legend("topright", inset = c(0.01, 0.02), bty = "n", cex = 0.82,
       legend = c("published parameterisation", "this study (posterior)"),
       fill = c("grey86", adjustcolor(MODEL_COL[["Yasso15"]], alpha.f = 0.9)),
       border = c("grey45", MODEL_COL[["Yasso15"]]))
mtext(expression("ratio = our median " %/% " published median"),
      side = 1, line = 2.4, cex = 0.76, col = "grey40")
mtext("Yasso07: line, not box \u2014 no published MCMC sample available",
      side = 1, line = 3.3, cex = 0.72, col = "grey45")

dev.off()
cat("\nwrote manuscript/figures/F12_mrt_yasso.png\n")
