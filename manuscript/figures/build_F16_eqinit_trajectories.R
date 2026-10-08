# =============================================================================
# build_F16_eqinit_trajectories.R   (2026-10-08)
#
# F16 -- THE COUNTERFACTUAL: transient vs equilibrium start, 1985-2084.
# Design: NEXT_RUN_equilibrium_init.md. Two outputs, same data:
#   F16_eqinit_trajectories.png      ONE panel, the ensemble (decided layout):
#                                    six thin lines per arm (model colour; solid =
#                                    transient, dashed = equilibrium), the ensemble
#                                    median bold, observed campaign means with 95% CI.
#                                    Lower strip: effective litter input per arm.
#   F16_eqinit_trajectories_six.png  fallback, one panel per model with 90% bands.
#
# Read with care: the projection (shaded, after 2024) holds litter at its 2024 value
# and recycles the 2005-2024 climate -- no climate change. Its slope is therefore
# each arm's remaining disequilibrium, which is exactly what differs between them.
# Ensemble = mean over the available models of each model's draw-median trajectory.
#
# Usage:  Rscript manuscript/figures/build_F16_eqinit_trajectories.R
# =============================================================================

setwd("/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling")
source("manuscript/figures/eqinit_common.R")
z  <- load_eqinit_comparison()
MS <- z$models
yr <- as.integer(colnames(z$arms[[MS[1]]]$prod$traj))
Y_PROJ <- 2024

med_traj <- function(m, arm) apply(z$arms[[m]][[arm]]$traj, 2, median)
band     <- function(m, arm) apply(z$arms[[m]][[arm]]$traj, 2, quantile, c(.05, .95))
ens      <- function(arm) rowMeans(sapply(MS, med_traj, arm = arm))

obs <- z$obs$camp
proj_shade <- function(ylim) {
  rect(Y_PROJ, ylim[1] - 1e3, 2100, ylim[2] + 1e3, col = "grey94", border = NA)
  abline(v = Y_PROJ, col = "grey55", lty = 3)
}
draw_obs <- function() {
  arrows(obs$year, obs$lo, obs$year, obs$hi, angle = 90, code = 3, length = 0.035,
         col = CAMPAIGN_COL, lwd = 1.6)
  points(obs$year, obs$m, pch = 21, bg = CAMPAIGN_COL, col = "grey15", cex = 1.4)
}
foot <- function() {
  rd <- function(r) paste(unique(substr(r, 1, 8)), collapse = "/")      # run dates; full ids in the cache
  mtext(sprintf("%s; transient = run %s, equilibrium = run %s%s", basis_note(seq_len(z$obs$n)),
                rd(z$stamp$prod), rd(z$stamp$eq),
                if (length(z$missing)) paste0("; missing: ", paste(z$missing, collapse = ", ")) else ""),
        side = 1, outer = TRUE, line = 0.4, cex = 0.62, col = "grey35", adj = 0)
}

# --- (1) ensemble, one panel ----------------------------------------------------
png("manuscript/figures/F16_eqinit_trajectories.png", width = 9, height = 7.2, units = "in", res = 200)
layout(matrix(1:2, 2), heights = c(3.1, 1.15))
par(oma = c(1.6, 0, 1.8, 0), mgp = c(2.5, 0.65, 0), las = 1)

all_med <- unlist(lapply(MS, function(m) c(med_traj(m, "prod"), med_traj(m, "eq"))))
ylim <- range(c(all_med, obs$lo, obs$hi)) + c(-2, 2)
par(mar = c(1.2, 4.6, 0.8, 1.2))
plot(NA, xlim = range(yr), ylim = ylim, xlab = "", xaxt = "n",
     ylab = expression("Soil carbon  (tC ha"^-1*")"))
proj_shade(ylim); axis(1, labels = FALSE)
for (m in MS) for (a in c("prod", "eq"))
  lines(yr, med_traj(m, a), col = adjustcolor(MODEL_COL[[m]], 0.85), lwd = 1.3, lty = ARM_LTY[[a]])
for (a in c("prod", "eq")) lines(yr, ens(a), col = "grey10", lwd = 3.2, lty = ARM_LTY[[a]])
draw_obs()
text(max(yr), ylim[1] + 0.3, "projection: litter held at 2024,\nrecent climate recycled",
     adj = c(1, 0), cex = 0.72, col = "grey35", font = 3)
legend("topleft", bty = "n", cex = 0.78, ncol = 2,
       legend = c(MS, "ensemble, transient start", "ensemble, equilibrium start", "observed (95% CI)"),
       col = c(MODEL_COL[MS], "grey10", "grey10", "grey15"),
       lty = c(rep(1, length(MS)), 1, 2, NA), lwd = c(rep(1.6, length(MS)), 3.2, 3.2, NA),
       pch = c(rep(NA, length(MS)), NA, NA, 21), pt.bg = c(rep(NA, length(MS) + 2), CAMPAIGN_COL[2]))
mtext("solid: transient start (production)    dashed: equilibrium start in 1985",
      side = 3, line = -1.1, adj = 0.99, cex = 0.72, col = "grey25")

# lower strip: effective litter input sigma_input x J, history then held at 2024
par(mar = c(3.2, 4.6, 0.4, 1.2))
flux_line <- function(m, a) {
  f <- z$arms[[m]][[a]]$flux
  x <- c(f$year, max(f$year) + seq_len(max(yr) - max(f$year)))
  list(x = x, y = c(f$med, rep(tail(f$med, 1), length(x) - nrow(f))))
}
fl <- unlist(lapply(MS, function(m) c(flux_line(m, "prod")$y, flux_line(m, "eq")$y)))
yl2 <- range(fl) * c(0.95, 1.05)
plot(NA, xlim = range(yr), ylim = yl2, xlab = "", ylab = "")
mtext(expression("Litter input (tC ha"^-1*" yr"^-1*")"), side = 2, line = 2.6, cex = 0.72, las = 0)
mtext("Year", side = 1, line = 2.2)
proj_shade(yl2)
for (m in MS) for (a in c("prod", "eq")) {
  l <- flux_line(m, a); lines(l$x, l$y, col = MODEL_COL[[m]], lwd = 1.3, lty = ARM_LTY[[a]]) }
mtext("Soil carbon under a transient vs an equilibrium start", outer = TRUE, line = 0.4, cex = 1.1, font = 2)
foot()
dev.off()

# --- (2) six panels -------------------------------------------------------------
png("manuscript/figures/F16_eqinit_trajectories_six.png", width = 11, height = 7, units = "in", res = 200)
par(mfrow = c(2, 3), mar = c(3.2, 4.2, 2.2, 0.8), oma = c(1.6, 0, 2, 0), mgp = c(2.3, 0.6, 0), las = 1)
yl <- range(c(unlist(lapply(MS, function(m) c(band(m, "prod"), band(m, "eq")))), obs$lo, obs$hi))
for (m in FIG_MODELS) {
  plot(NA, xlim = range(yr), ylim = yl, xlab = "Year", ylab = expression("SOC (tC ha"^-1*")"), main = m)
  if (!m %in% MS) { text(mean(range(yr)), mean(yl), "equilibrium arm\nnot available", col = "grey40"); next }
  proj_shade(yl)
  for (a in c("prod", "eq")) {
    b <- band(m, a)
    polygon(c(yr, rev(yr)), c(b[1, ], rev(b[2, ])), border = NA,
            col = adjustcolor(if (a == "prod") MODEL_COL[[m]] else "grey45", 0.25))
    lines(yr, med_traj(m, a), col = if (a == "prod") MODEL_COL[[m]] else "grey25", lwd = 2, lty = ARM_LTY[[a]])
  }
  draw_obs()
  if (m == FIG_MODELS[1])
    legend("topleft", bty = "n", cex = 0.8, legend = ARM_LAB, lty = ARM_LTY, lwd = 2,
           col = c(MODEL_COL[[m]], "grey25"))
}
mtext("Transient vs equilibrium start, per model (median and 90% posterior band of the national mean)",
      outer = TRUE, line = 0.4, cex = 1.0, font = 2)
foot()
dev.off()

# --- numbers ------------------------------------------------------------------
tab <- do.call(rbind, lapply(MS, function(m) do.call(rbind, lapply(c("prod", "eq"), function(a) {
  t <- z$arms[[m]][[a]]$traj
  data.frame(model = m, arm = a, C_1985 = median(t[, "1985"]), C_2024 = median(t[, "2024"]),
             C_2084 = median(t[, "2084"]))
}))))
cat("\nWrote F16_eqinit_trajectories.png and F16_eqinit_trajectories_six.png\n\n")
print(transform(tab, C_1985 = round(C_1985, 1), C_2024 = round(C_2024, 1), C_2084 = round(C_2084, 1)),
      row.names = FALSE)
