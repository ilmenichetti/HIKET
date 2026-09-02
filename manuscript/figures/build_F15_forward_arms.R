# =============================================================================
# build_F15_forward_arms.R   (2026-09-02)
#
# F15 -- the consequence figure. Two parameterisations of the same model that
# agree on everything observable, asked what the soil does under warming.
#
# Read left to right as an argument, because a divergence in projection is
# otherwise dismissed as "their models just differ":
#   (a) they agree on everything we can observe   [Yasso07, 1985-2084]
#   (b) they diverge once warming is applied      [response at 2084, all three]
#   (c) the SIZE of the disagreement, and it grows [ours - published, zero line]
#
# Yasso20 is the NEGATIVE CONTROL: its published parameterisation already sits
# where our calibration puts it (MRT x1.02), so it must not separate.
# Yasso07 is the positive control (MRT x0.73) and must.
#
# Data: doublechecks/f15_forward_experiment.rds
# Units: tC/ha AND % -- the absolute response carries kinetics + input, the
# relative response is EXACTLY sigma_input-invariant and isolates kinetics.
# =============================================================================

setwd("/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling")
source("manuscript/figures/model_palette.R")
source("manuscript/figures/obs_basis.R")
R <- readRDS("doublechecks/f15_forward_experiment.rds")
MODELS <- c("Yasso07","Yasso15","Yasso20")

ARM_COL <- c(published = "#4C72A8", ours = "#C26B51")
ARM_LAB <- c(published = "published parameters", ours = "our calibration")

# --- extract: response at the final projection year, per draw ----------------
resp <- function(M, arm, dT, rel = TRUE) {
  z <- R$out[[M]][[arm]]; if (!length(z)) return(numeric(0))
  vapply(z, function(d) {
    c0 <- tail(d$traj[["dT0"]], 1); cw <- tail(d$traj[[paste0("dT", dT)]], 1)
    if (rel) 100*(cw/c0 - 1) else cw - c0
  }, numeric(1))
}
mrt_of <- function(M, arm) vapply(R$out[[M]][[arm]], function(d) d$mrt, numeric(1))

png("manuscript/figures/F15_forward_arms.png", width = 13, height = 4.6,
    units = "in", res = 220)
layout(matrix(1:3, nrow = 1), widths = c(1.15, 1, 1))
par(mar = c(4.2, 4.4, 3.0, 1.0), mgp = c(2.6, 0.7, 0), las = 1, cex.axis = 0.95)

# ============================== PANEL (a) ====================================
M <- "Yasso07"
yrs <- 1985:(1985 + R$n_proj + 39)
med <- function(arm, dT) {
  z <- R$out[[M]][[arm]]
  apply(do.call(cbind, lapply(z, function(d) d$traj[[paste0("dT", dT)]])), 1, median)
}
a_ctrl <- lapply(c("published","ours"), med, dT = 0); names(a_ctrl) <- c("published","ours")
a_warm <- lapply(c("published","ours"), med, dT = 5); names(a_warm) <- c("published","ours")
yl <- range(unlist(c(a_ctrl, a_warm)), 60, 90)
plot(NA, xlim = range(yrs), ylim = yl, xlab = "", ylab = expression("SOC (tC ha"^-1*")"))
rect(2024.5, yl[1]-10, max(yrs)+5, yl[2]+10, col = "#F5F5F2", border = NA)
abline(v = 2024.5, col = "grey45", lwd = 1.2, lty = 2)
for (arm in c("published","ours")) {
  lines(yrs, a_ctrl[[arm]], col = ARM_COL[arm], lwd = 2.6)
  lines(yrs, a_warm[[arm]], col = ARM_COL[arm], lwd = 2.0, lty = 3)
}
om <- readRDS(sprintf("Data/model_inputs/%s_inputs_%s.rds", M, R$rid[[M]]))$obs_meta
oc <- obs_campaigns(om, balanced_plots(om))
arrows(oc$year, oc$lo, oc$year, oc$hi, angle = 90, code = 3, length = 0.03, lwd = 1.6, col = "#6E2414")
points(oc$year, oc$m, pch = 21, bg = "#C85A3C", col = "#6E2414", cex = 1.5)
text(2005, yl[2]-1, "observed window", cex = 0.85, col = "grey30")
text(2056, yl[2]-1, "projection", cex = 0.85, col = "grey30")
legend("bottomleft", bty = "n", cex = 0.85,
       legend = c(ARM_LAB, "control", "+5 °C", "observed"),
       col = c(ARM_COL, "grey30", "grey30", "#6E2414"),
       lty = c(1,1,1,3,NA), lwd = c(2.6,2.6,2,2,NA), pch = c(NA,NA,NA,NA,21), pt.bg = "#C85A3C")
mtext("(a) indistinguishable where the data are", side = 3, line = 0.9, adj = 0, font = 2, cex = 0.95)
d_pre <- max(abs(a_ctrl$published[1:40] - a_ctrl$ours[1:40]))
mtext(sprintf("Yasso07  ·  max arm gap before 2024: %.1f tC/ha  ·  misfit 0.020 vs 0.029",
              d_pre), side = 3, line = -0.1, adj = 0, cex = 0.72, col = "grey35")

# ============================== PANEL (b) ====================================
par(mar = c(4.2, 4.4, 3.0, 1.0))
xs <- seq_along(MODELS)
plot(NA, xlim = c(0.5, 3.5), ylim = c(-30, 0), xaxt = "n",
     xlab = "", ylab = "SOC change under warming (%)")
axis(1, at = xs, labels = MODELS)
abline(h = 0, col = "grey60")
off <- c(published = -0.13, ours = 0.13)
for (i in seq_along(MODELS)) for (arm in c("published","ours")) for (dT in c(2,5)) {
  v <- resp(MODELS[i], arm, dT); if (!length(v)) next
  x <- i + off[arm] + ifelse(dT == 2, -0.045, 0.045)
  q <- quantile(v, c(.05,.5,.95))
  segments(x, q[1], x, q[3], col = ARM_COL[arm], lwd = ifelse(dT==2, 2.0, 3.4))
  points(x, q[2], pch = ifelse(dT==2, 21, 23), bg = ARM_COL[arm], col = "white", cex = 1.25)
}
legend("bottomleft", bty = "n", cex = 0.85, legend = c(ARM_LAB, "+2 °C", "+5 °C"),
       col = c(ARM_COL, "grey30","grey30"), lwd = c(3,3,2,3.4),
       pch = c(NA,NA,21,23), pt.bg = "grey30")
mtext("(b) and they diverge once you ask about warming", side = 3, line = 0.9, adj = 0, font = 2, cex = 0.95)
mtext("median, 90% of draws  ·  Yasso07 published = point (no published posterior)",
      side = 3, line = -0.1, adj = 0, cex = 0.72, col = "grey35")

# ============================== PANEL (c) ====================================
par(mar = c(4.2, 4.6, 3.0, 1.0))
dif <- function(M, dT) {
  o <- resp(M, "ours", dT); p <- resp(M, "published", dT)
  if (!length(o) || !length(p)) return(c(NA,NA,NA))
  d <- outer(o, p, "-")
  quantile(as.vector(d), c(.05,.5,.95))
}
plot(NA, xlim = c(0.5, 3.5), ylim = c(-12, 5), xaxt = "n",
     xlab = "", ylab = "difference, ours − published (pp)")
axis(1, at = xs, labels = MODELS)
abline(h = 0, col = "grey30", lwd = 1.4)
for (i in seq_along(MODELS)) for (dT in c(2,5)) {
  q <- dif(MODELS[i], dT); if (any(is.na(q))) next
  x <- i + ifelse(dT == 2, -0.12, 0.12)
  segments(x, q[1], x, q[3], col = "grey25", lwd = ifelse(dT==2, 2.2, 3.6))
  points(x, q[2], pch = ifelse(dT==2, 21, 23), bg = MODEL_COL[MODELS[i]], col = "white", cex = 1.6)
}
mtext("(c) the size of the disagreement", side = 3, line = 0.9, adj = 0, font = 2, cex = 0.95)
mtext("intervals crossing zero = arms indistinguishable", side = 3, line = -0.1,
      adj = 0, cex = 0.72, col = "grey35")
text(3.0, -9.0, "Yasso20: no arm\ndisplacement\n(control)", cex = 0.78, col = "grey30", font = 3)
text(1.0, -10.8, "only Yasso07\nexcludes zero", cex = 0.78, col = "grey20", font = 2)
dev.off()
cat("wrote manuscript/figures/F15_forward_arms.png\n")

# --- numbers for the caption / text -----------------------------------------
cat("\nresponse at 2084 (%), median [90%]:\n")
for (M in MODELS) for (dT in c(2,5)) {
  for (arm in c("published","ours")) { v <- resp(M, arm, dT)
    if (length(v)) cat(sprintf("  %-8s +%dC %-10s %6.2f [%6.2f, %6.2f]  n=%d\n",
                               M, dT, arm, median(v), quantile(v,.05), quantile(v,.95), length(v))) }
  q <- dif(M, dT); cat(sprintf("  %-8s +%dC DIFFERENCE %6.2f [%6.2f, %6.2f] pp\n\n", M, dT, q[2], q[1], q[3]))
}
