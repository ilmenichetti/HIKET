# =============================================================================
# build_F15_forward_arms.R   (2026-09-02)
#
# F15 -- the consequence figure. Two parameterisations of the same model that
# agree on everything observable, asked what the soil does under warming.
#
#   (a) trajectories, climate held stationary   [3 models x 2 arms]
#   (b) trajectories, warmed to SSP2-4.5        [same design, shared y-axis]
#   (c) the warming response, by arm            [boxplots]
#   (d) the difference between arms             [boxplots, zero line]
#
# THE RESULT IS STRUCTURAL, and the wording matters: Yasso15/20 prove RESILIENT
# to the input-anchor choice, not "right". Their pool-specific climate modifiers
# leave the carbon-holding pools (H, N: ~60% of the stock) responding to warming
# almost identically in both arms, while the arms differ mainly in the AWE
# modifier -- and AWE holds ~13% of the stock. Yasso07 has ONE modifier for every
# pool, so displacing it displaces the sensitivity of all the carbon at once.
#
# ⚠ PANEL TITLES DESCRIBE, THEY DO NOT JUDGE. An earlier version titled panel (a)
# "indistinguishable where the data are", which asserts the reading rather than
# letting the reader make it. Titles now say what is plotted; the interpretation
# belongs in the caption and the text.
#
# SCENARIO: linear ramp to SSP2-4.5 (Ruosteenoja & Jylha 2021, Geophysica
# 56(1-2), 39-69; CMIP6, Finland, annual mean), baseline-adjusted for our
# 2005-2024 recycled climate. All scenarios stored; only the mid one is plotted.
# =============================================================================

setwd("/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling")
source("manuscript/figures/model_palette.R")
source("manuscript/figures/obs_basis.R")
R <- readRDS("doublechecks/f15_forward_experiment.rds")
MODELS <- c("Yasso07","Yasso15","Yasso20")
PLOT_S <- "ssp245"
# ⚠ lty=2 read as a near-solid line at this line width. A tight explicit dash
# pattern ("22" = 2 on, 2 off) is legible against the solid arm.
ARM_LTY <- c(published = "22", ours = "solid")

resp <- function(M, arm, sc = PLOT_S) {
  z <- R$out[[M]][[arm]]; if (!length(z)) return(numeric(0))
  vapply(z, function(d) { c0 <- tail(d[["traj"]][["control"]], 1)
                          cw <- tail(d[["traj"]][[sc]], 1); 100*(cw/c0 - 1) }, numeric(1))
}
dif <- function(M, sc = PLOT_S) {
  o <- resp(M,"ours",sc); p <- resp(M,"published",sc)
  if (!length(o) || !length(p)) return(numeric(0))
  v <- as.vector(outer(o, p, "-")); if (length(v) > 6000) sample(v, 6000) else v
}
med <- function(M, arm, sc) {
  z <- R$out[[M]][[arm]]
  apply(do.call(cbind, lapply(z, function(d) d[["traj"]][[sc]])), 1, median)
}
bx <- function(v, at, col, w) {
  if (!length(v)) return(invisible())
  if (length(v) == 1) { segments(at-w, v, at+w, v, col=col, lwd=3.2); return(invisible()) }
  b <- boxplot(v, plot = FALSE)
  rect(at-w, b$stats[2], at+w, b$stats[4], col=adjustcolor(col,.35), border=col, lwd=1.6)
  segments(at-w, b$stats[3], at+w, b$stats[3], col=col, lwd=3)
  segments(at, b$stats[1], at, b$stats[2], col=col, lwd=1.4)
  segments(at, b$stats[4], at, b$stats[5], col=col, lwd=1.4)
}

yrs <- 1985:(1985 + R$n_proj + 39)
om  <- readRDS(sprintf("Data/model_inputs/Yasso07_inputs_%s.rds", R$rid[["Yasso07"]]))$obs_meta
oc  <- obs_campaigns(om, balanced_plots(om))
YL  <- range(unlist(lapply(MODELS, function(M) lapply(c("published","ours"),
        function(a) c(med(M,a,"control"), med(M,a,PLOT_S))))), oc$lo, oc$hi) + c(-0.8, 0.8)

png("manuscript/figures/F15_forward_arms.png", width = 11.2, height = 8.6,
    units = "in", res = 220)
par(mfrow = c(2,2), mar = c(4.1, 4.6, 3.1, 0.9), mgp = c(2.7, 0.7, 0), las = 1,
    cex.axis = 0.95)

traj_panel <- function(sc, tag, title, sub) {
  plot(NA, xlim=range(yrs), ylim=YL, xlab="", ylab=expression("SOC (tC ha"^-1*")"))
  usr <- par("usr"); rect(2024.5, usr[3], usr[2], usr[4], col="#F6F5F1", border=NA); box()
  abline(v = 2024.5, col="grey45", lwd=1.2, lty=3)
  for (M in MODELS) for (arm in c("published","ours"))
    lines(yrs, med(M, arm, sc), col=MODEL_COL[M], lwd=2.4, lty=ARM_LTY[arm])
  arrows(oc$year, oc$lo, oc$year, oc$hi, angle=90, code=3, length=.03, lwd=1.7, col="#4A1A0E")
  points(oc$year, oc$m, pch=21, bg="#C85A3C", col="#4A1A0E", cex=1.55)
  text(2004, usr[4]-1.4, "observed", cex=.8, col="grey35")
  text(2056, usr[4]-1.4, "projection", cex=.8, col="grey35")
  mtext(sprintf("%s %s", tag, title), 3, line=1.0, adj=0, font=2, cex=.95)
  mtext(sub, 3, line=0.0, adj=0, cex=.72, col="grey35")
}

# (a) stationary climate
traj_panel("control", "(a)", "modelled SOC, climate held stationary",
  "posterior medians; climate recycled from 2005-2024, no further warming")
legend("bottomright", bty="n", cex=.85, ncol=2, seg.len=2.8,
  legend=c(MODELS, "observed mean", "our calibration", "published parameters"),
  col=c(MODEL_COL[MODELS], "#4A1A0E", "grey30", "grey30"),
  lty=c("solid","solid","solid",NA,"solid","22"), lwd=c(2.4,2.4,2.4,NA,2.4,2.4),
  pch=c(NA,NA,NA,21,NA,NA), pt.bg="#C85A3C")

# (b) warmed
traj_panel(PLOT_S, "(b)", "modelled SOC, warming to SSP2-4.5",
  "same models and axes; projection warmed on a linear ramp to +2.6 °C at 2084")

# (c) response
par(mar = c(4.1, 4.6, 3.1, 0.9))
rr <- unlist(lapply(MODELS, function(M) lapply(c("published","ours"), function(a) resp(M,a))))
plot(NA, xlim=c(.5,3.5), ylim=range(rr)+c(-.7,.7), xaxt="n", xlab="",
     ylab="SOC change under SSP2-4.5 (%)")
axis(1, at=1:3, labels=MODELS); abline(h=0, col="grey65")
for (i in seq_along(MODELS)) {
  bx(resp(MODELS[i],"published"), i-0.22, "#4C72A8",           w=0.16)
  bx(resp(MODELS[i],"ours"),      i+0.22, MODEL_COL[MODELS[i]], w=0.16)
}
legend("bottomleft", bty="n", cex=.83,
  fill=c(adjustcolor("#4C72A8",.35), adjustcolor(MODEL_COL[["Yasso15"]],.35)),
  border=c("#4C72A8", MODEL_COL[["Yasso15"]]),
  legend=c("published parameters", "our calibration (model colour)"))
mtext("(c) warming response at 2084, by arm", 3, line=1.0, adj=0, font=2, cex=.95)
mtext("boxes = interquartile range of 200 draws; Yasso07 published = point",
      3, line=0.0, adj=0, cex=.72, col="grey35")

# (d) difference
par(mar = c(4.1, 4.9, 3.1, 0.9))
dd <- lapply(MODELS, dif)
plot(NA, xlim=c(.5,3.5), ylim=range(unlist(dd))+c(-.4,.4), xaxt="n", xlab="",
     ylab="difference, ours − published (pp)")
axis(1, at=1:3, labels=MODELS); abline(h=0, col="grey25", lwd=1.5)
for (i in seq_along(MODELS)) bx(dd[[i]], i, MODEL_COL[MODELS[i]], w=0.28)
mtext("(d) difference between the two arms", 3, line=1.0, adj=0, font=2, cex=.95)
mtext("boxes crossing zero: the input anchor does not change the response",
      3, line=0.0, adj=0, cex=.72, col="grey35")
dev.off()
cat("wrote manuscript/figures/F15_forward_arms.png\n\n")

cat("SOC change at 2084 (%), median [90%] -- ALL SCENARIOS\n")
for (sc in setdiff(names(R$scen), "control")) {
  cat(sprintf("\n-- %s (ramp to +%.1f C) --\n", R$scen_lab[[sc]], R$scen[[sc]]))
  for (M in MODELS) {
    for (arm in c("published","ours")) { v <- resp(M, arm, sc); if (!length(v)) next
      cat(sprintf("   %-8s %-10s %6.2f [%6.2f, %6.2f]\n", M, arm, median(v),
                  quantile(v,.05), quantile(v,.95))) }
    d <- dif(M, sc)
    if (length(d)) cat(sprintf("   %-8s %-10s %6.2f [%6.2f, %6.2f] pp%s\n\n", M, "DIFF",
        median(d), quantile(d,.05), quantile(d,.95),
        if (quantile(d,.05) > 0 || quantile(d,.95) < 0) "   EXCLUDES ZERO" else ""))
  }
}
