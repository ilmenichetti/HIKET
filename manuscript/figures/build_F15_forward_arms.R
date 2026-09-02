# =============================================================================
# build_F15_forward_arms.R   (2026-09-02)
#
# F15 -- the consequence figure. Two parameterisations of the same model that
# agree on everything observable, asked what the soil does under warming.
#
#   (a) they agree where the data are        [all three Yassos, 1985-2084]
#   (b) the warming response, by arm         [boxplots]
#   (c) the SIZE of the disagreement         [boxplots, zero line]
#
# THE RESULT IS STRUCTURAL, and the wording matters: Yasso15/20 prove RESILIENT
# to the input-anchor choice, not "right". Their pool-specific climate modifiers
# mean the carbon-holding pools (H, N: ~60% of the stock) respond to warming
# almost identically in both arms, while the arms differ mainly in the AWE
# modifier -- and AWE holds ~13% of the stock. Yasso07 has ONE modifier for every
# pool, so displacing it displaces the sensitivity of all the carbon at once.
#
# SCENARIO: linear ramp to the SSP2-4.5 level (Ruosteenoja & Jylha 2021,
# Geophysica 56(1-2), 39-69; CMIP6, Finland, annual mean), baseline-adjusted for
# our 2005-2024 recycled climate. All scenarios are computed and stored; only the
# MID one is plotted -- the others are in the .rds and the appendix.
# =============================================================================

setwd("/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling")
source("manuscript/figures/model_palette.R")
source("manuscript/figures/obs_basis.R")
R <- readRDS("doublechecks/f15_forward_experiment.rds")
MODELS  <- c("Yasso07","Yasso15","Yasso20")
PLOT_S  <- "ssp245"                       # the one scenario shown in the MS figure
ARM_LTY <- c(published = 2, ours = 1)

resp <- function(M, arm, sc = PLOT_S) {
  z <- R$out[[M]][[arm]]; if (!length(z)) return(numeric(0))
  vapply(z, function(d) { c0 <- tail(d[["traj"]][["control"]], 1)
                          cw <- tail(d[["traj"]][[sc]], 1); 100*(cw/c0 - 1) }, numeric(1))
}
dif <- function(M, sc = PLOT_S) {
  o <- resp(M, "ours", sc); p <- resp(M, "published", sc)
  if (!length(o) || !length(p)) return(numeric(0))
  v <- as.vector(outer(o, p, "-")); if (length(v) > 6000) v <- sample(v, 6000) else v
}
bx <- function(v, at, col, w = 0.26) {
  if (!length(v)) return(invisible())
  if (length(v) == 1) { segments(at-w, v, at+w, v, col = col, lwd = 3); return(invisible()) }
  b <- boxplot(v, plot = FALSE)
  rect(at-w, b$stats[2], at+w, b$stats[4], col = adjustcolor(col, .35), border = col, lwd = 1.6)
  segments(at-w, b$stats[3], at+w, b$stats[3], col = col, lwd = 3)
  segments(at, b$stats[1], at, b$stats[2], col = col, lwd = 1.4)
  segments(at, b$stats[4], at, b$stats[5], col = col, lwd = 1.4)
}

png("manuscript/figures/F15_forward_arms.png", width = 13.4, height = 4.7,
    units = "in", res = 220)
layout(matrix(1:3, nrow = 1), widths = c(1.25, 1, 1))
par(mar = c(4.2, 4.5, 3.2, 0.8), mgp = c(2.6, 0.7, 0), las = 1, cex.axis = 0.95)

# ============================== PANEL (a) ====================================
yrs <- 1985:(1985 + R$n_proj + 39)
med <- function(M, arm, sc) {
  z <- R$out[[M]][[arm]]
  apply(do.call(cbind, lapply(z, function(d) d[["traj"]][[sc]])), 1, median)
}
# ⚠ CONTROL ONLY in (a). Adding the scenario lines makes 12 curves and the panel
# stops being readable. (a)'s job is LEVELS -- the arms agree where the data are,
# and drift apart afterwards; (b)/(c) carry the warming RESPONSE, which is a
# difference and nets the baseline out. Scenario trajectories go to the appendix.
allv <- unlist(lapply(MODELS, function(M) lapply(c("published","ours"),
          function(a) med(M, a, "control"))))
plot(NA, xlim = range(yrs), ylim = range(allv, 58, 82),
     xlab = "", ylab = expression("SOC (tC ha"^-1*")"))
usr <- par("usr"); rect(2024.5, usr[3], usr[2], usr[4], col = "#F6F5F1", border = NA); box()
abline(v = 2024.5, col = "grey45", lwd = 1.2, lty = 2)
for (M in MODELS) for (arm in c("published","ours"))
  lines(yrs, med(M, arm, "control"), col = MODEL_COL[M], lwd = 2.4, lty = ARM_LTY[arm])
om <- readRDS(sprintf("Data/model_inputs/Yasso07_inputs_%s.rds", R$rid[["Yasso07"]]))$obs_meta
oc <- obs_campaigns(om, balanced_plots(om))
arrows(oc$year, oc$lo, oc$year, oc$hi, angle=90, code=3, length=.03, lwd=1.7, col="#4A1A0E")
points(oc$year, oc$m, pch = 21, bg = "#C85A3C", col = "#4A1A0E", cex = 1.6)
text(2004, usr[4]-1.2, "observed window", cex=.82, col="grey30")
text(2056, usr[4]-1.2, "projection", cex=.82, col="grey30")
legend("bottomright", bty="n", cex=.85, ncol=2, seg.len=2.6,
  legend=c(MODELS, "observed", "our calibration", "published params"),
  col=c(MODEL_COL[MODELS], "#4A1A0E", "grey30", "grey30"),
  lty=c(1,1,1,NA,1,2), lwd=c(2.4,2.4,2.4,NA,2.4,2.4),
  pch=c(NA,NA,NA,21,NA,NA), pt.bg="#C85A3C")
mtext("(a) indistinguishable where the data are", 3, line=1.0, adj=0, font=2, cex=.95)
gap <- max(vapply(MODELS, function(M) max(abs(med(M,"published","control")[1:40] -
                                             med(M,"ours","control")[1:40])), numeric(1)))
mtext(sprintf("posterior medians, no further warming  ·  largest arm gap before 2024: %.1f tC/ha", gap),
      3, line=0.0, adj=0, cex=.72, col="grey35")

# ============================== PANEL (b) ====================================
par(mar = c(4.2, 4.5, 3.2, 0.8))
plot(NA, xlim=c(.5,3.5), ylim=range(unlist(lapply(MODELS,function(M)
       lapply(c("published","ours"), function(a) resp(M,a))))) + c(-.6,.6),
     xaxt="n", xlab="", ylab="SOC change under SSP2-4.5 (%)")
axis(1, at=1:3, labels=MODELS); abline(h=0, col="grey65")
for (i in seq_along(MODELS)) {
  bx(resp(MODELS[i],"published"), i-0.19, "#4C72A8")
  bx(resp(MODELS[i],"ours"),      i+0.19, MODEL_COL[MODELS[i]])
}
legend("bottomleft", bty="n", cex=.82,
       fill=c(adjustcolor("#4C72A8",.35), adjustcolor(MODEL_COL[["Yasso15"]],.35)),
       border=c("#4C72A8", MODEL_COL[["Yasso15"]]),
       legend=c("published parameters", "our calibration (model colour)"))
mtext("(b) the warming response", 3, line=1.0, adj=0, font=2, cex=.95)
mtext("200 draws per arm  ·  Yasso07 published = point (no published posterior)",
      3, line=0.0, adj=0, cex=.72, col="grey35")

# ============================== PANEL (c) ====================================
par(mar = c(4.2, 4.8, 3.2, 0.8))
dd <- lapply(MODELS, dif)
plot(NA, xlim=c(.5,3.5), ylim=range(unlist(dd))+c(-.4,.4), xaxt="n",
     xlab="", ylab="difference, ours − published (pp)")
axis(1, at=1:3, labels=MODELS); abline(h=0, col="grey25", lwd=1.5)
for (i in seq_along(MODELS)) bx(dd[[i]], i, MODEL_COL[MODELS[i]], w=0.3)
mtext("(c) the size of the disagreement", 3, line=1.0, adj=0, font=2, cex=.95)
mtext("boxes crossing zero = the anchor choice does not matter", 3, line=0.0,
      adj=0, cex=.72, col="grey35")
dev.off()
cat("wrote manuscript/figures/F15_forward_arms.png\n\n")

# --- ALL scenarios recorded, for the appendix and later use -------------------
cat("SOC change at 2084 (%), median [90%] -- ALL SCENARIOS\n")
for (sc in setdiff(names(R$scen), "control")) {
  cat(sprintf("\n-- %s (ramp to +%.1f C over the projection) --\n", R$scen_lab[[sc]], R$scen[[sc]]))
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
