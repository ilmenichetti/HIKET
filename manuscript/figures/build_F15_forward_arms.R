# =============================================================================
# build_F15_forward_arms.R   (2026-09-02)
#
# F15 -- the consequence figure. Two parameterisations of the same model that
# agree on everything observable, asked what the soil does under warming.
#
#   (a) trajectories, climate held stationary   [3 models x 2 arms]
#   (b) trajectories, warmed to SSP2-4.5        [same design and axes as (a)]
#   (c) (b) minus (a) at 2084, by arm           [boxplots]
#   (d) the difference between arms in that      [boxplots, zero line]
#
# THE RESULT IS STRUCTURAL, and the wording matters: Yasso15/20 prove RESILIENT
# to the input-anchor choice, not "right". Their pool-specific climate modifiers
# leave the carbon-holding pools (H, N: ~60% of the stock) responding to warming
# almost identically in both arms, while the arms differ mainly in the AWE
# modifier -- and AWE holds ~13% of the stock. Yasso07 has ONE modifier for every
# pool, so displacing it displaces the sensitivity of all the carbon at once.
#
# EVERYTHING IS IN tC/ha. (a) and (b) are the raw simulations; (c) and (d) are the
# difference between them -- warmed minus stationary -- so the running titles must
# SAY that, which is the whole fix for the confusion below.
#
# ⚠ (b) MINUS (a) IS NOT THE GAP YOU SEE IN (b). The visible arm gap in (a)/(b) is
# a LEVEL difference and it inverts against the response: Yasso15's arms differ
# +0.73 in level but -0.38 in response, Yasso20's -1.80 in level but +0.81. Do not
# try to read (d) off (b) as a single gap; it is the CHANGE in that gap between
# (a) and (b). The subtitles now state the quantity outright.
#
# ⚠ (c)/(d) ARE PAIRED WITHIN A DRAW (warmed and stationary share the parameter
# vector), which cancels most of the parameter uncertainty and is why the
# intervals are tight. The UNPAIRED alternative -- comparing predicted 2084 stocks
# directly -- gives medians of the same size but intervals that all cross zero,
# because absolute stock carries the full within-arm spread (our Yasso07 arm spans
# 60.6-76.2 at 2084). Both are legitimate answers to different questions; if the
# unpaired version is quoted anywhere, say so explicitly. Do not let the framing
# silently create or destroy significance.
#
# ⚠ The framing also changes which model looks resilient: on the RESPONSE both
# Yasso15 and Yasso20 are insensitive to the anchor; on the PREDICTED STOCK
# Yasso15 is robust and Yasso20 carries a real 1.8 tC/ha disagreement.
#
# ⚠⚠ THE FRAMING CHANGES WHICH MODEL LOOKS RESILIENT, and both readings are true
# of different things. On the RESPONSE, Yasso15 and Yasso20 are both insensitive
# to the anchor. On the PREDICTION, Yasso15 is by far the most robust (+0.7) while
# Yasso20 carries a real 1.8 tC/ha disagreement that the response framing hides.
# State whichever one the sentence in the text is actually about.
#
# ⚠ PANEL TITLES DESCRIBE, THEY DO NOT JUDGE. An earlier version titled panel (a)
# "indistinguishable where the data are", which asserts the reading rather than
# letting the reader make it. Titles now say what is plotted; the interpretation
# belongs in the caption and the text.
#
# ARM "eqinit" (2026-10-22, storyline R3b): the balanced-start calibration, the
# steady-state assumption applied in 1985. Drawn only when present in the .rds.
# (d) then carries TWO comparisons per model: ours - published (the input anchor)
# and ours - balanced (the starting state).
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
ARM_LTY <- c(published = "22", ours = "solid", eqinit = "13")
ARMS    <- intersect(c("published","ours","eqinit"), names(R$out[["Yasso15"]]))
EQ_COL  <- "#8C6D1F"                  # balanced-start boxes in (c)

# (b) minus (a) at 2084, tC/ha -- PAIRED within each draw
resp <- function(M, arm, sc = PLOT_S) {
  z <- R$out[[M]][[arm]]; if (!length(z)) return(numeric(0))
  vapply(z, function(d) tail(d[["traj"]][[sc]], 1) - tail(d[["traj"]][["control"]], 1), numeric(1))
}
resp_pct <- function(M, arm, sc = PLOT_S) {          # appendix / console only
  z <- R$out[[M]][[arm]]; if (!length(z)) return(numeric(0))
  vapply(z, function(d) { c0 <- tail(d[["traj"]][["control"]],1)
                          100*(tail(d[["traj"]][[sc]],1)/c0 - 1) }, numeric(1))
}
dif <- function(M, sc = PLOT_S, vs = "published") {   # ours - <vs> of the (b)-(a) change
  o <- resp(M,"ours",sc); p <- resp(M,vs,sc)
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
YL  <- range(unlist(lapply(MODELS, function(M) lapply(ARMS,
        function(a) c(med(M,a,"control"), med(M,a,PLOT_S))))), oc$lo, oc$hi) + c(-0.8, 0.8)

png("manuscript/figures/F15_forward_arms.png", width = 11.2, height = 8.6,
    units = "in", res = 220)
par(mfrow = c(2,2), mar = c(4.1, 4.6, 3.1, 0.9), mgp = c(2.7, 0.7, 0), las = 1,
    cex.axis = 0.95)

traj_panel <- function(sc, tag, title, sub) {
  plot(NA, xlim=range(yrs), ylim=YL, xlab="", ylab=expression("SOC (tC ha"^-1*")"))
  usr <- par("usr"); rect(2024.5, usr[3], usr[2], usr[4], col="#F6F5F1", border=NA); box()
  abline(v = 2024.5, col="grey45", lwd=1.2, lty=3)
  for (M in MODELS) for (arm in ARMS)
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
  "posterior medians; published arm: litter multiplier fitted to the three campaign means")
legend("bottomright", bty="n", cex=.85, ncol=2, seg.len=2.8,
  legend=c(MODELS, "observed mean", "our calibration", "published parameters",
           if ("eqinit" %in% ARMS) "balanced start"),
  col=c(MODEL_COL[MODELS], "#4A1A0E", "grey30", "grey30", if ("eqinit" %in% ARMS) "grey30"),
  lty=c("solid","solid","solid",NA,"solid","22", if ("eqinit" %in% ARMS) "13"),
  lwd=c(2.4,2.4,2.4,NA,2.4,2.4, if ("eqinit" %in% ARMS) 2.4),
  pch=c(NA,NA,NA,21,NA,NA, if ("eqinit" %in% ARMS) NA), pt.bg="#C85A3C")

# (b) warmed levels -- same design and axes as (a)
traj_panel(PLOT_S, "(b)", "modelled SOC, warming to SSP2-4.5",
  "same models and axes as (a); projection warmed on a linear ramp to +2.6 °C at 2084")

# (c) response
par(mar = c(4.1, 4.6, 3.1, 0.9))
rr <- unlist(lapply(MODELS, function(M) lapply(ARMS, function(a) resp(M,a))))
OFF <- if ("eqinit" %in% ARMS) c(published=-0.27, ours=0, eqinit=0.27) else c(published=-0.22, ours=0.22)
BW  <- if ("eqinit" %in% ARMS) 0.12 else 0.16
plot(NA, xlim=c(.5,3.5), ylim=range(rr)+c(-.7,.7), xaxt="n", xlab="",
     ylab=expression("SOC change, warmed - stationary (tC ha"^-1*")"))
axis(1, at=1:3, labels=MODELS); abline(h=0, col="grey65")
for (i in seq_along(MODELS)) {
  bx(resp(MODELS[i],"published"), i+OFF[["published"]], "#4C72A8",           w=BW)
  bx(resp(MODELS[i],"ours"),      i+OFF[["ours"]],      MODEL_COL[MODELS[i]], w=BW)
  if ("eqinit" %in% ARMS) bx(resp(MODELS[i],"eqinit"), i+OFF[["eqinit"]], EQ_COL, w=BW)
}
legend("bottomleft", bty="n", cex=.83,
  fill=adjustcolor(c("#4C72A8", MODEL_COL[["Yasso15"]], if ("eqinit" %in% ARMS) EQ_COL), .35),
  border=c("#4C72A8", MODEL_COL[["Yasso15"]], if ("eqinit" %in% ARMS) EQ_COL),
  legend=c("published parameters", "our calibration (model colour)",
           if ("eqinit" %in% ARMS) "balanced start"))
mtext("(c) SOC in 2084: warmed minus stationary, by arm", 3, line=1.0, adj=0, font=2, cex=.95)
mtext(sprintf("panel (b) minus panel (a) at 2084; boxes = IQR of %d draws; Yasso07 published = point", R$n_draw),
      3, line=0.0, adj=0, cex=.72, col="grey35")

# (d) difference
par(mar = c(4.1, 4.9, 3.1, 0.9))
dd <- lapply(MODELS, dif)
de <- if ("eqinit" %in% ARMS) lapply(MODELS, dif, vs = "eqinit") else NULL
plot(NA, xlim=c(.5,3.5), ylim=range(unlist(c(dd, de)))+c(-.4,.4), xaxt="n", xlab="",
     ylab=expression("difference in that change (tC ha"^-1*")"))
axis(1, at=1:3, labels=MODELS); abline(h=0, col="grey25", lwd=1.5)
for (i in seq_along(MODELS)) {
  if (is.null(de)) bx(dd[[i]], i, MODEL_COL[MODELS[i]], w=0.28) else {
    bx(dd[[i]], i-0.18, "#4C72A8", w=0.14); bx(de[[i]], i+0.18, EQ_COL, w=0.14) }
}
if (!is.null(de)) legend("bottomleft", bty="n", cex=.83,
  fill=adjustcolor(c("#4C72A8", EQ_COL), .35), border=c("#4C72A8", EQ_COL),
  legend=c("ours - published parameters", "ours - balanced start"))
mtext("(d) difference between arms in that change", 3, line=1.0, adj=0, font=2, cex=.95)
mtext("distance between (c)'s boxes; crossing zero = the choice changes nothing",
      3, line=0.0, adj=0, cex=.72, col="grey35")
dev.off()
cat("wrote manuscript/figures/F15_forward_arms.png\n\n")

cat("WARMED MINUS STATIONARY at 2084 -- ALL SCENARIOS (tC/ha)\n")
for (sc in setdiff(names(R$scen), "control")) {
  cat(sprintf("\n-- %s (ramp to +%.1f C) --\n", R$scen_lab[[sc]], R$scen[[sc]]))
  for (M in MODELS) {
    for (arm in ARMS) { v <- resp(M, arm, sc); if (!length(v)) next
      w <- resp_pct(M, arm, sc)
      cat(sprintf("   %-8s %-10s %6.2f [%6.2f, %6.2f] tC/ha   (%6.2f%%)\n",
                  M, arm, median(v), quantile(v,.05), quantile(v,.95), median(w))) }
    d <- dif(M, sc)
    if (length(d)) cat(sprintf("   %-8s %-10s %6.2f [%6.2f, %6.2f] tC/ha%s\n\n", M, "DIFF",
        median(d), quantile(d,.05), quantile(d,.95),
        if (quantile(d,.05) > 0 || quantile(d,.95) < 0) "   EXCLUDES ZERO" else ""))
    if ("eqinit" %in% ARMS) { d <- dif(M, sc, "eqinit")
      cat(sprintf("   %-8s %-10s %6.2f [%6.2f, %6.2f] tC/ha%s\n\n", M, "OURS-EQ",
        median(d), quantile(d,.05), quantile(d,.95),
        if (quantile(d,.05) > 0 || quantile(d,.95) < 0) "   EXCLUDES ZERO" else "")) }
  }
}
