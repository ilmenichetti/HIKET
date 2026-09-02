# =============================================================================
# build_A1_forward_appendix.R   (2026-09-02)
#
# Appendix companion to F15. Three panels plus a table:
#   (a) MECHANISM -- warming response against each draw's own intrinsic MRT.
#       ⚠ PREDICTION: the clouds should NOT collapse onto one curve. The response
#       is governed by the STOCK-WEIGHTED xi sensitivity, sum_i w_i dlog(xi_i)/dT,
#       not by bulk MRT: Yasso15/20 carry pool-specific modifiers, and the arms
#       differ mainly in the AWE one while H and N hold ~60% of the carbon. If the
#       clouds DID collapse, bulk MRT would be the right summary statistic and the
#       structural argument in the main text would be wrong.
#   (b) the same contrast as F15(c)/(d) in RELATIVE units (%), which are exactly
#       invariant to sigma_input and therefore isolate kinetics.
#   (c) all three SSP scenarios, not just the mid one plotted in F15.
#   T_A1: how the published arm was constructed.
#
# Data: doublechecks/f15_forward_experiment.rds, doublechecks/published_arm_dat.rds
# =============================================================================

setwd("/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling")
source("manuscript/figures/model_palette.R")
R   <- readRDS("doublechecks/f15_forward_experiment.rds")
ARM <- readRDS("doublechecks/published_arm_dat.rds")
MODELS <- c("Yasso07","Yasso15","Yasso20")

pull <- function(M, arm, sc) {
  z <- R$out[[M]][[arm]]; if (!length(z)) return(NULL)
  data.frame(
    mrt = vapply(z, function(d) d$mrt, numeric(1)),
    abs = vapply(z, function(d) tail(d$traj[[sc]],1) - tail(d$traj[["control"]],1), numeric(1)),
    rel = vapply(z, function(d) 100*(tail(d$traj[[sc]],1)/tail(d$traj[["control"]],1) - 1), numeric(1)),
    model = M, arm = arm, stringsAsFactors = FALSE)
}
bx <- function(v, at, col, w) {
  if (!length(v)) return(invisible())
  if (length(v) == 1) { segments(at-w, v, at+w, v, col=col, lwd=3.2); return(invisible()) }
  b <- boxplot(v, plot=FALSE)
  rect(at-w, b$stats[2], at+w, b$stats[4], col=adjustcolor(col,.35), border=col, lwd=1.5)
  segments(at-w, b$stats[3], at+w, b$stats[3], col=col, lwd=2.8)
  segments(at, b$stats[1], at, b$stats[2], col=col, lwd=1.3)
  segments(at, b$stats[4], at, b$stats[5], col=col, lwd=1.3)
}

png("manuscript/figures/A1_forward_appendix.png", width = 13.4, height = 4.7,
    units="in", res=220)
layout(matrix(1:3, nrow=1)); par(mar=c(4.3,4.7,3.2,0.9), mgp=c(2.7,0.7,0), las=1, cex.axis=.95)

# ---- (a) does the response collapse onto MRT? -------------------------------
D <- do.call(rbind, lapply(MODELS, function(M)
       do.call(rbind, lapply(c("published","ours"), function(a) pull(M,a,"ssp245")))))
D <- D[is.finite(D$mrt) & is.finite(D$rel), ]
plot(D$mrt, D$rel, type="n", xlab="intrinsic MRT of the draw (yr)",
     ylab="SOC change under SSP2-4.5 (%)")
for (M in MODELS) for (a in c("published","ours")) {
  k <- D$model==M & D$arm==a
  points(D$mrt[k], D$rel[k], pch=ifelse(a=="ours",16,1), cex=.55,
         col=adjustcolor(MODEL_COL[M], ifelse(a=="ours",.55,.75)))
}
legend("bottomright", bty="n", cex=.8, ncol=2,
  legend=c(MODELS,"our calibration","published"),
  col=c(MODEL_COL[MODELS],"grey30","grey30"), pch=c(16,16,16,16,1))
mtext("(a) response against bulk MRT", 3, line=1.0, adj=0, font=2, cex=.95)
mtext("one point per draw; the clouds do NOT lie on a common curve", 3, line=0.0,
      adj=0, cex=.72, col="grey35")

# ---- (b) the F15 contrast in relative units ---------------------------------
rr <- unlist(lapply(MODELS, function(M) lapply(c("published","ours"),
        function(a) pull(M,a,"ssp245")$rel)))
plot(NA, xlim=c(.5,3.5), ylim=range(rr)+c(-.4,.4), xaxt="n", xlab="",
     ylab="SOC change, warmed - stationary (%)")
axis(1, at=1:3, labels=MODELS); abline(h=0, col="grey65")
for (i in seq_along(MODELS)) {
  bx(pull(MODELS[i],"published","ssp245")$rel, i-0.22, "#4C72A8", w=0.16)
  bx(pull(MODELS[i],"ours","ssp245")$rel,      i+0.22, MODEL_COL[MODELS[i]], w=0.16)
}
legend("bottomleft", bty="n", cex=.8,
  fill=c(adjustcolor("#4C72A8",.35), adjustcolor(MODEL_COL[["Yasso15"]],.35)),
  border=c("#4C72A8", MODEL_COL[["Yasso15"]]),
  legend=c("published parameters","our calibration (model colour)"))
mtext("(b) the same contrast in relative units", 3, line=1.0, adj=0, font=2, cex=.95)
mtext("exactly invariant to sigma_input, so this isolates kinetics", 3, line=0.0,
      adj=0, cex=.72, col="grey35")

# ---- (c) all three scenarios ------------------------------------------------
SC <- setdiff(names(R$scen), "control")
dif <- function(M, sc) { o <- pull(M,"ours",sc)$abs; p <- pull(M,"published",sc)$abs
  v <- as.vector(outer(o,p,"-")); if (length(v)>4000) sample(v,4000) else v }
alld <- unlist(lapply(MODELS, function(M) lapply(SC, function(s) dif(M,s))))
plot(NA, xlim=c(.5,3.5), ylim=range(alld)+c(-.3,.3), xaxt="n", xlab="",
     ylab=expression("difference, ours - published (tC ha"^-1*")"))
axis(1, at=1:3, labels=MODELS); abline(h=0, col="grey25", lwd=1.4)
offs <- c(-0.26, 0, 0.26)
for (i in seq_along(MODELS)) for (j in seq_along(SC))
  bx(dif(MODELS[i], SC[j]), i+offs[j],
     adjustcolor(MODEL_COL[MODELS[i]], alpha.f = c(.45,.7,1)[j]), w=0.10)
legend("topleft", bty="n", cex=.8,
       fill=vapply(c(.45,.7,1), function(a) adjustcolor("grey35", a), character(1)),
       border="grey35", legend=unname(R$scen_lab[SC]))
mtext("(c) all three scenarios", 3, line=1.0, adj=0, font=2, cex=.95)
mtext("warming level changes the size, not the pattern", 3, line=0.0,
      adj=0, cex=.72, col="grey35")
dev.off()
cat("wrote manuscript/figures/A1_forward_appendix.png\n")

# ---- T_A1: how the published arm was built ----------------------------------
f <- "manuscript/figures/T_A1_published_arm.tex"; con <- file(f,"w")
wl <- function(...) cat(..., "\n", sep="", file=con)
wl("% auto-generated by build_A1_forward_appendix.R")
wl("\\begin{tabular}{lccccc}")
wl("\\hline")
wl("Model & source of the & MRT & $\\sigma_{\\mathrm{input}}$ & eff.\\ flux & misfit \\\\")
wl(" & published vector & (yr) & (refitted) & (tC\\,ha$^{-1}$yr$^{-1}$) & (rms log) \\\\")
wl("\\hline")
src <- c(Yasso07="own defaults (no \\texttt{.dat})",
         Yasso15="\\texttt{Yasso15.dat}", Yasso20="\\texttt{Yasso20.dat}")
for (M in MODELS) with(ARM[[M]], wl(sprintf(
  "%s & %s & %.2f & %.3f & %.2f & %.4f \\\\", M, src[[M]], mrt_pub, s_pub, s_pub*2.511, rms_pub)))
wl("\\hline")
wl("\\multicolumn{6}{l}{\\footnotesize Kinetics fixed at published values; only")
wl("$\\sigma_{\\mathrm{input}}$ refitted, to the three campaign means.} \\\\")
wl(sprintf("\\multicolumn{6}{l}{\\footnotesize $\\sigma_{\\mathrm{init}}$ pinned to our posterior in BOTH arms (%.2f/%.2f/%.2f).} \\\\",
           ARM$Yasso07$sigma_init, ARM$Yasso15$sigma_init, ARM$Yasso20$sigma_init))
wl("\\end{tabular}")
close(con); cat("wrote", f, "\n")

# ---- does it collapse? quantify ---------------------------------------------
cat("\nR^2 of a single common curve (response ~ log MRT) pooled over all models:\n")
m0 <- lm(rel ~ log(mrt), data = D)
cat(sprintf("   pooled                : R2 = %.3f\n", summary(m0)$r.squared))
for (M in MODELS) { k <- D$model==M
  cat(sprintf("   within %-8s       : R2 = %.3f\n", M, summary(lm(rel ~ log(mrt), D[k,]))$r.squared)) }
