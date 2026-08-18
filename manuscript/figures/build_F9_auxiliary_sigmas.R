source("manuscript/figures/run_ids.R")   # auto-selects current RUN_IDs
setwd("/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling")
# F9 (merged, was F9+F10) — the two auxiliary uncertainty parameters, one plate:
#   (a) sigma_input as EFFECTIVE LITTER FLUX vs the physical NPP envelope (bounded run only)
#   (b) sigma_init as the below-equilibrium 1917 start (post-exploitation recovery)
# Both are per-model forest plots sharing the same model rows, so they read as a pair.
# The unbounded "before" run is NOT shown (main story = bounded only). It used to have
# its own figure, build_F9_effective_flux.R -> F9_effective_flux_vs_ceiling.png in
# appendix_sigma_input.tex; both were DELETED 2026-08-17 (the builder's "old"
# posteriors had been cleaned up, so it could no longer be regenerated).
suppressMessages(library(BayesianTools))

new <- as.list(RID)
models <- names(new); n <- length(models); ypos <- rev(seq_len(n))

littercols <- function(df) grep("^(nwl_|fwl_|cwl_)", names(df), value=TRUE)
totlit <- function(df){ if("J_total" %in% names(df)) return(mean(df$J_total)); lc<-littercols(df); mean(rowSums(df[,lc,drop=FALSE])) }
post_of <- function(m) readRDS(sprintf("Calibration_real_data_transient/runs/%s_posterior_%s.rds", m, new[[m]]))

# --- panel A data: effective litter flux = sigma_input * mean litter -----------
effflux <- function(m){
  s  <- getSample(post_of(m))
  inp<- readRDS(sprintf("Data/model_inputs/%s_inputs_%s.rds", m, new[[m]]))
  si <- quantile(s[,"sigma_input"], c(.025,.5,.975))
  Jm <- median(sapply(inp$plots, function(id) totlit(inp$inputs_by_plot[[id]])))
  si * Jm
}
# --- panel B data: sigma_init -------------------------------------------------
siginit <- function(m) quantile(getSample(post_of(m))[,"sigma_init"], c(.025,.5,.975))

A  <- t(sapply(models, effflux))   # bounded effective flux
SI <- t(sapply(models, siginit))   # sigma_init
rawJ <- median(sapply(models, function(m){ inp<-readRDS(sprintf("Data/model_inputs/%s_inputs_%s.rds",m,new[[m]])); median(sapply(inp$plots,function(id) totlit(inp$inputs_by_plot[[id]]))) }))
# Two NPP envelopes, drawn as nested bands (Gower 2001, Class I evergreen maxima):
#   GLOBAL  [0.5, 8.7]   the original window -- so wide it never bound any model
#   NORDIC  [0.05, 4.62] the decided replacement (2026-08-14), region-matched
# ⚠ BOTH are MAXIMA, while sigma_input x Jbar is a national MEAN (sigma_input is one
# global scalar). Bounding a mean with a maximum is the category error that made the
# global window vacuous; the Nordic window narrows it but does not change its kind.
# Read an excursion as a flag, not as proof of physical impossibility.
LO_G <- 0.5;  HI_G <- 8.7        # global
LO_N <- 0.05; HI_N <- 4.62       # Nordic (decided)
LO <- LO_G; HI <- HI_G
source("manuscript/figures/model_palette.R")   # shared per-model palette
col <- MODEL_COL
darken <- function(c, f=0.65){ v <- col2rgb(c)/255; rgb(v[1]*f, v[2]*f, v[3]*f) }

png("manuscript/figures/F9_auxiliary_sigmas.png", width=11, height=4.8, units="in", res=200)
par(mfrow=c(1,2), mar=c(4.6,5.8,3.0,1.2), mgp=c(2.7,0.7,0), cex.axis=1.05, cex.lab=1.2, cex.main=1.15)

## ---- Panel A: sigma_input as effective litter flux -------------------------
plot(NA, xlim=c(0.5,12), ylim=c(0.4,n+0.6), log="x", xaxs="i", axes=FALSE, xlab="", ylab="")
rect(LO_G, 0.2, HI_G, n+0.8, col="#eef6ee", border=NA)             # global envelope (pale)
rect(0.5,  0.2, HI_N, n+0.8, col="#c9e2c9", border=NA)             # Nordic envelope (stronger)
abline(v=HI_G, col="#7fb07f", lwd=1.4, lty=2)
abline(v=HI_N, col="#1b5e20", lwd=2.0, lty=1)
abline(v=LO_G, col="#7fb07f", lwd=1.2, lty=2)
abline(v=rawJ, col="grey45", lwd=1.1, lty=3)                        # tree litter mean
xt <- c(0.5,1,2,5,10)
axis(1, at=xt, labels=xt); axis(2, at=ypos, labels=models, las=1, tick=FALSE); box()
mtext(expression("Effective litter flux  " * sigma[input] %*% bar(J) * "  (tC " * ha^-1 * " " * yr^-1 * ")"),
      side=1, line=2.9, cex=0.9)
for(i in seq_len(n)){
  y <- ypos[i]; cm <- col[models[i]]
  segments(A[i,1], y, A[i,3], y, col=cm, lwd=2.6)
  points(A[i,2], y, pch=21, bg=cm, col=darken(cm), cex=1.6, lwd=1.3)
}
text(0.55, 0.62, "NPP envelope (Gower 2001 Class I maxima)", col="#2e7d32", font=3, cex=0.68, adj=0)
text(0.55, 0.42, "dashed = global 8.7    solid = Nordic 4.62 (decided)", col="#2e7d32", font=3, cex=0.68, adj=0)
text(rawJ, n+0.55, expression(bar(J) %~~% "2.5"), col="grey35", pos=4, cex=0.72)
title(main=expression("A  "*sigma[input]*": litter input multiplier, as effective flux"),
      cex.main=0.98, font.main=1)

## ---- Panel B: sigma_init (below-equilibrium 1917 start) --------------------
plot(NA, xlim=c(0,1.15), ylim=c(0.4,n+0.6), xaxs="i", axes=FALSE, xlab="", ylab="")
rect(0, 0.2, 1, n+0.8, col="#f0f0f0", border=NA)                    # below-equilibrium region
abline(v=1, col="grey35", lwd=1.6, lty=1)                           # equilibrium (full) 1917 start
axis(1, at=seq(0,1,0.25)); axis(2, at=ypos, labels=models, las=1, tick=FALSE); box()
mtext(expression(sigma[init]*"  (initial carbon state, fraction of equilibrium)"),
      side=1, line=2.9, cex=0.9)
for(i in seq_len(n)){
  y <- ypos[i]; cm <- col[models[i]]
  segments(SI[i,1], y, SI[i,3], y, col=cm, lwd=2.6)
  points(SI[i,2], y, pch=21, bg=cm, col=darken(cm), cex=1.6, lwd=1.3)
}
text(1, 0.72, "equilibrium\nstart = 1", col="grey35", pos=2, cex=0.72, font=3)
text(0.5, n+0.55, "below-equilibrium (post-exploitation recovery)", col="grey35", font=3, cex=0.72)
title(main=expression("B  "*sigma[init]*": below-equilibrium 1917 start"),
      cex.main=0.98, font.main=1)

dev.off()
cat("Panel A effective flux (median):\n"); print(round(A[,2],1))
cat("Panel B sigma_init (median):\n"); print(round(SI[,2],2))
cat("Wrote manuscript/figures/F9_auxiliary_sigmas.png\n")
