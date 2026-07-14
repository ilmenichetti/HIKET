setwd("/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling")
suppressMessages(library(BayesianTools))
new <- list(SP1="20260710_104903", TP2="20260710_104904", TP3="20260710_104904",
            Yasso07="20260710_104902", Yasso15="20260710_104902", Yasso20="20260710_102431")
old <- list(SP1="20260608_015026", TP2="20260608_020212", TP3="20260630_090644",
            Yasso07="20260611_032825", Yasso15="20260611_032825", Yasso20="20260611_033140")
littercols <- function(df) grep("^(nwl_|fwl_|cwl_)", names(df), value=TRUE)
totlit <- function(df){ if("J_total" %in% names(df)) return(mean(df$J_total)); lc<-littercols(df); mean(rowSums(df[,lc,drop=FALSE])) }
effflux <- function(m, rid){
  post <- readRDS(sprintf("Calibration_real_data_transient/runs/%s_posterior_%s.rds", m, rid))
  inp  <- readRDS(sprintf("Data/model_inputs/%s_inputs_%s.rds", m, rid))
  si <- quantile(getSample(post)[,"sigma_input"], c(.025,.5,.975))
  Jm <- median(sapply(inp$plots, function(id) totlit(inp$inputs_by_plot[[id]])))
  si * Jm
}
models <- names(new)
B <- t(sapply(models, function(m) effflux(m, old[[m]])))   # before
A <- t(sapply(models, function(m) effflux(m, new[[m]])))   # after
rawJ <- median(sapply(models, function(m){ inp<-readRDS(sprintf("Data/model_inputs/%s_inputs_%s.rds",m,new[[m]])); median(sapply(inp$plots,function(id) totlit(inp$inputs_by_plot[[id]]))) }))
LO <- 0.5; HI <- 8.7

png("manuscript/figures/F9_effective_flux_vs_ceiling.png", width=9, height=5.2, units="in", res=200)
par(mar=c(4.6,6.2,2.4,1.4), xaxs="i")
n <- length(models); ypos <- rev(seq_len(n))
xlim <- c(0.5, 60)
plot(NA, xlim=xlim, ylim=c(0.4,n+0.6), log="x", axes=FALSE, xlab="", ylab="")
# physical envelope band
rect(LO, 0.2, HI, n+0.8, col="#e6f2e6", border=NA)
abline(v=HI, col="#2e7d32", lwd=1.6, lty=1)
abline(v=rawJ, col="grey45", lwd=1.1, lty=3)
# axes
xt <- c(0.5,1,2,5,10,20,50)
axis(1, at=xt, labels=xt); axis(2, at=ypos, labels=models, las=1, tick=FALSE)
mtext(expression("Effective litter flux  " * sigma[input] %*% bar(J) * "  (tC " * ha^-1 * " " * yr^-1 * ", log scale)"), side=1, line=2.9)
box()
off <- 0.16
for(i in seq_len(n)){
  y <- ypos[i]
  # before (unbounded) — grey open, upper offset
  segments(B[i,1], y+off, B[i,3], y+off, col="grey55", lwd=2)
  points(B[i,2], y+off, pch=21, bg="white", col="grey35", cex=1.3, lwd=1.6)
  # after (bounded) — filled blue, lower offset
  ab <- if(B[i,3] > HI) "#c62828" else "#1565c0"   # red if the before-run breached ceiling
  segments(A[i,1], y-off, A[i,3], y-off, col="#1565c0", lwd=2.4)
  points(A[i,2], y-off, pch=21, bg="#1565c0", col="#0d3c78", cex=1.4, lwd=1.2)
  # annotate the impossible before-values
  if(B[i,2] > HI) text(B[i,2], y+off, sprintf("%.0f", B[i,2]), pos=3, offset=0.35, col="#c62828", font=2, cex=0.8)
}
# labels for ceiling / rawJ
text(HI, n+0.55, "boreal NPP ceiling ≈ 8.7", col="#2e7d32", pos=2, cex=0.78, font=3)
text(rawJ, 0.55, expression("tree litter " * bar(J) %~~% "2.5"), col="grey35", pos=4, cex=0.78)
legend("bottomright", inset=c(0.01,0.06), bty="n", cex=0.85,
       legend=c("unbounded run (before)","flux-bounded run (after)","physical envelope [0.5, 8.7]"),
       pch=c(21,21,22), pt.bg=c("white","#1565c0","#e6f2e6"), col=c("grey35","#0d3c78","#e6f2e6"), pt.cex=c(1.3,1.4,2))
title(main="Effective litter flux against the physical NPP envelope: bounding closes the escape hatch", cex.main=0.98)
dev.off()
cat("BEFORE (median effJ):\n"); print(round(B[,2],1))
cat("AFTER  (median effJ):\n"); print(round(A[,2],1))
cat("Wrote manuscript/figures/F9_effective_flux_vs_ceiling.png\n")
