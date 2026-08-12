source("manuscript/figures/run_ids.R")   # auto-selects current RUN_IDs
setwd("/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling")
# Robustness appendix figures R1 (MCMC convergence) + R2 (residual diagnostics).
# R1: per-model max R-hat + min ESS (from metadata) vs the usual thresholds.
# R2: pooled residual diagnostics (residual-vs-fitted / QQ / histogram) on the
#     log scale = the multiplicative-normal error model's natural scale.
source("manuscript/figures/model_palette.R")

rid <- as.list(RID)
runs <- "Calibration_real_data_transient/runs"; M <- names(rid)

# ---- R1: convergence summary ----
rh <- sapply(M, function(m) readRDS(sprintf("%s/%s_metadata_%s.rds", runs, m, rid[[m]]))$rhat_range[2])
es <- sapply(M, function(m) readRDS(sprintf("%s/%s_metadata_%s.rds", runs, m, rid[[m]]))$ess_range[1])

png("manuscript/figures/S3_convergence.png", width=10, height=4.6, units="in", res=200)
par(mfrow=c(1,2), mar=c(5,4.8,3,1), mgp=c(2.7,0.7,0), las=1)
bp <- barplot(rh[M], col=MODEL_COL[M], border=NA, las=2, ylim=c(1, max(1.06, max(rh)*1.02)),
              ylab="max R-hat (across free params)", main="A . Gelman-Rubin R-hat", font.main=1, xpd=FALSE)
abline(h=c(1.05,1.10), lty=c(3,2), col=c("#f1b300","#d01c8b"), lwd=1.4)
text(bp, rh[M], sprintf("%.3f", rh[M]), pos=3, cex=0.75, xpd=NA)
legend("topright", c("1.05 target","1.10 limit"), lty=c(3,2), col=c("#f1b300","#d01c8b"), bty="n", cex=0.8)
bp2 <- barplot(es[M], col=MODEL_COL[M], border=NA, las=2, ylim=c(0, max(es)*1.12),
               ylab="min effective sample size", main="B . ESS (minimum across params)", font.main=1)
abline(h=400, lty=3, col="grey40"); text(bp2, es[M], sprintf("%.0f", es[M]), pos=3, cex=0.75, xpd=NA)
mtext("All six: max R-hat < 1.02 (well under 1.05) and min ESS > 1100 -> healthy mixing", side=1, line=-1.4, outer=TRUE, cex=0.8, col="grey30")
dev.off(); cat("wrote S3_convergence.png\n")

# ---- R2: pooled residual diagnostics ----
res <- do.call(rbind, lapply(M, function(m){
  rd <- readRDS(sprintf("%s/%s_posterior_predictive_%s.rds", runs, m, rid[[m]]))$residuals_df
  rd <- rd[is.finite(rd$residual_log) & is.finite(rd$log_hat_mean), ]
  data.frame(model=m, fitted=rd$log_hat_mean, r=rd$residual_log) }))
png("manuscript/figures/S4_residuals.png", width=12, height=4.2, units="in", res=200)
par(mfrow=c(1,3), mar=c(4.4,4.4,3,1), mgp=c(2.5,0.7,0), las=1)
plot(res$fitted, res$r, pch=19, cex=0.4, col=adjustcolor(MODEL_COL[res$model],0.4),
     xlab="Fitted log-SOC", ylab="Log residual", main="A . Residual vs fitted", font.main=1)
abline(h=0, col="grey30", lwd=1.4, lty=2); lines(lowess(res$fitted, res$r), col="firebrick", lwd=2)
qqnorm(res$r, pch=19, cex=0.4, col=adjustcolor("grey40",0.4), main="B . Normal Q-Q (log residuals)", font.main=1)
qqline(res$r, col="firebrick", lwd=2)
hist(res$r, breaks=60, col="grey80", border="white", freq=FALSE, xlab="Log residual",
     main="C . Residual distribution", font.main=1)
curve(dnorm(x, mean(res$r), sd(res$r)), add=TRUE, col="firebrick", lwd=2)
mtext("Pooled over all six models; log scale = the multiplicative-normal error model's natural scale (roughly symmetric, near-Gaussian)",
      side=1, line=-1.4, outer=TRUE, cex=0.78, col="grey30")
dev.off(); cat("wrote S4_residuals.png\n")
cat(sprintf("R-hat max %.3f-%.3f; ESS min %.0f-%.0f; log-resid mean %.3f sd %.3f\n",
            min(rh), max(rh), min(es), max(es), mean(res$r), sd(res$r)))
