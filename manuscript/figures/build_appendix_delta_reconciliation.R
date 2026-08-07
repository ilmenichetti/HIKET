source("manuscript/figures/run_ids.R")   # auto-selects current RUN_IDs
setwd("/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling")
# Appendix reconciliation: the OLD cumulative-average-rate plot and the NEW absolute-stock
# F4 are the SAME six trajectories under different transforms -- no contradiction.
#   (a) absolute mean SOC          -- rises toward saturation (F3 / F4)
#   (b) cumulative-average rate     -- (SOC(t)-SOC(t0))/(t-t0)  [the OLD multimodel_delta_soc]
#   (c) instantaneous increment     -- dSOC/dt (5-yr mean)
# (b)'s decline is (a)'s DECELERATION, not a loss of SOC; all three stay >0 early.

rid <- as.list(RID)
source("manuscript/figures/model_palette.R")   # shared per-model palette
col <- MODEL_COL
roll <- function(x,k=5){ n<-length(x); s<-rep(NA,n); h<-(k-1)/2
  for(i in seq_len(n)){ s[i]<-mean(x[max(1,i-h):min(n,i+h)]) }; s }
mt <- list()
for(m in names(rid)){
  b <- readRDS(sprintf("Calibration_real_data_transient/runs/%s_posterior_predictive_%s.rds", m, rid[[m]]))
  ps <- b$posterior_summary; v <- tapply(ps$soc_mean, ps$year, mean, na.rm=TRUE)
  mt[[m]] <- data.frame(yr=as.integer(names(v)), soc=as.numeric(v))
}
obs_rate <- (105.3-63.3)/(2024-1985)   # observed mean annual rate (VMI8 -> Komeetta)

png("manuscript/figures/appendix_delta_reconciliation.png", width = 12, height = 4.3, units = "in", res = 200)
par(mfrow = c(1,3), mar = c(4.0, 4.5, 3.0, 1.0), mgp = c(2.5, 0.7, 0), las = 1)

## (a) absolute stock
plot(NA, xlim=c(1985,2024), ylim=c(84,112), xlab="Year", ylab="Mean SOC (tC/ha)",
     main="(a)  Absolute stock: rises to saturation")
for(m in names(rid)) lines(mt[[m]]$yr, mt[[m]]$soc, col=col[m], lwd=2)
legend("bottomright", bty="n", cex=0.8, lwd=2, col=col, legend=names(col))

## (b) cumulative-average rate = the OLD plot
plot(NA, xlim=c(1985,2024), ylim=c(0,2.1), xlab="Year",
     ylab=expression("(SOC(t)-SOC(t"[0]*")) / (t-t"[0]*")  (tC/ha/yr)"),
     main="(b)  Cumulative-average rate  [the OLD figure]")
abline(h=obs_rate, col="grey30", lwd=1.6, lty=2); text(2005, obs_rate, "observed", pos=3, cex=0.75, col="grey30")
for(m in names(rid)){ d<-mt[[m]]; r<-(d$soc-d$soc[1])/(d$yr-d$yr[1]); lines(d$yr[-1], r[-1], col=col[m], lwd=2) }
mtext("declining = DECELERATION, not SOC loss; below obs = 1985 over-prediction", side=3, line=-1.1, cex=0.62, col="grey35")

## (c) instantaneous increment
plot(NA, xlim=c(1986,2024), ylim=c(-0.5,1.7), xlab="Year",
     ylab="dSOC/dt (tC/ha/yr, 5-yr mean)", main="(c)  Instantaneous increment")
abline(h=0, col="grey70", lty=3)
for(m in names(rid)){ d<-mt[[m]]; inc<-roll(diff(d$soc)); lines(d$yr[-1], inc, col=col[m], lwd=2) }
dev.off()
cat("Wrote manuscript/figures/appendix_delta_reconciliation.png\n")
