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
PS <- list(); obs <- NULL
for(m in names(rid)){
  b <- readRDS(sprintf("Calibration_real_data_transient/runs/%s_posterior_predictive_%s.rds", m, rid[[m]]))
  PS[[m]] <- b$posterior_summary
  if (is.null(obs)) obs <- PS[[m]][is.finite(PS[[m]]$soc_obs_tCha), c("plot_id","year","soc_obs_tCha")]
}

# --- PLOT-SET ALIGNMENT (added 2026-08-12) -----------------------------------
# The model curves used to be averaged over ALL 512 plots while the observed
# reference rate comes from the 316 plots measured at BOTH ends. Those are
# different populations, and the mismatch is not cosmetic: for Yasso15 the model
# rate is +0.422 over all plots and +0.379 over the 316, and the 1985 level moves
# 63.3 -> 65.1. Comparing the two inflated the apparent model-observation gap by
# roughly 40%. Both sides are now restricted to the same plots.
o85 <- setNames(obs$soc_obs_tCha[obs$year == 1985L], obs$plot_id[obs$year == 1985L])
o24 <- setNames(obs$soc_obs_tCha[obs$year == 2024L], obs$plot_id[obs$year == 2024L])
both <- as.integer(intersect(names(o85), names(o24)))

mt <- lapply(PS, function(ps) {
  d <- ps[ps$plot_id %in% both, ]
  v <- tapply(d$soc_mean, d$year, mean, na.rm = TRUE)
  data.frame(yr = as.integer(names(v)), soc = as.numeric(v))
})

# --- observed reference rate -------------------------------------------------
# FIXED 2026-08-12. This was hardcoded as (105.3-63.3)/39 = 1.077 tC/ha/yr, from
# the SOC series RETIRED on 2026-08-04 (the pre-homogenization 63/102/105
# medians, which over-counted mineral C ~1.6x for want of a stoniness
# correction). The true value on the live target is ~0.36 -- so the "observed"
# line in panel (b) sat 3x too high, i.e. every model was being shown as failing
# far worse than it is. Now COMPUTED from the bundle's own soc_obs_tCha so it
# cannot go stale again.
# Restricted to plots observed at BOTH endpoints: the campaigns cover different
# plot subsets (404/456/409), and differencing means over different subsets is
# not an estimate of change.
obs_rate <- (mean(o24[as.character(both)]) - mean(o85[as.character(both)])) / (2024 - 1985)
message(sprintf("observed reference rate: %+.3f tC/ha/yr  (n=%d plots with both endpoints)",
                obs_rate, length(both)))
message(sprintf("model rates on the SAME plots: %s",
        paste(sprintf("%s %+.3f", names(mt),
              sapply(mt, function(d) (tail(d$soc,1)-d$soc[1])/39)), collapse="  ")))
# NB two corrections are pending and will both move this number: adding the
# missing 1985 litter (LM) layer, and using the true 1986-1995 sampling years
# (the "1985" campaign is really ~1989, so the interval is ~34.7 yr, not 39).
# The 39 here is deliberate: it must match the model trajectories' own x-axis.

png("manuscript/figures/appendix_delta_reconciliation.png", width = 12, height = 4.3, units = "in", res = 200)
par(mfrow = c(1,3), mar = c(4.0, 4.5, 3.0, 1.0), mgp = c(2.5, 0.7, 0), las = 1)

# axis limits from the data, not hardcoded (they were also stale: ylim=c(84,112)
# against an actual model range of 63-81, so every trajectory drew off-scale)
soc_rng <- range(sapply(mt, function(d) range(d$soc)))
rate_all <- unlist(lapply(mt, function(d) (d$soc - d$soc[1])[-1] / (d$yr - d$yr[1])[-1]))
inc_all  <- unlist(lapply(mt, function(d) roll(diff(d$soc))))

## (a) absolute stock
plot(NA, xlim=c(1985,2024), ylim=soc_rng + c(-1,1)*diff(soc_rng)*0.06,
     xlab="Year", ylab="Mean SOC (tC/ha)",
     main="(a)  Absolute stock: rises to saturation")
for(m in names(rid)) lines(mt[[m]]$yr, mt[[m]]$soc, col=col[m], lwd=2)
legend("bottomright", bty="n", cex=0.8, lwd=2, col=col, legend=names(col))

## (b) cumulative-average rate = the OLD plot
plot(NA, xlim=c(1985,2024), ylim=range(0, rate_all, obs_rate) * c(1, 1.12), xlab="Year",
     ylab=expression("(SOC(t)-SOC(t"[0]*")) / (t-t"[0]*")  (tC/ha/yr)"),
     main="(b)  Cumulative-average rate  [the OLD figure]")
abline(h=obs_rate, col="grey30", lwd=1.6, lty=2); text(2005, obs_rate, "observed", pos=3, cex=0.75, col="grey30")
for(m in names(rid)){ d<-mt[[m]]; r<-(d$soc-d$soc[1])/(d$yr-d$yr[1]); lines(d$yr[-1], r[-1], col=col[m], lwd=2) }
mtext("declining = DECELERATION, not SOC loss; below obs = 1985 over-prediction", side=3, line=-1.1, cex=0.62, col="grey35")

## (c) instantaneous increment
plot(NA, xlim=c(1986,2024), ylim=range(0, inc_all, na.rm=TRUE) * c(1.15, 1.15), xlab="Year",
     ylab="dSOC/dt (tC/ha/yr, 5-yr mean)", main="(c)  Instantaneous increment")
abline(h=0, col="grey70", lty=3)
for(m in names(rid)){ d<-mt[[m]]; inc<-roll(diff(d$soc)); lines(d$yr[-1], inc, col=col[m], lwd=2) }
dev.off()
cat("Wrote manuscript/figures/appendix_delta_reconciliation.png\n")
