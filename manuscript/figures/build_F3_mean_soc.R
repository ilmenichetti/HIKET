source("manuscript/figures/run_ids.R")   # auto-selects current RUN_IDs
setwd("/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling")
# F3 (rebuild) -- headline: cross-plot MEAN SOC trajectory 1985-2024, all six models,
# tracking the campaign means. Uniform model palette (Temperature Diverging), larger text.
# Reuses the F4 spin-up/trajectory cache (stored = posterior mean + 95% band per year).
source("manuscript/figures/model_palette.R")
# Must match the RUN_ID-keyed cache written by build_F4_initialization.R (run that
# first). Keyed so a re-calibration cannot leave this figure silently stale.
.f4_cache <- sprintf("manuscript/figures/F4_cache_%s.rds",
                     substr(paste(RID[FIG_MODELS], collapse = "-"), 1, 120))
if (!file.exists(.f4_cache))
  stop("F4 cache for the current RUN_IDs is missing -- run build_F4_initialization.R first:\n  ",
       .f4_cache, call. = FALSE)
cache <- readRDS(.f4_cache)$stored

# observed campaign means +/- 95% CI (same source as F4)
om <- readRDS(sprintf("Data/model_inputs/Yasso20_inputs_%s.rds", RID[["Yasso20"]]))$obs_meta
obs <- do.call(rbind, lapply(names(om), function(p){ z<-om[[p]]; if(!length(z$soc_obs)) return(NULL)
  data.frame(year=1984L+z$idx, soc=z$soc_obs) }))
cm <- aggregate(soc~year, obs, function(x) c(m=mean(x), lo=mean(x)-1.96*sd(x)/sqrt(length(x)),
                                             hi=mean(x)+1.96*sd(x)/sqrt(length(x))))
cm <- data.frame(year=cm$year, m=cm$soc[,"m"], lo=cm$soc[,"lo"], hi=cm$soc[,"hi"])

YR <- 1985:2024
png("manuscript/figures/F3_mean_soc.png", width=9.5, height=6, units="in", res=200)
par(mar=c(4.6,5.0,3.2,1.2), mgp=c(3.0,0.8,0), las=1, cex.axis=1.15, cex.lab=1.3)
yl <- range(sapply(MODEL_ORDER, function(m){ d<-cache[[m]]; i<-d$year %in% YR; c(d$lo[i],d$hi[i]) }), cm$lo, cm$hi)
plot(NA, xlim=range(YR), ylim=yl, xlab="Year", ylab="Mean SOC across plots (tC/ha)",
     main="F3 . Accumulation, then incipient saturation --- six models track the campaigns")
for (m in MODEL_ORDER) { d<-cache[[m]]; i<-d$year %in% YR
  polygon(c(d$year[i],rev(d$year[i])), c(d$lo[i],rev(d$hi[i])), col=adjustcolor(MODEL_COL[m],0.13), border=NA)
  lines(d$year[i], d$m[i], col=MODEL_COL[m], lwd=2.6) }
arrows(cm$year, cm$lo, cm$year, cm$hi, angle=90, code=3, length=0.05, col="black", lwd=2)
points(cm$year, cm$m, pch=21, bg="firebrick", col="black", cex=1.7, lwd=1.3)
text(cm$year, cm$hi, c("VMI8\n1985","Biosoil\n2006","Komeetta\n2024"), pos=c(4,3,2), offset=0.7, cex=0.85, font=2)
legend("bottomright", legend=MODEL_ORDER, col=MODEL_COL[MODEL_ORDER], lwd=2.8, bty="n",
       cex=1.05, ncol=2, title="model means (simple -> complex)")
dev.off(); cat("wrote F3_mean_soc.png\n")
