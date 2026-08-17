setwd("/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling")
# F6b — skill (R^2) does NOT climb the structural-complexity ladder, on EITHER axis.
# Source: multimodel_metrics(_holdout)_20260710_* (flux-bounded run). Values hardcoded
# (they are the headline R^2 already reported in the metrics tables).

models <- c("SP1","TP2","TP3","Yasso07","Yasso15","Yasso20")
pools  <- c(1, 2, 3, 5, 5, 5)                       # structural pool count
calib  <- c(0.062, 0.070, 0.077, 0.067, 0.060, 0.047)
hold   <- c(0.026, 0.027, 0.037, 0.027, 0.026, 0.015)
source("manuscript/figures/model_palette.R")   # shared per-model palette
col_m  <- MODEL_COL

png("manuscript/figures/F6b_skill_vs_complexity.png", width=10, height=4.8, units="in", res=200)
par(mfrow=c(1,2), mar=c(4.4,4.6,3.0,1.0), mgp=c(2.6,0.7,0), cex.axis=1.05, cex.main=1.15)
yl <- c(0, 0.085)

## ---- Panel A: pool-complexity axis ----------------------------------------
xj <- pools + c(0,0,0,-0.16,0,0.16)                 # jitter the three 5-pool models
plot(NA, xlim=c(0.6,5.6), ylim=yl, axes=FALSE, xlab="", ylab="")
axis(1, at=c(1,2,3,5), labels=c("1\n(SP1)","2\n(TP2)","3\n(TP3)","5\n(Yasso 07/15/20)"),
     padj=0.5, cex.axis=0.85)
axis(2, las=1); box()
mtext("Number of pools (structural complexity)", side=1, line=2.9, cex=0.95)
mtext(expression("Predictive "*R^2), side=2, line=2.9, cex=0.95)
# flat OLS reference through calibration R^2 (shows no upward climb)
abline(lm(calib ~ pools), col="grey55", lwd=1.4, lty=2)
# holdout (open) then calibration (filled)
points(xj, hold,  pch=21, bg="white",    col=col_m, cex=1.5, lwd=1.8)
points(xj, calib, pch=21, bg=col_m,       col=col_m, cex=1.7, lwd=1.2)
segments(xj, hold, xj, calib, col=adjustcolor(col_m,0.5), lwd=1.2)
text(xj, calib, models, pos=3, offset=0.55, cex=0.68, col=col_m, font=2)
legend("topright", bty="n", cex=0.8, pch=21,
       pt.bg=c("black","white"), col="black",
       legend=c("calibration", "independent holdout"))
title(main="A · Pool-complexity axis", cex.main=0.95, font.main=1)

## ---- Panel B: climate-integration sub-axis (Yasso only) -------------------
yx <- 1:3; yi <- 4:6
plot(NA, xlim=c(0.7,3.3), ylim=yl, axes=FALSE, xlab="", ylab="")
axis(1, at=yx, labels=c("Yasso07","Yasso15","Yasso20"), cex.axis=0.9)
axis(2, las=1); box()
mtext("Climate-integration structure (all 5 pools)", side=1, line=2.9, cex=0.95)
mtext(expression("Predictive "*R^2), side=2, line=2.9, cex=0.95)
lines(yx, calib[yi], col="#e65100", lwd=1.6)
points(yx, hold[yi],  pch=21, bg="white",   col="#e65100", cex=1.6, lwd=1.8)
points(yx, calib[yi], pch=21, bg="#e65100",  col="#e65100", cex=1.8, lwd=1.2)
segments(yx, hold[yi], yx, calib[yi], col=adjustcolor("#e65100",0.5), lwd=1.2)
text(yx, calib[yi], sprintf("%.3f", calib[yi]), pos=3, offset=0.6, cex=0.72, col="#e65100")
title(main="B · Climate axis", cex.main=0.95, font.main=1)

dev.off()
cat("Wrote manuscript/figures/F6b_skill_vs_complexity.png\n")
cat(sprintf("Calib slope vs pools: %+.5f R^2/pool\n", coef(lm(calib~pools))[2]))
