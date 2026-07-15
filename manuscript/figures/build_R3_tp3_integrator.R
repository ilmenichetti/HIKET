setwd("/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling")
# Robustness R3 -- TP3 integrator check. Two panels:
#  (A) constant forcing: the OLD explicit-Euler cascade RINGS when alpha_A*xi>1
#      (numerical overshoot), while the exact matrix-exp step is monotone -> the
#      artefact that was removed.
#  (B) the real calibrated TP3 mean trajectory (exact integrator) still shows a
#      BOUNDED year-to-year oscillation -> that residual is PHYSICAL (fast pool
#      tracking climate), not the Euler artefact. Documents [[tp3-forcing-oscillation]].
suppressMessages({
  source("Calibration_real_data_transient/calibration_engine_transient.R")
  source("Model_functions_real_data_transient/Decomposition_functions/SimpleModels/tp3_wrapper_transient.R")
})
RID <- "20260710_104904"
post <- readRDS(sprintf("Calibration_real_data_transient/runs/TP3_posterior_%s.rds", RID))
mp <- apply(BayesianTools::getSample(post), 2, median)
aA <- unname(mp["alpha_A"]); aS <- unname(mp["alpha_S"]); aH <- unname(mp["alpha_H"])
pS <- unname(mp["p_S"]); pH <- unname(mp["p_H"])

# (A) constant forcing, xi set so alpha_A*xi = 1.6 (the posterior-median regime, >1)
xi <- 1.6 / aA; J <- 2.5; n <- 40L
kA <- aA*xi; kS <- aS*xi; kH <- aH
C0 <- c(A=0, S=0, H=0)   # start depleted: relaxation to steady state reveals Euler ringing
ex <- eu <- matrix(NA, n, 3)
cA<-cS<-cH<-NULL
Ce <- C0
Cu <- C0
for (t in 1:n) {
  Ce <- .tp3_step(Ce["A"], Ce["S"], Ce["H"], kA, kS, kH, pS, pH, J); ex[t,] <- Ce
  nA <- Cu["A"] + J - kA*Cu["A"]; nS <- Cu["S"] + pS*kA*Cu["A"] - kS*Cu["S"]
  nH <- Cu["H"] + pH*kS*Cu["S"] - kH*Cu["H"]; Cu <- c(A=unname(nA),S=unname(nS),H=unname(nH)); eu[t,] <- Cu
}
ex_tot <- rowSums(ex); eu_tot <- rowSums(eu)

# (B) real calibrated mean trajectory (exact integrator, from the predictive bundle)
pp <- readRDS(sprintf("Calibration_real_data_transient/runs/TP3_posterior_predictive_%s.rds", RID))$posterior_summary
tr <- aggregate(soc_mean ~ year, pp, mean)

png("manuscript/figures/R3_tp3_integrator.png", width=11, height=4.6, units="in", res=200)
par(mfrow=c(1,2), mar=c(4.4,4.6,3,1), mgp=c(2.6,0.7,0), las=1)
plot(1:n, eu_tot, type="l", col="#d01c8b", lwd=2.2, ylim=range(ex_tot,eu_tot),
     xlab="Year (constant forcing)", ylab="Total SOC (tC/ha)",
     main=sprintf("A . Constant forcing: Euler rings, exact is monotone\n(alpha_A*xi = 1.6)"), font.main=1, cex.main=0.95)
lines(1:n, ex_tot, col="#1a9850", lwd=2.4)
legend("topright", c("explicit Euler (old)","exact matrix-exp (now)"), col=c("#d01c8b","#1a9850"), lwd=2.3, bty="n", cex=0.85)
plot(tr$year, tr$soc_mean, type="l", col="#C26B51", lwd=2, xlab="Year", ylab="Mean SOC across plots (tC/ha)",
     main="B . Real calibrated TP3 (exact): bounded oscillation survives", font.main=1, cex.main=0.95)
points(tr$year, tr$soc_mean, pch=19, cex=0.5, col="#C26B51")
dev.off()
cat(sprintf("Euler roughness %.3f vs exact %.4f (constant forcing); real-traj roughness %.3f\n",
            mean(abs(diff(diff(eu_tot)))), mean(abs(diff(diff(ex_tot)))), mean(abs(diff(diff(tr$soc_mean))))))
cat("wrote R3_tp3_integrator.png\n")
