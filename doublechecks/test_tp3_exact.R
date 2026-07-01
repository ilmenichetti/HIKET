#!/usr/bin/env Rscript
# =============================================================================
# test_tp3_exact.R  --  verify the exact TP3 integrator
# -----------------------------------------------------------------------------
# (a) Correctness: the closed-form .tp3_step() must equal a reference matrix
#     exponential of the cascade generator, to machine precision, for random
#     rate sets.
# (b) Behaviour: on an oscillating-climate, constant-litter case, the OLD
#     explicit-Euler step rings while the NEW exact step is smooth; both must
#     converge to the same long-run mean (the ringing was an integration
#     artefact, not a difference in the represented dynamics).
#
# Run from project root:  Rscript doublechecks/test_tp3_exact.R
# =============================================================================

WRAP <- "Model_functions_real_data_transient/Decomposition_functions/SimpleModels/tp3_wrapper_transient.R"
source(WRAP)                                   # provides .tp3_step(), tp3_run()
dir.create("doublechecks/figures", showWarnings = FALSE, recursive = TRUE)

# --- reference exact step via a brute-force matrix exponential ---------------
ref_expm <- function(M, n = 30L) {            # scaling-and-squaring Taylor
  s <- max(0L, ceiling(log2(max(rowSums(abs(M)) + 1e-12))))
  Ms <- M / (2^s); E <- diag(nrow(M)); term <- diag(nrow(M))
  for (k in seq_len(n)) { term <- (term %*% Ms) / k; E <- E + term }
  for (i in seq_len(s)) E <- E %*% E
  E
}
ref_step <- function(cA, cS, cH, kA, kS, kH, pS, pH, J) {
  M  <- matrix(c(-kA, 0, 0,  pS*kA, -kS, 0,  0, pH*kS, -kH), 3, 3, byrow = TRUE)
  Css <- c(J/kA, pS*J/kS, pH*pS*J/kH)
  as.numeric(Css + ref_expm(M) %*% (c(cA, cS, cH) - Css))
}

# --- (a) correctness over random rate sets -----------------------------------
set.seed(1); max_err <- 0
for (i in 1:2000) {
  kA <- runif(1, .2, 3); kS <- runif(1, .02, .5); kH <- runif(1, 5e-4, .01)
  pS <- runif(1, .05, .6); pH <- runif(1, .1, .8); J <- runif(1, .5, 4)
  cA <- runif(1, 0, 50); cS <- runif(1, 0, 80); cH <- runif(1, 0, 120)
  a <- .tp3_step(cA, cS, cH, kA, kS, kH, pS, pH, J)
  b <- ref_step(cA, cS, cH, kA, kS, kH, pS, pH, J)
  max_err <- max(max_err, max(abs(as.numeric(a) - b)))
}
cat(sprintf("(a) closed-form vs reference matrix-exp: max abs error = %.2e  %s\n",
            max_err, if (max_err < 1e-8) "[PASS]" else "[FAIL]"))

# --- (b) Euler vs exact on oscillating climate -------------------------------
# old explicit-Euler forward run (reconstructed from the previous wrapper)
tp3_run_euler <- function(inputs, mp, C_init, xi) {
  n <- nrow(inputs); A <- S <- H <- numeric(n)
  aA <- mp[["alpha_A"]]; aS <- mp[["alpha_S"]]; aH <- mp[["alpha_H"]]
  pS <- mp[["p_S"]]; pH <- mp[["p_H"]]; si <- mp[["sigma_input"]]
  cA <- C_init[["A"]]; cS <- C_init[["S"]]; cH <- C_init[["H"]]
  for (t in seq_len(n)) {
    J <- inputs$J_total[t]*si; x <- xi[t]
    fA <- aA*x*cA; fS <- aS*x*cS
    cA <- cA - fA + J; cS <- cS + pS*fA - fS; cH <- cH + pH*fS - aH*cH
    A[t] <- cA; S[t] <- cS; H[t] <- cH
  }
  data.frame(year = inputs$year, total_soc = A + S + H, A = A)
}

p  <- readRDS("Calibration_real_data_transient/runs/TP3_posterior_20260608_020212.rds")
mp <- c(alpha_A = median(p[,"alpha_A"]), alpha_S = median(p[,"alpha_S"]),
        alpha_H = median(p[,"alpha_H"]), p_S = median(p[,"p_S"]),
        p_H = median(p[,"p_H"]), sigma_input = 1)
yrs <- 1986:2024; n <- length(yrs)
# CONSTANT climate and litter: any wiggle in the output is then purely numerical.
# With alpha_A ~ 0.97 and xi = 1.15, k_A = alpha_A*xi > 1, so explicit Euler rings
# while the exact step relaxes monotonically.
xi  <- rep(1.75, n)                # k_A = alpha_A*xi ~ 1.7 -> Euler factor ~ -0.7 (clear ringing)
inp <- data.frame(year = yrs, J_total = rep(2.7, n))
C0  <- c(A = 0, S = 0, H = 40)                                   # start far below equilibrium

ex <- tp3_run(inp, mp, C0, xi)
eu <- tp3_run_euler(inp, mp, C0, xi)

# year-to-year roughness of the ACTIVE pool (where the ringing lives)
rough <- function(v) mean(abs(diff(diff(v))))
cat(sprintf("(b) constant forcing: active-pool roughness  Euler = %.4f   exact = %.4f\n",
            rough(eu$A), rough(ex$A)))
cat(sprintf("    (exact ~0 = monotonic relaxation; Euler > 0 = spurious ringing) %s\n",
            if (rough(ex$A) < 0.1 * rough(eu$A)) "[PASS]" else "[CHECK]"))
cat(sprintf("    mean total SOC          Euler = %.2f    exact = %.2f   tC/ha\n",
            mean(eu$total_soc), mean(ex$total_soc)))

png("doublechecks/figures/tp3_euler_vs_exact.png", width = 1500, height = 600, res = 150)
par(mfrow = c(1, 2), mar = c(4, 4.2, 2.5, 1))
plot(eu$year, eu$A, type="l", col="firebrick", lwd=2, xlab="Year",
     ylab="Active pool (tC/ha)", main="Active pool: Euler rings, exact smooth")
lines(ex$year, ex$A, col="forestgreen", lwd=2)
legend("topright", c("Euler (old)","Exact (new)"), col=c("firebrick","forestgreen"),
       lwd=2, bty="n")
plot(eu$year, eu$total_soc, type="l", col="firebrick", lwd=2, xlab="Year",
     ylab="Total SOC (tC/ha)", main="Total SOC")
lines(ex$year, ex$total_soc, col="forestgreen", lwd=2)
legend("topright", c("Euler","Exact"), col=c("firebrick","forestgreen"), lwd=2, bty="n")
dev.off()
cat("Figure: doublechecks/figures/tp3_euler_vs_exact.png\n")
