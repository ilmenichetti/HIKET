# =============================================================================
# test_tp3_c2_xi_on_H.R  (C2 verification)
#
# After C2 (k_H = alpha_H * xi on the humus pool), confirm the closed-form
# .tp3_step still equals a reference matrix-exponential step, and that the
# near-coincident-eigenvalue guard is safe now that C1 puts alpha_S and alpha_H
# both at k2-scale (kS ~ kH).
# =============================================================================

source("Model_functions_real_data_transient/Decomposition_functions/SimpleModels/tp3_wrapper_transient.R")

# Reference: exact one-year step of dC/dt = M C + b via augmented-matrix expm.
# State [A,S,H]; b = [J,0,0]. M lower-triangular cascade generator.
ref_step <- function(cA, cS, cH, kA, kS, kH, pS, pH, J) {
  M <- matrix(c(-kA,    0,    0,
                 pS*kA, -kS,   0,
                 0,      pH*kS, -kH), nrow = 3, byrow = TRUE)
  b <- c(J, 0, 0)
  Aug <- rbind(cbind(M, b), c(0,0,0,0))          # 4x4 augmented
  E <- as.matrix(Matrix::expm(Aug))
  step <- E[1:3, 1:3] %*% c(cA, cS, cH) + E[1:3, 4]
  c(A = step[1], S = step[2], H = step[3])
}

check <- function(label, kA, kS, kH, pS, pH, J, C0 = c(A=10, S=20, H=40)) {
  cf  <- .tp3_step(C0["A"], C0["S"], C0["H"], kA, kS, kH, pS, pH, J)
  rf  <- ref_step(C0["A"], C0["S"], C0["H"], kA, kS, kH, pS, pH, J)
  err <- max(abs(unname(unlist(cf)) - unname(rf)))
  cat(sprintf("%-34s closed-form vs matrix-exp max abs err = %.2e  %s\n",
              label, err, if (err < 1e-9) "OK" else "*** FAIL ***"))
  invisible(err)
}

# --- C1-era boreal params, xi ~ 0.9 (H now carries xi) ---
xi <- 0.90
aA <- 0.851; aS <- 0.0074; aH <- 0.00644; pS <- 0.13; pH <- 0.13; J <- 2.5
cat("=== C2: kH = alpha_H * xi ===\n")
check("C1 centres, xi=0.90",      aA*xi, aS*xi, aH*xi, pS, pH, J)
check("cold xi=0.55",             aA*0.55, aS*0.55, aH*0.55, pS, pH, J)
check("kS ~ kH (aS=aH=0.0065)",   aA*xi, 0.0065*xi, 0.0065*xi, pS, pH, J)   # near-coincident
check("kS == kH exactly",         aA*xi, 0.0065*xi, 0.0065*xi + 0, pS, pH, J)
check("high humification pH=0.4", aA*xi, aS*xi, aH*xi, 0.4, 0.4, J)

# --- Steady-state convergence: iterate under constant forcing, compare to the
#     analytical H_ss = pH*pS*J/(alpha_H*xi) (now with xi on H). ---
cat("\n=== Steady-state convergence (constant forcing, 4000 yr) ===\n")
kA <- aA*xi; kS <- aS*xi; kH <- aH*xi
C  <- c(A=0, S=0, H=0)
for (i in 1:4000) C <- .tp3_step(C["A"], C["S"], C["H"], kA, kS, kH, pS, pH, J)
Ass <- J/kA; Sss <- pS*J/kS; Hss <- pH*pS*J/kH   # analytical with xi on H
cat(sprintf("  A: iter %.4f vs analytic %.4f\n", C["A"], Ass))
cat(sprintf("  S: iter %.4f vs analytic %.4f\n", C["S"], Sss))
cat(sprintf("  H: iter %.4f vs analytic %.4f  (xi-dependent now)\n", C["H"], Hss))
cat(sprintf("  max rel err = %.2e\n",
            max(abs(c(C["A"]-Ass, C["S"]-Sss, C["H"]-Hss)/c(Ass,Sss,Hss)))))
