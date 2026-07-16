# =============================================================================
# icbm_anchor_sanity.R  (C1 pre-flight)
#
# Question (user's caveat): do the ICBM-anchored rates carry REAL, boreal-
# sensible information about alpha, or are we just swapping one weak point
# estimate for another? Test BEFORE editing any prior / spending a recalibration.
#
# Checks (local, no HPC), all at real Finnish plot climate + litter (site_raw):
#   1. Reference & effective bulk MRT — is boreal turnover in the right range?
#   2. Steady-state stock at real J with sigma_input = 1 — does it land in the
#      observed ~55-90 tC/ha window WITHOUT an inflated input multiplier?
#      (If stocks come out far too low at sigma_input=1, the slow rate is too
#       fast and the ICBM anchor is suspect for boreal forest.)
#   3. Implied sigma_input to match observed SOC (should be ~1 if anchor good).
#   4. k1 fast-pool cross-check vs half-life expectation.
# =============================================================================

# --- xi function (Yasso07 form; identical to the wrappers) -------------------
compute_xi <- function(temp_mean, temp_amp, precip, beta1, beta2, gamma) {
  sq2 <- sqrt(2)
  T1 <- temp_mean + 4 * temp_amp / pi * (1 / sq2 - 1)
  T2 <- temp_mean - 4 * temp_amp / (sq2 * pi)
  T3 <- temp_mean + 4 * temp_amp / pi * (1 - 1 / sq2)
  T4 <- temp_mean + 4 * temp_amp / (sq2 * pi)
  temp_mod <- (exp(beta1*T1 + beta2*T1^2) + exp(beta1*T2 + beta2*T2^2) +
               exp(beta1*T3 + beta2*T3^2) + exp(beta1*T4 + beta2*T4^2)) / 4
  temp_mod * (1 - exp(gamma * precip / 1000))
}

# --- Constants ---------------------------------------------------------------
# ICBM (Andren & Katterer 1997, Ultuna bare fallow, r = 1 = central Sweden)
k1 <- 0.8;  k2 <- 0.00605;  h <- 0.13
# Simple-model climate centres (SP1/TP2/TP3 *_FREE_DEFAULTS)
b1 <- 0.095;  b2 <- -0.00014;  g <- -1.21
# Ultuna reference climate (half-range T_amp convention, matching the plots)
ULT_T <- 5.4;  ULT_P <- 520
# Humification centre (re-centred, kept FREE in calibration)
pH <- 0.13;  pS <- 0.13
# TP3 slow rate: ICBM is 2-timescale (fast k1, slow k2) with NO intermediate.
# S and H TOGETHER = ICBM's "old" subsystem, so alpha_S must be SLOW (k2 scale),
# not intermediate. Test around k2; include the plan's subsystem identity
# alpha_S = k2/(1 - p_H) = 0.00605/0.87 = 0.00695 (S+H subsystem MRT pinned 165 yr).
alphaS_grid <- c(k2/(1 - pH), 0.006, 0.010, 0.020)

# --- xi_Ultuna: sensitivity to the T_amp convention --------------------------
cat("=== xi_Ultuna (T=5.4, P=520) vs T_amp convention ===\n")
for (amp in c(9, 10, 11, 12)) {
  cat(sprintf("  T_amp=%2d -> xi_Ultuna = %.4f\n", amp, compute_xi(ULT_T, amp, ULT_P, b1, b2, g)))
}
xi_U <- compute_xi(ULT_T, 10, ULT_P, b1, b2, g)   # adopt half-range ~10
cat(sprintf("  adopted xi_Ultuna = %.4f (T_amp=10)\n\n", xi_U))

# --- Fixed rate constants (alpha = k_ICBM / xi_Ultuna) -----------------------
alpha_SP1 <- (1/(1/k1 + h/k2)) / xi_U     # SP1 single pool = ICBM bulk MRT
alpha_A   <- k1 / xi_U
alpha_H   <- k2 / xi_U
cat("=== Fixed rate constants (reference-climate, /xi_Ultuna) ===\n")
cat(sprintf("  SP1 alpha   = %.4f  (bulk MRT ref = %.1f yr)\n", alpha_SP1, 1/k1 + h/k2))
cat(sprintf("  TP2/3 alpha_A = %.4f, alpha_H = %.5f\n\n", alpha_A, alpha_H))

# --- Real plot climate + litter ----------------------------------------------
d  <- read.csv("Data/model_inputs/site_raw.csv")
d  <- d[d$calib_ready & !is.na(d$mean_soc_Mgha) & !is.na(d$mean_litter), ]
Tm <- d$mean_temp
Ta <- (d$warmest_month_T - d$coldest_month_T) / 2
Pr <- d$mean_precip_annual
J  <- d$mean_litter                       # tC/ha/yr, raw (sigma_input = 1)
SOCobs <- d$mean_soc_Mgha
xi <- compute_xi(Tm, Ta, Pr, b1, b2, g)
cat(sprintf("Plots: %d | xi_FI median %.3f (IQR %.3f-%.3f) | xi_FI/xi_U median %.3f\n",
            nrow(d), median(xi), quantile(xi,.25), quantile(xi,.75), median(xi)/xi_U))
cat(sprintf("Observed SOC median %.1f (IQR %.1f-%.1f) | J median %.2f tC/ha/yr\n\n",
            median(SOCobs), quantile(SOCobs,.25), quantile(SOCobs,.75), median(J)))

# --- Steady-state stock at sigma_input = 1, per plot -------------------------
# effective rate k_X = alpha_X * xi_FI.  Report with AND without xi on the
# humus pool (current convention = no xi on H; C2 would add it).
report <- function(name, MRT_vec) {
  stock <- J * MRT_vec                       # C_ss = J * bulk_MRT
  simp  <- SOCobs / stock                     # sigma_input needed to match obs
  cat(sprintf("%-22s stock med %6.1f (IQR %5.1f-%5.1f) | MRT med %5.1f yr | sigma_input med %.2f\n",
              name, median(stock), quantile(stock,.25), quantile(stock,.75),
              median(MRT_vec), median(simp)))
}
cat("=== Steady-state SOC at sigma_input = 1 (vs observed median", round(median(SOCobs),1), ") ===\n")

kA <- alpha_A * xi;  kH_x <- alpha_H * xi;  kH_nox <- alpha_H  # with / without xi on H
kSP1 <- alpha_SP1 * xi

report("SP1", 1/kSP1)
report("TP2 (H has xi)",   1/kA + pH/kH_x)
report("TP2 (H no xi)",    1/kA + pH/kH_nox)
for (aS in alphaS_grid) {
  kS <- aS/xi_U * xi
  report(sprintf("TP3 aS=%.4f (H xi)",  aS), 1/kA + pS/kS + pS*pH/kH_x)
  report(sprintf("TP3 aS=%.4f (H noxi)", aS), 1/kA + pS/kS + pS*pH/kH_nox)
}

cat(sprintf("\n=== Fast-pool cross-check: k1 = %.2f -> active half-life %.2f yr ===\n",
            k1, log(2)/k1))
