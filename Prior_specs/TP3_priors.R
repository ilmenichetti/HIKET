# TP3 prior specification — HIKET pipeline
# Hand-set defaults: no published global calibration exists.
# Extends TP2 with a Slow pool; alpha_S intermediate between alpha_A and
# alpha_H; p_S matches TP2's p_H; p_H neutral at 0.50.
# Climate widths on the Yasso07 empirical scale (Tuomi 2009 T3, 1σ); p_S/p_H on
# the Tier-2 logit SD 0.4 (see PRIOR_HOMOGENIZATION_PLAN.md).
# All other widths weakly informative — Finnish data dominates.

TP3_FREE_DEFAULTS <- c(
  alpha_A     = 0.73,
  alpha_S     = 0.10,
  alpha_H     = 0.0015,
  p_S         = 0.028,
  p_H         = 0.50,
  beta1       = 0.095,
  beta2       = -0.00014,
  gamma       = -1.21,
  sigma_init  = 1.00,
  sigma_input = 1.00
)

# sigma_ppm: prior SDs in unconstrained (transformed) space.
TP3_SIGMA_PPM <- c(
  alpha_A     = 0.50,
  alpha_S     = 0.50,
  alpha_H     = 0.50,
  p_S         = 0.40,     # Tier-2 logit (was 1.00)
  p_H         = 0.40,     # Tier-2 logit (was 1.00)
  beta1       = 0.26,     # Yasso07 scale (was 0.20)
  beta2       = 0.00065,  # Yasso07 scale (was 0.05 — explosive)
  gamma       = 0.20,     # Yasso07 scale (was 0.30)
  sigma_init  = 0.50,
  sigma_input = 0.50
)
