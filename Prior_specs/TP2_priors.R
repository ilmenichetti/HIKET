# TP2 prior specification — HIKET pipeline
# Hand-set defaults: no published global calibration exists.
# ICBM-style starting values (Andren & Katterer 1997); climate copied from
# Yasso07 MAP. Climate widths on the Yasso07 empirical scale (Tuomi 2009 T3, 1σ);
# p_H on the Tier-2 logit SD 0.4 (see PRIOR_HOMOGENIZATION_PLAN.md).
# All other widths weakly informative — Finnish data dominates.

TP2_FREE_DEFAULTS <- c(
  alpha_A     = 0.73,
  alpha_H     = 0.0015,
  p_H         = 0.028,
  beta1       = 0.095,
  beta2       = -0.00014,
  gamma       = -1.21,
  sigma_init  = 1.00,
  sigma_input = 1.00
)

# sigma_ppm: prior SDs in unconstrained (transformed) space.
TP2_SIGMA_PPM <- c(
  alpha_A     = 0.50,
  alpha_H     = 0.50,
  p_H         = 0.40,     # Tier-2 logit (was 1.00)
  beta1       = 0.26,     # Yasso07 scale (was 0.20)
  beta2       = 0.00065,  # Yasso07 scale (was 0.05 — explosive)
  gamma       = 0.20,     # Yasso07 scale (was 0.30)
  sigma_init  = 0.50,
  sigma_input = 0.50
)
