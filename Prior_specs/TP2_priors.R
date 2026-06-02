# TP2 prior specification — HIKET pipeline
# Hand-set defaults: no published global calibration exists.
# ICBM-style starting values (Andren & Katterer 1997); climate copied from
# Yasso07 MAP. All widths weakly informative — Finnish data dominates.

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
  p_H         = 1.00,
  beta1       = 0.20,
  beta2       = 0.05,
  gamma       = 0.30,
  sigma_init  = 0.50,
  sigma_input = 0.50
)
