# SP1 prior specification — HIKET pipeline
# Hand-set defaults: no published global calibration exists.
# Climate parameters (beta1/beta2/gamma) copied from Yasso07 MAP as
# biologically plausible starting values for the Finnish boreal context.
# All widths weakly informative — Finnish data dominates.

SP1_FREE_DEFAULTS <- c(
  alpha       = 0.09,
  beta1       = 0.095,
  beta2       = -0.00014,
  gamma       = -1.21,
  sigma_init  = 1.00,
  sigma_input = 1.00
)

# sigma_ppm: prior SDs in unconstrained (transformed) space.
SP1_SIGMA_PPM <- c(
  alpha       = 0.50,
  beta1       = 0.20,
  beta2       = 0.05,
  gamma       = 0.30,
  sigma_init  = 0.50,
  sigma_input = 0.50
)
