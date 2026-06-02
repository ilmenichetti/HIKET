# Yasso15 prior specification — HIKET pipeline
# To regenerate: set WRITE_PRIORS <- TRUE and run Priors_model_matching.R.
# Prior centres (transfer fractions): YASSO15_DEFAULT_PARAMS, FMI Ryassofortran
# Prior centres (climate/size): posterior means, Yasso15.dat, FMI Ryassofortran
# Prior widths: posterior SDs, Yasso15.dat; gammaH capped at 1.5
# Transfer fraction widths: 1.0 (weakly informative — Finnish data drives structure)

# Complete physical-space prior centres for all free parameters.
YASSO15_FREE_DEFAULTS <- c(
  # Transfer fractions (12)
  p_WA = 0.43628932,
  p_EA = 0.24997402,
  p_NA = 0.91512685,
  p_AW = 0.99258227,
  p_EW = 0.083853738,
  p_NW = 0.011476783,
  p_AE = 6.08e-04,
  p_WE = 4.76e-04,
  p_NE = 0.066037729,
  p_AN = 7.71e-04,
  p_WN = 0.10401742,
  p_EN = 0.64880756,
  # AWE climate response
  beta1  = 0.09062,
  beta2  = -0.000215,
  gamma  = -1.80897,
  # N climate response
  betaN1 = 0.04878,
  betaN2 = -0.0000792,
  gammaN = -1.17294,
  # H climate response
  betaH1 = 0.03518,
  betaH2 = -0.000208,
  gammaH = -12.54094,
  # Woody size modifier
  delta1 = -0.43883,
  delta2 = 1.26838,
  r      = 0.25687,
  # Auxiliary uncertainty parameters
  sigma_init  = 1.00,
  sigma_input = 1.00
)

# sigma_ppm overrides: only non-1.0 values listed.
# Transfer fractions default to 1.0 (weakly informative).
# gammaH capped at 1.5: near-unidentifiable at Finnish precipitation levels.
YASSO15_SIGMA_PPM <- c(
  beta1       = 0.04593,
  beta2       = 0.00014,
  gamma       = 0.07022,
  betaN1      = 0.10502,
  betaN2      = 0.00008,
  gammaN      = 0.15543,
  betaH1      = 0.13477,
  betaH2      = 0.00016,
  gammaH      = 1.50000,
  delta1      = 0.32810,
  delta2      = 0.28643,
  r           = 0.05504,
  sigma_init  = 0.50,
  sigma_input = 0.50
)
