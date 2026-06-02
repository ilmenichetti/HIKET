# Yasso20 prior specification — HIKET pipeline
# To regenerate: set WRITE_PRIORS <- TRUE and run Priors_model_matching.R.
# Prior centres (transfer fractions): same as Yasso15 (YASSO15_DEFAULT_PARAMS)
# Prior centres (climate/size): Yasso20_sample_parameters.rda MAP, FMI Ryassofortran
# Prior widths: posterior SDs, Yasso20.dat, FMI Ryassofortran
# Transfer fraction widths: 1.0 (weakly informative — Finnish data drives structure)

# Complete physical-space prior centres for all free parameters.
YASSO20_FREE_DEFAULTS <- c(
  # Transfer fractions (12) — same defaults as Yasso15
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
  beta1  = 0.158000,
  beta2  = -0.002000,
  gamma  = -1.440000,
  # N climate response
  betaN1 = 0.170000,
  betaN2 = -0.005000,
  gammaN = -2.000000,
  # H climate response
  betaH1 = 0.067000,
  betaH2 =  0.000000,
  gammaH = -6.900000,
  # Woody size modifier
  delta1 = -2.550000,
  delta2 =  1.240000,
  r      =  0.250000,
  # Auxiliary uncertainty parameters
  sigma_init  = 1.00,
  sigma_input = 1.00
)

# sigma_ppm overrides: only non-1.0 values listed.
# Transfer fractions default to 1.0 (weakly informative).
YASSO20_SIGMA_PPM <- c(
  beta1       = 0.10355,
  beta2       = 0.00054,
  gamma       = 0.14583,
  betaN1      = 0.06246,
  betaN2      = 0.00041,
  gammaN      = 0.03532,
  betaH1      = 0.09118,
  betaH2      = 0.00009,
  gammaH      = 1.39056,
  delta1      = 0.43530,
  delta2      = 0.21098,
  r           = 0.05446,
  sigma_init  = 0.50,
  sigma_input = 0.50
)
