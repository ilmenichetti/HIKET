# Yasso15 prior specification — HIKET pipeline
# To regenerate: set WRITE_PRIORS <- TRUE and run Priors_model_matching.R.
# Prior centres (transfer fractions): YASSO15_DEFAULT_PARAMS, FMI Ryassofortran
# Prior centres (climate/size): posterior means, Yasso15.dat, FMI Ryassofortran
# Prior widths: posterior SDs, Yasso15.dat; gammaH capped at 1.5
# Transfer fraction widths: logit SD 0.4 (Tier-2 common weak prior — Finnish
#   data drives structure; see Prior_specs/PRIOR_HOMOGENIZATION_PLAN.md §4.2)

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
  # P3 (2026-08-07): 0.90 = J_1917/J_1985 from the NFI growing-stock record with the
  # fitted litter-growing-stock elasticity (eps 0.45-0.66 over 1986-2023 => litter
  # scales ~ sqrt(growing stock); GS_1917/GS_1985 = 0.789 => R = 0.789^eps ~ 0.90).
  # Meaningful only together with P1 (common J_t0 anchor in the wrappers): before
  # P1, sigma_init = R * 0.818 and a centre of 1.00 silently asserted R = 1.22.
  sigma_init  = 0.90,
  sigma_input = 1.30   # C4b: re-centred >1 for missing (understorey-dominated) litter (D2)
)

# sigma_ppm in unconstrained (transformed) space. All free params listed
# explicitly (decision #5: explicit per-fraction listing for traceability).
# gammaH capped at 1.5: near-unidentifiable at Finnish precipitation levels.
YASSO15_SIGMA_PPM <- c(
  # Transfer fractions (12) — Tier-2 common logit SD 0.4
  p_WA = 0.4, p_EA = 0.4, p_NA = 0.4, p_AW = 0.4,
  p_EW = 0.4, p_NW = 0.4, p_AE = 0.4, p_WE = 0.4,
  p_NE = 0.4, p_AN = 0.4, p_WN = 0.4, p_EN = 0.4,
  # Climate & size (posterior SDs, Yasso15.dat — Tier-1 unchanged)
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

# Physical litter-flux envelope (tC/ha/yr), homogeneous across all six models.
# Bounds the EFFECTIVE flux sigma_input*J (and the 1917 flux) to the boreal NPP
# range (Gower et al. 2001); floor relaxed below the Gower min as a small,
# strictly-positive, non-binding guardrail. Enforced via the `flux_pair`
# transform in calibration_engine.R. See sigma_input_physical_bounds_note.
YASSO15_INPUT_FLUX_WINDOW <- c(0.05, 8.7)
