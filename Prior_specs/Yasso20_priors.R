# Yasso20 prior specification — HIKET pipeline
# To regenerate: set WRITE_PRIORS <- TRUE and run Priors_model_matching.R.
# Prior centres (transfer fractions): same as Yasso15 (YASSO15_DEFAULT_PARAMS)
# Prior centres (climate/size): Yasso20_sample_parameters.rda MAP, FMI Ryassofortran
# Prior widths: posterior SDs, Yasso20.dat, FMI Ryassofortran
# Transfer fraction widths: logit SD 0.4 (Tier-2 common weak prior — Finnish
#   data drives structure; see Prior_specs/PRIOR_HOMOGENIZATION_PLAN.md §4.2)

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
  # P3 (2026-08-07): 0.90 = J_1917/J_1985 from the NFI growing-stock record with the
  # fitted litter-growing-stock elasticity (eps 0.45-0.66 over 1986-2023 => litter
  # scales ~ sqrt(growing stock); GS_1917/GS_1985 = 0.789 => R = 0.789^eps ~ 0.90).
  # Meaningful only together with P1 (common J_t0 anchor in the wrappers): before
  # P1, sigma_init = R * 0.818 and a centre of 1.00 silently asserted R = 1.22.
  sigma_init  = 0.90,
  # C4b/2026-08-17: centre DERIVED, not chosen. Tupek J_bar = 2.511 tC/ha/yr is TREE
  # litter only (understorey excluded); Lehtonen & Heikkinen 2015 give TOTAL litter
  # (tree + understorey) = 2.70 [2.43, 2.97]. Adding the missing understorey is exactly
  # this parameter's job => 2.70/2.511 = 1.08. The old 1.30 was directionally right but
  # had no source. Coupled to sigma_init through J_1917 -- do not move one alone.
  #
  # RECENTRED 1.08 -> 1.27, 2026-08-20 (decision: Lorenzo), on Liski et al. 2006
  # (Ann. For. Sci. 63:687-697), which is the closer benchmark: same country, same
  # period, same model (Yasso), same NFI source -- and A. Lehtonen and M. Peltoniemi
  # are coauthors. Their Fig. 7 (1990s, kg C/m2/yr -> tC/ha/yr) splits as tree litter
  # 1.58 + ground vegetation 0.61 + harvest residues 0.63 + natural mortality 0.06,
  # using their statement (S4.2) that ground vegetation is 28% of the litter production
  # of living vegetation. The subtotal COMPARABLE to our J_bar (tree + residues +
  # mortality, understorey excluded) is 2.27; their TOTAL input to soil is 2.88.
  #   => understorey correction = 2.88/2.27 = 1.27
  # ⚠ TWO defensible recentrings; we take the RATIO, not the absolute. Matching their
  # absolute total instead would give 2.88/2.511 = 1.15. They differ by the 11% gap
  # between the two litter products (our J_bar 2.511 vs their 2.27). The RATIO is the
  # right invariant here because sigma_input exists to correct OUR J_bar for the
  # understorey gap -- but it does accept our J_bar's level as given.
  # ⚠ Liski (1.27) and Lehtonen & Heikkinen (1.08) DISAGREE on the understorey share
  # (27% vs 8%). Candidate reason: whether L&H's "total litter" carries harvest
  # residues, which our J_bar does. Unresolved -- ask A. Lehtonen (author of both).
  # Independent support for the larger share: Liski S4.2 warns that ignoring ground
  # vegetation underestimates "not only these parameters but also the soil carbon
  # stock and sink" -- exactly this parameter's failure mode.
  # SENSITIVITY ARM B (2026-08-31): Lehtonen & Heikkinen anchor, 2.70/2.511 = 1.08.
  # Arm A (Liski, 1.27) is run 20260820_1554*. The two anchors disagree by 18% on the
  # understorey share and the disagreement is UNRESOLVED (see above) -- we report both
  # rather than choose. Revert to 1.27 to reproduce arm A.
  sigma_input = 1.08
)

# Published Viskari 2022 parameter vector — used by run_yasso20_baseline.R only.
# Fractions: from Yasso20_sample_parameters.rda MAP (Ryassofortran / FMI).
#   Structural zeros (SD=0 in posterior): p_NW, p_AE, p_WE, p_NE, p_AN, p_EN.
#   p_EA = 0 at MAP: algebraically derived fraction whose constraint evaluates to 0.
# Climate/size: identical to YASSO20_FREE_DEFAULTS (same rda source).
# Do NOT use for calibration — YASSO20_FREE_DEFAULTS is the correct prior centre.
YASSO20_PUBLISHED_DEFAULTS <- c(
  # Transfer fractions — rda MAP values
  p_WA = 0.500,
  p_EA = 0.000,
  p_NA = 1.000,
  p_AW = 1.000,
  p_EW = 0.990,
  p_NW = 0.000,
  p_AE = 0.000,
  p_WE = 0.000,
  p_NE = 0.000,
  p_AN = 0.000,
  p_WN = 0.163,
  p_EN = 0.000,
  # Climate/size — rda MAP (same as YASSO20_FREE_DEFAULTS)
  beta1  = 0.158000,
  beta2  = -0.002000,
  gamma  = -1.440000,
  betaN1 = 0.170000,
  betaN2 = -0.005000,
  gammaN = -2.000000,
  betaH1 = 0.067000,
  betaH2 =  0.000000,
  gammaH = -6.900000,
  delta1 = -2.550000,
  delta2 =  1.240000,
  r      =  0.250000,
  # Auxiliary uncertainty parameters
  sigma_init  = 1.00,
  sigma_input = 1.00
)

# sigma_ppm in unconstrained (transformed) space. All free params listed
# explicitly (decision #5: explicit per-fraction listing for traceability).
YASSO20_SIGMA_PPM <- c(
  # Transfer fractions (12) — Tier-2 common logit SD 0.4
  p_WA = 0.4, p_EA = 0.4, p_NA = 0.4, p_AW = 0.4,
  p_EW = 0.4, p_NW = 0.4, p_AE = 0.4, p_WE = 0.4,
  p_NE = 0.4, p_AN = 0.4, p_WN = 0.4, p_EN = 0.4,
  # Climate & size (posterior SDs, Yasso20.dat — Tier-1 unchanged)
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
  # Tier-3 width HALVED 0.50 -> 0.25 (2026-08-17), same value for both, preserving the
  # Tier-3 principle that the two auxiliary sigmas share a width.
  # WHY: the sigma_init CENTRE (0.90) was already well derived -- Korhonen et al. 2024
  # growing stock V(1917)/V(1985) = 0.789 raised to the litter-growing-stock elasticity
  # eps = 0.43 [0.23, 0.63], fitted on our own litter record => 0.789^0.43 = 0.90. The
  # 0.50 width was what let posteriors sit 2.3-5.0x BELOW that centre (Yasso15 0.182,
  # i.e. a 1917 soil equilibrated to 18% of 1985 litter, against a forest record saying
  # 79%). Width rebuilt from its parts: elasticity CI 0.03, volume-record uncertainty
  # 0.07, structural litter-model error 0.09 => 0.12 in quadrature, prudential x2 = 0.25.
  # Deliberately looser than the evidence supports so the sampler is PULLED, not WALLED.
  # Check with doublechecks/sigma_init_vs_growing_stock.R.
  sigma_init  = 0.25,
  sigma_input = 0.25
)

# --- Fraction-prior tightening switch (2026-08-10) ----------------------------
# HIKET_PRIOR_TIGHTEN=<f> multiplies the TRANSFER-FRACTION prior SDs by f, and
# nothing else. The Tier-2 logit SD of 0.4 is the ONLY width in the whole scheme
# chosen by us rather than derived from a source, so it is the only one we are
# entitled to narrow arbitrarily. Justification: the fractions are weakly
# identified (all within 1.6 sigma of prior) yet high-leverage -- a ~1 sigma move
# in p_WA/p_WN halves bulk MRT.
#
# Climate and woody-size widths are NOT scaled here: they come from sources and
# were corrected directly above (Tuomi tables) or are genuine posterior SDs
# (Yasso15/20 .dat). Rate anchors and the two auxiliary sigmas are untouched.
#
# Unset or 1 => production behaviour unchanged.
.hiket_tighten <- suppressWarnings(as.numeric(Sys.getenv("HIKET_PRIOR_TIGHTEN", "1")))
if (!is.finite(.hiket_tighten) || .hiket_tighten <= 0)
  stop("HIKET_PRIOR_TIGHTEN must be a positive number")
if (!isTRUE(all.equal(.hiket_tighten, 1))) {
  .fr <- grep("^p_", names(YASSO20_SIGMA_PPM), value = TRUE)
  if (length(.fr)) {
    YASSO20_SIGMA_PPM[.fr] <- YASSO20_SIGMA_PPM[.fr] * .hiket_tighten
    message(sprintf("[PRIOR TIGHTENING] Yasso20 fraction SDs x %.3f (%d fractions)",
                    .hiket_tighten, length(.fr)))
  } else {
    message("[PRIOR TIGHTENING] Yasso20 has no transfer fractions -- no effect")
  }
}

# Physical litter-flux envelope (tC/ha/yr), homogeneous across all six models.
# Bounds the EFFECTIVE flux sigma_input*J (and the 1917 flux) to the boreal NPP
# range (Gower et al. 2001); floor relaxed below the Gower min as a small,
# strictly-positive, non-binding guardrail. Enforced via the `flux_pair`
# transform in calibration_engine.R. See sigma_input_physical_bounds_note.
YASSO20_INPUT_FLUX_WINDOW <- c(0.05, 8.7)
