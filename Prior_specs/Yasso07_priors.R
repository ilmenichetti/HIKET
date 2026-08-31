# Yasso07 prior specification — HIKET pipeline
# To regenerate: run Priors_model_matching.R interactively (Yasso07 is not
# auto-written; sigma_ppm values are from published posterior limits, not a
# posterior file).
# Prior centres (transfer fractions): YASSO07_DEFAULT_PARAMS (Tuomi et al. 2009)
# Prior centres (climate/size): GUI MAP, y07par_gui.csv (Tuomi et al. 2011 EMS)
# Prior widths (climate): Tuomi et al. 2009 (Ecol. Modelling) Table 3
# Prior widths (woody size): Tuomi et al. 2011 (Ecol. Modelling) Table 4
#   Published "±" limits read as 1σ (conservative — not divided by 1.96); for
#   log-transform params the relative SD = (1σ limit)/(paper MAP) via delta method.
# Transfer fraction widths: logit SD 0.4 (Tier-2 common weak prior — Finnish
#   data drives structure; see Prior_specs/PRIOR_HOMOGENIZATION_PLAN.md §4.2)

# Complete physical-space prior centres for all free parameters.
# r stored as positive (abs() applied; Tuomi reference stores it negative).
YASSO07_FREE_DEFAULTS <- c(
  # Transfer fractions (12)
  p_WA = 0.4888527989,
  p_EA = 0.01905768365,
  p_NA = 0.9696374536,
  p_AW = 0.9872559905,
  p_EW = 0.002843263559,
  p_NW = 0.003396461252,
  p_AE = 1.399370376e-05,
  p_WE = 1.796692413e-05,
  p_NE = 0.01218125969,
  p_AN = 0.002777846763,
  p_WN = 0.01269555371,
  p_EN = 0.9713827968,
  # Climate response
  beta1  = 0.09873183817,
  beta2  = -0.001571640489,
  gamma  = -1.271691799,
  # Woody size modifier
  delta1 = -1.708411336,
  delta2 = 0.8585553765,
  r      = 0.3068014085,
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

# sigma_ppm in unconstrained (transformed) space. All free params listed
# explicitly (decision #5: explicit per-fraction listing for traceability).
# CORRECTED 2026-08-10: both tables state '95% confidence limits', but the
# locked convention read '+-' as 1 sigma -- ~1.96x too wide. Widths below are
# the published half-widths divided by 1.959964. Centres unchanged.
# Climate/woody widths from Tuomi 2009 T3 / 2011 T4, "±" read as 1σ (§4.1).
YASSO07_SIGMA_PPM <- c(
  # Transfer fractions (12) — Tier-2 common logit SD 0.4
  p_WA = 0.4, p_EA = 0.4, p_NA = 0.4, p_AW = 0.4,
  p_EW = 0.4, p_NW = 0.4, p_AE = 0.4, p_WE = 0.4,
  p_NE = 0.4, p_AN = 0.4, p_WN = 0.4, p_EN = 0.4,
  # Climate response (Tuomi 2009 Table 3)
  beta1       = 0.1326555,     # log; rel SD 0.020/0.076 (paper MAP)
  beta2       = 0.0003316,  # unconstrained; T3 −8.9 ±6.5 ×10⁻⁴
  gamma       = 0.10204270,     # unconstrained; T3 −1.27 ±0.20
  # Woody size modifier (Tuomi 2011 Table 4)
  delta1      = 0.0816342,     # unconstrained; T4 −1.71 ±0.16
  delta2      = 0.0612256,     # log; rel SD 0.10/0.86
  r           = 0.021429,    # log; rel SD 0.013/0.306
  # Auxiliary uncertainty
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
  .fr <- grep("^p_", names(YASSO07_SIGMA_PPM), value = TRUE)
  if (length(.fr)) {
    YASSO07_SIGMA_PPM[.fr] <- YASSO07_SIGMA_PPM[.fr] * .hiket_tighten
    message(sprintf("[PRIOR TIGHTENING] Yasso07 fraction SDs x %.3f (%d fractions)",
                    .hiket_tighten, length(.fr)))
  } else {
    message("[PRIOR TIGHTENING] Yasso07 has no transfer fractions -- no effect")
  }
}

# Physical litter-flux envelope (tC/ha/yr), homogeneous across all six models.
# Bounds the EFFECTIVE flux sigma_input*J (and the 1917 flux) to the boreal NPP
# range (Gower et al. 2001); floor relaxed below the Gower min as a small,
# strictly-positive, non-binding guardrail. Enforced via the `flux_pair`
# transform in calibration_engine.R. See sigma_input_physical_bounds_note.
YASSO07_INPUT_FLUX_WINDOW <- c(0.05, 8.7)
