# TP3 prior specification — HIKET pipeline
# Hand-set defaults: no published global calibration exists.
# Extends TP2 with a Slow pool. NOTE: alpha_S is NOT intermediate between alpha_A
# and alpha_H, and p_H is NOT neutral at 0.50 — both were superseded by C1 below
# (an intermediate alpha_S starves the cascade; see the C1 note).
# Climate widths on the Yasso07 empirical scale (Tuomi 2009 T3, 1σ); p_S/p_H on
# the Tier-2 logit SD 0.4 (see PRIOR_HOMOGENIZATION_PLAN.md).
# All other widths weakly informative — Finnish data dominates.

# C1 (2026-07-16): ICBM kinetic anchor, homogeneous with Yasso (fixes rates, frees
# fractions). alpha_A (fast) is FIXED — NOT in this free vector; injected as a
# constant in assemble_model_params (= k1/xi_Ultuna = 0.851). alpha_S and alpha_H
# are BOTH slow (k2-scale = ICBM "old" subsystem: S+H together represent ICBM's old
# pool) and stay free with VERY-INFORMATIVE priors. CRITICAL: alpha_S must be
# k2-scale (~0.0074), NOT intermediate — an intermediate value starves the cascade
# (verified, icbm_anchor_sanity.R). p_S, p_H stay FREE (Tier-2 logit), re-centred at
# ICBM h=0.13. ICBM: Andren & Katterer 1997. See NEXT_SESSION.md §5 (C1).
TP3_FREE_DEFAULTS <- c(
  alpha_S     = 0.0074,   # k2/(1-p_H)/xi_Ultuna, k2-scale (was 0.10 = cascade-starving)
  alpha_H     = 0.00644,  # k2/xi_Ultuna (was 0.0015)
  p_S         = 0.13,     # ICBM humification h (was 0.028); FREE, logit SD 0.4
  p_H         = 0.13,     # ICBM humification h (was 0.50); FREE, logit SD 0.4
  beta1       = 0.095,
  beta2       = -0.00014,
  gamma       = -1.21,
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

# sigma_ppm: prior SDs in unconstrained (transformed) space.
# alpha_A omitted — FIXED constant (injected in assemble_model_params).
TP3_SIGMA_PPM <- c(
  alpha_S     = 0.15,     # C1: VERY-INFORMATIVE (was 0.50)
  alpha_H     = 0.15,     # C1: VERY-INFORMATIVE (was 0.50)
  p_S         = 0.40,     # Tier-2 logit (unchanged)
  p_H         = 0.40,     # Tier-2 logit (unchanged)
  beta1       = 0.1326555,     # Yasso07 scale (was 0.20)
  beta2       = 0.0003316,  # Yasso07 scale (was 0.05 — explosive)
  gamma       = 0.10204270,     # Yasso07 scale (was 0.30)
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
  .fr <- grep("^p_", names(TP3_SIGMA_PPM), value = TRUE)
  if (length(.fr)) {
    TP3_SIGMA_PPM[.fr] <- TP3_SIGMA_PPM[.fr] * .hiket_tighten
    message(sprintf("[PRIOR TIGHTENING] TP3 fraction SDs x %.3f (%d fractions)",
                    .hiket_tighten, length(.fr)))
  } else {
    message("[PRIOR TIGHTENING] TP3 has no transfer fractions -- no effect")
  }
}

# FIXED fast rate (injected in the run script's assemble_model_params, like Yasso's
# fixed a-vector). Value = k1/xi_Ultuna at the beta centre.
TP3_ALPHA_A_FIXED <- 0.8 / 0.9397   # = 0.851

# Physical litter-flux envelope (tC/ha/yr), homogeneous across all six models.
# Bounds the EFFECTIVE flux sigma_input*J (and the 1917 flux) to the boreal NPP
# range (Gower et al. 2001); floor relaxed below the Gower min as a small,
# strictly-positive, non-binding guardrail. Enforced via the `flux_pair`
# transform in calibration_engine.R. See sigma_input_physical_bounds_note.
TP3_INPUT_FLUX_WINDOW <- c(0.05, 8.7)
