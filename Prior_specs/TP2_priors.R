# TP2 prior specification — HIKET pipeline
# Hand-set defaults: no published global calibration exists.
# ICBM-style starting values (Andren & Katterer 1997); climate copied from
# Yasso07 MAP. Climate widths on the Yasso07 empirical scale (Tuomi 2009 T3, 1σ);
# p_H on the Tier-2 logit SD 0.4 (see PRIOR_HOMOGENIZATION_PLAN.md).
# All other widths weakly informative — Finnish data dominates.

# C1 (2026-07-16): ICBM kinetic anchor, homogeneous with Yasso (which FIXES its
# a-vector and frees its lateral fractions). alpha_A (fast) is FIXED — it is NOT in
# this free vector; it is injected as a constant in assemble_model_params (run
# script), = k1/xi_Ultuna = 0.8/0.9397 = 0.851. alpha_H (slow) stays free with a
# VERY-INFORMATIVE prior centred at k2/xi_Ultuna. p_H (humification) stays FREE
# under the Tier-2 logit prior, re-centred at ICBM h=0.13. ICBM: Andren & Katterer
# 1997 (k1=0.8, k2=0.00605, h=0.13). See NEXT_SESSION.md §5 (C1).
TP2_FREE_DEFAULTS <- c(
  alpha_H     = 0.00644,  # k2/xi_Ultuna (was 0.0015); very-informative prior
  p_H         = 0.13,     # ICBM humification h (was 0.028); FREE, logit SD 0.4
  beta1       = 0.095,
  beta2       = -0.00014,
  gamma       = -1.21,
  # P3 (2026-08-07): 0.90 = J_1917/J_1985 from the NFI growing-stock record with the
  # fitted litter-growing-stock elasticity (eps 0.45-0.66 over 1986-2023 => litter
  # scales ~ sqrt(growing stock); GS_1917/GS_1985 = 0.789 => R = 0.789^eps ~ 0.90).
  # Meaningful only together with P1 (common J_t0 anchor in the wrappers): before
  # P1, sigma_init = R * 0.818 and a centre of 1.00 silently asserted R = 1.22.
  sigma_init  = 0.90,
  sigma_input = 1.30   # C4b: re-centred >1 for missing (understorey-dominated) litter (D2)
)

# sigma_ppm: prior SDs in unconstrained (transformed) space.
# alpha_A omitted — FIXED constant (injected in assemble_model_params).
TP2_SIGMA_PPM <- c(
  alpha_H     = 0.15,     # C1: VERY-INFORMATIVE (was 0.50); k2 pinned but nudgeable
  p_H         = 0.40,     # Tier-2 logit (unchanged)
  beta1       = 0.1326555,     # Yasso07 scale (was 0.20)
  beta2       = 0.0003316,  # Yasso07 scale (was 0.05 — explosive)
  gamma       = 0.10204270,     # Yasso07 scale (was 0.30)
  sigma_init  = 0.50,
  sigma_input = 0.50
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
  .fr <- grep("^p_", names(TP2_SIGMA_PPM), value = TRUE)
  if (length(.fr)) {
    TP2_SIGMA_PPM[.fr] <- TP2_SIGMA_PPM[.fr] * .hiket_tighten
    message(sprintf("[PRIOR TIGHTENING] TP2 fraction SDs x %.3f (%d fractions)",
                    .hiket_tighten, length(.fr)))
  } else {
    message("[PRIOR TIGHTENING] TP2 has no transfer fractions -- no effect")
  }
}

# FIXED fast rate (injected in the run script's assemble_model_params, like Yasso's
# fixed a-vector). Value = k1/xi_Ultuna at the beta centre.
TP2_ALPHA_A_FIXED <- 0.8 / 0.9397   # = 0.851

# Physical litter-flux envelope (tC/ha/yr), homogeneous across all six models.
# Bounds the EFFECTIVE flux sigma_input*J (and the 1917 flux) to the boreal NPP
# range (Gower et al. 2001); floor relaxed below the Gower min as a small,
# strictly-positive, non-binding guardrail. Enforced via the `flux_pair`
# transform in calibration_engine.R. See sigma_input_physical_bounds_note.
TP2_INPUT_FLUX_WINDOW <- c(0.05, 8.7)
