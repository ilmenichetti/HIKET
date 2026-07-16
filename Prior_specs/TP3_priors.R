# TP3 prior specification — HIKET pipeline
# Hand-set defaults: no published global calibration exists.
# Extends TP2 with a Slow pool; alpha_S intermediate between alpha_A and
# alpha_H; p_S matches TP2's p_H; p_H neutral at 0.50.
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
# ICBM h=0.13. ICBM: Andren & Katterer 1997. See REVISION_PLAN.md §C1.
TP3_FREE_DEFAULTS <- c(
  alpha_S     = 0.0074,   # k2/(1-p_H)/xi_Ultuna, k2-scale (was 0.10 = cascade-starving)
  alpha_H     = 0.00644,  # k2/xi_Ultuna (was 0.0015)
  p_S         = 0.13,     # ICBM humification h (was 0.028); FREE, logit SD 0.4
  p_H         = 0.13,     # ICBM humification h (was 0.50); FREE, logit SD 0.4
  beta1       = 0.095,
  beta2       = -0.00014,
  gamma       = -1.21,
  sigma_init  = 1.00,
  sigma_input = 1.30   # C4b: re-centred >1 for missing (understorey-dominated) litter (D2)
)

# sigma_ppm: prior SDs in unconstrained (transformed) space.
# alpha_A omitted — FIXED constant (injected in assemble_model_params).
TP3_SIGMA_PPM <- c(
  alpha_S     = 0.15,     # C1: VERY-INFORMATIVE (was 0.50)
  alpha_H     = 0.15,     # C1: VERY-INFORMATIVE (was 0.50)
  p_S         = 0.40,     # Tier-2 logit (unchanged)
  p_H         = 0.40,     # Tier-2 logit (unchanged)
  beta1       = 0.26,     # Yasso07 scale (was 0.20)
  beta2       = 0.00065,  # Yasso07 scale (was 0.05 — explosive)
  gamma       = 0.20,     # Yasso07 scale (was 0.30)
  sigma_init  = 0.50,
  sigma_input = 0.50
)

# FIXED fast rate (injected in the run script's assemble_model_params, like Yasso's
# fixed a-vector). Value = k1/xi_Ultuna at the beta centre.
TP3_ALPHA_A_FIXED <- 0.8 / 0.9397   # = 0.851

# Physical litter-flux envelope (tC/ha/yr), homogeneous across all six models.
# Bounds the EFFECTIVE flux sigma_input*J (and the 1917 flux) to the boreal NPP
# range (Gower et al. 2001); floor relaxed below the Gower min as a small,
# strictly-positive, non-binding guardrail. Enforced via the `flux_pair`
# transform in calibration_engine.R. See sigma_input_physical_bounds_note.
TP3_INPUT_FLUX_WINDOW <- c(0.05, 8.7)
