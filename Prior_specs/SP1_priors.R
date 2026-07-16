# SP1 prior specification — HIKET pipeline
# Hand-set defaults: no published global calibration exists.
# Climate parameters (beta1/beta2/gamma) copied from Yasso07 MAP as
# biologically plausible starting values for the Finnish boreal context.
# Climate widths on the Yasso07 empirical scale (Tuomi 2009 T3, 1σ) — the old
# hand-set beta2=0.05 detonated xi via the T² lever (see PRIOR_HOMOGENIZATION_PLAN.md).
# All other widths weakly informative — Finnish data dominates.

# C1 (2026-07-16): the single rate is ICBM-anchored via a VERY-INFORMATIVE prior
# (it stays free but tightly pinned). SP1 collapses the two ICBM timescales into a
# bulk residence time (1/k1 + h/k2 = 22.7 yr), so the rate is slow-dominated and
# carries the k2 (arable->forest) uncertainty -> pinned, not fixed. Centre placed on
# our xi via the one-time /xi_Ultuna offset:
#   alpha = (1 / (1/k1 + h/k2)) / xi_Ultuna = 0.0440 / 0.9397 = 0.0468
# (ICBM Andren & Katterer 1997 k1=0.8, k2=0.00605, h=0.13; xi_Ultuna at the beta
# centre, Ultuna T=5.4, T_amp=10 half-range, P=520). See REVISION_PLAN.md §C1.
SP1_FREE_DEFAULTS <- c(
  alpha       = 0.0468,   # ICBM bulk-MRT anchor (was 0.09, ~2x too fast)
  beta1       = 0.095,
  beta2       = -0.00014,
  gamma       = -1.21,
  sigma_init  = 1.00,
  sigma_input = 1.30   # C4b: re-centred >1 for missing (understorey-dominated) litter (D2)
)

# sigma_ppm: prior SDs in unconstrained (transformed) space.
# sigma_init/sigma_input widths are in the bounded flux_pair coordinate (a scaled
# logit onto the flux window), not log space — 0.50 remains weakly informative
# (~ +/-33% on the multiplier near J_bar).
SP1_SIGMA_PPM <- c(
  alpha       = 0.15,     # C1: VERY-INFORMATIVE (was 0.50) — pins the ICBM bulk
                          # rate but lets SOC data nudge the k2-dominated centre;
                          # posterior width = arable->forest transferability signal
  beta1       = 0.26,     # Yasso07 scale (was 0.20)
  beta2       = 0.00065,  # Yasso07 scale (was 0.05 — explosive)
  gamma       = 0.20,     # Yasso07 scale (was 0.30)
  sigma_init  = 0.50,
  sigma_input = 0.50
)

# Physical litter-flux envelope (tC/ha/yr), homogeneous across all six models.
# Bounds the EFFECTIVE flux sigma_input*J (and the 1917 flux) to the boreal NPP
# range (Gower et al. 2001); floor relaxed below the Gower min as a small,
# strictly-positive, non-binding guardrail. Enforced via the `flux_pair`
# transform in calibration_engine.R. See sigma_input_physical_bounds_note.
SP1_INPUT_FLUX_WINDOW <- c(0.05, 8.7)
