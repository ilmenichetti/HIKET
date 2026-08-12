# run_config.R — shared MCMC run settings for all HIKET models.
# Edit here to change settings across all five models simultaneously.

# --- Production defaults -----------------------------------------------------
N_PLOTS_TEST <- NA    # NA = full dataset; set to e.g. 20L for quick tests
N_CHAINS     <- 5L
N_ITER       <- 50000L
N_BURNIN     <- 5000L
N_LOG        <- 200L

# --- Ablation overrides (env vars; ALL default to the values above) ----------
# Set by doublechecks/run_ablation.R so a config can be varied WITHOUT editing
# this file. Unset => production behaviour is bit-for-bit unchanged.
.envnum <- function(nm, default) {
  v <- Sys.getenv(nm, NA_character_)
  if (is.na(v) || !nzchar(v)) default else as.numeric(v)
}
N_CHAINS <- as.integer(.envnum("HIKET_N_CHAINS", N_CHAINS))
N_ITER   <- as.integer(.envnum("HIKET_N_ITER",   N_ITER))
N_BURNIN <- as.integer(.envnum("HIKET_N_BURNIN", N_BURNIN))

# C5 (2026-07-16) -- WITHDRAWN 2026-08-05. The 1985 (VMI8) campaign is NOT down-weighted.
#
# C5 inflated the observation SD of the 1985 campaign because the VMI8 mineral stocks looked
# systematically low (the 1985->2006 jump of +61%). That premise no longer holds: the SOC
# homogenization (2026-08-04) traced most of that jump to a missing coarse-fragment correction
# and to cross-campaign processing drift, leaving a physical +12%. Ablations then showed C5 to
# be inert for prediction and consequential in the wrong way:
#   - turning it off changes trusted-campaign (2006+2024) RMSE by 0.45% (TP2) / 0.25% (SP1),
#     against a ~2% noise floor;
#   - but it moves sigma_init from 0.199 (C5=3) to 0.780 (1985 withheld) -- a factor of ~4 on
#     a headline parameter. A knob that does not improve prediction while substantially moving
#     what we report is a researcher degree of freedom.
#   - with 1985 fully trusted the model over-predicts it by only +3.2 tC/ha relative to its
#     general bias, and by +6.3 when it has never seen it: no anomaly to correct.
#   - the mechanism never matched the claim. Inflating INDEPENDENT per-observation sigma
#     averages away over 441 plots (SE ~ sigma/sqrt(n)), whereas a protocol difference is a
#     SHARED offset whose uncertainty does not shrink with n at all.
#   - the faithful version (a campaign-level offset delta) IS identifiable, contrary to the
#     revision plan -- but it is informed only by the 1985 residual itself, requires
#     sigma_init ~0.98 (80% of draws past the pre-run inversion threshold), and implies no
#     accumulation ever occurred. It behaves as a misfit sink, not a measurement diagnostic.
# The residual concern -- 2006/2024 were validated against official LUKE stocks and 1985 could
# not be, and the protocols differ (VMI8 0-5/5-20 cm vs 0-10/10-20 cm) -- is real, and belongs
# in the limitations as a stated sensitivity rather than in a tuning constant.
# Full record: manuscript/HIKET_data_and_ablation_tests.pdf
#
# Kept as a switch rather than deleted, so the sensitivity stays reproducible: the ablation
# harness sets HIKET_SIGMA_1985_INFL to sweep it. Default 1.0 = no down-weighting.
# -----------------------------------------------------------------------------
# REINSTATED 2026-08-12 at f = 2.0, on different grounds from the withdrawn C5.
#
# C5 was withdrawn because its premise ("the VMI8 mineral stocks look low") was
# informed only by the 1985 residual itself -- i.e. it down-weighted a campaign
# because the model fitted it poorly, which is circular. What has changed is that
# there are now DOCUMENTED, MODEL-INDEPENDENT reasons to distrust the first
# campaign, established from the source workbook and the field protocol:
#   * ~75% of its mineral carbon concentrations were PREDICTED from loss-on-
#     ignition by regression, against only the 20-40 cm layer in 2006, and on a
#     different instrument (Kramarenko 2012 sec. 4.4 calls this a possible
#     campaign-level "tasoero");
#   * its subplots sat OUTSIDE the plot at 11 m, while 2006 sampled INSIDE at 9 m
#     -- literally different soil;
#   * its 1985->2006 mineral gain is DEPTH-INVERTED (subsoil +19% vs topsoil
#     +8.5%), and a real input-driven gain must be surface-weighted.
# The two defects that could be repaired have been (missing LM litter layer; the
# 1986-1995 sampling dates). This term covers only what remains.
#
# f = 2.0 is a JUDGEMENT, not an estimate. It is PRE-REGISTERED: fixed before the
# corrected data were run, and it must NOT be revisited on the basis of results --
# that is the entire defence against the circularity that sank C5. Weight scales
# as 1/f^2, so the first campaign carries a quarter of the weight of the others.
#
# ⚠ Not a "cautionary" or neutral choice: the C5 ablation moved the 1985-2024
# change from +0.638 to -0.184 as the weighting weakened. It can flip the sign of
# the sink, so NEVER report a single f. Sensitivity: short-chain Yasso15 runs at
# f = 1 / 1.5 / 3 via doublechecks/run_ablation.R (A1/A2/A3); A1 doubles as the
# attribution arm separating "the data fixes moved it" from "f = 2 moved it".
# Full record: NEXT_SESSION.md sec. 0c, manuscript/HIKET_discussion_memo.tex.
# -----------------------------------------------------------------------------
SIGMA_1985_INFL <- .envnum("HIKET_SIGMA_1985_INFL", 2.0)   # 2.0 = ON, pre-registered 2026-08-12