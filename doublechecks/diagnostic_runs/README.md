# Diagnostic runs — 2026-08-07

Short local calibrations (3 chains x 6000) run while diagnosing the error model and the
initialisation. **Not production runs.** They live here rather than in
`Calibration_real_data_transient/runs/` because `run_*_predictive.R` and the figure builders
select a posterior by sorting filenames and taking the newest, with no notion of what kind of
run produced it — so leaving them in `runs/` silently hijacks every downstream stage.

| run | configuration | key result |
|---|---|---|
| `TP2_*_20260807_142117` | log-normal likelihood | bias **+12.7 → −0.9** tC/ha; median log-residual −0.205 → −0.023. Level fixed. Trend *worsened* 0.168 → 0.116. |
| `TP2_*_20260807_145728` | log-normal + slow-rate prior widened to SD 0.50 | `alpha_H` drifted **down** 0.0053 → 0.0031 (τ 210 → 365 yr); trend moved only 0.116 → 0.130 against an observed 0.309. Freeing the rate does **not** fix the trend. |
| `Yasso20_*_20260807_145739` | log-normal likelihood | bias +11.2 → −1.9; median log-residual −0.170 → −0.001. Level fix **generalises** to a fixed-rate model. Trend −0.114 → −0.064: still a source, so the **inverted sink is not a likelihood artefact**. σ_init 1.236 → 1.042, still above the (then) 0.818 inversion threshold. |

## Important: these predate the P1 anchor change

All three were fitted with the **old** pre-run anchor (`J_full_mean` for 1917). Commit `44c2093`
changed the wrappers to use `J_t0_mean` at both ends, so **these posteriors cannot be reproduced
by current code**, and evaluating them with current wrappers gives predictions ~25% low on the
1917 flux. To re-analyse them, restore the pre-P1 wrappers first:

    git show 44c2093^:Model_functions_real_data_transient/Decomposition_functions/SimpleModels/tp2_wrapper_transient.R

(and equivalently for the others). This bit us once during the session — a comparison run with
mismatched wrappers gave a wrong answer that looked plausible.

Full analysis: `manuscript/M&M_parameterization_working_document.pdf`.
