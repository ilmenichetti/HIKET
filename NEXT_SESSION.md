# NEXT SESSION — start here

**Rewritten 2026-08-10 (evening).** Supersedes the 2026-08-07 version.

> **Companion documents:**
> 1. `manuscript/HIKET_next_session.pdf` (4 pp) — the level/trend discrepancy and its six
>    hypotheses. Still valid EXCEPT H2, now closed (see §4).
> 2. `manuscript/M&M_parameterization_working_document.pdf` — the parameterisation assumptions and
>    the four defects. **Defect #4 (likelihood conflates observation with model error) is the
>    subject of the run now in flight.**

---

## 0. RUN IN FLIGHT — Roihu jobs 563524–563529, launched 2026-08-10 ~15:30

| job | model | | job | model |
|---|---|---|---|---|
| 563524 | SP1 | | 563527 | Yasso07 |
| 563525 | TP2 | | 563528 | Yasso15 |
| 563526 | TP3 | | 563529 | Yasso20 |

Commit `3b0d533`. Verified at launch in all six: `total sigma OVERRIDDEN: 0.720`,
`LOG-NORMAL likelihood (default)`, `Forward-run sanity: PASS`, `5 chains x 50000`,
`Cores per chain: 40`. ~13–19 h ⇒ results morning of **2026-08-11**.

**ONE FACTOR CHANGED ON PURPOSE.** The fraction prior stays at 0.4 so the error model is tested
alone. Everything in this run is a **defect fix**, not a design choice:

1. **`HIKET_SIGMA_TOTAL=0.72`.** `sigma_obs_fixed` (0.442) is the *measurement* CV but was used as
   the *total* error; measured log-residual spread is 0.708–0.735 in all six models, implying
   model error 0.55–0.59 — larger than the observation error. Because the level penalty goes as
   1/σ², this amplified stock-level pressure ~2.6× relative to the priors, a candidate driver of
   the short bulk MRT. Set via `SINGULARITYENV_HIKET_SIGMA_TOTAL` in the SLURM scripts — **a bare
   export never reaches R inside the r-env container.**
2. **Double precision** in `yasso15.f90` + wrapper (Yasso15/20 were at ~7 significant digits while
   SP1/TP2/Yasso07 were at ~16). Verified to change no result.
3. **Prior widths corrected to source.** Tuomi 2009 T3 and 2011 T4 both state *95% confidence
   limits*; the locked convention read "±" as 1σ, ~1.96× too wide. Corrected for Yasso07 and for
   SP1/TP2/TP3 (which inherit "Yasso07 scale"). Yasso15/20 untouched — genuine posterior SDs.

---

## 1. POST-RUN CHECKLIST (in order)

```bash
rsync -av roihu:/scratch/project_2019134/HIKET/Calibration_real_data_transient/runs/ \
          ./Calibration_real_data_transient/runs/
rsync -av roihu:/scratch/project_2019134/HIKET/Calibration_real_data_transient/diagnostics/ \
          ./Calibration_real_data_transient/diagnostics/
rsync -av roihu:/scratch/project_2019134/HIKET/Data/model_inputs/ ./Data/model_inputs/
Rscript --no-save Calibration_real_data_transient/run_hiket_pipeline.R --skip-calibration
```

Use the `roihu:` alias, not the raw hostname. Sync `Data/model_inputs/` too — the predictive stage
hard-loads the bundle keyed to each RUN_ID.

Then, in order:

1. **σ self-consistency.** Recompute the log-residual spread. **≈0.72 ⇒ the plug-in choice is
   validated retrospectively**; otherwise iterate once at the new value. This is what makes 0.72
   defensible rather than circular — at convergence it *is* the free-σ answer.
2. **TP3 (and TP2) R-hat.** The specific test of whether the weaker likelihood suppresses the
   degenerate mode. **Fraction tightening did NOT** — see §3.
3. **`Rscript doublechecks/intrinsic_mrt.R`** → rebuild F12 (`manuscript/figures/build_F12_mrt_yasso.R`).
   Baseline to beat: ours **11.3 / 17.7 / 14.7** vs published **33.4 / 30.5 / 25.0** (Yasso07/15/20).
4. **F3 shape, and BIAS separately from RMSE.** A wider σ penalises the systematic +6–11% level
   offset less, so bias may persist or grow even as the trajectory shape improves.

---

## 2. ⚠ ALL MRT NUMBERS BEFORE 2026-08-10 ARE SUPERSEDED

The engine binding named `steady_state` is **not** a steady state — for Yasso it is
`*_transient_init` (1917 equilibrium + 68-year ramp to 1985). It is contaminated by `sigma_init`
(25.05 at 0.90 vs 33.84 at 0.35) and diverges for near-conservative draws via `model_step`'s Euler
fallback. That divergence is the failure Lorenzo originally reported to FMI, and it is **not** a
property of the published parameters.

**Use `doublechecks/intrinsic_mrt.R` only**: unit litter input at a fixed reference (dataset-mean
climate and AWEN × size composition), pure steady-state routine. Independent of `sigma_input`,
`sigma_init` and the SOC data; verified invariant to `sigma_input` over a 12× range.

⚠ **And MRT is conditional on the `sigma_input` prior.** MRT as a *function of parameters* is
invariant to `sigma_input`, but as an *inferred quantity* it is not: the likelihood mainly pins the
product `MRT × σ_input × J_raw` (= observed stock), so the σ_input prior decides the split.
Yasso15 published kinetics need σ_input 1.305; ours sit at 2.432 — **ratios inverse to within 8%**.
So *"MRT too short"* and *"σ_input too high"* are one finding stated from two ends, and F12 must be
read against the σ_input posterior with the caption saying so. Not purely prior-driven, though:
published kinetics with both auxiliaries optimised still lose 93 nats, i.e. there is information
beyond the level (plot-distribution shape, temporal profile) that a global multiplier cannot absorb.

---

## 3. TP3 DOES NOT CONVERGE — AND THAT IS THE RESULT

Its two modes fit within **4.2 nats** (the collapsed mode slightly *better*), posterior mass split
60/40. So the data cannot distinguish a three-pool cascade from one effective pool with a flat
climate response: **the third pool is not identifiable from two SOC observations per plot.**
R-hat 18 is the correct output for a genuinely bimodal posterior, not a sampler failure.

Structurally: TP2's bracket `[1/a_A + p_H/a_H]` has 2 free parameters (1-D ridge); TP3's
`[1/a_A + p_S/a_S + p_S·p_H/a_H]` has 4 (3-D manifold). And `p_S → 0` empties S *and* H, freeing
`a_S`, `a_H` and `p_H` at once — a 3-D flat region versus TP2's 1-D.

The fraction tightening could not fix it because the degeneracy spans the **rates** too, and those
were deliberately left at SD 0.15 as the ICBM transferability diagnostic. **You can have that
diagnostic or TP3's convergence, not obviously both.**

Report the non-identifiability rather than fixing it: it gives the complexity thread a mechanism
(SP1 and TP2 identifiable, TP3 not, boundary at the third pool) instead of only a skill comparison.

---

## 4. LITTER — H2 CLOSED

**The 2006 peak is real and the series stands as published** (confirmed with B. Tupek). A modest
post-2006 decline in modelled SOC is therefore *expected and defensible*, not an artefact to
remove — possibly weaker than currently simulated. Do not re-open.

Context: the product rises +51.7% over 1986–2006 while NFI growing stock rises +22.4%, then falls
−12.8% over 2006–2021 while growing stock rises +16.0%. Net 1986–2021: litter ×1.32, stock ×1.42.

---

## 5. NEXT LEVERS, IF THE ERROR MODEL IS NOT ENOUGH

In order:

1. **`HIKET_PRIOR_TIGHTEN=0.5`** — Tier-2 fraction SDs only (0.4 → 0.2). Already committed and
   inert by default. The one width in the whole scheme *we* chose rather than sourced, so the one
   we are entitled to narrow. Expect it to matter most for Yasso15/20 (fractions carry 58% and 94%
   of their MRT gap) and not at all for Yasso07 (whose gap is 111% climate).
2. **Free `sigma_model`** (`σ_total² = σ_obs² + σ_model²`). ~25 edits across 13 files; the two
   model families build `sigma_ppm` differently. Buys a *reportable* quantity — σ_model ≈ 0.55,
   i.e. structural error exceeds measurement error — and fixes the predictive-coverage caveat.
   ⚠ A free per-model σ reintroduces variance inflation across models, so pair it with a common
   fixed σ for cross-model comparison.
3. **Emergency plan: uniform *relative* prior width ~0.05** on all structural parameters, centres
   at published. ⚠ A uniform **0.2 in the transformed space is not feasible** — `beta2` would be
   370–2500× looser and would reproduce the historical `beta2` detonation.

---

## 6. STILL OPEN

- **Level offset** +6–11%, direction unchanged since the log-normal switch.
- **Predictive coverage** still a parameter-CI, not a posterior-predictive interval.
- **Benchmark trio unusable** until TP2/TP3 are settled (§3), so the complexity thread has no
  evidence yet.
- **NextGenC bundle, `HIKET_calibration.Rmd`, and the manuscript figure set** are all keyed to
  `20260807_1655*` and become stale the moment 563524–563529 land.
