# NEXT SESSION — start here

## ⭐ 0. STATE AT 2026-08-12 — a run is IN FLIGHT, and next session is DATA WORK

**Launched 2026-08-12**, commit `3e0c807`, error model **scale 0.80 + Student-t ν = 6**.
Jobs **597028** SP1, **597029** TP2, **597030** TP3, **597031** Yasso07, **597032** Yasso15,
**597033** Yasso20. Verified at launch, 6× each: `total sigma OVERRIDDEN: 0.800`,
`LOG-NORMAL likelihood`, `STUDENT-t tails: df = 6`, `Cores per chain: 40`,
`5 chains x 50000`. Results ~2026-08-13. ⚠ Restrict log globs to `*_597*` or you read the
previous run's logs.
**Next session's work is the SOC data, especially 1985** — not this run. See §0c.

### 0a. Why this run

The posteriors are **overconfident** — see CLAUDE.md §"THE POSTERIORS ARE OVERCONFIDENT" and
`doublechecks/ridge_test.R`, which reproduces every number. In one line: Yasso15's `sigma_input`
moved 2 posterior SD from a single defensible choice, so between-run sensitivity exceeds within-run
uncertainty.

This run fixes the two *simplest* defects only, both measured, neither assumed:

- **σ 0.72 → 0.80.** 0.72 was the in-sample spread; held-out plots give 0.784–0.799, and the
  measured spatial structure adds a little more (√(0.79² + 0.13²) ≈ 0.80). **This is the largest
  value the data support** — beyond it you contradict the residuals and forfeit the
  self-consistency argument.
- **Student-t, ν = 6.** Residual kurtosis is 6.9–7.1 in all six against a Gaussian's 3.0.

⚠ **Expect this to CONFIRM the overconfidence, not fix it** — worth ~1.3–1.4× broadening where
Yasso15 needs 2.5×. The reason to run it is that heavy tails may move parameter **locations**: under
a Gaussian a handful of badly-missed plots steer the fit. Campaign-specific σ was deliberately
**excluded** so it cannot confound the trend.

### 0b. Post-run checklist

1. **Self-consistency**, twice over: residual sd ≈ 0.80 *and* kurtosis ≈ 6.
2. **`Rscript doublechecks/ridge_test.R`** — did corr(log MRT, log σ_input) move off −0.3? Did
   `sd(log product)/sd(log MRT)` drop below 1?
3. **Did locations move?** The real question. Compare against `20260810_1529*`.
4. `doublechecks/intrinsic_mrt.R` (update the RIDs) → rebuild F12.

⚠ Log-likelihoods are **not** comparable to any previous run — both family and scale changed.

### 0c. Then: the SOC data, 1985 first

Evidence that 1985 is the weak link, none of it chosen to suit the answer:

- residual spread **0.97** vs 0.54 for 2024, near-identical across all six models (0.967–0.978), so
  it is a property of the campaign, not of any model;
- differencing does **not** clean it — Δ(2006→2024) drops to 0.37–0.45 while Δ(1985→2024) stays at
  0.77–0.93, the signature of *independent* error rather than a shared offset;
- all six models over-predict 1985, i.e. independently say the true value is **higher** than
  recorded, by ~3–6 tC/ha beyond the general bias.

**Calibrated expectations.** A 10% upward revision of 1985 raises the level ~3% ⇒ Yasso15 MRT
21.9 → ~22.6 against the 30.3 needed. **1985 barely touches the MRT tension** — it constrains
`sigma_init`, not MRT. It would however reshape the **trend**: observed 1985→2024 is +0.312, and a
5 tC/ha higher 1985 takes it to ~+0.08, *below* several models, flipping the sign of the
discrepancy rather than closing it.

⚠ **And the 2006→2024 window contains no 1985 data at all**, yet every model produces a source
(−0.02 to −0.43) against an observed +0.209. **No 1985 revision can fix that.** Don't let a 1985
finding be mistaken for a resolution of the trend problem.

---


**Rewritten 2026-08-10 (evening).** Supersedes the 2026-08-07 version.

> **Companion documents:**
> 1. `manuscript/HIKET_next_session.pdf` (4 pp) — the level/trend discrepancy and its six
>    hypotheses. Still valid EXCEPT H2, now closed (see §4).
> 2. `manuscript/M&M_parameterization_working_document.pdf` — the parameterisation assumptions and
>    the four defects. **Defect #4 (likelihood conflates observation with model error) is the
>    subject of the run now in flight.**

---

## 0a. RESULTS — Roihu jobs 563524–563529 LANDED 2026-08-11

RUN_IDs `20260810_152914` (SP1/TP2/TP3), `_152915` (Yasso07), `_152916` (Yasso15), `_152917`
(Yasso20). Stages 2–4 and the full figure set regenerated locally 2026-08-11.

**1. σ self-consistency — VALIDATED.** sd(log residual) = 0.712 / 0.735 / 0.737 / 0.737 / 0.744 /
0.732 (SP1…Yasso20) against the assumed 0.720. It is a fixed point: widening the likelihood 1.63×
moved the residual spread ~2%. **This result is bankable and needs no revisiting.**

**2. ALL SIX CONVERGED — including TP3.** Multivariate psrf 1.001–1.024, no parameter > 1.05.
TP2 was 15.6 and TP3 18.7 with ESS 4–5. **⚠ This falsifies §3 below** — see the rewritten §3.

**3. MRT rose 18–34%** (`doublechecks/intrinsic_mrt.R`, run IDs updated in place):

| | before | now | published | ratio |
|---|---|---|---|---|
| Yasso07 | 11.3 | **15.17** | 33.39 | 0.34 → 0.45 |
| Yasso15 | 17.7 | **21.91** | 30.27 | 0.58 → 0.72 |
| Yasso20 | 14.7 | **17.37** | 25.02 | 0.59 → 0.69 |

Attribution is clean for Yasso15/20: their climate priors were NOT touched (genuine posterior SDs),
so their +24%/+18% is the σ change alone. Yasso07 got both corrections and moved most.
**So the error-model scale WAS part of the MRT story — but not all of it.**

**4. σ_input fell 7–16% everywhere** — SP1 2.06→1.84, TP2 2.55→2.38, TP3 2.59→2.39,
Yasso07 2.92→2.49, Yasso15 2.43→2.04, Yasso20 2.67→2.42. The level ridge holds exactly:
Yasso15's MRT×σ_input goes 43.0 → 44.7, constant to 4%.

**5. Still unfixed.** Bias GREW (+3.7 to +6.8 tC/ha, Yasso20 worst) — expected, a wider σ penalises
the level offset less. R² still ~0 (0.004–0.019). And **2006–2024 is a source in all six**
(−0.018 to −0.433) against observed +0.209; full-window positive in five of six (Yasso15 +0.379 vs
obs +0.312) but by cancellation.

**Next lever (NOT launched — pending discussion 2026-08-11).** Tighten `sigma_input` log-SD
0.50 → 0.20, centre 1.30 (Lehtonen & Heikkinen 2015). The ridge makes the outcome predictable:
σ_input forced to 1.30 implies Yasso15 MRT ≈ 44.7/1.30 ≈ **34 yr** vs published 30.3. Two caveats
to state out loud: (a) it partly IMPOSES the MRT answer — defensible because the width comes from
independent NFI sampling error rather than the SOC being fitted, but it must be declared in the
methods; (b) the note to apply it to the CONTEMPORARY flux only is NOT implemented, so the global
version risks σ_init absorbing the slack (watch it — currently 0.53–1.11).

**Known gap:** the multimodel projection plot is skipped — `nfi_region` is absent from
`site_raw.csv`; run `assign_nfi_regions.R` first.

---

## 0b. (historical) RUN IN FLIGHT — Roihu jobs 563524–563529, launched 2026-08-10 ~15:30

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

## 3. ⚠ RETRACTED 2026-08-11 — TP3 *DOES* CONVERGE

**Everything below this box is superseded. Do not report it.** In run `20260810_152914` TP3 has
all 9 parameters R-hat < 1.05 and multivariate psrf **1.0053**. The collapsed mode is gone, not
outvoted: `p_S` is cleanly unimodal at 0.455 (90% 0.354–0.568), where the old quantiles ran
0.004 / 0.006 / 0.583 / 0.635 — the bimodal fingerprint — and `gamma` no longer reaches the
climate-off value (old 97.5% −0.14, new −1.19).

**The bimodality was an artefact of the over-wide climate priors** (the Tuomi 95%-as-1σ error:
`beta1` 0.26→0.133, `gamma` 0.20→0.102), not structural non-identifiability of the third pool. The
degenerate mode required climate switched off; corrected to source, that region is unaffordable.

**Attribution is not fully settled.** This run moved climate widths *and* σ together. The
`HIKET_PRIOR_TIGHTEN=0.5` test cited below looks like it separates them, but its chains file is
1.8 MB against production's 20.9 MB (~11× fewer iterations), so its R-hat 1.3–18.4 may be
incomplete mixing rather than genuine bimodality. **Do not lean on it.** A clean separation needs a
full-length run at σ=0.442 with corrected climate priors — cheap, and worth doing if the
prior-width mechanism is to be claimed in print.

**Consequence for the manuscript:** the complexity thread loses its mechanism and reverts to a
skill comparison. The mode-anatomy numbers below stay valid as a description of the OLD posterior
and are still usable as a methods cautionary point about prior width.

<details><summary>Superseded text</summary>

### TP3 DOES NOT CONVERGE — AND THAT IS THE RESULT

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

</details>

---

## 4. LITTER — H2 CLOSED

**The 2006 peak is real and the series stands as published** (confirmed with B. Tupek). A modest
post-2006 decline in modelled SOC is therefore *expected and defensible*, not an artefact to
remove — possibly weaker than currently simulated. Do not re-open.

Context: the product rises +51.7% over 1986–2006 while NFI growing stock rises +22.4%, then falls
−12.8% over 2006–2021 while growing stock rises +16.0%. Net 1986–2021: litter ×1.32, stock ×1.42.

---

## 5. THE IDENTIFICATION PROBLEM — how to choose a lever

**The role assignment is right and is empirically supported:**

| knob | what it sets | evidence |
|---|---|---|
| `sigma_init` | the **trend** (how far below equilibrium 1985 sits) | orders the paired stock change monotonically across four structurally different models, crossing zero at σ_init ≈ 0.9 |
| `sigma_input × MRT` | the **level** | the two ratios are inverse to within 8% (Yasso15: published 30.45×1.305 vs ours 17.68×2.432) |
| `MRT` alone | the **responsiveness** (trajectory curvature) | fast models track the 2006 litter peak; slow ones smooth it |

Three features — level, slope, curvature — against three knobs, with three campaigns to see them.
**So it is identified in principle.** It fails in practice because the likelihood sums 1269
plot-year observations at equal weight with residual spread 0.72: the level is worth ~93 nats, the
whole 1985–2024 trend ~3.4, and curvature less still. **The knobs are right; the likelihood cannot
see the features that separate them.**

### ⚠ A standing claim that does NOT survive this

"Models need 2–3× the measured litter to fit" is **not a finding at current values.** The test is
whether *any* physically reasonable point on the ridge fits:

- **σ_input 13–20× (pre-`flux_pair`, TP2/TP3):** no — closing it would need MRT ≈ 5 yr, absurd for
  a whole profile. Every point on that ridge was impossible ⇒ **genuine finding, still stands**,
  and `F9_effective_flux_vs_ceiling` keeps its meaning because it is pinned to that comparison.
- **σ_input 2.4× (now):** yes — MRT ≈ 30 with σ_input ≈ 1.3 fits the level (scenario test). So the
  model needs nothing; our prior placed us at one end of an unconstrained ridge. Reporting it as a
  finding would attribute to Yasso what our own σ_input prior did.

What may survive is much narrower: published kinetics with *both* auxiliaries optimised still lose
93 nats (~35 at the corrected σ). That is a claim about plot-level distribution shape and timing,
**not** about litter magnitude. Check whether it survives run 563524 at all.

### Closing the ridge: two ends, not equally defensible

- **(a) Narrow `sigma_input`** ⇒ the level then identifies MRT.
- **(b) Narrow the kinetics** ⇒ the level then identifies σ_input.

For **σ_input** we have a physical argument with bounds (the Tupek product excludes understorey;
understorey is a reasonably constrained fraction of boreal litter) — independent of the SOC being
fitted. For the **kinetics** the external information is weaker and its transferability is exactly
what this study questions: the ICBM anchor is arable (Ultuna), and Yasso's values come from
litterbag experiments elsewhere — Toni's own reason that local recalibration is right.

**The current setup has this backwards.** The kinetics are already tightly constrained (Yasso rates
fixed, fractions prior-pinned, ICBM anchors) while `sigma_input` sits at log SD **0.50** (±65%) —
we constrain what we know less about and leave loose what we know better. That argues for (a)
independently of which answer it produces.

### (c) Add information instead of constraining

MRT governs responsiveness, σ_input does not, so any observable measuring how sharply the system
tracks forcing separates them. Candidates: **radiocarbon** (constrains turnover directly,
independent of input magnitude), or the **vertical distribution** — the homogenised target already
carries organic / mineral 0–40 / deep tail, and pool-resolved information would identify routing
without touching the input. Layers are not pools so the mapping is not clean, but this is the
principled route and belongs in the discussion even if not attempted.

### Decision order

1. See what the error model alone did (run 563524).
2. If the ridge is still open, close it at the **σ_input** end and **measure what it costs**.
3. If constraining σ_input to its physical range degrades the fit badly, that is a genuine conflict
   between the litter estimate, the models and the SOC target — a **finding**, not a knob to turn.

---

## 6. NEXT ACTION — TIGHTEN THE `sigma_input` PRIOR (Lehtonen & Heikkinen 2015)

**Do this regardless of how the error-model run comes out.** `sigma_input`'s prior is log SD
**0.50** (±65%) on a quantity we can bound from the literature — that is *unused information*, the
same class of defect as the Tuomi widths, not a lever pulled to get an answer.

**Source** (`literature/lehtonen-heikkinen-2015-uncertainty-of-upland-soil-carbon-sink-estimate-for-finland.pdf`,
Aleksi's pointer; Figs 2–3):

- Total litter input (tree + understorey), 1990–2013: **~3.15 tC/ha/yr south, ~2.15 north**,
  95% CI ≈ ±15–20% ⇒ **σ ≈ 8–10%**. *(Fig. 2 read by eye — ±0.2 on the levels.)*
- **The spread is predominantly EMPIRICAL.** Their own decomposition: *"uncertainty about the
  volumes of living trees due to sampling error in NFI was clearly more influential than the other
  components … the effects of uncertainty in logging volumes, biomass models, litter rates, and
  understorey litter were relatively small."* The stipulated CVs (5% logging, 10% understorey) are
  the *small* ones. ⚠ That decomposition is reported for stock CHANGE, not for the litter band
  itself — litter derives from tree volumes so the same term should dominate, but it is inference.
- ⚠ It is a **floor**: leaf-to-fine-root ratios and AWEN proportions were excluded.
- ⚠ The **year-to-year variation** that dominates their sensitivity analysis is assumed (5/10/20%
  tested) and is **Fig. 3, not the Fig. 2 band** — do not import it as level uncertainty.

**USE THE SPREAD, KEEP OUR CENTRE.** Their *total* (~2.7 national, area-weighted ~60/40) is
essentially identical to our Tupek **tree-only** 1990–2013 mean of **2.697** — which is not
agreement but a hidden discrepancy (Tupek's implied total exceeds theirs by the understorey). It
cannot be resolved from the paper: the text points to Appendix A Table A1 for the understorey
litter, but that table is *stem volumes*. So do **not** adopt their level. Keeping our centre works
anyway:

    sigma_input centre 1.30 x J_raw 2.697 = 3.51 tC/ha/yr
    published kinetics matched the level at sigma_input ~1.305 -> 3.52

**The centre was never the problem — the width was.** At 0.50 the posterior ran to 2.43
(flux 6.56, ~2.4x Lehtonen's total).

### Construction

| | |
|---|---|
| quantity | the **contemporary** effective flux `sigma_input x J_bar` |
| centre | **unchanged at 1.30** |
| width | **0.50 → 0.15–0.20** (0.25 = exactly halving, the conservative choice) |
| form | **two-sided informative prior**, not a ceiling — a conflict should show as a reportable displacement, not be clipped at a boundary that also mixes badly |
| untouched | `sigma_init` and the **1917** flux |
| backstop | NPP ceiling, now inactive |

⚠ **Implementation trap:** `flux_pair` currently bounds BOTH `sigma_input x J_bar` AND
`sigma_init x sigma_input x J_full` (the 1917 flux). Applying this naively would tighten the
historical flux too, push `sigma_init`, and **leak straight into the trend**. Apply to the
contemporary side only, and verify against the transform code first.

### Expected effect — helps, but not decisive alone

Static balance for Yasso15 at the current position: tightening to 0.20 raises the `sigma_input`
prior penalty 0.8 → 4.9, so total prior cost ~21 → ~25, against a likelihood advantage of ~35 nats
once the error model is in (93 at the old sigma). **Near balance, not past it.** The error model
does the heavy lifting; this narrows the remainder. If BOTH together still fail to move it, that is
the informative outcome — the data demanding the fast solution against two independent external
constraints — and it should be reported as a conflict rather than met with a third lever.

---

## 7. FRACTION-TIGHTENING TEST (local, 2026-08-10) — WHAT IT DID AND DIDN'T SHOW

`HIKET_PRIOR_TIGHTEN=0.5`, 3 chains x 6000, all six, on the double-precision build with the
source-corrected climate widths. Outputs quarantined to `doublechecks/prior_tighten_test/`.

| model | R-hat | ESS | |
|---|---|---|---|
| SP1 | 1.003–1.013 | 331–632 | clean |
| TP2 | 1.007–1.099 | 15–124 | OK |
| **TP3** | **1.329–18.372** | **2–34** | **broken** |
| Yasso15 | 1.002–2.326 (20 warn) | 9–136 | poor |
| Yasso07 | 1.022–3.066 (18 warn) | 6–61 | poor |
| Yasso20 | 1.004–2.135 (16 warn) | 13–153 | poor |

**ESTABLISHED:** TP3 is *not* fixed by fraction tightening — confirmed now at two run lengths
(5x50000 on 2026-08-07 and 3x6000 here), so it is not a sampling fluke. Consistent with §3: the
degeneracy spans the RATES, which the fraction prior does not touch. SP1 and TP2 unaffected.

**NOT ESTABLISHED — do not over-read the Yasso rows.** At 6000 iterations with 20–26 free
parameters and ESS 6–153, R-hat is unreliable (it is inflated at ESS ~10). This is almost certainly
insufficient sampling, **not** evidence that tightening harmed the Yasso models. The two cannot be
separated without production-length chains.

**CONSEQUENCE:** the test did **not** deliver the fraction-only MRT effect for the Yasso family —
those posteriors are too poorly mixed to quote an MRT from. So there is no clean single-factor
fraction result to pair with the error-model run. Obtaining one needs a production-length Roihu
job, which is only worth spending if §6 (the `sigma_input` tightening) leaves a gap.

---

## 8. FURTHER LEVERS

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

## 9. STILL OPEN

- **Level offset** +6–11%, direction unchanged since the log-normal switch.
- **Predictive coverage** still a parameter-CI, not a posterior-predictive interval.
- **Benchmark trio unusable** until TP2/TP3 are settled (§3), so the complexity thread has no
  evidence yet.
- **NextGenC bundle, `HIKET_calibration.Rmd`, and the manuscript figure set** are all keyed to
  `20260807_1655*` and become stale the moment 563524–563529 land.
