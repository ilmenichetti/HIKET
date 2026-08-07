# NEXT SESSION — start here

**Rewritten 2026-08-07 (evening).** Supersedes the earlier 2026-08-07 version and
`manuscript/REVISION_PLAN.md`, which has been removed (recoverable from git commit `5324383`;
its C1–C5 decisions are summarised in §5).

> **Read `manuscript/M&M_parameterization_working_document.pdf` first** (15 pp). It states every
> parameterisation assumption, the four defects, the proposals and the evidence for each.

---

## 1. What this run is

The first calibration with a **correctly specified error model**, plus the initialisation defect
fixed. Three changes:

| change | why | evidence |
|---|---|---|
| **log-normal likelihood** | the multiplicative normal sets the variance from the parameter being estimated, so a prediction widens its own error bar and the location estimate is contaminated by dispersion | TP2 refit: bias **+12.7 → −0.9** tC/ha; median log-residual **−0.205 → −0.023** |
| **common pre-run anchor** | the two ends of the pre-run used different litter aggregates (`J_full` vs `J_t0`), so σ_init was the 1917/1985 ratio × 0.75 — a centre of 1.00 silently asserted **33% more** litter in 1917 than 1985 | plain defect; wrappers are shared, so calibration and prediction follow together |
| **ratio prior centre 0.90** | NFI growing stock with the fitted litter–growing-stock elasticity (ε ≈ 0.45–0.66 ⇒ litter ~ √GS; GS₁₉₁₇/GS₁₉₈₅ = 0.789 ⇒ R ≈ 0.90) | meaningful only *together with* the common anchor |

**Implementation state.** The log-normal is a switch (`HIKET_LOGNORMAL_LIK=1`, default off) — it
must be set explicitly in the SLURM scripts. The common anchor and the 0.90 centre are **direct
edits and live by default**. The `20260805` posteriors are therefore no longer reproducible from
current code; wrapper backups are in the session scratchpad.

### Deliberately not included

- **Plain priors** (drop `flux_pair` for two lognormals) — tested safe (identical pre-flight
  blow-ups in all six models, 0% prior mass above the NPP ceiling) but it changes the prior
  *shape* as well as its centre, which would make this run unattributable. Next.
- **Free slow rate** — tested, and it does **not** fix the trend (§3). Left at C1.
- **Ratio coordinate** — superseded; the common anchor achieves the same thing more simply.

---

## 2. What to check when it lands

| | expectation |
|---|---|
| **Level** | bias near zero at all three campaigns (TP2 achieved −0.9) |
| **Yasso15/20 sink sign** | should turn positive; they were −0.063 and −0.155 against observed +0.312 |
| **Yasso trajectory shape** | rise then flatten (saturation), as in `20260710` — *not* monotonic decline |
| **F9, Yasso** | effective flux well inside the NPP envelope — **non-negotiable** |
| **F9, SP1/TP2/TP3** | high flux tolerable; report as a finding |
| **Trend** | within roughly a third of the observed **+0.312 ± 0.035**; report +0.209 (LUKE-validated campaigns only) as sensitivity |
| **σ_init** | inversion threshold is now **1.0**, not 0.818 — the common anchor removed the conversion |

**Before launching**: run `preflight_prior_pushforward.R` at full N. The common anchor lowers
every plot's 1917 flux by ~25%, so the initial state changes everywhere. This is the gate that
caught the `beta2` detonation.

---

## 3. What changed the framing

**The level bias was ours, not the models'.** ~20%, uniform across all six, caused by the error
model. Every level-based number from `20260805` carries it.

**The trend is invisible to the likelihood.** A *perfect* trend fit is worth **3.4 nats**; the
parameter move delivering it costs **8.9** even at a widened prior. The calibration correctly
ignores it — confirmed empirically: freeing the slow rate sent `alpha_H` *down* (0.0053 → 0.0031)
and moved the trend only 0.116 → 0.130 against an observed 0.309. **No reparameterisation fixes
the trend.** Fitting it needs a change to what the likelihood targets — paired per-plot
differences, or the aggregate stock change as an explicit observation with its own (0.035)
uncertainty.

**The trend deficit is not new.** In `20260710` the models captured 25–50% of the observed rate;
now 42%. What changed is that the old inflated target coincided with the models' levels at
2006/2024, so F4 *looked* right. Yasso15/20 already had negative 2006→2024 trends in July.

**The two families fail differently.** Yasso's effective response time (17–38 yr) is adequate but
it starts at or above its own equilibrium; TP2/TP3 start sensibly but are far too slow
(τ ≈ 195 yr). Only the Yasso failure is fixable by initialisation — and Yasso's fix is
**flux-neutral**, whereas the cascades' would cost input.

**Simple structures *can* fit.** TP2 reaches the observed level and trend with `alpha_H` ≈ 0.025,
`p_H` ≈ 0.30, σ_input ≈ 2.3 — flux 5.8 against a ceiling of 8.7. The solution exists; the
likelihood just doesn't reward it. The ICBM anchor wasn't making the comparison fair, it was
pinning a parameter the data has almost no opinion about.

**σ_init is conditionally identified.** The data narrows it 2–7× and shifts the median up to 44%
off the prior — it is *not* blind. But the information sits almost entirely in the 1985 campaign,
so the inference is conditional: *given that VMI8 is taken at face value…*. C5 moved it 4× while
trusted-campaign RMSE moved 0.45%. A scope condition, not a defect.

**Retracted today**: the historical-depletion hypothesis (litter raking, slash-and-burn). It came
from a fixed-equilibrium fit that ignored the rising litter. Litter rises ×1.35 over the window,
and a soil merely *tracking* that rise reproduces the observed accumulation with no
below-equilibrium start. The hypothesis may be true; the SOC record doesn't demand it.

---

## 4. Next session

1. **Analyse the run** against §2.
2. **Plain priors** — drop `flux_pair`, two lognormals plus a flux guard. Tested; deferred only
   for attribution.
3. **Predictive-stage defects**, cheap and independent:
   - the coverage column is a parameter-CI mislabelled as posterior-predictive (omits observation
     error — the "Param cov ≈ 0.05" artefact)
   - no achieved-plot-count assertion (this is what let the `alpha_A` bug run silently)
   - Yasso20 lacks the `tryCatch` guards Yasso15 has
4. **NextGenC** — still on stale `20260710` RUN_IDs in `build_soc_matrices.R`,
   `build_soc_maps.R`, `build_vulnerability_map.R`, `diagnostic_semivariogram.R`. They run
   *successfully* on stale data because those posteriors were deliberately kept.
5. **For collaborators**: B. Tupek — is the litter product derived from biomass? If so the
   elasticity behind R = 0.90 is internal to it rather than independent evidence; this is the
   load-bearing caveat. A. Lehtonen — stock-QC corroboration and understorey.

---

## 5. C1–C5, for the record (`REVISION_PLAN.md` removed)

| | what | status |
|---|---|---|
| C1 | ICBM-anchored kinetics for SP1/TP2/TP3 (fast rate fixed, slow rate pinned at SD 0.15) | in. Its transferability diagnostic **fired** — posteriors sit 1.6–2.2σ off the anchor. But relaxing it does **not** fix the trend |
| C2 | TP3 climate response on all pools | in |
| C3 | growing-stock *shape* for the pre-run ramp | in. Its *level* was never taken from the same record — that is what the common anchor and ratio prior now fix |
| C4 | σ_input reframe; C4b re-centred at 1.30 | in. Note the σ_input prior is nearly irrelevant: the likelihood decides it by ~112:1 |
| C5 | down-weight the 1985 campaign | **withdrawn**. See the 2026-08-07 addendum in `HIKET_data_and_ablation_tests.pdf` — it was adjudicated on stock RMSE, which is blind to the stock *change* it actually controls |

---

## 6. Housekeeping

- **The 2026-08-07 diagnostic posteriors are in `runs/`** (`TP2_posterior_20260807_142117`,
  `_145728`, plus a Yasso20). They sort newest and **will be auto-selected by every downstream
  script** — the `run_ids.R` resolver already refused once because of it. Quarantine before
  running anything downstream.
- **Uncommitted code from this session**: TP2/TP3 predictive `alpha_A` fix; two auto-selecting
  doublechecks; the error-model and rate-prior switches; the common anchor; the 0.90 centre;
  `preflight_naive_priors.R`, `prerun_direction.R`, `ablation_stock_change.R`,
  `production_fit_by_campaign.R`, `build_F5_stock_change.R`, `run_ids.R`, and the figure-builder
  auto-selection pass.
