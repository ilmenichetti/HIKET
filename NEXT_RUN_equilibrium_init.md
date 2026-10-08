# The equilibrium-init counterfactual — run design

**Status (2026-10-08): DESIGN SETTLED, CODE WRITTEN + LOCALLY TESTED, NOT RUN.** Rewritten from the 2026-09-04
draft after Samuli's storyline comment and the 2026-10-08 discussion. Lorenzo: design and figures
first; launch later.

## Why

The paper argues that a transient start is needed to reproduce the observed trajectory, but nothing
in it **shows** that. No experiment so far has an equilibrium-initialised arm. A reviewer will ask
for one, and Samuli already has.

**A forward-only swap is not the test** (Lorenzo, 2026-09-04): posteriors calibrated with a
transient start would fail under an equilibrium start simply because they were never refitted. The
equilibrium arm must be **calibrated** with the same freedom, so that it shows the best an
equilibrium start can do.

## The arm, in one sentence

**The production calibration with `σ_init ≡ 1`.** Nothing else changes.

Why that is exactly an equilibrium start: the pre-run flux is
`J(i) = J_1917 + (J_1985 − J_1917)·shape[i]`, with `J_1917 = J_t0_mean·σ_init·σ_input` and
`J_1985 = J_t0_mean·σ_input` (`tp2_wrapper_transient.R:362`, `yasso15_wrapper_transient.R:290`,
same in all six). With `σ_init = 1` the bracket is zero, so the ramp is flat **whatever the shape**
(Liski or growing stock), and the pre-run returns the steady state it started from:

- **litter:** mean of 1985–89 (`J_t0_mean`, `N_PREINIT_SMOOTH = 5`), with that window's AWEN and size
  composition for Yasso;
- **climate:** mean of 1985–2004 (`STEADY_STATE_YEARS = 20`; the climate series starts in 1985). The
  window reaches past t0, but the transient arm has the same look-ahead.

So the two arms differ in one thing only: **whether the soil carries a deficit into 1985.**

| | transient (production) | equilibrium arm |
|---|---|---|
| state in 1985 | below equilibrium, still accumulating | at steady state for 1985–89 litter and 1985–2004 climate |
| `σ_init` | free (prior centre 0.90, log-SD 0.25) | **fixed at 1**, removed from the free set |
| `σ_input` | free, arm B | free, **identical effective prior** (see trap 1) |
| fractions, climate, size, pinned rates | free / pinned as now | identical |
| target, plots, holdout split, correlated likelihood, σ_total 0.800 | production | byte-identical |

**Relation to the GHG inventory.** Same kind of start, different date and windows. The NID
(`literature/10_national_soil_terms/FI_NID_2025.txt:16476`) equilibrates Yasso07 to NFI6 litter
(1971–74 South, 1975–76 North) and the **1960–1990** mean climate, then runs forward with 30-year
moving-average climate and MAP parameters. Our litter starts in 1985, so we cannot reproduce that
exactly. In the paper, call the arm **"the inventory's convention applied at our start date"**,
never "the inventory's procedure". A 1960–1990-climate variant is possible (the climate file
starts in 1961) but is a second factor. **Not in this run.**

**Rejected designs:** fixing `σ_input` so that the level must come from turnover alone (not
interesting, and two factors); the Liski-prescribed transient arm (`σ_init` = 0.861 from his own
per-hectare reconstruction, raw unclamped shape). The latter is parked, not dead.

## Implementation — DONE 2026-10-08 (uncommitted until reviewed)

- **`Calibration_real_data_transient/equilibrium_init.R`** (sourced by `calibration_engine_transient.R`,
  so by all twelve scripts): reads `HIKET_EQUILIBRIUM_INIT`; `hiket_eq_tag()` (MODEL_NAME →
  `<model>_eqinit`), `hiket_eq_param_spec()` (flux_pair → flux_now), `hiket_eq_assemble()`
  (forces `sigma_init = 1`, also over the 0.90 in `free_defaults` used by the sanity run),
  `hiket_eq_names()` (drops `sigma_init` from plot lists). Inert when the variable is unset.
- **`calibration_engine.R`:** new transform type **`flux_now`** (see trap 1).
- **Six calibration scripts:** four one-line hooks (tag, param_spec, assemble, plot lists) + the
  `J_bar` injection now matches `flux_now` too. `make_likelihood(…, transient_init = TRUE)` untouched.
- **Six predictive scripts:** tag, tagged file patterns, assemble hook. `forward_scenarios.R`
  needs nothing (it already sets `sigma_init = 1` for the equilibrium stock).
- **Outputs:** `<MODEL>_eqinit_posterior_<RUN_ID>.rds`, `<MODEL>_eqinit_inputs_<RUN_ID>.rds`,
  `diagnostics/<MODEL>_eqinit/`. `run_ids.R` never matches them. ⬜ `run_ids_eqinit.R` for the figures.
- **Launch:** `Calibration_real_data_transient/submit_eqinit.sh` submits the UNCHANGED production
  `hiket_<model>.sh` with `--export=ALL,SINGULARITYENV_HIKET_EQUILIBRIUM_INIT=1` and `eqinit_*` logs,
  so the configuration is production's by construction.
- **Predictive stage:** run each `run_<M>_transient_predictive.R` with `HIKET_EQUILIBRIUM_INIT=1`
  (it auto-detects the newest `_eqinit` posterior). No residual/multimodel stage needed.
- Pre-run speed-up (skipping the 68 flat steps) deliberately NOT done: the nesting test is exact
  because both arms run the same code path.

### ⚠ Trap 1 — the `σ_input` prior is NOT the first coordinate alone

The prior is Gaussian in sampling space, and `log_jacobian()` is added to the log-likelihood on top
of it. So the effective prior is `N(x)·|J(x)|`. For `flux_pair` the Jacobian is
`g'(x1)·g'(x2)/(J̄·F_now)`. Integrating out `x2` leaves, for `σ_input`,

```
N(x1) · g'(x1) / F_now        (g = bounded logistic onto the flux window)
```

The `1/F_now` comes from `σ_init = F_1917/F_now` and **stays in the marginal**. A naive
one-parameter bounded group would give `N(x1)·g'(x1)` and tilt the equilibrium arm's prior
**towards larger `σ_input`**, by exactly the factor the arm is expected to exploit. **Fix:** in the
equilibrium arm, the `σ_input` Jacobian term must be `log_jac_bounded(F_now) − log(F_now)`.
**Test:** `doublechecks/eqinit_prior_check.R` — ✅ PASS 2026-10-08, all six: total-variation
distance 5e-16. The naive version would have moved the prior median of `σ_input` 1.066 → 1.097.
(Side question for later, NOT for this run: whether adding a Jacobian to a prior already defined in
sampling space is what the engine intends. It affects every run equally, so it does not bias this
comparison.)

### Checks before any launch (all local)

1. **Nesting test** `doublechecks/eqinit_nesting_test.R`: at θ with `σ_init = 1`, the data ll
   (ll − log-Jacobian) of the two arms must agree, all six models, real scripts, production
   likelihood config.
   ✅ **PASS 2026-10-08, all six** (4 θ each): max |Δ| = 0 (SP1, TP2, TP3, Yasso20), 2e-13
   (Yasso07), 1e-13 (Yasso15). Free parameters 6→5, 7→6, 9→8, 20→19, 26→25, 26→25.
   It caught one bug first: the Yasso scripts filled `sigma_ppm` by name from the prior file,
   which still lists `sigma_init` (fixed: only the free names are assigned; no-op in production).
2. **Prior test:** `σ_input` prior marginals identical in the two arms (trap 1).
3. **Pre-flight pushforward** (`preflight_prior_pushforward.R`) with the switch on: blow-up rates,
   no boundary pinning. ⬜ Optional, not run: the arm changes no prior except removing one
   coordinate, and the forward sanity passed in all six.
4. **Launch log** must echo `EQUILIBRIUM INIT: ON` and a free-parameter count one shorter than
   production.

### Prerequisites on Roihu

`git checkout main && git pull` (the clone is on the old branch), and the Fortran recompile (guard
0.9999 → 1.0), with one `R CMD SHLIB` per `.f90`. Then six jobs, same footprint as production.

## What the comparison can and cannot show

- **Level:** both arms will reach it. `σ_input` scales the equilibrium stock directly. If the
  equilibrium arm cannot reach the level, the arm is broken, not the convention.
- **⚠ The equilibrium arm will NOT be flat.** Litter rises ~52% over 1986–2006 and falls ~13% after,
  so an equilibrium start still rises and then bends. **The slowdown after 2006 therefore does not
  discriminate between the arms on its own.** What discriminates is the **size** of the rise:
  whether litter and climate after 1985 can produce +0.399 (1985→2006) and +0.259 (1985→2024)
  without an inherited deficit.
- **It can partly compensate through turnover.** A shorter transit time follows the litter rise
  faster. Expect the arm to move along the input–turnover trade-off. **That movement is a result,**
  and it is the title's claim tested directly.
- **The log-likelihoods are comparable** (same data, error model and σ), and the arms are nested, so
  production ll ≥ equilibrium ll. One parameter is cheap to charge for, but information criteria
  inherit the documented overconfidence (n_eff ≈ 221). **Decide on the trajectory and the rates,
  not on ΔLL alone.**

## Pre-registered (write into the paper's Methods before the results exist)

1. Both arms fit the campaign levels.
2. The equilibrium arm under-produces the 1985→2006 and 1985→2024 rates in **most** models.
3. The equilibrium arm moves along the trade-off: **`σ_input` higher and/or transit time shorter**
   than production, per model.
4. Yasso07, whose single climate modifier rescales time (`MTT = MTT_ref/ξ`), has the most room to
   compensate. If any model closes the gap, expect it to be Yasso07.

**Falsifier:** the equilibrium arm reproduces the observed rates within their intervals, at no
material fit cost, in most models ⇒ the transient start is not needed, and the central claim must
be restated.

## Figures

Observed basis everywhere: the **310 balanced plots, unweighted, whole profile, true observation
years** (`manuscript/figures/obs_basis.R`). Colours: one per model as elsewhere; **solid =
transient, dashed = equilibrium**.

**F16 — the counterfactual (the thesis figure).** DECIDED: try the **ensemble in ONE panel** first
(six thin lines per arm, ensemble median bold); fall back to six panels only if unreadable.
1985–2024. Each panel shows the mean trajectory over the balanced plots for both arms (median plus
90% band from the predictive draws) and the three observed campaign means with CI at their mean
sampling years. A seventh panel, or an inset, shows the litter input (`σ_input·J`) for both arms,
so the reader sees that the equilibrium arm's rise is input-driven.

**F17 — the rates.** DECIDED: build it; whether it is used is open. Dot plot with three columns (1985→2006, 2006→2024, 1985→2024) and
one row per model. Each row shows the transient dot, the equilibrium dot (posterior median + 90%
interval) and the observed rate as a vertical band. This is the quantitative core that F16 shows
qualitatively.

**F18 — displacement along the trade-off.** DECIDED: supplement for now. The plane of `σ_input` (or the
effective flux) against intrinsic transit time, both log scales. One posterior cloud or ellipse per
model and arm, with an arrow from transient to equilibrium, plus the published transit times as
reference ticks. Uses `doublechecks/intrinsic_mrt.R` (⚠ check the echoed RUN_IDs).

**S-figures / tables:**
- **T-arm:** per model, ll at the posterior median and maximum, ΔLL, `σ_input`, transit time, and
  calibration and holdout RMSE for both arms.
- **Campaign residuals:** mean log residual per campaign and arm. Shows where the equilibrium arm
  misses (expected: 1985 too high, or 2024 too low).
- **Projection consequence:** the 2024–2084 sink and equilibrium headroom for both arms. This is
  the inventory-relevant point: does the start-date assumption change the projected sink?

**Optional, decide later:** a decomposition of the production rise into the inherited-deficit part
(production posterior forward with forcing held at 1985) and the forcing part. This is cheap and
needs no new calibration.

## Figure layer — BUILT 2026-10-08, tested on production-as-both-arms

All in `manuscript/figures/`; one command after the rsync:
`bash manuscript/figures/run_eqinit_figures.sh` (equilibrium-arm predictive stage for all six,
then everything below; `SKIP_PREDICTIVE=1` for figures only).

| file | what |
|---|---|
| `run_ids_eqinit.R` | `RID` (production, via run_ids.R) + `RID_EQ` (newest `_eqinit`); partial landing allowed (`EQ_MODELS`) |
| `eqinit_common.R` | extraction + cache `eqinit_comparison.rds` (stamped with both arms' RUN_IDs) |
| `eqinit_draws.R` | per-draw DATA ll (Jacobian removed — it differs between arms) + transit time → `eqinit_draws.rds` |
| `build_F16_eqinit_trajectories.R` | F16 ensemble panel 1985–2084 + litter strip; `_six.png` fallback |
| `build_F17_eqinit_rates.R` | F17 rates, three intervals, posterior 90% |
| `build_F19_eqinit_forecast.R` | **F19 (new) forecast**: change since 2024, sink 2025–44, headroom |
| `build_F18_eqinit_tradeoff.R` | F18 supplement: input vs transit time, arrow per model |
| `build_T_eqinit.R` | `T_eqinit.csv`: ll max/median, ΔLL, RMSE, σ's, MTT, rates, sink, headroom |

Test (production linked as both arms) reproduced the recorded numbers: observed rates
+0.399/+0.117/+0.259; Yasso07/15 1985→2024 +0.260/+0.245 (F5: +0.258/+0.250); headroom 26–42%,
SP1 ≈ 0; Yasso20 MTT 22.1. ⚠ On the production arm itself, 2006→2024 is already negative for
SP1/Yasso15/Yasso20 — so that interval discriminates less than pre-registration item 2 assumed;
read it as a CHANGE between arms, not against the observation alone.
⬜ Posterior comparison (appendix: parameters the equilibrium arm pushes to unreasonable values) —
to design once the calibrations land (Lorenzo).

## Decisions and sequencing (2026-10-08)

- F16 ensemble panel; F17 built, use undecided; F18 supplement.
- **The storyline and manuscript are NOT to be rewritten for this yet** (Lorenzo). Order: run the
  analysis → look at the results → discuss the story → only then touch the text. The comparison will
  probably enter the story, but how is decided after the results.
- Projection comparison: main text or supplement, decided with the story.
