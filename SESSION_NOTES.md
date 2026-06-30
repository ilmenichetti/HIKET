# HIKET Session Notes

Running log of development sessions — append a new section each time.

---

## 2026-06-02 — Transient init implementation, Prior_specs refactor, bug marathon

### Context
Working toward: (1) proper transient initialization for all 6 models, (2) documentation
update, (3) uncalibrated Yasso baseline runs. Session covered 1 and a lot of debugging.

### What was in the Rmd before this session
`Calibration_real_data_transient/documentation/HIKET_calibration.Rmd` existed but
was stale in several ways:
- Said "five models" throughout — TP3 was missing entirely
- §5.4 described steady-state initialization as universal — wrong; had never been a real
  pre-run, and TP3's loop was a no-op
- "nuisance parameters" terminology (CLAUDE.md says use "auxiliary uncertainty parameters")
- No data section
- No uncalibrated baseline section
- Parameter table missing TP3, Komeetta 2024 data listed as pending
The Rmd has NOT been edited yet. All documentation work is still pending.

### Major code changes

#### 1. Transient initialization — full hybrid pre-run implemented
**What it is:** Each model now runs a 68-year pre-run (1917→1985) before the calibration
period. Litter linearly interpolates from a 1917 historical anchor to the 1985 observed
series start. Climate is constant at xi_mean (no pre-1985 climate data exists).

**Endpoint formula (same for all 6 models):**
- `J_1917 = J_full_mean * sigma_init * sigma_input` (historical anchor)
- `J_1985 = J_t0_mean  * sigma_input`              (observed-period start)
- `frac = (i-1)/(68-1)` for i = 1..68
- `J(t) = J_1917 + (J_1985 - J_1917) * frac`

**sigma_init interpretation:** ratio of 1917 litter to contemporary mean. σ_init > 1 =
historically more productive (old-growth pre-intensification); σ_init < 1 = historically
less (post-disturbance recovery). σ_input scales the entire contemporary time series and
is factored out of σ_init cleanly — they have non-overlapping roles.

**Implementation per model:**
- SP1/TP2: inline 68-year R loop using exact per-step solver
- TP3: same, but Euler steps (matches tp3_run's own stepping method)
- Yasso07/15/20: build 68-row AWEN input_df → call `yasso*_run` → extract
  `attr(pre_out, "C_final")` (15-element per-cohort terminal state)
  NOTE: Yasso15/20 Fortran is **single precision** despite `REAL(dp)` declarations —
  `dp = SELECTED_REAL_KIND(6, 37)` = REAL(4) in yasso15_mod. Wrappers correctly
  use `as.single()`.

**Files changed:**
- `Model_functions_real_data_transient/Decomposition_functions/SimpleModels/sp1_wrapper_transient.R`
- `Model_functions_real_data_transient/Decomposition_functions/SimpleModels/tp2_wrapper_transient.R`
- `Model_functions_real_data_transient/Decomposition_functions/SimpleModels/tp3_wrapper_transient.R`
- `Model_functions_real_data_transient/Decomposition_functions/Yasso/yasso07_wrapper_transient.R`
- `Model_functions_real_data_transient/Decomposition_functions/Yasso/yasso15_wrapper_transient.R`
- `Model_functions_real_data_transient/Decomposition_functions/Yasso/yasso20_wrapper_transient.R`
- `Calibration_real_data_transient/run_Yasso15_transient_calibration.R` (litter_means)
- `Calibration_real_data_transient/run_Yasso20_transient_calibration.R` (litter_means)

**litter_means additions for Yasso15/20** (were missing, all others already had them):
- `nwl/fwl/cwl_full_mean`: mean over full observed series → 1917 anchor
- `nwl/fwl/cwl_t0_mean`: mean over first N_PREINIT_SMOOTH=5 years → 1985 endpoint

#### 2. Prior_specs refactor
Priors were hardcoded in each calibration script and `Priors_model_matching.R` was
sourced at runtime, causing shadowing bugs (see below).

**New structure:** `Prior_specs/` directory at project root:
```
Prior_specs/
  SP1_priors.R       # hand-set, static
  TP2_priors.R       # hand-set, static
  TP3_priors.R       # hand-set, static
  Yasso07_priors.R   # literature CIs, static (not auto-generated)
  Yasso15_priors.R   # FMI posterior, regenerate with WRITE_PRIORS=TRUE
  Yasso20_priors.R   # FMI posterior, regenerate with WRITE_PRIORS=TRUE
```

Each file defines `MODEL_FREE_DEFAULTS` (physical-space prior centres) and
`MODEL_SIGMA_PPM` (prior SDs in transformed space; overrides only — fractions default 1.0).

Calibration scripts source the file at top, then:
- `free_defaults <- MODEL_FREE_DEFAULTS`
- `sigma_ppm <- setNames(rep(1.0,N_FREE), FREE_NAMES)` + override for Yasso
- `sigma_ppm <- setNames(MODEL_SIGMA_PPM[FREE_NAMES], FREE_NAMES)` for SP1/TP2/TP3

`Priors_model_matching.R` is **no longer sourced at runtime**. To regenerate Yasso15/20
prior files from FMI data: `WRITE_PRIORS <- TRUE; source("Priors_model_matching.R")`.

#### 3. Bug fixes (all required to get Yasso15/20 running)

| Bug | Root cause | Fix |
|---|---|---|
| TP3 pre-run no-op | Loop started at steady state, used constant J_pre → never moved | Made J time-varying (was the whole point) |
| Yasso15/20 segfault | `Priors_model_matching.R` sourced at runtime, redefining `yasso15_run` with legacy version lacking `C_final` → `attr(pre_out,"C_final")` = NULL → NULL C_init passed to Fortran | Remove Priors from runtime path (Prior_specs refactor) |
| Yasso15/20 segfault (earlier) | `dyn.load` pointed at non-transient `yasso15.so` (Feb 17 compile, no `C_final` in `yasso15_run_r`) | Change to transient `yasso15.so` (May 25 compile, has `C_final`) |
| Wrong numeric output | Wrapper used `as.double()` but Fortran is `REAL(dp)=REAL(4)` (single) | Revert wrapper to `as.single()` |

**Key permanent gotcha for Yasso15/20:**
`yasso15.f90` defines `dp = SELECTED_REAL_KIND(6, 37)` = single precision.
Despite `REAL(dp)` declarations everywhere, all Fortran arguments are 4-byte.
The R wrapper must use `as.single()` / `single()`. Yasso07 is genuine double
(`dp = SELECTED_REAL_KIND(15, 307)`).

### What still needs to happen

#### Baseline runs (do before or alongside documentation)
Uncalibrated Yasso07/15/20 runs at published defaults, to serve as reference
baselines for the intercomparison. Script `run_published_yasso07_baseline.R`
exists in `Calibration_real_data_transient/` but was written for Yasso07 only.
Need:
1. Equivalent scripts for Yasso15 and Yasso20
2. Equivalent scripts for SP1/TP2/TP3 at their prior centres (optional but useful)
These should run BEFORE calibration results in the document narrative.

#### Documentation — HIKET_calibration.Rmd
File: `Calibration_real_data_transient/documentation/HIKET_calibration.Rmd`
Target audience: generalist scientists + Yasso development team.

**Structural changes needed:**
1. Update model count throughout: 5 → 6 (add TP3)
2. Add TP3 to model table: 3-pool ASH cascade, 10 free params, same Yasso07 xi form
3. Add complexity-ladder narrative update (SP1→TP2→TP3→Yasso07 for pool count axis)
4. Add Data section (currently missing): Finnish NFI plots ~350-400, VMI8/Biosoil/
   Komeetta SOC campaigns, Tupek litter dataset, gridded climate 1961-2025
5. Add Uncalibrated baseline section (before calibration section): what published
   defaults predict on Finnish data; what systematic biases exist
6. Rewrite §5 Transient initialization:
   - The approach is now a TRUE hybrid pre-run for all models (not steady-state)
   - Litter linearly interpolated 1917→1985
   - sigma_init as historical productivity ratio (document clearly for FMI audience)
   - sigma_input as global litter scaling (separate non-overlapping role)
   - Climate constant at xi_mean during pre-run (data constraint, not approximation)
   - Yasso07/15/20 use C_final from chained Fortran run (not analytical steady state)
   - All 6 models now comparable — same initialization protocol
7. Fix "nuisance parameters" → "auxiliary uncertainty parameters"
8. Check Komeetta 2024 data: confirm whether it's actually in the current calibration
   data or still pending

**Tone/style for documentation:**
- Keep math but add plain-language orientation paragraphs before each section
- "Transient" should be explained in plain terms early on
- AWEN fractionation needs a brief explanation
- The complexity ladder narrative (why these 6 models?) should be front-loaded

### Run status as of end of session
| Model | Posterior exists | Notes |
|---|---|---|
| SP1 | Yes (2026-05-19) | Pre-runs with new transient init |
| TP2 | Yes (2026-05-19) | Same |
| TP3 | No | Scripts complete, no run yet |
| Yasso07 | Yes (2026-05-19) | |
| Yasso15 | Yes (2026-05-19) | Crash fixed, reruns needed with new transient init |
| Yasso20 | Yes (2026-05-19, 2026-05-26) | Same |

All existing posteriors were generated WITHOUT the proper pre-run (steady-state only or
no-op). They need to be re-run on Puhti to incorporate the hybrid transient initialization.

### Recommended next session agenda
Option A (baseline first): write Yasso15/Yasso20 baseline scripts, run all 3 Yasso
baselines, then start documentation with the full picture available.
Option B (documentation first): write the methods section while code changes are fresh,
then add baseline results when available.
Recommendation: Option A — the baselines are needed for a complete document, and
writing the uncalibrated baseline section without those results means writing placeholders.

---

## 2026-06-03 — Baseline scripts, Yasso20 mass-balance discovery, full documentation rewrite

### Context
Calibration jobs submitted on Puhti (status unknown — queued/running). This session:
built the Yasso baseline scripts, discovered a mass-balance issue in the published
Yasso20 parameters, and rewrote both documentation files. Decision taken: baselines
first (Option A), but most of the session went to documentation + the Yasso20 finding.

### Baseline scripts (DONE)
- `run_yasso07_baseline.R` — full rewrite of the legacy `run_published_yasso07_baseline.R`
  (legacy file left untouched on disk; can be deleted when ready).
- `run_yasso15_baseline.R`, `run_yasso20_baseline.R` — new.
- `hiket_baselines.sh` — single SLURM job running all three sequentially (4h, 40 cpus).
- Output dirs renamed `*_published/` → `*_baseline/` (both transient and non-transient
  diagnostics trees; folders were empty).
- New plots vs legacy: **mean absolute SOC trajectory** (predicted ±1 SE + observed
  campaign means) and **mean annual ΔSOC** (tC/ha/yr, matches calibrated predictive
  form). Colour scheme switched KA → `kasvup_tyyppi` to match calibrated outputs.
- Init: steady-state (at published defaults σ_init=σ_input=1, so SS ≡ pre-run endpoint).
- SP1/TP2/TP3 baselines NOT written (no published global params; would use prior centres).

### Yasso20 mass-balance discovery (IMPORTANT)
Investigating why the Yasso20 baseline "explodes" at the tail, we traced it to the
**published Viskari (2022) Yasso20 parameters violating pool-level mass conservation**:
- Donor pools A and N each redistribute `Σp_ij + p_H = 1.0042` (> 1) → respired
  fraction −0.0042. Root cause: `pAW` and `pNA` are fixed at exactly 1.0, and
  `p_H = 0.0042` is applied on top.
- Verified using **FMI sources only**: `ParY20.dat` ≡ `Yasso20_sample_parameters.rda`
  (identical 35-vectors), matrix code from FMI `Yasso20.r` / `yassofortran.f90`.
- Our HIKET Fortran has a diagonal-dominance guard (a 2025 addition, not in FMI code)
  that flags this and falls back to Euler+clamp, which then accumulates error over the
  ~20-yr run → the "explosion". FMI code just runs `expm` and tolerates it; the
  violation is LOCAL (per-pool) and masked system-wide by W's strong respiration, so a
  total-stock run still shows carbon DECREASING — which is why it was never noticed.
- Report written for FMI: `documentation/Yasso20_mass_balance_note.md` — bombproof,
  every input linked to an FMI source, neutral tone, ends asking the dev team whether
  `pAW`/`pNA` should be re-normalised to `1 − p_H`. **Caveat to check before sending:**
  confirm with Toni that `ParY20.dat`/`.rda` IS the parameter set from the 2022 GMD paper.

### Yasso20 baseline parameter decision (OPEN)
- `Prior_specs/Yasso20_priors.R`: added `YASSO20_PUBLISHED_DEFAULTS` (Viskari fractions
  incl. the 6 structural zeros + p_EA=0) — kept for reference but **reverted the baseline
  script to use `YASSO20_FREE_DEFAULTS`** (Yasso15 fractions) so it runs stably.
- User wants to decide later whether the Yasso20 baseline should instead reproduce the
  original FMI result using the unmodified FMI code+params. Left provisional; documented
  as such in both the data... no, in `HIKET_calibration.Rmd` baseline section ("Status").

### Documentation rewrite (DONE)
**`HIKET_calibration.Rmd`** — full rewrite based on the code:
- 5→6 models (TP3 added everywhere: 3-pool A→S→H sequential cascade, Euler stepping,
  10 free params, Yasso07 xi form).
- New: Data section, Uncalibrated baselines section, Yasso-family overview with a
  **TikZ pool diagram** (all 12 p_ij arrows + humification + respiration + litter) and a
  parameter summary table.
- Rewrote §Transient initialisation as the central story: 68-yr pre-run 1917→1985,
  σ_init = historical productivity ratio, why steady-state init fails the observed SOC
  trend. Storyline confirmed with user (see below).
- New §"Fortran implementation: differences from the FMI reference code" — documents all
  7 HIKET modifications (diagonal-dominance guard, matrixexp scaling cap, Taylor 10→20,
  C_final output, xi precompute in R, single vs double precision, param ordering, r sign).
- "nuisance" → "auxiliary uncertainty parameters" throughout.

**Code-vs-docs review findings (corrected in the Rmd):**
- **Free-param counts were wrong** (stick_break gives K free params per column, not K−1;
  sum strictly < B). Yasso07 18→**20**, Yasso15/20 28→**26**. Fixed in all locations.
  NB: the Yasso15 calibration SCRIPT's own header comment ("28 total") is also wrong —
  left for the user to fix in code (out of scope this session).
- Stick-breaking description corrected (was K−1 "third determined").
- **Holdout split** (present in all 6 calibration scripts, excluded from likelihood,
  separate predictive metrics) was undocumented — added a subsection.

**`HIKET_data_preparation.Rmd`** — updated stale parts from `Data/Data_work.R`:
- Komeetta 2024 now INTEGRATED (was "pending"): section 1.0b details — layer mapping,
  `KOM_ORGANIC_ALL_SUBLAYERS=TRUE`, 58 Mg/ha QC check, >150 Mg/ha per-layer filter
  (plots 55652/3314/71471), input_monthly extended to 2024 with copied 2023 litter.
- Model list → all six; SOC "two campaigns" → three; month-count 39→40 yr.
- Left PLACEHOLDERs for counts needing a fresh `Data_work.R` run (input_raw_monthly
  rows, site_raw total/calib_ready, which OFH sublayers 1985/2006 measured).

Both Rmd files render to PDF cleanly (tikz fixed via header-includes load order).

### Storyline (confirmed with user)
Finnish forest SOC accumulates 1985→2006 and Komeetta 2024 overall still saturating
(layer-level nuances aside); the saturation is real and politically relevant (Sweden etc.
don't report it). Transient init is what lets the models track the initial trend rather
than predicting flat SOC from a wrong steady-state start. Conclusions come from the
ENSEMBLE; Yasso20 is less reliable and not over-weighted. Plot-level litter (Tupek) is
fairly flat, unlike older NFI inputs which rose — to be discussed in the data doc.

### Run status as of end of session
Unchanged from 2026-06-02 (existing posteriors predate the proper pre-run; calibration
re-runs submitted on Puhti, outcome pending).

### Recommended next session agenda
1. Check Puhti calibration job status (`squeue -u menichet`); pull results if done.
2. Launch the baseline job on Puhti: `sbatch Calibration_real_data_transient/hiket_baselines.sh`
   (after `git pull` + Fortran recompile if any .f90 changed — none did this session).
3. Then either: debug whatever the runs surface, OR start assembling a draft report.
4. Decide the Yasso20 baseline parameterisation (Yasso15 fractions vs original FMI code).
5. Optional: send `Yasso20_mass_balance_note.md` to Toni (after confirming the param-file
   provenance caveat).

---

## 2026-06-12 — Production-run inspection (all 6), local downstream pipeline, coverage finding

### Context
The full post-homogenization production calibration finished: SP1/TP2/TP3 `20260608_*`
(Jun 8) and Yasso07/15/20 re-run `20260611_*` (Jun 11, after the stale-`.so` fix). User
flagged that R-hat looked "surprisingly good." Session: verify that's real not a trap,
run Stages 2–4, look at headline results.

### Convergence inspection — all 6
- **Yasso07/15/20:** all params R-hat < 1.05, ESS 1265–2369, inf_rate ~0%, forward sanity
  PASS, distinct non-degenerate ll-at-defaults (−4868/−5083/−7497) → NOT the stale-`.so`
  trap. The 12 fractions converged in R-hat BUT stay non-identified by the right measures:
  prior-pinned (13–17/26 within 0.5σ of prior centre), ridge anti-correlations ≈ −0.98
  (`p_EN↔p_EA`, `p_NE↔p_NA`), KL < 1 nat. **Key correction:** the old guardrail "fraction
  R-hat must stay poor" conflated non-identifiability with poor *mixing*; the old bad R-hat
  was the `beta2` detonation freezing the sampler. Fixed → good R-hat + intact
  non-identifiability. CLAUDE.md §"Convergence expectations" rewritten.
- **SP1/TP2/TP3:** TP2 clean (max R-hat 1.006). TP3 R-hat clean but ENTIRE 3-pool kinetics
  gain ~0 nats — `sigma_input` (KL ≈ 21) absorbs everything; fits by rescaling inputs, not
  learning rates (8–17% forward-blowup, documented, tolerable). SP1 `alpha` R-hat 1.09 /
  ESS 208 = the lone genuine R-hat signal, but it's the MOST identified param (KL ≈ 19),
  just heavy-right-tailed (97.5% ~6× median) → hard to sample, benign.
- **Yasso20 chain-2 `final_ll` WARN (945):** single-last-draw artifact. Loaded posterior
  (chain-major rbind, 5×45006), per-chain param means agree to <2.5% → not displaced. (The
  saved posterior is a params-only matrix; no `LL` column, so a fully ll-based check would
  need re-evaluation — not worth it given param agreement.)

### Downstream pipeline — RAN LOCALLY, end-to-end, zero errors (~12 min)
- Fixed two known rough edges first: added `rownames(inputs_proj/clim_proj) <- NULL` to
  SP1/TP2/Yasso07/Yasso15 predictive scripts (TP3/Yasso20 already had it); added
  `nrow/ncol(imp_mat) < 2` guard around the multimodel RF heatmap. All parse-checked.
- `Rscript run_hiket_pipeline.R --skip-calibration` auto-picked correct RUN_IDs, ran
  predictive (100 draws/model) + residuals (with 897/224 calib/holdout split, RF OOB R²
  ~0.42) + multimodel. **Local is the right venue** — SLURM queue would dwarf the compute.
- **Sync gotcha hit:** predictive hard-loads `Data/model_inputs/<MODEL>_inputs_<RUN_ID>.rds`
  (no fallback). Had to sync `Data/model_inputs/` from Puhti in addition to `runs/`+
  `diagnostics/`.

### Headline results (447 plots)
All models R² 0.05–0.11 (calib) / 0.02–0.05 (holdout), RMSE ~50–62, slight + bias.
**TP3 best, Yasso20 worst, spread small** — on-thesis (structure → modest skill diffs vs
generally low SOC predictability).

### PAUSED — open methods decision
"95% cov" ≈0.05 is a methods artifact: predictive interval = 100-draw spread of the model
**mean**, omits the multiplicative-normal obs error (`sigma_obs` ≈ 0.4723, loaded but never
injected). NOT under-dispersion; R²/RMSE/bias/rankings unaffected. **Recommended (option 1,
not applied):** add true posterior-predictive coverage to all 6 predictive scripts —
`total_soc * exp(rnorm(0, sigma_obs))` before quantiles, keep both numbers. See memory
[[predictive-coverage-artifact]].

### Uncommitted
All this session's edits (CLAUDE.md, 4 predictive scripts, multimodel script,
SESSION_NOTES.md) are uncommitted, on top of pre-existing staged changes. Nothing committed.

### Next steps
1. Decide on the coverage fix (option 1 vs relabel), then rerun pipeline.
2. Commit the session's code + doc changes (CLAUDE.md doc update can be its own commit).
3. Inspect multimodel figures (obs-vs-pred, residuals, ΔSOC, RF heatmaps) for structure.
4. Surface the `sigma_input` cross-model contrast (TP3 ≈21 nats; Yasso20 inputs ↑3×; SP1
   inputs ↓3×) in the writeup — it's a structural-intercomparison result.

---

## 2026-06-16 — Results+Discussion writing; doc figures/sections; revision workflow

### Discussion findings (now in `HIKET_calibration.Rmd`)
- **Over-prediction = likelihood artifact.** All 6 models over-predict slightly because
  the multiplicative-normal error scales SD on the prediction (penalty on fractional error
  δ ~ δ²/(1+δ)², asymmetric → hedges upward). Log-scale re-scoring of the `*_residuals_*.csv`:
  mean log-residual −0.10..−0.16 (≈10–16% multiplicative offset), ~identical for mean vs
  median predictor, tracks R² inversely (TP3/TP2 least, Yasso20 most). Only the bias is
  likelihood-sensitive; rankings + KL/identifiability invariant. → Rmd §`sec:biasdisc`.
  Memory [[over-prediction-likelihood-artifact]].
- **TP3 forcing oscillation.** TP3 (mildly TP2) rings year-to-year in mean-SOC (fig:meansoc);
  source is climate not litter — the new driver plot shows T/precip oscillatory, litter
  smooth. Low-pass-filter argument: fast prior-set pools + inflated σ_input (≈13) → cutoff
  above annual forcing → transmits ξ. One non-identifiability, three faces (KL / ringing /
  invalid-proposal 8–17%). → Rmd §`sec:tp3osc`. Memory [[tp3-forcing-oscillation]].
- **Residual RF.** `basal_area_85` is the #1 residual driver in 5/6 models (#2 TP3); soil
  chem + climate next; diffuse (<~5–9% each), shared ordering across models; OOB R²
  0.34–0.52; `dev_class_85` weak; no direct age var. Stand-age reading OPEN. → Rmd
  §`sec:residuals` + fig:rfheat + OOB table. Memory [[residual-rf-basal-area]].

### Doc/code changes
- New **3-panel driver figure** (T, precip, total litter; cross-plot mean ±95% CI, 1985
  dropped). Wired into `run_multimodel_comparison.R` §4b "Plot B2" → regenerates every run
  (`multimodel_drivers_<COMPID>.png`); chunk `fig:drivers` added to the Rmd.
- Reworded the Yasso pre-run "genuine Fortran run, not an analytical approximation" sentence
  — it's the SAME analytical (matrix) solver as runtime; the real point is `C_final`
  per-cohort bookkeeping. (User caught this; their mental model "analytical SS then transient"
  is correct.)
- PDF re-rendered cleanly (exit 0, no undefined refs).
- **Placeholders left in the Rmd to decide together:** stand-age reading; TP3 genuine-vs-
  Euler-artifact; residual Yasso over-accumulation (lognormal refit).

### Revision workflow established (next step)
- User does a full read-pass on a COPY: `documentation/HIKET_calibration_annotated.pdf`,
  highlight+Note in Preview. Extract with **PyMuPDF (installed this session)** —
  `page.annots()`, highlight rect → anchor text, nearest same-page Text annot → comment.
  Verified on 12 notes. Never render to `_annotated.pdf`.
- **Open contradiction flagged by the user (note #8):** the doc says litter is
  "approximately flat over the observation period", but fig:drivers panel (c) shows it
  rising ~2.1→3.0 to a mid-2000s peak then plateauing. Must reconcile — load-bearing for the
  flat-vs-rising tension in [[transient-phase-open-question]].
- Memory [[revision-pass-workflow]] holds the 12-note snapshot + themes.

### Next steps
1. **New session:** re-extract all annotations from `_annotated.pdf`, build a numbered
   agenda, work through them one at a time (edit → confirm → re-render at breakpoints).
2. Resolve the flat-vs-rising litter contradiction.
3. Still pending from prior session: coverage fix decision; commit the uncommitted edits.

---

## 2026-06-30 — NextGenC gridded SOC maps (kriging) deliverable

### Context
New deliverable for the NextGenC report: interpolated Finland SOC maps (one per year)
from the per-plot posterior-predictive output, destined for a Zenodo deposit. Spun off
from [[nextgenc-soc-report]] (which had the per-plot matrices + trajectory plot).

### What was built (all in `Reporting/NextgenC_report/`)
- **`build_soc_maps.R`** — ordinary kriging in log space of per-plot SOC → 2 km grid
  (EPSG:3067), masked to Finland. Per model + 6-model ensemble; multiband GeoTIFF
  (`_mean.tif` / `_sd.tif`, 40 bands = years 1985–2024). One pooled variogram/model.
  Per-model SD = kriging var + interpolated model SD (option C); ensemble SD = within +
  between-model (law of total variance). 14 GeoTIFFs (~99 MB), ~25 min runtime.
- **`make_thumbnails.R`** — 14 PNG previews (time-averaged layers, 2040×2000 px, 300 dpi,
  ~320 KB each), `SOC_maps/thumbnails/`.
- **`SOC_maps/`** = self-contained Zenodo deposit: GeoTIFFs + copied tabular matrices
  (CSV/ODS) + thumbnails + `README.md` + `LICENSE` (**CC-BY-4.0**).
- **`SOC_maps_README.{tex,pdf}`** — 2-page methods note (procedural).

### Key finding (documented in memory + READMEs)
Temporal-mean log-SOC is **~85–90 % nugget** at the ~27.5 km plot spacing → SOC is
spatially weakly autocorrelated at the network scale. Maps are plot-anchored "bullseyes"
reverting to the national mean, NOT smooth regional gradients — the spatial echo of HIKET's
low point-predictability (R² ~0.05–0.11). Nugget ~identical across models → property of
sites + single-pit sampling, not the model. SD layer honestly inflates between plots.
Discussed (not implemented): this is NOT detectable "spatial chaos" (low-dim determinism);
the defensible framing is high-dimensional/path-dependent determinism, or multifractal
heterogeneity below the sampling scale. Upgrade path if wanted: external-drift kriging on
gridded climate. Also: between-model structural spread is SMALL vs within-model SD for SOC
*stock level* (model means 90.6–95.4 tC/ha) — models agree on level, differ on dynamics.

### Git / repo hygiene
`.gitignore` now excludes `Reporting/NextgenC_report/SOC_maps/*.{tif,csv,ods}` (large +
regenerable + Zenodo-bound; canonical CSV/ODS kept in the parent folder). Committed:
scripts, docs, README, LICENSE, thumbnails, parent matrices.

### Caveat carried forward
TP3 maps use the pre-fix (Euler) posterior; regenerate all 3 TP3 products
(`build_soc_matrices.R` → `build_soc_maps.R` → `make_thumbnails.R`) when the exact-
integrator run syncs from Puhti. SOC levels change negligibly. See [[tp3-forcing-oscillation]].

### Extension (same day, commit dfffae2): peat/water masking + delta map
- **Peat mask** (`Data/Peat_extension/`): Luke MS-NFI `paatyyppi` 16 m, INSPIRE-Atom
  fetch. 2 km peat fraction = peat / **whole cell** (water in denom). *Bug fixed:*
  dividing by classified-land-only made Saimaa lakeland read as ~100% peat (water in
  `paatyyppi` is **NaN**, not 32766). Mask >50% peat (~9%); drop Tu/Ve input plots.
- **Water mask** (`Data/Water_extension/`): SYKE Ranta10 lakes (direct URL), 2 km water
  fraction, mask >50% (~3%). MS-NFI has no water; Natural Earth far too coarse.
- **Rendering:** thumbnails show **peat black, water/sea white**. SD palette → Purples.
- **Delta map** (`make_delta_map.R`): ensemble SOC change 2015–2024 (OLS slope), RdYlBu.
  *Bug fixed:* extreme cells rendered white (terra clips at `range=` edge) → looked like
  holes; fix = clamp + pad range. ~76% gaining, SW loss hotspot.
- **Diagnostic:** `diagnostic_semivariogram.R` (outside deposit) shows the ~88% nugget.
- **Zenodo:** zip `HIKET_SOC_maps_Finland_v1.zip` (~112 MB, gitignored); description
  drafted (credits HIKET + NextGenC; CC-BY; attributes Luke + SYKE). Upload pending user.
