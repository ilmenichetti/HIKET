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
