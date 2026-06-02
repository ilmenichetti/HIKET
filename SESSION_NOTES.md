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
