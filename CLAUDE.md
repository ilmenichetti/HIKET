# HIKET Project — Claude Code Context

## Project overview

**HIKET** is a Bayesian calibration and structural intercomparison of five soil
organic carbon (SOC) decomposition models applied to Finnish forest ecosystems.
The scientific question: do divergent model projections of Finnish forest carbon
sink saturation reflect genuine structural differences, or calibration artefacts?

**Models in scope:** SP1, TP2, TP3, Yasso07, Yasso15, Yasso20  
**Application:** Finnish National Forest Inventory (~350–400 permanent plots)  
**End use:** Finnish greenhouse gas inventory

---

## Repository layout

```
HIKET/
├── Calibration_real_data_transient/
│   ├── calib_config.R                        # shared MCMC settings (edit here)
│   ├── calibration_engine.R                  # base engine (non-transient)
│   ├── calibration_engine_transient.R        # transient extension (sources base)
│   ├── run_hiket_pipeline.R                  # top-level orchestrator
│   ├── run_SP1_transient_calibration.R       # per-model calibration scripts
│   ├── run_TP2_transient_calibration.R
│   ├── run_TP3_transient_calibration.R
│   ├── run_Yasso07_transient_calibration.R
│   ├── run_Yasso15_transient_calibration.R
│   ├── run_Yasso20_transient_calibration.R
│   ├── run_SP1_transient_predictive.R        # per-model predictive scripts
│   ├── run_TP2_transient_predictive.R
│   ├── run_TP3_transient_predictive.R
│   ├── run_Yasso07_transient_predictive.R
│   ├── run_Yasso15_transient_predictive.R
│   ├── run_Yasso20_transient_predictive.R
│   ├── run_residual_analysis.R               # model-agnostic residual stage
│   ├── run_multimodel_comparison.R           # cross-model comparison & figures
│   ├── hiket_sp1.sh / hiket_tp2.sh / ...    # SLURM job scripts
│   ├── runs/                                 # posterior RDS files (gitignored)
│   ├── diagnostics/<MODEL>/                  # per-model PNGs, CSVs, RDS
│   └── diagnostics/multimodel/              # cross-model figures & metrics
├── Data/
│   ├── model_inputs/                         # input bundles (gitignored)
│   └── ...
├── Model_functions_real_data/
│   └── Decomposition_functions/Yasso/        # Fortran .f90 and compiled .so
├── Reporting/
│   └── NextgenC_report/                      # NextGenC deliverables
│       ├── build_soc_matrices.R              # per-plot SOC mean/sd matrices (CSV+ODS)
│       ├── build_soc_maps.R                  # kriged 2 km SOC maps, EPSG:3035 EEA-snapped (GeoTIFF, +ensemble)
│       ├── make_thumbnails.R                 # time-averaged PNG previews
│       ├── ZENODO_description.md             # paste-ready Zenodo deposit description
│       ├── SOC_maps_README.{tex,pdf}         # maps methods note
│       └── SOC_maps/                          # Zenodo deposit (maps gitignored; CC-BY-4.0)
├── Data_work.R                               # upstream data preparation
├── Additional_layers.R
└── HIKET_calibration.Rmd                     # methods documentation
```

**Gitignored:** `Data/`, `runs/`, `diagnostics/**/*.rds`, `diagnostics/**/*.csv`,
Fortran `.so` binaries, `*.rds`/`*.csv` in diagnostics,
`Reporting/NextgenC_report/SOC_maps/*.{tif,csv,ods}` (large/regenerable; Zenodo-bound).

---

## MCMC settings (`calib_config.R`)

| Setting | Value | Notes |
|---|---|---|
| `N_CHAINS` | 5 | DEzs sampler (BayesianTools) |
| `N_ITER` | 50 000 | per chain |
| `N_BURNIN` | 5 000 | |
| `N_LOG` | 200 | progress log interval |
| `N_PLOTS_TEST` | `NA` | `NA` = full dataset; set e.g. `20L` for quick tests |

---

## Environments

### Mac (development)
- Architecture: arm64, R 4.3
- Run scripts interactively or via `Rscript`
- Fortran `.so` files compiled locally for arm64 — **do not push to git**

### CSC Puhti HPC (production)
- Path: `/scratch/project_2019134/HIKET/`
- Module: `module load r-env`
- Fortran compilation: `apptainer_wrapper exec R CMD SHLIB <file>.f90`
- SLURM user: `menichet`
- **⚠ DECOMMISSION ~end July 2026.** Plan: finish the in-flight TP3 exact-
  integrator run + downstream + docs on Puhti and **close out HIKET here first**;
  do NOT start new development on Puhti. Further work migrates to **Roihu** (CSC's
  Puhti successor) — re-clone the repo, recompile the Fortran `.so` there (each
  `.f90` in its own `SHLIB` call, see Fortran section), and update paths/module/
  SLURM details below once the Roihu environment is set up.

### Sync workflow
```bash
# Code → Puhti
git add ... && git commit -m "..." && git push
# [on Puhti]
cd /scratch/project_2019134/HIKET/ && git pull
# Recompile .so if any .f90 changed (see Fortran section below)

# Results → Mac
rsync -av menichet@puhti.csc.fi:/scratch/project_2019134/HIKET/Calibration_real_data_transient/runs/ ./Calibration_real_data_transient/runs/
rsync -av menichet@puhti.csc.fi:/scratch/project_2019134/HIKET/Calibration_real_data_transient/diagnostics/ ./Calibration_real_data_transient/diagnostics/
```

---

## Running the pipeline

### Stage 1 — Calibration (SLURM, one job per model)
```bash
# From /scratch/project_2019134/HIKET/
sbatch Calibration_real_data_transient/hiket_sp1.sh
sbatch Calibration_real_data_transient/hiket_tp2.sh
sbatch Calibration_real_data_transient/hiket_tp3.sh
sbatch Calibration_real_data_transient/hiket_yasso07.sh
sbatch Calibration_real_data_transient/hiket_yasso15.sh
sbatch Calibration_real_data_transient/hiket_yasso20.sh

squeue -u menichet    # monitor
```

### Stages 2–4 — Predictive + Residuals + Comparison (after calibration)
```bash
module load r-env
apptainer_wrapper exec Rscript --no-save \
  Calibration_real_data_transient/run_hiket_pipeline.R --skip-calibration
```

### Fortran recompilation (after any `.f90` change or fresh clone)
**Compile each `.f90` in its OWN `SHLIB` call.** A combined
`R CMD SHLIB yasso07.f90 yasso15.f90` links both objects into a single
`yasso07.so` and never creates `yasso15.so` — which breaks Yasso15 AND Yasso20
(Yasso20 loads `yasso15.so`; it shares Yasso15's Fortran), and risks symbol
shadowing. The `.so` files are gitignored, so this MUST be redone on Puhti after
any `.f90` change — a stale `.so` silently returns garbage (e.g. zero `C_init`
from the transient init → see the Yasso07 non-convergence post-mortem, 2026-06).
```bash
cd Model_functions_real_data_transient/Decomposition_functions/Yasso/
rm -f yasso07.so yasso07.o yasso07_mod.mod yasso15.so yasso15.o yasso15_mod.mod
apptainer_wrapper exec R CMD SHLIB yasso07.f90    # separate call
apptainer_wrapper exec R CMD SHLIB yasso15.f90    # separate call
```
Quick check a `.so` is current: `yasso07_transient_init` at
`YASSO07_DEFAULT_PARAMS` should give ~69 tC/ha; ~0 means the binary is stale.

---

## Key design decisions

### Intercomparison design principle
All models must have **equal degrees of freedom** in flow fractions. Any
structural asymmetry (fixed vs. free parameters) confounds climate computation
differences with calibration asymmetry and makes results uninterpretable.

### Yasso20 structural fix
All 12 inter-pool transfer fractions are **free parameters** (previously 6 were
fixed as structural zeros and 3 algebraically derived following Viskari 2022).
Prior centres use Yasso15 published defaults. The original Viskari (2022)
constraints reflected a warm-site global dataset, not Finnish boreal data.

### Transient initialisation (`transient_init = TRUE`)
Each MCMC draw initialises carbon pools via a pre-run from `PREINIT_YEAR`
(1917) to VMI8 (1985). `sigma_init` is a **proportional SD on the initial
carbon state**, not a likelihood patch. Parameter uncertainty propagates into
both the initial state and forward dynamics — this legitimately widens
posterior predictive uncertainty; it is not a bug.

### Jensen's inequality correction
`xi(mean_climate) ≠ mean(xi_annual)` for nonlinear climate response functions.
`compute_xi_mean()` must be used instead of `mean(xi_array)`. This corrects
5–15% bias in boreal settings and is applied in all model engines.

### Uniform exact integration — TP3 fixed (2026-06-17)
All six models now integrate their linear pool system **exactly**: SP1/TP2
closed-form per-step; **TP3** closed-form matrix exponential of its lower-
triangular cascade generator (divided differences of eigenvalues
`-aA*xi,-aS*xi,-aH`); Yasso07/15/20 Fortran matrix exponential. **TP3 was the
sole approximate integrator** — it used explicit forward Euler, which rings
when the active pool turns over faster than the 1-yr step (`alpha_A*xi>1`;
median ~1.1, warm years ~1.6, 95th pct ~3.3 → Euler factor `(1-k)<0` →
overshoot). That produced the interannual oscillation previously attributed to
a "low-pass climate response"; it was a **numerical artefact**. Fixed in
`Model_functions_real_data_transient/.../SimpleModels/tp3_wrapper_transient.R`
(`.tp3_step`). Verified by `doublechecks/test_tp3_exact.R`: closed-form ==
reference matrix-exp to 1.6e-12; under constant forcing Euler rings (roughness
0.28) while exact is monotonic (0.006), mean SOC identical (40.90 vs 40.88).
Diagnostic clincher: TP2's active pool is *faster* (`alpha_A*xi~1.5`) yet smooth
because it integrates *exactly*. **Consequence: TP3 must be re-calibrated** for
the doc's TP3 posterior/figures (RUN 20260608) to reflect the exact integrator;
stocks/skill change negligibly, the ringing disappears. The only remaining non-
exact path is the Yasso diagonal-dominance **Euler fallback** (a guard that the
stick-breaking budget keeps from firing).

**Deployment correction (2026-06-30):** the exact integrator was NOT actually
on Puhti for the "2026-06-17" re-run — the `.tp3_step` rewrite sat *uncommitted*
in the Mac working tree, so Puhti (commit 92da86f) re-ran **Euler**; the
20260617_085813 posterior is Euler-based. The local exact code also carried a
bug: `.tp3_step` called with named scalars (`tp3_transient_init` passes `C["A"]`)
produced compound names (`A.A`,…) so `C["A"]` next step → `NA`, `C_init` all NA,
every plot dropped (predictive stage failed). Fixed by `unname()`-ing the scalar
inputs at entry; committed + pushed as **ff216ce** (exact integrator + bug fix).
**The genuine exact-integrator TP3 calibration is RUNNING on Puhti as of
2026-06-30** (job 35323262, from ff216ce). When it finishes follow the TP3
re-calibration steps in "Known outstanding items".

### Likelihood / error model
Multiplicative-normal (proportional) error model. Asymmetric penalisation of
underprediction is a known consequence — noted in methods.

### Convergence expectations
With ~19 free parameters and only 2 SOC observations per plot, the 12 flow
fractions are **structurally non-identifiable**. The correct signature of this is
**not** poor R-hat — that conflates non-identifiability with poor *mixing*. Once
mixing is healthy (the `beta2` detonation fixed; inf_rate ~0%), DEzs traverses
the non-identified ridge cleanly, so fraction R-hat can be **good (≈1.0) even
though the fractions are not data-identified**. The genuine non-identifiability
signature is: (1) fractions **prior-pinned** (most within 0.5σ of prior centre;
posterior ≈ prior), (2) **near-perfect anti-correlation ridges** between fractions
out of the same pool (e.g. `p_EN↔p_EA`, `p_NE↔p_NA` ≈ −0.98 — data constrains the
net flux out of a pool, not the split), and (3) **low KL(posterior‖prior) < 1 nat**
for fractions. Confirmed in the 2026-06-11 production run (Yasso07/15/20): all
params R-hat < 1.05, yet fractions show all three non-identifiability signatures —
good convergence is the *expected* result of fixing mixing, not a red flag.

### Prior homogenization across models — ✅ IMPLEMENTED (2026-06-06)
A three-tier prior scheme, now **applied to all six `Prior_specs/*_priors.R`** —
see `Prior_specs/PRIOR_HOMOGENIZATION_PLAN.md` for full tables, rationale, and
decision record. Summary:
- **Tier 1 (climate & size):** per-model centre & width, each from that model's
  *own* published calibration (Yasso07 ← Tuomi 2009 Table 3 [climate] + Tuomi
  2011 Table 4 [woody]; Yasso15/20 ← `Yasso15.dat`/`Yasso20.dat`). Homogeneous
  in *method + scale*, not in raw number.
- **Tier 2 (12 transfer fractions):** common weakly-informative logit prior,
  SD **0.4**, identical across models; per-model published centres. Wired by
  **explicit per-fraction listing** in each `*_SIGMA_PPM` (decision #5).
  Guardrail: fractions must stay **non-identified** (prior-pinned + ridge-
  correlated + low KL) — *not* "R-hat must stay poor" (that wrongly conflates
  non-identifiability with poor mixing; see Convergence expectations above).
- **Tier 3 (`sigma_init`, `sigma_input`):** identical log SD 0.50 everywhere.

Concept: priors are homogeneous in *construction method, scale, constraints,
and degrees of freedom* — **not** identical numbers (identical numbers is what
broke Yasso07). Any residual divergence is then attributable to structure/data,
not prior asymmetry.

**Decisions (locked 2026-06-05; #5 settled 2026-06-06):** published "±" limits
read as **1σ** (conservative); **GUI centres** kept; woody widths from Tuomi
2011 Table 4; fraction logit SD **0.4**; fraction-width mechanism = **explicit
per-fraction listing**.

**Bug fixed:** Yasso07 (and SP1/TP2/TP3) carried a hand-set `beta2` prior SD of
**0.05**, ~100–360× looser than the empirical Yasso15/20 widths. Since `beta2`
multiplies T² (~600 in boreal climate), this detonated the climate modifier
`xi`, giving 25–46 % `-Inf` proposal rates and R-hat up to 38 in run `Yasso07
20260603_053323`. Applied fix (1σ reading): `beta2` 0.05→**0.00065**, `beta1`
0.20→**0.26**, `gamma` 0.30→**0.20**, `delta1` 0.15→**0.16**, `delta2`
0.10→**0.12**, `r` 0.015→**0.042**.

**Local de-risk complete (2026-06-06):** a faithful prior-pushforward pre-flight
(`Calibration_real_data_transient/preflight_prior_pushforward.R`) sources each
model's real calibration script up to the MCMC launch and pushes prior draws
through the genuine `ll_fn` at full N — no Puhti, no fairshare. Result: the
`beta2` detonation is eliminated in all six models. A `beta2` sweep through the
real `ll_fn` shows the failure is *prior-driven, not engine-driven* (positive
`beta2` excursions give ll ≤ −10⁹, → −Inf only past +0.12; the matrix-exp cap
catches only the extreme tail). OLD `beta2~N(−0.0016,0.05)` puts ~40 % of draws
in fatal territory (reproduces the 25–46 % failure); NEW `~N(−0.0016,0.00065)`
≈0 %. NEW-prior forward-blowup rates: SP1/TP2/Yasso15/Yasso20 **0 %**, Yasso07
~2–3 % (`delta2`/`sigma` baseline, not `beta2`), TP3 ~10–12 % — and the TP3
residual is **genuine under-constraint** (`alpha` log SD 0.50, ~5 %; climate
~3 %; fractions & `sigma` **0 %**), a finding left untouched, confirmed to run
acceptably on Puhti. The methods-doc parameter-class × prior-criteria table is
done (`…/documentation/HIKET_calibration.Rmd`, §"Priors").

**Production run done (2026-06-11; Yasso07/15/20 inspected 2026-06-12):** all
params R-hat < 1.05, ESS 1265–2369, inf_rate ~0%, forward sanity PASS, distinct
non-degenerate ll-at-defaults (not the stale-`.so` trap). Fractions converged in
R-hat **as expected** but remain non-identified by every other measure (prior-
pinned, ridge-correlated, KL < 1 nat). The earlier expectation "fraction R-hat
stays poor" was wrong — see revised Convergence expectations. Minor: report
`final_ll` spread `[WARN]` (Yasso20 chain 2) is a single-last-draw artifact, not a
displaced chain (R-hat/ESS/traces/marginals all agree the 5 chains share one
posterior).

**SP1/TP2/TP3 cross-checked (2026-06-08 batch, inspected 2026-06-12):** same
unified picture — where structure exceeds what 2 SOC obs/plot resolve, surplus DoF
go to prior-pinning + correlations + low KL + input-multiplier compensation, never
to bad R-hat. **TP2** textbook clean (max R-hat 1.006). **TP3** R-hat clean but the
**entire 3-pool kinetics (`alpha_A/S/H`, `p_S/p_H`) gain ~0 nats — data informs
none of them; `sigma_input` absorbs everything (KL ≈ 21)**. The model fits by
rescaling inputs, not learning rates (8–17% forward-blowup `-Inf`, documented,
tolerable). **SP1** is the *only* genuine R-hat signal: `alpha` R-hat 1.092, ESS
208 — but it's the most strongly *identified* param in the study (KL ≈ 19), just
heavy-right-tailed (97.5% ~6× median), so hard to sample (chain-5 12% `-Inf` from
tail draws). Identified-but-heavy-tailed ≠ non-identified. Interpretation flag for
the writeup: `sigma_input` does very different work across models (TP3 KL ≈ 21;
Yasso20 median ≈ 3.0, inputs ↑3×; SP1 median ≈ 0.31, inputs ↓3×) — a structural-
intercomparison result in its own right.

---

## Parameter conventions

- `sigma_init` — proportional SD on initial carbon state (log-transformed, prior centre 1.0)
- `sigma_input` — multiplier on litter inputs (conflates tree litter model uncertainty with understorey; not a pure correction)
- `MODEL_FREE_NAMES` — vector of free parameter names, used instead of hardcoded index ranges
- Individual `param_spec` entries per parameter (never grouped)
- Split transforms: `beta1`/`betaN1`/`betaH1`/`delta2`/`r` → `log`; flow fractions → `logit`; others unconstrained

---

## Data

| Dataset | Description |
|---|---|
| NFI/Biosoil/MUSTIKKA/Komeetta plots | ~350–400 calibration-ready Finnish plots |
| Litter inputs | Tupek et al., Zenodo DOI: 10.5281/zenodo.19736499 |
| Climate | `nfi_plot_weather_data_1961_2025.nc` (gridded daily) |
| SOC campaigns | VMI8 (1985–86), Biosoil (2006), Komeetta (pending) |

**Litter units:** `input_raw_monthly.csv` is already in **tC/ha/yr** — no
multiplier needed in pipeline scripts (fix applied upstream in `Data_work.R`).

**Excluded plots:** zero-litter plots, OFH-absent plots, MRT > 100 years.

---

## Output file conventions

| Location | Contents |
|---|---|
| `runs/<MODEL>_posterior_<YYYYMMDD_HHMMSS>.rds` | MCMC posterior |
| `runs/<MODEL>_posterior_predictive_<RUN_ID>.rds` | Predictive bundle |
| `diagnostics/<MODEL>/` | Per-model PNGs, CSVs, residual RDS |
| `diagnostics/multimodel/` | Cross-model figures and metrics |
| `Data/model_inputs/` | Input bundles (gitignored) |

- PNG only (no PDFs)
- `set.seed(2025)` standardised across all scripts
- `RUN_ID` = `format(Sys.time(), "%Y%m%d_%H%M%S")`

---

## Known outstanding items

- **⚠ TP3 exact-integrator re-calibration RUNNING on Puhti (2026-06-30, job
  35323262, commit ff216ce), awaiting results — THIS is the real one.**
  History: the "2026-06-17" launch actually ran **Euler** (the `.tp3_step`
  rewrite was uncommitted on the Mac; Puhti was at 92da86f) → the
  20260617_085813 posterior is Euler-based, discard it. The local exact code
  also had a named-scalar bug in `.tp3_step` (NA `C_init`); fixed + committed +
  pushed as **ff216ce**, pulled on Puhti, re-launched 2026-06-30. **When it
  finishes:** (1) `rsync` back `runs/`, `diagnostics/`, AND `Data/model_inputs/`
  (the predictive stage hard-loads the input bundle keyed to the new TP3 RUN_ID);
  (2) re-run stages 2–4 locally (`run_hiket_pipeline.R --skip-calibration`) to
  refresh TP3 predictive + multimodel figures; (3) refresh the **NextGenC report**
  (`Reporting/NextgenC_report/`: `build_soc_matrices.R` → update the TP3 RUN_ID,
  then `plot_soc_trajectories.R`) — TP3 panel should go smooth; (4) re-render
  `HIKET_calibration.Rmd` and drop the "figures predate the fix" caveat in
  §sec:tp3osc. Expected: ringing gone, fit metrics ~unchanged. This run also
  produces the **restyled marginal plots** (filled polygons, class colours) — the
  `plot_one_marginal_honest` change in `calibration_engine.R` is already in code
  (previewed via `doublechecks/preview_marginal_style.R`); the §15 appendix
  auto-picks them up via `latest()`.
  **Then: close out HIKET on Puhti before the Roihu migration (Puhti decommissions
  ~end July 2026). Get this TP3 run + downstream + docs fully complete on Puhti
  first; further development continues on Roihu (CSC's Puhti successor) — see
  Environments.**
- **Prior homogenization applied + locally de-risked (2026-06-06)** — all six
  `*_priors.R` edited; methods-doc table done; the `preflight_prior_pushforward.R`
  pre-flight confirms the `beta2` detonation is gone in all six models at full N
  (see Prior-homogenization design note above). **One follow-up remains:** the
  production run on Puhti — **DONE & all six inspected 2026-06-12**. Yasso07/15/20
  (06-11) clean: all R-hat < 1.05; fractions non-identified by prior-pinning +
  ridge correlations + KL, not by poor R-hat. SP1/TP2/TP3 (06-08) cross-checked:
  TP2 clean; TP3 kinetics fully prior-dominated (`sigma_input` absorbs all); SP1
  `alpha` R-hat 1.09 = identified-but-heavy-tailed, benign. See Prior-
  homogenization note for the unified picture.
- **Downstream pipeline run locally — DONE 2026-06-12.** Stages 2–4
  (`run_hiket_pipeline.R --skip-calibration`) ran end-to-end on the Mac for all
  six models, zero errors, ~12 min total. Auto-selected RUN_IDs SP1/TP2/TP3
  `20260608_*`, Yasso07/15/20 `20260611_*`. **Local is the right venue** for this
  stage (predictive = `N_PP_DRAWS = 100` draws/model, trivial vs calibration; SLURM
  queue wait dwarfs the compute). Outputs in `diagnostics/multimodel/` + per-model.
  **Headline predictive metrics (447 plots): all models R² 0.05–0.11 (calib) /
  0.02–0.05 (holdout), RMSE ~50–62, slight + bias. TP3 best, Yasso20 worst, spread
  small — on-thesis (structure gives modest skill differences against generally low
  SOC predictability).**
- **⚠ Predictive-coverage caveat (PENDING decision, paused 2026-06-12).** The "95%
  cov" column reads ≈0.05–0.09, *not* a bug: the predictive interval
  (`quantile(total_soc, .025/.975)` in each `run_*_predictive.R`) is built from the
  100-draw spread of the model **mean** and **omits the multiplicative-normal
  observation error** (`sigma_obs_fixed` ≈ 0.4723, loaded but never injected). So it
  is a *parameter CI on the mean*, mislabeled as a posterior-predictive interval. R²/
  RMSE/bias and the rankings are unaffected. **Recommended fix (option 1, not yet
  applied):** in all six predictive scripts add a true posterior-predictive coverage
  alongside the current one — push each draw's `total_soc` through
  `* exp(rnorm(0, sigma_obs))` before quantiles (should land ≈0.95 if the error model
  is adequate); keep both numbers. Minimum alternative: relabel column to "param-CI
  coverage". ~5-min change + rerun.
- Hard dependency in `calibration_engine_transient.R` on original non-transient
  engine file existing in project
- ✅ DONE 2026-06-12: `rownames(inputs_proj/clim_proj) <- NULL` now applied in all
  six predictive scripts (was missing in SP1/TP2/Yasso07/Yasso15); heatmap
  `nrow/ncol(imp_mat) < 2` guard added in `run_multimodel_comparison.R`. (CI `mtext`
  overflow not hit in the 2026-06-12 run; revisit only if it recurs.)
- **Predictive stage needs the input bundle keyed to the posterior RUN_ID** —
  `run_*_predictive.R` hard-loads `Data/model_inputs/<MODEL>_inputs_<RUN_ID>.rds`
  (no fallback). After syncing `runs/`+`diagnostics/` from Puhti, ALSO sync
  `Data/model_inputs/` or Stage 2 fails with "cannot open … _inputs_<RUN_ID>.rds".
- Holdout validation infrastructure (`prepare_holdout_split.R`,
  `run_holdout_validation.R`) not yet implemented
- Six ambiguous BIOSOIL→MUSTIKKA plot mappings deferred
- Ground vegetation data for residual analysis requires LUKE MUSTIKKA database
  request (absent locally)

---

## Code style

- Concise, well-commented R with consistent section headers (`# === ... ===`)
- Short comments, not verbose prose
- Explicit rationale for design decisions in headers
- Logical placement descriptions preferred over line number references
- Documentation (`HIKET_calibration.Rmd`): rationale over implementation;
  "auxiliary uncertainty parameters" not "nuisance parameters";
  Puhti infrastructure details excluded
- **FUNDAMENTAL doc element:** the general calibration documentation
  (`HIKET_calibration.Rmd`) **must** contain a *parameter-class × prior-criteria*
  table — for each class of parameter (climate, size, transfer fractions,
  auxiliary), the criteria used to set its prior (centre source, width source,
  transform/constraint, and whether it was consistently calibrated across the
  original models). This table is the backbone of the homogenization argument;
  the worked version lives in `Prior_specs/PRIOR_HOMOGENIZATION_PLAN.md` §2 and
  must be mirrored (rationale form) into the Rmd.

---

## Key collaborators

| Person | Role |
|---|---|
| Toni Viskari | Yasso20 model author; structural clarifications |
| Boris Tupek | Litter dataset author |
| Jani Anttila | LUKE, GHG portal |

**GitHub:** `ilmenichetti/HIKET`  
**Prior sources:** FMI Ryassofortran repository (`https://github.com/YASSOmodel/Ryassofortran/tree/master/data`); `Priors_model_matching.R` derives prior specs.
