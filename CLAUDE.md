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
├── Data_work.R                               # upstream data preparation
├── Additional_layers.R
└── HIKET_calibration.Rmd                     # methods documentation
```

**Gitignored:** `Data/`, `runs/`, `diagnostics/**/*.rds`, `diagnostics/**/*.csv`,
Fortran `.so` binaries, `*.rds`/`*.csv` in diagnostics.

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
```bash
cd Model_functions_real_data/Decomposition_functions/Yasso/
rm -f yasso07.so yasso07.o yasso07_mod.mod yasso15.so yasso15.o yasso15_mod.mod
apptainer_wrapper exec R CMD SHLIB yasso07.f90 yasso15.f90
```

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

### Likelihood / error model
Multiplicative-normal (proportional) error model. Asymmetric penalisation of
underprediction is a known consequence — noted in methods.

### Convergence expectations
With ~19 free parameters and only 2 SOC observations per plot, poor R-hat for
flow fractions is **structural non-identifiability**, not a bug.

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
  Guardrail: fraction R-hat must stay poor — we are not forcing convergence.
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

**Still pending:** one production run on Puhti (all 6, full N, full iter, in
parallel) — the only thing local testing cannot show is the MCMC convergence
behaviour (climate R-hat should *improve* while fraction R-hat *stays poor*).

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

- **Prior homogenization applied + locally de-risked (2026-06-06)** — all six
  `*_priors.R` edited; methods-doc table done; the `preflight_prior_pushforward.R`
  pre-flight confirms the `beta2` detonation is gone in all six models at full N
  (see Prior-homogenization design note above). **One follow-up remains:** the
  production run on Puhti (all 6, full N, full iter, parallel) to confirm MCMC
  convergence behaviour (climate R-hat improves; fraction R-hat stays poor).
- Hard dependency in `calibration_engine_transient.R` on original non-transient
  engine file existing in project
- `rownames(inputs_proj) <- NULL` / `rownames(clim_proj) <- NULL` not yet
  applied in SP1/TP2/Yasso07/Yasso15 predictive scripts
- CI max `mtext` overflow and heatmap `ncol(imp_mat) < 2` guard in multimodel
  comparison script
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
