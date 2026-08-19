# HIKET Project — Claude Code Context

## Project overview

**HIKET** is a Bayesian calibration and structural intercomparison of five soil
organic carbon (SOC) decomposition models applied to Finnish forest ecosystems.
The scientific question: do divergent model projections of Finnish forest carbon
sink saturation reflect genuine structural differences, or calibration artefacts?

**Models in scope:** SP1, TP2, TP3, Yasso07, Yasso15, Yasso20  
**Application:** Finnish National Forest Inventory (512 calibration-ready permanent plots)  
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

**Environment switches** (all inert unless set; on Roihu they need the `SINGULARITYENV_` prefix
to reach R inside the r-env container):

| variable | effect | default |
|---|---|---|
| `HIKET_LOGNORMAL_LIK` | `0` reverts to multiplicative normal | log-normal |
| `HIKET_SIGMA_TOTAL` | overrides the likelihood error scale (total **SD**) | `sigma_obs_fixed` (0.442) |
| `HIKET_LIK_DF` | Student-t degrees of freedom on the log scale; unset/`Inf` = Gaussian | Gaussian |
| `HIKET_PRIOR_TIGHTEN` | multiplies Tier-2 **fraction** SDs only | 1 (0.4 unchanged) |
| `HIKET_N_CHAINS` / `_N_ITER` / `_N_BURNIN` | short test runs | 5 / 50000 / 5000 |

---

## Environments

### Mac (development)
- Architecture: arm64, R 4.3
- Run scripts interactively or via `Rscript`
- Fortran `.so` files compiled locally for arm64 — **do not push to git**

### CSC Roihu HPC (production, since 2026-07-02)
- Path: `/scratch/project_2019134/HIKET/` (same project number as Puhti; **180-day
  scratch cleanup** — working area only, keep code in git + final outputs on Zenodo)
- Login: `ssh roihu` — **cert-based auth, re-sign the SSH certificate daily** via
  MyCSC → Profile → SSH PUBLIC KEYS → ⋮ → "Sign and download SSH certificate" (24 h
  validity). It downloads as `cert.pub`; save to `~/.ssh/cert.pub` (the `roihu` alias
  in `~/.ssh/config` points `CertificateFile` there). Host key ED25519 SHA256
  `YNdesHbXhxN0hKD4mWvYGQONebjRqY+CGXDqPiZyByQ`. **The user handles all SSH cert/key
  operations and Roihu logins MANUALLY — do not touch `~/.ssh` files or `ssh roihu`
  yourself unless explicitly asked (see memory `ssh-cert-handled-manually`); give the
  user the commands to run instead.**
- Module: `module load r-env`
- **No explicit `apptainer_wrapper` call needed** — you invoke `srun Rscript
  --no-save …` directly (batch) or `start-r`/`R --no-save` (interactive). BUT the
  `r-env` module's `Rscript` shim still `exec`s a **singularity container** under the
  hood (confirmed in job 147311/147312 `.err`: `/usr/bin/singularity … exec …`). This
  matters for parallelism: the container does NOT see `SLURM_CPUS_PER_TASK`, so
  `parallelly::availableCores()` returns the **full node (383)**, ignoring
  `--cpus-per-task` → 383 mclapply forks → **OOM** (jobs 147307-312, 2026-07-06). Fix
  applied in all six `run_*_transient_calibration.R` (cap cores by the SLURM alloc) +
  `hiket_*.sh` (`export SINGULARITYENV_SLURM_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK`
  passes the alloc into the container). Watch the `Cores per chain:` log line — it must
  read **40**, not 383.
- Fortran compilation: `R CMD SHLIB <file>.f90` (native, no wrapper; AMD Zen 5 x86-64)
- Partitions: `small` (72 h, 384 cores/node), `test` (15 min), `medium`/`large`
  (36 h), `longrun` (10 d). Existing `--partition=small --time=36:00:00` is valid.
- SLURM user: `menichet`

### CSC Puhti HPC (legacy — DECOMMISSIONS ~end July 2026)
- HIKET close-out on Puhti is **DONE** (2026-07-01). Do NOT start new work here;
  all further development is on Roihu. Retained only as a source to `rsync` any
  remaining `runs/`/`diagnostics/`/`Data/` from before it goes offline.
- `ssh menichet@puhti.csc.fi` (no certificate needed, unlike Roihu).

### Sync workflow
```bash
# Code → Roihu
git add ... && git commit -m "..." && git push
# [on Roihu]
cd /scratch/project_2019134/HIKET/ && git pull
# Recompile .so if any .f90 changed (see Fortran section below)

# Data (gitignored) → Roihu: from Mac (or one-time rsync from Puhti while it lives)
# NOTE: use the `roihu:` SSH alias, NOT `menichet@roihu-cpu.csc.fi` — only the alias
# carries `CertificateFile ~/.ssh/cert.pub`; the raw hostname fails `Permission denied
# (publickey)` because it never presents the daily-signed cert.
rsync -av "<mac-repo>/Data/" roihu:/scratch/project_2019134/HIKET/Data/

# Results → Mac (also sync Data/model_inputs/ — the predictive stage hard-loads the
# input bundle keyed to each posterior's RUN_ID)
rsync -av roihu:/scratch/project_2019134/HIKET/Calibration_real_data_transient/runs/ ./Calibration_real_data_transient/runs/
rsync -av roihu:/scratch/project_2019134/HIKET/Calibration_real_data_transient/diagnostics/ ./Calibration_real_data_transient/diagnostics/
rsync -av roihu:/scratch/project_2019134/HIKET/Data/model_inputs/ ./Data/model_inputs/
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
Rscript --no-save \
  Calibration_real_data_transient/run_hiket_pipeline.R --skip-calibration
```

### Fortran recompilation (after any `.f90` change or fresh clone)
**Compile each `.f90` in its OWN `SHLIB` call.** A combined
`R CMD SHLIB yasso07.f90 yasso15.f90` links both objects into a single
`yasso07.so` and never creates `yasso15.so` — which breaks Yasso15 AND Yasso20
(Yasso20 loads `yasso15.so`; it shares Yasso15's Fortran), and risks symbol
shadowing. The `.so` files are gitignored, so this MUST be redone on Roihu after
a fresh clone or any `.f90` change — a stale `.so` silently returns garbage (e.g.
zero `C_init` from the transient init → see the Yasso07 non-convergence post-
mortem, 2026-06).
```bash
module load r-env
cd Model_functions_real_data_transient/Decomposition_functions/Yasso/
rm -f yasso07.so yasso07.o yasso07_mod.mod yasso15.so yasso15.o yasso15_mod.mod
R CMD SHLIB yasso07.f90    # native on Roihu (no apptainer_wrapper); separate call
R CMD SHLIB yasso15.f90    # separate call
```
Quick check a `.so` is current: `yasso07_transient_init` at
`YASSO07_DEFAULT_PARAMS` should give ~69 tC/ha; ~0 means the binary is stale.

---

## Key design decisions

### Intercomparison design principle
All models must have **equal degrees of freedom** in flow fractions. Any
structural asymmetry (fixed vs. free parameters) confounds climate computation
differences with calibration asymmetry and makes results uninterpretable.

### Round-1 kinetic homogenization — ICBM anchor (finalized + validated 2026-07-16)
Refines the principle above to **equal external *information***, not just equal DoF,
by parameterizing the simple models (SP1/TP2/TP3) the way Yasso already is. **Yasso
does not fit its decomposition rates** — the `a`-vector is a fixed constant
(`run_Yasso07_*`: `FIXED_RATE_NAMES <- c("alpha_A","alpha_W","alpha_E","alpha_N",
"p_H","alpha_H")`), injected in `assemble_model_params`; what Yasso *calibrates* is the
12 lateral transfer fractions, climate (β), woody size, and the two σ's. So the simple
models split kinetics into two classes:
- **Intrinsic rates → externally anchored, NOT fit to SOC data** (ICBM, Andrén &
  Kätterer 1997: `k1=0.8`, `k2=0.00605`, `h=0.13`; Ultuna `r=1`). *Fast* rate
  (`alpha_A`/SP1 fast) **fixed** (litterbag-grounded, transferable); *slow* rate
  (`alpha_H`, TP3 `alpha_S`, SP1's slow-dominated single rate) gets a **very-informative
  prior** — free but tightly pinned at ICBM, so the data can nudge the one weak link
  (`k2`, arable→forest) and the posterior width is a transferability diagnostic.
- **Humification fractions → FREE** (analog of Yasso's free lateral fractions): `p_H`
  (TP2), `p_S`+`p_H` (TP3) keep the Tier-2 logit prior (SD 0.4), **re-centred at
  `h≈0.13`**. Fixing them would make the simple models *more* constrained than Yasso.

Free set is then **identical across SP1/TP2/TP3** (`β1,β2,γ,σ_input,σ_init`); complexity
adds only free *partitioning* splits (0→1→2→12). Rates are fixed/pinned **constants** set
**once** at the β centre with a one-time `/xi_Ultuna` (~6%) offset — **no per-draw
derivation, no TP3 reparam** (both superseded). **TP3 critical detail:** `alpha_S` must be
**k2-scale** (`≈k2/(1−p_H)≈0.0070`), NOT intermediate — an intermediate value starves the
cascade (two 0.13 humification steps → MRT ~6 yr, σ_input 4–5×). **Validated** by
`doublechecks/icbm_anchor_sanity.R` (447 plots, σ_input=1): SP1/TP2/TP3 all reach observed
SOC (median 71) at bulk MRT ~24 yr and **σ_input 1.05–1.36** — physical, inside the
understorey window, closing the σ_input escape hatch via kinetics. Implements C1 of
`NEXT_SESSION.md` §5 (REVISION_PLAN.md removed 2026-08-07); edits land in `Prior_specs/{SP1,TP2,TP3}_priors.R` +
each `assemble_model_params`. Memory: [[icbm-anchor-validated]], [[revision-round1-plan]].

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
**Consequence: TP3 must be re-calibrated** for the doc's TP3 posterior/figures
(RUN 20260608) to reflect the exact integrator; stocks/skill change negligibly.
The only remaining non-exact path is the Yasso diagonal-dominance **Euler
fallback** (a guard that the stick-breaking budget keeps from firing).

**IMPORTANT correction (2026-07-01, verified on the real exact run
`TP3_posterior_20260630_090644`):** the exact integrator does **NOT** make TP3
smooth — the earlier "the ringing disappears" expectation was WRONG. Two effects
were conflated. (1) Euler added a huge *numerical* overshoot (posterior-median
trajectory roughness up to ~10³–10⁴ tC/ha on some plots); exact removes it, cutting
roughness ~5–700×. (2) A **bounded, genuinely physical** interannual oscillation
*survives* exact integration: at the posterior median `alpha_A*xi` has median
**~1.4–2.0** (not ~1.1) and exceeds 1 in **~95–100%** of years, so the fast active
pool re-equilibrates to each year's climate (`A_ss=J/(alpha_A*xi)`) → total SOC
genuinely tracks interannual climate. **Sub-annual integration does NOT fix it**:
re-integrating the same posterior with 12 monthly sub-steps (seasonal T from
`temp_mean`+amplitude) cuts roughness only ~10–20% and shifts mean SOC <2 tC/ha —
the seasonal cycle is already inside `xi`, the residual is year-to-year. **The
TP2-vs-TP3 contrast is STRUCTURAL, not the integrator** (both exact): governed by
active-pool turnover × that pool's SOC share. Posterior-median clincher (30 plots):
SP1 `alpha*xi`≈0.003, share 100%, roughness 0.08 (single slow pool, still rising);
TP2 `alpha_A*xi`≈1.61, active share ~23%, roughness ~3.7 (fast pool but mostly
buffered in humus); TP3 `alpha_A*xi`≈1.86, active share ~26%, roughness ~9.3 (fast
pool + 3-pool cascade propagates the swing + higher stocks). Yasso parks carbon in
slow humus → never sub-annual → smooth. **The oscillation is a correct symptom of
TP3's fast-pool/inflated-input (`sigma_input`≈13–20) under-identification, not a
defect.** Docs reframed: `HIKET_calibration.Rmd` §sec:tp3osc + the two forward refs.
Reproducible checks in `doublechecks/`: `tp3_exact_vs_euler_realforcing.R`,
`tp3_subannual_forcing_test.R`, `xmodel_turnover_share.R`.

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

### 🚩 THE POSTERIORS ARE OVERCONFIDENT (established 2026-08-11/12)

Measured on run `20260810_1529*`. **The operational test is met:** defensible analysis choices move
the estimate further than the posterior's own uncertainty. Yasso15's `sigma_input` moved 0.175 on
the log scale (σ 0.442→0.72 alone; its climate priors were untouched) against a posterior SD of
**0.087** — two standard deviations from one choice.

**The mechanism.** The likelihood treats 1269 plot-years as independent. They are not:

| shared component | measured | effect |
|---|---|---|
| plot-persistent | ICC **0.57–0.71**; r(2006,2024) = 0.66–0.81 | effective n 1269 → **~630** |
| spatial (latitude bands) | band means vary **5–6×** more than iid allows; regional sd 0.12–0.14 | SE of the national mean **2.3×** larger |
| campaign | means +0.131/−0.045/−0.031; spreads **0.97 / 0.62 / 0.54** | reweights the trend |
| heavy tails | kurtosis **6.9–7.1** (Gaussian = 3) | a few plots steer the fit |

**Why this specifically breaks MRT.** All of it attacks the *aggregate level*, and the level is what
pins `MRT × σ_input`. Yet across runs the product is constant to 4% (Yasso15 43.0 → 44.7) while the
split slides. Reproduce with `doublechecks/ridge_test.R`.

⚠ **CORRECTED 2026-08-13 — the ridge SPLITS BY FAMILY; the old blanket claim was mis-scoped.**
`ridge_test.R` §1 computes MRT in closed form and covers **SP1/TP2/TP3 only** (its header says so:
`intrinsic_mrt.R`'s per-draw values were not paired to `sigma_input` draws). Its −0.30/−0.37 was
therefore written up as "the posterior does not explore that ridge" in a paragraph whose evidence is
*Yasso15* — a model it never measured. `build_F14_mrt_ridge.R` closes that gap by pairing per-draw
intrinsic MRT with per-draw `sigma_input`, and the answer is **the opposite for Yasso**:

| | corr(log MRT, log σ_input) | sd(log product)/sd(log MRT) |
|---|---|---|
| TP2 / TP3 | −0.263 / −0.313 | 0.98 / 0.97 → pinned independently |
| Yasso07 / 15 / 20 | **−0.787 / −0.725 / −0.572** | **0.62 / 0.73 / 0.91** → a real ridge |

⚠ **Updated on run `20260817_12*` (2026-08-18):** Yasso07 **−0.733**, Yasso15 **−0.675**, Yasso20
**−0.637** (ratios 0.68 / 0.77 / 0.81). The family is now much more HOMOGENEOUS — Yasso20 no longer
stands apart as the intermediate case, which removes the tension noted below. S13 (the benchmark
trio, ξ included) reads SP1 **−0.902**, TP2 **−0.812**, TP3 **−0.805**: the simple models still ride
a STRONGER ridge than any Yasso, so "the simple models show no such ridge" remains falsified.

Not a run effect: TP2/TP3 give −0.34/−0.37 on `20260810_1529*` and −0.26/−0.31 on `20260812_0809*`.
So the overconfidence argument stands for the simple models but must NOT be stated for the Yasso
family, where the two genuinely trade off. Yasso20 is intermediate (0.92) — and is also the model
whose published MRT costs almost nothing in fit (below); those two facts are not yet reconciled.

**Why scale alone cannot fix it.** Uniform σ inflation broadens every direction by √k, so reaching
the 2.5× Yasso15 needs would require σ ≈ 1.8 against a measured residual spread of 0.79 — it would
invalidate the self-consistency argument that justifies σ in the first place. **Correlated-error
structure broadens the level direction at fixed total variance; nothing else does.**

**Options, in increasing cost:** (1) scale 0.72→0.80 + Student-t tails — **the run launched
2026-08-12**, worth ~1.3–1.4×, *not* a fix, but may move locations; (2) campaign-specific σ — cheap,
targets the **trend**, and note it *sharpens* the posterior ~10% (1985's weight 32%→14%, 2024's
32%→46%); (3) fit regional/campaign aggregates instead of plot-years — a data decision, not a
likelihood one, and defensible because there is no plot-level signal to lose (R² 0.004–0.019);
(4) marginalised compound-symmetry likelihood — O(n) rank-1 updates, no new parameters.

⚠ **Campaign-specific σ is safe; a campaign-specific BIAS is not** — the trend *is* the difference
between campaign levels, so a free `c_j` eats the signal (already demonstrated: the δ-offset test
gave δ = −0.137, requiring σ_init ≈ 0.98, implying no accumulation ever happened).

Memory: [[likelihood-overconfident-ridge]].

### 🚩 "MRT TOO SHORT" IS THREE DIFFERENT FINDINGS, NOT ONE (established 2026-08-13)

Run `20260812_0809*`. Two results, both reproducible.

**1. The cost of the published MRT is wildly uneven.** Best attainable fit near the published value vs
best anywhere (max-over-draws lower bound on the profile, from `F14_mrt_ridge.rds`):

| | P(MRT > published) | cost in ll units | draws in the band |
|---|---|---|---|
| Yasso07 | **0%** (max MRT 25.4 vs 33.5) | **NEVER REACHED — no estimate** | 0 |
| Yasso15 | 0.011% (24) | **≤ 14.7** | 741 |
| Yasso20 | **18.1%** (40 689) | **≤ 1.6** | 121 141 |

⚠ These are **over**estimates of the cost, worst where draws are thinnest. Only a profile likelihood
settles it — and for Yasso07 nothing at all can be read off the samples.

⚠⚠ **CORRECTED 2026-08-13 — the first version of this table was BUILT ON BURN-IN ARTEFACTS.**
`getSample()` on the saved chains returns the **first retained iteration of each internal DEzs chain**
(3 per sampler × 5 samplers = **15 rows per model**), sitting 100–200 ll below the bulk with a clean
gap. Those 15 rows were the ONLY draws near the published MRT for Yasso07/15, so the original
≤137.7 / ≤14.3 described initialisation, not the posterior. **Fix: extract per sampler with
`start = 2`** (also drops the 1-in-3 thinning `getSample` applies to the list ⇒ **225 015** draws, not
75 015). Yasso07's apparent reach 34.1 → **25.4 yr**; ll range 209 → 33.

⚠⚠ **TWO DIFFERENT "PUBLISHED MRT" NUMBERS EXIST — DO NOT MIX THEM.** `intrinsic_mrt.R` computes
both: the published **POINT** (`to_original(best_x)`) and the published **POSTERIOR** median (MRT over
the FMI `.dat` sample, Yasso15/20 only). MRT is a nonlinear many-to-one map, so `median(MRT) ≠
MRT(median)`. For Yasso20 the point is **19.03** and the posterior median **25.02** — and the cost of
reaching it goes **1.4 → 9.8** accordingly (29 draws in the band vs 40 369). The table above and F14's
dashed line use the **POINT** (33.47 / 30.38 / 19.03). CLAUDE.md's older line "published 33.39 / 30.27
/ 25.02" and memories `mrt-too-short` / `run-563524-error-model` mix the two bases. **Resolve which
comparator the paper uses before quoting any cost.**

**2. Yasso07's gap is ENTIRELY the climate modifier, and that is STRUCTURAL — not a prior artefact.**
Yasso07 applies **one** ξ to every pool ⇒ ξ is a pure rescaling of time and `MRT = MRT_ref/ξ` exactly.
Calibration moves ξ 0.855 → **1.822** (×2.13) at the Finnish mean climate, which alone predicts
**15.7 yr against an actual 15.2** — the whole gap, nothing left for the fractions. Yasso15/20 carry
**three** pool-specific modifiers (ξ_AWE, ξ_N, ξ_H); MRT is set by the slow humus pool, which has its
own `betaH1`, and that moved only ×1.04 / ×1.06. In Yasso20 the components oppose each other
(ξ_AWE ×0.76 vs ξ_N ×1.16) and largely cancel. **The leverage Yasso07 gives climate is structurally
unavailable to its successors** — so the family's MRT gaps must NOT be discussed as one phenomenon.

Corollaries: Yasso07's *published* ξ is **0.855 < 1** (its parameterisation says Finland decomposes
slower than its reference; Yasso15/20 already say faster, 1.16–1.66), so it starts furthest away and
is the only one able to travel the distance in one parameter. And **calibration inverts the family
ordering**: published 33.5 / 30.4 / 19.0 (07 slowest) → ours 15.2 / 22.1 / 17.5 (07 **fastest**).

⚠ OPEN: is ξ = 1.82 defensible? It needs `beta1` 0.0987 → **0.1578** (+60%), ~4.4 prior σ off centre
under the corrected Tuomi width. Suggestive but not proof: that value is essentially Yasso20's
*published* `beta1` (0.1580) — though Yasso07's β1 scales all pools while Yasso20's scales AWE only,
so they are not strictly the same quantity. Pair with the FMI warming-rate check.

Scripts: `doublechecks/xi_published_vs_ours.R`, `manuscript/figures/build_F14_mrt_ridge.R`.
Memory: [[mrt-climate-leverage-structural]].

### Likelihood / error model
**LOG-NORMAL is the DEFAULT since 2026-08-07** (`HIKET_LOGNORMAL_LIK=0` reverts to the old
multiplicative normal). The multiplicative normal let a prediction widen its own error bar, so the
penalty for a badly-missed plot plateaued and the fit bought tolerance by inflating everything.

**⚠ The error SCALE was wrong until 2026-08-10.** `sigma_obs_fixed` (0.442) is the *measurement*
CV from the SOC homogenisation, but it was used as the *total* error. Measured log-residual spread
is **0.708–0.735** across all six models (ratio 1.60–1.66), implying model error **0.55–0.59** —
larger than the observation error. Three consequences: posteriors overconfident; log-likelihood
differences inflated ~2.6×; and, because the level penalty goes as 1/σ², stock-LEVEL pressure
amplified ~2.6× relative to the priors, a candidate driver of the short bulk MRT.
`HIKET_SIGMA_TOTAL=<s>` overrides the scale (run 563524 uses **0.72**, the plug-in MLE of the
total error; self-consistency is testable by recomputing the residual spread afterwards).
The principled upgrade — free `sigma_model` with `σ_total² = σ_obs² + σ_model²` — is deferred:
~25 edits across 13 files, and a free per-model σ would reintroduce variance inflation across
models, so it needs pairing with a common fixed σ for cross-model comparison.

⚠ **Log-likelihoods are not comparable across different σ** (the normalising constant changes).

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
- **Tier 3 (`sigma_init`, `sigma_input`):** identical log SD everywhere — **0.25 since 2026-08-17**
  (was 0.50). Centres: `sigma_init` **0.90** (Korhonen growing stock ^ fitted elasticity 0.43),
  `sigma_input` **1.08** (Lehtonen & Heikkinen total litter 2.70 ÷ Tupek tree-only J̄ 2.511 — the
  understorey gap this parameter exists to close; was 1.30, unsourced). The two share a width by
  design; keep it that way. Width derivation: elasticity CI 0.03 + volume-record uncertainty 0.07 +
  structural litter-model error 0.09 → 0.12 in quadrature, prudential ×2. ⚠ Deliberately looser than
  the evidence supports so the sampler is pulled, not walled — and **never re-tune it to land MRT on
  a published value.**

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
| NFI/Biosoil/MUSTIKKA/Komeetta plots | **512** calibration-ready Finnish plots |
| Litter inputs | Tupek et al., Zenodo DOI: 10.5281/zenodo.19736499 |
| Climate | `nfi_plot_weather_data_1961_2025.nc` (gridded daily) |
| SOC campaigns | VMI8 (1985–86), Biosoil (2006), Komeetta (2024) — all three wired |

**Litter units:** `input_raw_monthly.csv` is already in **tC/ha/yr** — no
multiplier needed in pipeline scripts (fix applied upstream in `Data_work.R`).

**⚠ 1985 litter is RECONSTRUCTED, not observed (fixed 2026-08-04).** The Tupek product's
**first year is unusable**: source median 0.060 tC/ha/yr in 1985 vs 1.398 in 1986, with an
*identical* record count (2805) in both — not missing data, the values collapse. Signature
of litter derived from between-inventory biomass increments (first year has no predecessor
to difference against; **confirm with B. Tupek** before this wording goes in the paper).
Critical because 1985 is **t0** — it contaminated *four* things at once: `J_t0_mean` (pre-run
ENDPOINT, 1/5 of its window), `J_full_mean` (1917 anchor), `J_total_mean` → **`J_bar`** (the
flux_pair units bridge, so the σ_input physical window shifted too), and the forward run's
own first year, which lands exactly on the VMI8 observation. Uncorrected, `J_t0_mean` was
**1.543 vs 1.863** → the spin-up ended ~21% too low, biasing `C_init` low in all six models.
**Fix:** per plot × AWEN component, fit 1986–1990 linear trend and **backcast one year**
(`Data_work.R` §2, before the monthly expansion). Chosen over a flat 1986–1990 mean because
litter *rises* through that window — a mean would put 1985 *above* 1986 and contradict the
growing-stock history used for the C3 pre-run shape. Backcast puts it just below
(1.754 vs 1.838). All three options (carry-back 1.890 / mean 1.918 / backcast 1.863) agree
within ~3% vs the 21% error corrected, so the choice doesn't matter — using the artefact did.
Dropping 1985 outright was rejected: same answer, but moves t0 to 1986 and forces all
**441 VMI8 observations** to be re-mapped.

### SOC calibration target — homogenized baseline (✅ wired 2026-08-04)

`soc_obs_tCha` comes from `Data/SOC_homogeneized/` (built by
`build_soc_homogenized.R` from H. Ilvesniemi's LUKE workbook), **not** from the
retired `Data/SOC/soilC1985_2006.csv` + separate Komeetta ingest. Why it changed:
the old assembly applied **no coarse-fragment (stoniness) correction** to the
mineral layers, over-counting mineral C by a near-constant **~1.6×** (median
stoniness 43%; 1/(1−0.43)=1.76), and the three campaigns ran through separate,
drifted processing paths (the 2024 ingest was coded against the *1985* layer
protocol). The baseline reproduces LUKE's official national stocks **exactly**
(org+0–40, weighted: 2006 = 59.1, 2024 = 61.0 Mg/ha, n=446).

- **Target = `soc_profile_Mgha`**: measured organic + measured mineral 0–40 cm +
  modelled deep tail, integrated only to `z_cap = min(100 cm, depth augering
  reached)` — a per-plot **variable** depth, not a fixed 1 m. The old fixed-1 m
  rule invented ~14 tC/ha below bedrock on 38 refusal-at-10/20/40 cm plots.
  Extrapolation validated against the **measured** 40–80 cm Biosoil layer
  (Krs 204, 501 plots, held out): median pred/obs **0.96**, bias −2.5 Mg/ha.
- **Campaign medians now 59.4 / 66.7 / 67.4 tC/ha** (was 63 / 102 / 105 — the
  62% 1985→2006 jump was a cross-campaign artefact, and the old 2006 value
  contradicted the official same-year figure by 1.7×).
- Observation CV (fixed likelihood σ_obs) drops **0.472 → 0.440**.
- `Data_work.R` §1.0b (Komeetta append) and the in-script §1.3b depth fit are
  **removed**; §1.3b now merges the baseline. Backup of the pre-swap script:
  `Data/Data_work_pre_SOC_swap_20260804.R`. Full M&M rationale lives in the
  METHODS & MATERIALS block at the end of `Data_work.R` (canonical text) and in
  §"SOC stocks: a single homogenized basis" of `manuscript/HIKET_main_manuscript.tex`.

### ✅ The 1985 target comparability fix — APPLIED 2026-08-12 (commit `8ef5c62`)

Two defects, both established from the workbook's own structure (see `NEXT_SESSION.md` §0c for the
implementation, `manuscript/HIKET_discussion_memo.tex` for the findings):

1. **The 1985 organic layer is OFH only.** LM (litter+moss) is a separately coded layer (`Massat_2006`
   Krs 100 = LM vs 101 = OFH); it was measured in 2006/2024 and folded into their organic layer, and
   is **absent from 1985**. Worth **3.37 / 4.27 Mg C/ha** — the README's "minor, sub-Mg ha⁻¹" is wrong
   by 3–4×. It supplies **89%** of the apparent 1985→2006 organic gain and **47%** of the 2006→2024
   change. Proof: BiSo lays 1985/2006/2024 out in parallel blocks with identical columns; the 1985
   block is byte-identical to `Data_1985` (414/414); and the 2006 `C_kgha` equals `1e4*C_kgm2`
   exactly, i.e. litter-free. Corroborated by dry mass (1985 43.7 vs 2006 OFH 50.0 vs OFH+LM 56.7).
   **Fix = treatment C:** add each plot's own 2006 LM to 1985.
2. **"1985" is 1986–1995.** BiSo col 140 (`MAANAYTE`) and sheet `kiv_maat85_95` agree exactly
   (387 overlapping, 0 disagreements). 69/76/103/54/**72** plots for 1986/87/88/89/**1995** — so
   **19.5% are 10 years mis-dated**, the campaign mean represents ~**1989**, and the true mean
   interval to 2024 is **34.7 yr, not 39** (12% error in every rate denominator).

**Applied:** treatment C adds each plot's own 2006 LM to its 1985 organic layer (478 of 488 rows;
10 skipped where OFH is zero/absent, so `organic_zero` keeps firing). `samp_year` → `obs_year` →
`soc_obs_year` carries the true date and all six calibration scripts index on it; 82 undated 1985
plot-years dropped (71 South / 9 North). `SIGMA_1985_INFL` = **2.0**, pre-registered.
Switch: `HIKET_ADD_1985_LM=0` reverts treatment C.

⚠ **`year` was NOT renamed** — it doubles as the campaign key, so `sigma_infl` and
`HIKET_DROP_CAMPAIGN` would have broken *silently*. `obs_year` sits alongside it.

⚠ **The baseline now CONTAINS the imputed LM.** Any script treating the 1985 organic layer as OFH
must subtract `lm_added_1985` (kg/ha, 0 for 2006/2024) or it will count litter as humus — S8/S9/S10
were caught by exactly this.

**Converged target:** 1408 plot-years; weighted profile means 66.309 / 69.971 / 72.444; official
validation still exact; 456 calib-ready.

**What is NOT wrong with 1985**, despite intuition: its depth distribution (λ 0.038, *between*
2006's 0.031 and 2024's 0.041), stoniness (ρ=−0.08 n.s.), the organic/mineral boundary (r=−0.04
n.s. for this pair), and its noise level (subsoil repeatability 0.63 vs 0.65). Whatever is wrong is
a **level**, not scatter, and proportional bias appears in *every* campaign pair. The surviving
1985-specific mechanisms are the LOI-derived mineral C%, relocated subplots, and the physically
implausible **depth-inverted** subsoil gain (subsoil +19% vs topsoil +8.5% over 1985→2006 — a real
input-driven gain must be surface-weighted).

### 🔒 The data layer must not depend on the calibration layer (loop broken 2026-08-12, `44215b0`)

`build_soc_homogenized.R` used to read `site_raw.csv`, which is built from `plot_data` — whose row
set is `merge(avg_inputs, avg_SOC)` on `common_plots`, an **inner join on the SOC data the builder
itself produces**. The GTK soil-class extraction also ran only over `plot_data`'s plots, so
`soil_code` (→ per-class λ → deep tail, all three campaigns) was SOC-gated at source. The cycle did
**not** reach a fixed point in one pass (1411 vs 1408 plot-years, 38 differing λ, means ~0.05 Mg/ha
apart), which made the target non-reproducible from a clean checkout.

**Now:** `Data_work.R` extracts soil types over the SOC-independent universe (`plots_sf`, from
litter-input coordinates) and writes **`Data/model_inputs/site_attributes.csv`** (2719 plots) before
anything SOC-derived. `build_soc_homogenized.R` takes every target-affecting lookup from there;
`site_raw.csv` is retained only for descriptive covariates that cannot reach `soc_profile`.

> **INVARIANT.** Nothing SOC-dependent may be added to `site_attributes.csv`, and the **ROW SET**
> matters as much as the column list — gating those rows on SOC restores the loop invisibly.
> Verify with the fixed-point test: `Data_work → build → Data_work → build` must be byte-identical.

**Excluded plots:** zero-litter, OFH-absent (`organic_missing`/`organic_zero`),
peatland (Cajander KA 11–13), MRT > 100 years, **`const_litter`** (litter identical to
rel. SD < 1e-6 over 1986–2024 — a fixed repeated value, not a measured series; 8 plots,
detected by rule not hard-coded; borderline plots 39251/67631 deliberately kept), and
**`soc_outlier`** (whole-profile
stock > 250 Mg/ha in any campaign; 4 plots: 29232, 31751, 33631, 49571 — applied and
documented in `Data_work.R`, not silently upstream). The companion flag `high_change`
(|rate| > 3 tC/ha/yr, 23 plots) is **recorded but NOT excluded** — read as resampling
noise the error model should absorb, not as data error.

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

- **⭐⭐ NEXT ACTION — the CORRELATED LIKELIHOOD; design SETTLED 2026-08-18, not yet implemented.**
  `log y_ij = log f_ij(θ) + u^R_r(i) + u^P_i + u^C_j + e_ij`. **Nothing is estimated** — all four
  variances FIXED, offsets marginalised: **τ_R = 0.117** (latitude bands), **τ_P = 0.396**
  (2006–2024 pair covariance), **τ_C = 0.06 (1985) / 0.03 (2006, 2024)** prescribed, **σ_e = 0.685**
  as the remainder of an UNCHANGED total 0.800 (split, never added). Student-t **dropped**,
  `HIKET_SIGMA_1985_INFL` → **1.0** — ⚠ **REPLACED, not removed: 1985 keeps 2× weight via τ_C=0.06.**
  The old switch inflated independent per-observation noise; τ_C is a shared LEVEL offset, which is
  what the defect is. Stacking both would give an effective τ_C of 0.085 — the broad scheme by
  accident. **ONE launch, SIX jobs.** Full write-up (7 pp):
  `manuscript/HIKET_correlated_likelihood_proposal.tex`; memory [[correlated-likelihood-proposal]].
  Recovery point: `snapshots/20260818_pre_correlated_likelihood/`.
  **Design principle (Lorenzo): residuals may CHECK a design value, never SET one** — every τ comes
  from the observations with campaign levels removed.
  ⚠ **TRAP:** the observations' ICC is 0.799, but transplanting that RATIO onto a total of 0.800
  gives σ_e = 0.359, *below* the models' within-plot error — recreating overconfidence in the TREND.
  **Carry absolute SDs, never ratios.**
  ⚠ **The plot term ALONE is not worth doing** (precision ÷1.4, estimates move 1–2%). The REGION
  term is where the effect is: together they take n_eff 1205 → **260**, SE of the national mean ×2.15.
  ⚠ **Marginalise ≠ fit** — the "campaign bias forbidden" rule applies to the FITTED version
  (δ = −0.137). τ_C = 0 is also a choice, and the least defensible one.
  ⚠ Implementation is trivial: all variances fixed ⇒ Σ is CONSTANT ⇒ factorise once (~12 MB), then a
  ~1 ms solve per evaluation. **Woodbury/Sherman–Morrison NOT needed.** But the campaign term couples
  everything, so the per-plot `mclapply` sum becomes one global solve.
  ⚠ Build order = **local unit tests, not runs**: all τ=0 must reproduce the current ll EXACTLY.
  ⚠ **The non-negotiable check: τ = 0 must reproduce the current log-likelihood EXACTLY.**
  ⚠ Student-t does **not** decompose into shared + independent Gaussian, so the ν=6 tails and this
  feature do not stack — decide deliberately, do not leave both on by accident.
  ⚠ Log-likelihoods are **not comparable** across the change (the normalising constant moves) —
  compare on RMSE distributions (`S14_rmse_posterior`) instead.
  **What it should move, corrected 2026-08-18:** the correction removes LEVEL information; within-plot
  contrasts survive differencing untouched. So it loosens what the level pins — `sigma_input` and the
  `MRT × sigma_input` ridge — while `sigma_init` (pinned largely by trajectory SHAPE) moves only as a
  knock-on through their −0.29…−0.45 correlation. An earlier reading of the marginal prior/posterior
  widths said the opposite; that inference was wrong.
  ⚠ The per-plot ICC is the SMALLER of the two documented non-independence effects on the national
  level; the regional/spatial term (SE of the national mean ×2.3) is larger and is deliberately a
  SECOND step. Expect the first step to under-deliver on MRT.
  ⚠ **Submit with** `ssh roihu 'bash -ic "module load r-env && cd … && sbatch …"'` — `MODULEPATH`
  is interactive-only; a bare `ssh … sbatch` dies in 1 s.

- **✅ LANDED 2026-08-18 — jobs 695142–695147, the auxiliary-sigma-prior run.** RUN_IDs
  `20260817_120140` / `_120433` / `_120621` / `_120827` / `_120828` / `_120829`. One factor vs the
  previous run: sigma priors (widths 0.50→0.25, `sigma_input` centre 1.30→1.08); everything else
  byte-identical. All six COMPLETED, R-hat ≤1.008, ESS ≥1450, inf_rate 0.00%. **Snapshotted whole
  (825 MB) to `snapshots/20260818_pre_correlated_likelihood/` with a MANIFEST.**
  ✅ **The rc5113 co-location survived** — all three Yassos on the 2026-08-13 failure node, peaks
  13.4–15.6 GB against 160 GB. Node separation was never the whole story; nothing was killed.
  **Results:** (1) **MRT ROSE — and this was pre-registered as NOT an MRT run**: 21.68 / 25.13 /
  21.40 (Yasso07/15/20), i.e. +24 / +9 / +13% on this run alone. **Yasso20 now EXCEEDS its published
  POINT** (19.03); Yasso15's gap to published fell 24%→17%. (2) **It was nearly free in fit** —
  R² 0.011–0.024 → 0.010–0.022. ⚠ **This contradicts "MRT is likelihood-limited, not prior-limited"**;
  a prior-only change moved MRT further than the error-model lever ever did. (3) The pre-registered
  cost DID appear, but in RMSE not R²: calibration RMSE degraded in five of six (+0.25 to +1.06),
  holdout moved ≤+0.28, so the **calibration–holdout gap narrowed in ALL SIX** (Yasso15's went
  negative) — the constraint was paid for in fit that was not generalising. (4) `sigma_init` rose in
  five of six (Yasso15 0.180→0.298) but **the growing-stock bound is still violated in all six**
  (0.30–0.63 vs ~0.90). (5) `sigma_input` fell 15–25% but **tracked its prior centre nearly 1:1**,
  which sits awkwardly with the claim that the likelihood pins it ~3.6× more sharply than the prior.
  (6) Effective flux 4.25–5.28 ⇒ **five of six now exceed the decided Nordic window 4.62**.

- **⚠ ALL MRT NUMBERS RECORDED BEFORE 2026-08-18 MAY BE MIS-MEASURED.** `doublechecks/intrinsic_mrt.R`
  carried **hard-coded RUN_IDs** (line 46, last set 2026-08-13) and was never repointed, so the
  "15.23 / 22.10 / 17.70, unmoved" recorded for the 20260814 corrected-target run is actually the
  **20260812** posteriors — that run's MRT was never measured. True value **17.49 / 23.15 / 18.96**
  (+15 / +5 / +7%), so the corrected target DID move MRT. Fixed: the script defaults to the current
  run, **echoes its RUN_IDs at startup**, stamps them into `intrinsic_mrt.rds`, and accepts
  `HIKET_MRT_RID` for deliberate comparisons; `build_F12_mrt_yasso.R` now sources `run_ids.R` and
  **refuses to plot a mismatched .rds**. `manuscript/figures/run_ids.R` gained `HIKET_FIG_RID` for the
  same purpose. **Always check the echoed RUN_IDs before quoting any MRT.**

- **✅ LANDED 2026-08-15/17 — the corrected-target run, all six on one footing.** SP1/TP2/TP3
  `20260813_0803*`, Yasso07/15/20 `20260814_105614` / `_105757` / `_105907`. 1205 obs, σ=0.800,
  log-normal + Student-t ν=6, `SIGMA_1985_INFL` 2.00, R-hat ≤1.009, ESS ≥1727. **Snapshotted whole
  (823 MB) to `snapshots/20260817_pre_sigma_tightening/` with a MANIFEST** — `runs/` and
  `diagnostics/` are gitignored and Roihu scratch is purged at 180 days, so this is the only
  comparison basis for the next run. ⚠ `manuscript/figures/run_ids.R` auto-selects the NEWEST
  posterior, so the next run silently overwrites every figure — that is why the snapshot exists.
  **Results:** (1) ⚠ intrinsic MRT was recorded as **15.23 / 22.10 / 17.70**, "unmoved as
  pre-registered" — **that was a stale-RUN_ID artefact**; the true value is **17.49 / 23.15 / 18.96**
  and MRT DID move (+15 / +5 / +7%). See the correction item above; (2) **bias FLIPPED SIGN** — the uniform over-prediction is gone, now −0.6 to
  −2.5 (calib), −1.2 to −4.0 (holdout), so the "level offset" thread needs re-reading; (3) the
  **"2006→2024 is a source in all six" headline is DEAD** — it survives only in SP1 (−0.362) and
  Yasso20 (−0.181); TP2/TP3/Yasso07/Yasso15 bracket the observed +0.117 well; (4) over the full span
  every model is a sink but **1.6–2.5× too strong**, and over 1985→2006 **2.1–3.1×** — the misfit is
  concentrated EARLY; (5) skill unchanged, R² 0.011–0.024 calib / 0.000–0.004 holdout.

- **⭐⭐ THE OBSERVED-RATE REPORTING BASIS IS NOW DECIDED (Lorenzo, 2026-08-17): the baseline is the
  SAME DATA the models are calibrated on** — balanced plot set (observed in all three campaigns,
  **n=310**), unweighted, whole profile, **true observation years** from `obs_meta` (`year = 1984+idx`).
  This is what `manuscript/figures/obs_basis.R` already implements. ⚠ **The choice is worth ~1.8× on
  the headline rate:** 2006→2024 observed is **+0.117** on the balanced set but **+0.209** pairwise
  (n=409). **Every "+0.209" in this file and in older memories is the PAIRWISE figure** and must be
  restated before quoting. Full basis: 1985→2024 **+0.259** (Δt 35.1 yr), 1985→2006 **+0.399**.

- **✅ THE CHRONIC OOM WAS BROKEN 2026-08-15 — BY NODE SEPARATION, NOT MEMORY.** Jobs
  654390/654391/654392 completed (R-hat ≤1.009) with `cgroup_memlog` peaks of **14.5 / 12.0 / 14.0 GB
  against the 160 GB requested** and `events:max = 0` throughout — i.e. the same 12–16 GB every dying
  job ever reached, nowhere near any ceiling. The 16→40→80→160 GB ladder was treating the wrong
  variable. ⚠ **Do NOT "save" memory by lowering the request:** a smaller ask lets SLURM pack more
  jobs per node, which raises the node-level pressure that actually does the killing. The large
  request functions as a de-facto node reservation; `--exclusive` is the explicit form.
  Instrumentation retained: `Calibration_real_data_transient/cgroup_memlog.sh` logs `memory.peak` and
  `memory.events:max` at **both** job and step cgroup — the step's `memory.max` is `UNLIMITED`, so
  its breach counter alone reads 0 forever. **Verdict: `events:max > 0` ⇒ our own ceiling (and
  `peak` sizes it); `max = 0` with an `oom_kill` ⇒ the node did it and more memory buys nothing.**
  Deferred deliberately: per-chain checkpointing, and the fork-cluster rewrite (worth only ~3% —
  dispatch is 8.3 ms of a 227 ms evaluation). Memory: [[roihu-oom-instrumentation]].

- **⭐⭐ σ_input FLUX WINDOW IS ANCHORED ON THE WRONG STATISTIC (2026-08-14) — decided, NOT
  implemented.** `σ_input` is one global scalar ⇒ `σ_input × J̄` is a **national mean**; the 8.7
  ceiling is Gower 2001's **global** Class I evergreen **maximum**. **Bounding a mean with a
  maximum** is why it has never bound (posterior effective fluxes 4.99–6.56, all inside).
  **Decided: `[0.05, 4.62]`** (Gower **Nordic** Class I max) ⇒ σ_input ≤ 1.84, binds all six.
  ⚠ **UNIT TRAP: Gower 2001 = gC, Zheng 2004 = DRY MATTER** (Zheng cites Gower's world TNPP as
  109–1827 where Gower's own table gives 218–912 gC — ratio exactly 2.00–2.10); Zheng's 563 is
  **2.81**, not 5.63 tC/ha/yr. ⚠ Narrowing also **tightens the prior** (scaled logit onto the
  window). ⚠ Three claims died and must not be revived — "physically impossible", "LUKE is biased
  high", "Finnish forests are 2× Gower's" (that last compared *current annual increment* with
  Gower's *MAI*). **Also corrected: `J̄` INCLUDES harvest residues and natural mortality** (both
  previously assumed missing) and excludes understorey; Zenodo DOI 10.5281/zenodo.19736499 is now
  **live**. Write-up: M&M §"Prior specification: the litter-input flux window", `NEXT_SESSION.md`
  §0-quater. Memory: [[sigma-input-physical-bounds]].

- **✅ RESOLVED 2026-08-14 — NEXT SESSION STARTS BY IDENTIFYING THE ROIHU RUN (2026-08-13).** A corrected-target
  production run is believed to have been launched in an earlier session, but that is **unverified**
  and `Data/` may not have been synced first. **Decisive check in the launch log: 1205 observations /
  456 calib-ready = the CORRECTED target; 1269 = the OLD one** (discard and relaunch after rsync).
  Pre-registered expectation if corrected: 1985–2024 observations move to ≈+0.182, **2006–2024
  untouched at +0.107 against a source in all six**, MRT unmoved. **Then** the next launch is the
  `sigma_input` prior tightening (log-SD 0.50 → 0.15–0.20 from Lehtonen & Heikkinen 2015, ⚠ never
  tuned to land MRT on the published value), with F14's pre-registered fit costs as the test:
  Yasso20 ≈1.6 ll, Yasso15 ≈14.7, Yasso07 unreachable. **Watch R² and log-likelihood, not MRT.**
  Full protocol: `NEXT_SESSION.md` §0 and §0-ter.

- **✅ C5 WITHDRAWN 2026-08-05** — the 1985 (VMI8) campaign is NOT down-weighted;
  `SIGMA_1985_INFL` default is now **1.0**. Ablations: turning C5 off changes trusted-campaign
  RMSE by 0.45% (TP2) / 0.25% (SP1) vs a ~2% noise floor, while it moves σ_init by ~4×. Its
  premise (VMI8 reads low) is contradicted once the SOC homogenization is in: with 1985 fully
  trusted the model over-predicts it by only +3.2 tC/ha over its general bias, +6.3 when it has
  never seen it. The faithful δ-offset version IS identifiable (the revision plan's
  "non-identifiable" claim was too strong) but behaves as a misfit sink — δ = −0.137 requires
  σ_init ≈ 0.98 (80% of draws past the pre-run inversion threshold) and implies no accumulation
  ever happened. Residual concern (1985 unverifiable vs official LUKE stocks; different layer
  protocol) → **limitations, as a stated sensitivity**. Record:
  `manuscript/HIKET_data_and_ablation_tests.pdf`.
  **⚠ CORRECTED 2026-08-07 (addendum in the same PDF):** the "changes RMSE by <1%" ground was
  measured on **stock** RMSE only. σ_init barely moves the stock; it moves the **stock change**,
  which C5 controls strongly enough to **flip its sign**. Paired 1985–2024 change (obs +0.309)
  runs monotonically +0.638 → −0.184 (SP1) and +0.441 → −0.103 (Yasso20) as C5 weakens. C5 stays
  withdrawn (choosing the 1985 weighting to make the sink come out right is circular), but the
  **live question is now σ_init**, not the 1985 weighting. See "σ_init pre-run inversion" below.

- **✅ SETTLED 2026-08-19 — the σ_init growing-stock bound is RETIRED (provisionally).** Decision
  (Lorenzo): set it aside now, run the correlated likelihood without it, judge the result on
  BIOLOGICAL PLAUSIBILITY of the derived state; it may return. **Why it was the wrong target:**
  `0.903 = (1400/1775)^0.43` constrains the 1917 **litter flux**, but σ_init *also* equilibrates the
  1917 soil at that flux, and the NFI is silent on whether 1917 soils had caught up — after a century
  of slash-and-burn, raking and heavy cutting they plausibly had not, and the one-parameter transient
  init **cannot express "high flux, disequilibrated soil"**. ⚠ **The decider:** judging our own
  initial state by an *equilibrium-init* benchmark re-imports the exact convention the manuscript
  argues against (Lehtonen 2016, Palosuo 2008, Peltoniemi 2004). What survives is a criticism of the
  **parameterisation**, not of the posterior: the data want a low 1917 *stock*, the only dial that
  delivers one is σ_init, and it drags the *flux* down with it — on-thesis, and stronger than "five
  of six violate a bound".
  **Replacement test: `doublechecks/init_state_plausibility.R`** (25 draws, per-draw evaluation;
  C_1917 obtained by setting `lm$preinit_shape <- rep(0, 68)`, which holds the flux at J_1917 so the
  initialiser returns its own starting state — works for both families). Run `20260817_12*`:

  | | σ_init | J_1917 | J_1985 | **C_1917** | C_1985 | **pre-run rate** |
  |---|---|---|---|---|---|---|
  | SP1 | 0.630 | 2.32 | 3.68 | **39.5** | 49.1 | **+0.148** |
  | TP2 | 0.451 | 1.90 | 4.36 | **32.8** | 44.8 | **+0.181** |
  | TP3 | 0.425 | 1.73 | 4.14 | **32.6** | 44.6 | **+0.177** |
  | Yasso07 | 0.320 | 1.50 | 4.67 | **30.8** | 45.2 | **+0.217** |
  | Yasso15 | 0.294 | 1.27 | 4.22 | **30.1** | 47.0 | **+0.242** |
  | Yasso20 | 0.332 | 1.47 | 4.56 | **31.3** | 47.9 | **+0.235** |

  ⚠⚠ **TWO NUMBERS RECORDED 2026-08-18 WERE WRONG.** (1) The implied pre-run rate is **+0.148…+0.242,
  BELOW the observed +0.259** — not the "+0.37…+0.40, 1.5× observed" recorded. That criterion is
  currently **SATISFIED in all six**, and the "soil accumulates faster while its driver grew 4× slower"
  argument does **not** apply to this run. The old figure paired C_1917 with a 1985 endpoint of 57–63;
  the initialiser's actual endpoint is **44.6–49.1**. (2) 1917 stocks were mis-recorded (SP1 44.4,
  TP3 45.0); true 39.5 / 32.6. **So the STOCK criterion is the binding one, alone:** all six sit at
  **30–40 tC/ha** in 1917, and 20–30 is too little for a forest.
  ⬜ **TODO — source the floor.** Screen the boreal SOC literature and take the **MINIMUM**, so a
  failing model fails conservatively. ⚠ Comparable basis only: whole profile to per-plot `z_cap`
  (organic + mineral + modelled deep tail), **not** 0–30 cm, not mineral-only. `STOCK_FLOOR = 40`
  in the script is a **placeholder**. Peltoniemi 2004 is the natural anchor.
  ⚠ **Do NOT touch the σ_init prior before the correlated-likelihood run** — that launch is one-factor
  by design, and is not expected to relieve σ_init anyway. The τ_C revisit trigger stays **evidence**
  (post-run campaign residual spreads), never the σ_init value — that is the C5 circularity again.
  ⚠ Two extraction traps were hit building this and are live: the posterior `.rds` is **PHYSICAL**
  while `*_chains_*.rds` is **SAMPLING** space (running `to_original` on the posterior gave SP1
  `beta1` 0.099 → 1.104, ξ → 1.2e7), and `getSample()` on the saved list returns each internal DEzs
  chain's **first retained iteration** — extract per sampler with `start = 2`.
  `doublechecks/sigma_init_vs_growing_stock.R` is marked RETIRED in its header; keep it only as the
  record of where 0.90 came from, and **do not quote its verdict column**.
  Memory: [[sigma-init-growing-stock-bound]]. Superseded framing follows:

- **⭐⭐ ~~σ_init IS A MEASURABLE RATIO, AND FIVE OF SIX MODELS VIOLATE IT~~ (2026-08-17, SEE ABOVE).**
  The engine builds `J_1917 = J_t0_mean × σ_init × σ_input` while the forward run carries `σ_input`
  only, so σ_input cancels and **σ_init = J₁₉₁₇ / J₁₉₈₅**. The NFI speaks to that directly: fitting
  the litter–growing-stock elasticity on our own record gives **eps = 0.43 [0.23, 0.63]**, and with
  Korhonen 2024 `V(1917)/V(1985) = 1400/1775 = 0.789` that gives **σ_init ≈ 0.90** — reproducing the
  existing prior centre exactly. So the CENTRE was already right; the **0.50 WIDTH** was the defect,
  letting posteriors sit **2.3–5.0× below** it (Yasso15 **0.182** = a 1917 soil equilibrated to 18%
  of 1985 litter against a forest record saying 79%). That is how the models manufacture the sink:
  start implausibly depleted, still be climbing out. Same symptom as overshooting 1985→2006 by 2–3×.
  ⚠ The comparison is CONSERVATIVE — younger, more heavily cut 1917 stands shed *more* litter per m³,
  which pushes the implied value ABOVE 0.90. The gap is a floor. Check:
  `doublechecks/sigma_init_vs_growing_stock.R`. Fix wired 2026-08-17 (width → 0.25); memory
  [[sigma-init-growing-stock-bound]], [[sigma-priors-tightened-20260817]].

- **⚠ SUPERSEDED — the "0.818 inversion threshold" is OBSOLETE.** It came from
  `σ_init > J_t0_mean/J_full_mean`, but **P1 (2026-08-07) moved the 1917 anchor to `J_t0_mean`**
  (`tp2_wrapper_transient.R` ~354, `yasso15_wrapper_transient.R` ~282). Both ends now use the same
  flux, so the pre-run builds iff **σ_init < 1**, full stop. Any "[0.64, 0.82] window" is wrong on
  both ends. ⚠ `build_F4_initialization.R` still used `J_full_mean` until 2026-08-17 and misled an
  entire analysis — **read the wrappers, not the figure scripts, for what a parameter means.**
  Superseded text follows:

- **🚩 ~~σ_init PRE-RUN INVERSION~~ (SUPERSEDED, see above).**
  The pre-run **declines** exactly when `σ_init > J_t0_mean/J_full_mean` (σ_input cancels; the
  pre-run starts *at* equilibrium and the flux moves monotonically). Per-plot ratio median
  **0.818** over 512 plots — this reproduces the 0.826 from the ablations and is a property of
  the **litter record, not of any model**. Production run 20260805: Yasso15 σ_init **1.007**,
  Yasso20 **1.243** → pre-run declines for **78% / 92%** of plots, and both report a **source**
  (−0.063, −0.155) where observations show a **sink** (+0.312). Survives dropping 1985: on
  2006–2024 alone obs +0.209 vs Yasso15 −0.183, Yasso20 −0.220. Stock change crosses zero at
  σ_init ≈ 0.8 in **all four** models tested (SP1 single-pool through Yasso20). **Proposed fix:**
  truncated/informative prior at the *national aggregate* ratio `σ_init ≲ 0.82` — a per-plot form
  is impossible (σ_init is one global scalar; the ratio spans 0.009–5.983). Binds only
  Yasso15/20. Caveats: it partly *imposes* the direction reported (defensible — the constraint is
  independent NFI growing-stock data, not the SOC being fitted — but must be stated); check the
  combined feasible region against `flux_pair` (which bounds the flux *level*, not its
  *direction*); watch ESS on σ_init. Scripts: `doublechecks/prerun_direction.R`,
  `ablation_stock_change.R`, `production_fit_by_campaign.R`.

- **✅ LANDED 2026-08-11 — jobs 563524–563529.** RUN_IDs `20260810_152914` (SP1/TP2/TP3),
  `_152915` (Yasso07), `_152916` (Yasso15), `_152917` (Yasso20). Stages 2–4 and 24/25 manuscript
  figures regenerated locally the same day. **Results:** (1) σ=0.72 **self-consistent** — residual
  sd 0.712–0.744 in all six, so the plug-in choice is validated, not circular; (2) **all six
  converged**, including TP3 and TP2 (were psrf 18.7 / 15.6) ⇒ TP3 retraction above; (3) intrinsic
  MRT rose 18–34% to **15.17 / 21.91 / 17.37** vs published 33.39 / 30.27 / 25.02 — for Yasso15/20
  that rise is attributable to σ **alone** (their climate priors were untouched); (4) σ_input fell
  only 7–16% (now 1.84–2.49); (5) **not** fixed — bias grew to +3.7…+6.8 tC/ha, R² still ~0, and
  2006–2024 is a source in all six vs observed +0.209. Level ridge confirmed: Yasso15's
  MRT×σ_input 43.0→44.7, constant to 4%. **Next lever (NOT launched, pending coauthors):**
  σ_input log-SD 0.50→0.20 ⇒ predicted Yasso15 MRT ≈ 34 yr. Superseded launch note follows:

- **🚀 ~~ROIHU RUN IN FLIGHT~~ — jobs 563524–563529, 2026-08-10 ~15:30** (SP1 563524, TP2 563525,
  TP3 563526, Yasso07 563527, Yasso15 563528, Yasso20 563529). Commit `3b0d533`. ~13–19 h ⇒
  results 2026-08-11. **Tests ONE factor: the error-model scale.** `HIKET_SIGMA_TOTAL=0.72`
  replaces `sigma_obs_fixed` (0.442) — the latter is the MEASUREMENT CV but was used as the TOTAL
  error, while measured log-residual spread is 0.708–0.735 in all six models (implied model error
  0.55–0.59, larger than observation error). Since the level penalty goes as 1/σ², that amplified
  stock-level pressure ~2.6× relative to the priors — a candidate driver of the short bulk MRT.
  Also in this run: double precision (item below) and Tuomi 95%→1σ prior widths (item below).
  Fraction prior deliberately left at 0.4 so the error model is tested alone.
  ⚠ Set via `SINGULARITYENV_HIKET_SIGMA_TOTAL` in the SLURM scripts — **a bare export never
  reaches R inside the r-env singularity container.** Verify with
  `grep -H -E "ERROR MODEL|Cores per chain|chains x" *_5635*.err` — expect TWO ERROR MODEL lines.
  **→ `NEXT_SESSION.md` has the post-run checklist.**

- **⚠ ALL MRT NUMBERS BEFORE 2026-08-10 ARE SUPERSEDED.** The engine binding named `steady_state`
  is NOT a steady state — for Yasso it is `*_transient_init` (1917 equilibrium + 68-yr ramp to
  1985). It is contaminated by `sigma_init` (25.05 at 0.90 vs 33.84 at 0.35) and diverges for
  near-conservative draws through `model_step`'s Euler fallback. **Use
  `doublechecks/intrinsic_mrt.R`**: unit litter input at a fixed reference (dataset-mean climate +
  AWEN × size composition), pure steady-state routine ⇒ a property of the generator, independent
  of `sigma_input`/`sigma_init` and of the SOC data. Correct basis: ours **11.3 / 17.7 / 14.7** vs
  published **33.4 / 30.5 / 25.0** (Yasso07/15/20). The displacement is *demanded* (147.7 / 93.4 /
  42.6 nats at σ=0.442), not merely permitted by loose priors.

- **✅ FORTRAN PRECISION FIXED 2026-08-10.** `yasso15.f90` ran SINGLE precision while
  `yasso07.f90` ran double, so Yasso15/20 were computed at ~7 significant digits and
  SP1/TP2/Yasso07 at ~16 — a confound for a structural intercomparison. Switched to double, with
  `yasso15_wrapper_transient.R` converted in the same commit: **both sides must move together
  (`as.single`/`single()` ↔ `as.double`/`double()`) or the ABI mismatches.** VERIFIED to change no
  result (posterior MRT bit-identical; ll at defaults −1695.1308 → −1695.1312). The `.so` is
  gitignored ⇒ **must be recompiled on Roihu**, each `.f90` in its own `SHLIB` call.

- **✅ PRIOR WIDTHS CORRECTED TO SOURCE 2026-08-10.** Tuomi 2009 Table 3 and Tuomi 2011 Table 4
  both state *"95% confidence limits"*, but the locked convention (decision, 2026-06-05) read "±"
  as **1σ** — ~1.96× too wide. Corrected for Yasso07 and, since they inherit "Yasso07 scale", for
  SP1/TP2/TP3: `beta1` 0.26→0.1327, `gamma` 0.20→0.1020, `beta2`, and Yasso07's
  `delta1`/`delta2`/`r`. **Centres unchanged.** Yasso15/20 NOT touched — their widths are genuine
  posterior SDs from the `.dat` samples. Halves the ensemble's climate asymmetry (simple models vs
  Yasso15: 5.7× → 2.9×).

- **⚠ RETRACTED 2026-08-11 — TP3 *DOES* CONVERGE.** In run `20260810_152914` all 9 parameters have
  R-hat < 1.05 (multivariate psrf 1.0053); `p_S` is cleanly unimodal at 0.455 (90% 0.354–0.568)
  where the old quantiles ran 0.004/0.006/0.583/0.635, and `gamma` no longer reaches the
  climate-off value. **The bimodality was an artefact of the over-wide climate priors** (Tuomi
  95%-as-1σ), not structural non-identifiability of the third pool. Attribution is not fully
  settled — this run moved climate widths and σ together, and the local test that appears to
  separate them ran 3×6000 vs production 5×50000, so its R-hat may be incomplete mixing. The
  complexity thread loses this mechanism. Superseded text follows:

- **🚩 ~~TP3 DOES NOT CONVERGE — AND THAT IS THE RESULT.~~ (SUPERSEDED, see above)** Its two modes fit within **4.2 nats**
  (the collapsed `p_S`→0 mode slightly BETTER), posterior mass 60/40. The data cannot distinguish
  a 3-pool cascade from one effective pool with a flat climate response ⇒ **the third pool is not
  identifiable from 2 SOC obs/plot**, and R-hat 18 is the correct output for a bimodal posterior,
  not a sampler failure. Structural: TP2's `[1/a_A + p_H/a_H]` is a 1-D ridge; TP3's
  `[1/a_A + p_S/a_S + p_S·p_H/a_H]` is a 3-D manifold, and `p_S`→0 frees `a_S`, `a_H` AND `p_H`
  at once. Fraction tightening cannot fix it — the degeneracy spans the RATES, deliberately left
  at SD 0.15 as the ICBM transferability diagnostic. **Report it rather than fix it**: it gives
  the complexity thread a mechanism (SP1/TP2 identifiable, TP3 not) instead of only a skill
  comparison.

- **✅ LITTER 2006 PEAK IS REAL — H2 CLOSED 2026-08-10** (confirmed with B. Tupek). The series
  stands as published; a modest post-2006 decline in modelled SOC is EXPECTED and defensible, not
  an artefact to remove. Context: litter rises +51.7% over 1986–2006 vs growing stock +22.4%, then
  falls −12.8% over 2006–2021 vs growing stock +16.0% (net ×1.32 vs ×1.42). Do not re-open.

- **✅ SUPERSEDED — Roihu recalibration jobs 474800–474805, launched 2026-08-05 13:24.**
  First run to include C1–C4 *and* both data corrections. ~36 h → expect ~01:00 on 2026-08-07.
  Verified at launch: `Cores per chain: 40` (OOM trap avoided), `5 chains x 50000` (no ablation
  env vars leaked), all six reached `Chain 1/5` → the input guard, forward sanity and the Yasso
  `.so` all passed. **NB guard output goes to `.err`, not `.out`** (R `message()` → stderr), and
  restrict globs to `*_4748*` or you read July's logs.
  **→ `NEXT_SESSION.md` has the full check / sync / refresh sequence.** Everything downstream
  (figures F1–F14, S1–S7, T1/T2, NextGenC bundle, `HIKET_calibration.Rmd`) is stale until it lands.

- **✅ SUPERSEDED — six-model Roihu re-calibration (was: next action).**
  The SOC baseline swap is **DONE and pre-flighted** (2026-08-04, see Data §"SOC
  calibration target"). Everything downstream now predates it. State:
  - ✅ `Data_work.R` rewired + re-run; 66/66 sanity checks pass; 520 calib-ready
    (416/104); `soc_obs_tCha` == `soc_profile_Mgha` exactly; 0 raw-sum fallbacks.
  - ✅ Pre-flight gate cleared: `preflight_prior_pushforward.R` (K=300, full N) on the
    NEW target — forward sanity **PASS in all six**; blow-up rate **0%** for
    SP1/TP2/TP3/Yasso15/Yasso20 and **0.6%** for Yasso07 (1 of 156 constraint-passing
    draws; historically ~2–3%, a `delta2`/σ baseline effect, not `beta2`). Notably
    **TP3 went 10–12% → 0%**: the ICBM rate anchors and `flux_pair` σ_input bounds
    were tuned on the old inflated target, and they survive the gentler one. The
    feared σ_input lower-bound clamp did **not** materialise (window [0.021, 3.647],
    J̄ = 2.386 tC/ha/yr). (The 144/86/86 "prior_rejects" in the Yasso models are
    stick-breaking simplex rejections — legitimate, never visited by the sampler.)
  - ⬜ **TODO:** sync `Data/` + code to Roihu, launch all six calibrations, then
    stages 2–4 locally. Expect: σ_input to fall further (target no longer inflated);
    a possible R² shift; the 1985 over-prediction gap to change (1985 is now
    59.4 vs 2006 66.7, not 63 vs 102) — **the "26 tC/ha above observed 1985" figure
    and every level quoted against a campaign mean must be restated.**
  - ⬜ **Then refresh:** all figures F1–F14 + S1–S7, T1/T2, the NextGenC report
    bundle, `HIKET_calibration.Rmd`, and drop the provisional-levels reading note in
    `manuscript/HIKET_storyline_note.tex`.

- **🚩 ROIHU-PHASE PRIORITY — physically-bounded input priors (re-calibration).**
  *Decided 2026-07-01. A SessionStart hook (`.claude/settings.json`) auto-surfaces
  this on the first Roihu session; do it before any new calibration.* The problem:
  `sigma_input` is homogenized as a weak Tier-3 nuisance (log-transform, centre 1,
  log-SD 0.5) — but litter inputs have a **physical referent** (Tupek litter model +
  NPP ceiling) that this prior throws away. Result: the posterior blows through by
  5–6σ. Effective litter flux = `sigma_input × raw_J` (raw ≈ 2.5 tC/ha/yr, shared by
  all six models): **SP1 0.8, Yasso07 3.5, Yasso15 4.8, Yasso20 7.5** (all physical,
  ≤ ~9 NPP ceiling) but **TP2 34, TP3 50 tC/ha/yr — 13–20×, above boreal NPP =
  physically impossible.** TP2/TP3 buy their (best) R² by manufacturing carbon; same
  root cause as the TP3 oscillation (cascade-starved fast pool force-fed). **The fix:**
  put a PHYSICAL bound on the effective *flux*, not the abstract multiplier, homogeneous
  across models — `sigma_input × J_raw ∈ ~[0.5, 9]` tC/ha/yr AND `sigma_init ×
  sigma_input × J_full ∈ ~[0.5, 9]` (the two sigmas couple via the 1917 pre-run flux;
  `J_1917 = J_full × sigma_init × sigma_input`). Implement as a truncated-lognormal on
  the fluxes (or hard cap `sigma_input ≲ 3.5`, `sigma_init ∈ [0.1, 1.5]`). **Note:**
  `sigma_init` is ALREADY physically behaved (all posteriors 0.19–0.72, < 1 = the
  post-exploitation recovery narrative) — it needs only a light guardrail; `sigma_input`
  is the one to fix. **Consequence:** removes TP2/TP3's escape hatch → forces misfit
  onto the (now prior-pinned) rates or into a VISIBLE R² hit = exposes structural
  inadequacy instead of hiding it (sharpens the thesis). **Caveat:** full six-model
  re-calibration + mixing check (leveraged param clamped near a boundary may mix worse).
  Verification scripts: `doublechecks/sigma_input_spread.R` (per-model posterior +
  σ-distance) and `doublechecks/effective_litter_flux_vs_physical.R` (effective flux
  vs bounds). Docs: manuscript outline §"multiplier crosses a
  physical line" + §"Physically-bounded input priors"; memory [[sigma-input-physical-bounds]].
  To silence the hook once done: `touch .claude/.roihu_input_bounds_done`.

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
| Aleksi Lehtonen | LUKE; **coauthor**; lead of the national soil-C inventory (Lehtonen et al. 2016, GMD — the steady-state approach HIKET extends). Frame that lineage collegially (evolving the group's own method), not as a critique. Corroborator for the stock QC (D3) & natural owner of understorey litter (D2). **NB the stock QC is OURS, not delegated** — we did it 2026-08-04 by moving to a standardized source (see Data §SOC target) and we keep verifying it; Aleksi's role is to review and confirm, not to perform it. On holiday until ~mid-Aug 2026. |
| Samuli Launiainen | **coauthor**; gave the round-1 review comments (annotated `manuscript/revisions/HIKET_storyline_note_sl.pdf`, incl. "ridiculously flawed"). Straightforward, clear-thinking; friends with Aleksi. His candid margin notes are internal — keep the manuscript prose measured. |
| Mikko Peltoniemi | LUKE; author of the keystone Peltoniemi et al. 2004 (the direct antecedent — identified the non-equilibrium init problem). Close collaborator (~2 yr), on good terms — approachable for the litter/understorey/stock questions while Aleksi is away. |

**GitHub:** `ilmenichetti/HIKET`  
**Prior sources:** FMI Ryassofortran repository (`https://github.com/YASSOmodel/Ryassofortran/tree/master/data`); `Priors_model_matching.R` derives prior specs.
