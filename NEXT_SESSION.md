# NEXT SESSION — start here

## 🧪 0-NOW. ARM B PRE-REGISTERED, NOT YET LAUNCHED (2026-08-31)

> **What it is.** One factor vs the Liski-input run (`20260820_1554*`, arm A): the `sigma_input`
> prior CENTRE, 1.27 → **1.08**. Width stays 0.25; the Liski per-hectare pre-1985 ramp stays ON;
> everything else byte-identical. Wired in all six `Prior_specs/*_priors.R`.
>
> ### ⚠ This is a SENSITIVITY ARM, not a correction
> 1.08 is NOT "the right value we drifted from". The two anchors genuinely disagree on the
> understorey share — **Liski 2006 says 27% (⇒ 1.27), Lehtonen & Heikkinen 2015 says 8% (⇒ 1.08)** —
> and the prior files record the disagreement as UNRESOLVED, pending **A. Lehtonen, who authored
> both**. So: run both, report both, let Aleksi settle the anchor. Framing it as a revert would be
> the tuning trap CLAUDE.md warns about in three places, because 1.08 also happens to move MRT
> toward the published values.
>
> **Second, independent reason to run it:** arm A moved TWO factors (ramp + centre), so its MRT drop
> is attributed by inference from the near-conserved ridge, not demonstrated. Arm B isolates the ramp.
>
> ### Pre-registered expectations (extrapolated from arm A's response; NOT intervals)
> Arm A moved the centre +17.6% and the posteriors followed +6.3…+9.1%, ridge product near-conserved.
> Running that backwards:
> - `sigma_input` posteriors **−6 to −9%**
> - intrinsic MRT **+7 to +9%** ⇒ ≈ 24.0 / 27.6 / 22.5 (Yasso07/15/20), restoring Yasso15's overlap
>   with the published value
> - effective flux ⇒ ≈ SP1 3.14 · Yasso15 3.76 · Yasso07 3.98 · Yasso20 4.11 · TP2 4.36 · TP3 4.38,
>   i.e. all six back under the Gower Nordic max 4.62
> - **C_1917 unchanged**, `sigma_init` unchanged (arm A moved it −0.3…+5.1%)
> - RMSE within ±0.3 either way
>
> **If MRT moves much MORE than +9%**, the ridge is not as conserved as arm A suggested and the
> two-factor attribution needs redoing. **If C_1917 moves**, something other than the centre changed.
>
> ### ⚠ The flux window stays WIDE — decided 2026-08-31 (Lorenzo)
> The hard window `[0.05, 4.62]` is **deliberately NOT implemented**. Reason: on arm A the absurdity
> it existed to kill is gone (TP2/TP3 are at 4.71/4.73, not the old 34/50), and bounding would
> IMPOSE the input conclusion instead of finding it. **Consequence for the manuscript:** the
> Discussion's §"multiplier crosses a physical line" counterfactual (bounding costs skill only for
> TP2/TP3) describes an experiment that was never run and numbers that no longer exist — rewrite it
> as a straight comparison against the literature anchors.
>
> ### ⚠ Do NOT overstate the input finding
> Decided 2026-08-31: "previous calibrations underestimate inputs" is defensible and stays, but is
> NOT to be strengthened. It and "MRT too short" are the same ridge from two ends — arm B weakens
> the first and helps the second by construction. Do not bank both.
>
> ### Launch checklist
> - ✅ Snapshot of arm A: `snapshots/20260831_pre_sigma_input_armB/`
> - ⬜ commit + push the prior change; `git pull` on Roihu; confirm clean tree
> - ⬜ `ssh roihu 'bash -ic "module load r-env && cd /scratch/project_2019134/HIKET && sbatch ..."'`
>       — MODULEPATH is interactive-only; a bare `ssh … sbatch` dies in 1 s
> - ⬜ verify in the `.err` logs: `sigma_input` centre **1.08**, `Cores per chain: 40`,
>       `5 chains x 50000`, `liski | mean 0.450` (ramp still on), 1205 obs / 456 calib-ready
>
> ### Still open, independent of this run
> - ⬜ `metrics_calib` fix in `run_multimodel_comparison.R` §2 (quantified 2026-08-31: true calib
>       R² **0.0150–0.0206**, holdout R² **1.2e-07 to 1.1e-03** — holdout skill is ZERO, not "low")
> - ⬜ F3/T1 cache-key fix is FIXED but UNCOMMITTED
> - ⬜ T1 caption still prints the superseded `[0.5, 8.7]` window
> - ⬜ Discussion §"over-prediction is a property of the likelihood" — premise FLIPPED SIGN, dead as written
> - ⬜ manuscript is a skeleton where it matters: Results unwritten, 0 `\cite`, 0 `\includegraphics`,
>       no `.bib`, title/abstract/coauthors all `[PLACEHOLDER]`

---


## 🚀 0-NOW. THE LISKI-INPUT RUN IS ANALYSED — NEXT MOVE DECIDED (2026-08-24)

> All six of jobs **749703–749708** COMPLETED (Yasso20 took 30 h 25 min, clearing its wall by ~5.5 h).
> R-hat 1.000–1.009, ESS 1752–11002. Synced, stages 2–4 run, all 30 manuscript figure builders
> re-run. Provenance verified: **Roihu ran `d4a223e` with a clean tree**, local HEAD differs only in
> this file, so the priors the KL/marginal figures compare against are the ones that ran.
>
> ### Results, against the pre-registration
> - ⭐ **C_1917 HELD at 43.4–55.9 tC/ha** (prev 44–56) — the pre-registered failure mode did NOT fire.
>   σ_init barely moved (−0.3% to +5.1%).
> - **σ_input rose only ~half the prior's step**: centre +17.6%, posteriors **+6.3 to +9.1%**
>   (SP1 1.350 · Yasso15 1.619 · Yasso07 1.714 · Yasso20 1.770 · TP2 1.876 · TP3 1.885).
>   **Breaks the 20260817 1:1-tracking pattern.**
> - ⚠ **MRT FELL**: 25.66→**22.13** / 27.35→**25.60** / 22.42→**20.82** (Yasso07/15/20).
>   Yasso15 **lost** its overlap with the published posterior median.
> - **Ridge product back to near-conserved** (Yasso15 −0.5%, Yasso20 −0.6%, Yasso07 −6.3%). So the
>   MRT loss is the mirror of the σ_input rise ⇒ traces to factor 2. **Strong inference, not proof:
>   two factors moved.**
> - RMSE +0.28…+0.81 calib, +0.31…+0.84 holdout. Stock change: **1985→2006 best yet at 0.99–1.37×**;
>   full span 0.32–1.23×; 2006→2024 a SOURCE in 4 of 6.
> - ⚠ **TP2 (1.876) / TP3 (1.885) now exceed the σ_input ≤ 1.84 implied by the decided flux window
>   `[0.05, 4.62]`** (effective flux 4.64 / 4.66 vs a 4.62 ceiling). `T1_model_summary.tex` still
>   prints the superseded `[0.5, 8.7]` and claims "all within" — caption to update.
>
> ### ⭐⭐ NEXT ACTIONS, decided by Lorenzo 2026-08-24, in order
> 1. **Re-run the calibration with the `sigma_input` prior centre back to 1.08** — ONE factor,
>    isolating whether the Liski ramp alone is a gain. Only that prior line changes (the ramp has
>    `HIKET_PREINIT_SHAPE`). ⚠ Snapshot first — `run_ids.R` auto-selects the newest posterior.
> 2. **Then the forward exploration** — see `memory/equifinality-forward-exploration-design.md`:
>    LL-threshold cut on the posterior, **uniform resampling** in the retained region, compare two
>    MRT/σ_input regions of near-equal fitness for their **extrapolation**. Profiling was proposed and
>    **REJECTED** (parameter interactions; the optimiser's bias manufactures discrimination).
> 3. **Deferred a few days — Lorenzo is consulting FMI first.**
>
> ### ⭐ The threshold blocker is largely solved, at zero machine cost
> `manuscript/figures/F14_mrt_ridge.rds` already holds per-draw `mrt`/`si`/`ll` for **225,015 draws**.
> Δll at the published MRT POINT = **13.70 / 5.03 / 1.13** (Yasso07/15/20) against each posterior's
> own 95% LL spread of **9.9–12.1**. ⚠⚠ **"Yasso07 unreachable" is DEAD** — its MRT now reaches
> **35.3 yr** (above the published 33.47). **Δ = 6 is the sweet spot** (2.6–9.3% retained, spans both
> regions); Δ ≤ 3 never reaches the published values; Δ = 10 retains 30–65% and the contrast
> dissolves. **Better: stratify along MRT instead of picking one global Δ** — then Δ is an output, and
> machine time is set by (bins × draws per bin) with no wasted volume.
>
> ### 🚫 Blocking prerequisite for the climate half
> **There is NO climate change in any projection** — all six recycle the last `RECYCLE_YEARS = 20`
> observed years cyclically over `PROJ_YEARS = 60`. Scenario forcing must be added first, through
> `compute_xi_mean()` (never `mean(xi_array)`).
>
> ### ⚠ Two defects found 2026-08-24
> - **NOT FIXED:** the multimodel "calibration" metrics read `metrics` (all 1205 obs, holdout
>   included), not `metrics_calib`. True calibration R² is **0.015–0.021**, not 0.010–0.015; any
>   calib–holdout gap claim is mechanically shrunk. `T1` is correct ⇒ two builders disagree.
> - **FIXED but UNCOMMITTED:** `build_F3_mean_soc.R` and `build_T1_model_summary.R` read the old
>   `F4_cache_bal_` key after F4 was renamed to `F4_cache_bal_c3_`, so both had been erroring since
>   2026-08-20 and the F3/T1 in the manuscript came from a **pre-C3 linear-ramp cache**.
> - Also unwired: `doublechecks/paired_stock_change.R` **ignores `HIKET_FIG_RID`**.

---

## 📦 SUPERSEDED — 0-NOW as of 2026-08-21 16:15

> **✅ Jobs 749703–749708** (SP1 / TP2 / TP3 / Yasso07 / Yasso15 / Yasso20), submitted
> 2026-08-20 ~14:40 from commit `d4a223e`. **They started 2026-08-20T15:54 — 1.2 h after
> submission, NOT the 2026-08-23 SLURM estimated.** Ignore `squeue --start` estimates in this
> project; they have been pessimistic by days.
>
> **FIVE COMPLETED CLEAN** (exit 0:0, all chains 100% finite, 0 R-hat warnings):
>
> | job | model | RUN_ID | wallclock | R-hat | ESS |
> |---|---|---|---|---|---|
> | 749703 | SP1 | `20260820_155430` | 888 min | 1.000–1.001 | 9911–11002 |
> | 749704 | TP2 | `20260820_155431` | 949 min | 1.000–1.002 | 6999–7837 |
> | 749705 | TP3 | `20260820_155430` | 1046 min | 1.000–1.004 | 5730–6541 |
> | 749706 | Yasso07 | `20260820_155430` | 1002 min | 1.001–1.005 | 2192–3643 |
> | 749707 | Yasso15 | `20260820_155430` | 934 min | 1.001–1.008 | 1817–3175 |
> | 749708 | Yasso20 | `20260820_155433` | chain 4/5 at 16:15 | — | — |
>
> **Yasso20 was projected to finish ≈23:45 on 2026-08-21**, ~4 h inside its 36 h wall
> (2026-08-22T03:54). It was the wall risk of this run — chains ran 287 / 457 / 461 min under
> heavy `rc5120` contention, then chain 4 *sped up* to 2.37 eval/s while the node got busier
> (another entry for [[roihu-runtime-is-node-contention]]). **VERIFY it completed before using
> anything**: `sacct -j 749708 --format=State,Elapsed,ExitCode`. If it hit the wall it is a total
> loss (no per-chain checkpointing) and needs a solo relaunch.
>
> ### ✅ FIRST CHECK — PASSED IN ALL SIX, nothing to redo
> ```
> [PREINIT SHAPE] liski | mean 0.450 | 1950 0.760 | 1970 0.645
> [ERROR MODEL] CORRELATED likelihood: n = 959 obs, 365 plots, 8 bands, 3 campaigns
> Cores per chain: 40        5 chains x 50000
> ```
> `sigma_input = 1.27` confirmed in `Prior_specs/*_priors.R` at `d4a223e`. So this IS the clean
> two-factor step; the likelihood is unchanged from the 08-19 correlated run.
> Logs live in `Calibration_real_data_transient/progress_logs/*_7497*.{err,out}` (NOT the repo root).
>
> ### NEXT ACTIONS
> 1. `rsync` back `runs/`, `diagnostics/` **and** `Data/model_inputs/` — the predictive stage
>    hard-loads `<MODEL>_inputs_<RUN_ID>.rds` with no fallback.
> 2. Stages 2–4 locally: `Rscript --no-save …/run_hiket_pipeline.R --skip-calibration`.
> 3. Read the results against the pre-registration below.
>
> **TWO FACTORS, deliberately** (decision Lorenzo 2026-08-20; every run since 08-17 was
> one-factor, this one is not):
> 1. **Pre-run ramp**: national-total growing stock → the Liski et al. 2006 reconstructed
>    **input to soil, per hectare, tree basis**. Normalised shape mean **0.194 → 0.450**.
> 2. **`sigma_input` prior centre 1.08 → 1.27** (Liski's understorey correction, 2.88/2.27).
>
> Comparison basis: `snapshots/20260820_pre_liski_recentring/` (862 MB + MANIFEST).
>
> ### PRE-REGISTERED — read the results against this
> - **`sigma_input`** should rise toward ~1.27 (effective flux → Liski's 2.88). Lorenzo expects
>   little change. If it tracks the centre 1:1 *again*, that repeats the 20260817 pattern and
>   undercuts "the likelihood pins σ_input".
> - ⭐ **WATCH THE 1917 STOCK, NOT `sigma_init`.** (Correction, Lorenzo 2026-08-20: my
>   pre-registration had this wrong.) `sigma_init` is a flux ratio and `C_1917 ∝ sigma_init ×
>   sigma_input`, so σ_init falling while σ_input rises leaves the stock untouched. **σ_init
>   falling is not a problem; C_1917 falling a lot is.** Run
>   `doublechecks/init_state_plausibility.R` — this run's C_1917 was **44–56 tC/ha**, up from
>   30–40. Losing that is the failure mode.
> - **MRT** should fall slightly (MRT × σ_input has been near-conserved).
> - **The litter-history disagreement will NOT close**: models infer +66–97% per-hectare litter
>   rise 1917→1985, Liski gives **+14.6%**, NFI+elasticity **+10.7%**. That lives in σ_init's
>   amplitude, which neither factor touches.
> - Compare on **RMSE distributions**, never log-likelihood.
>
> ### Reverting without a code change
> `HIKET_PREINIT_SHAPE=growing_stock` restores the old ramp **bit-identically** (verified,
> max |diff| = 0); `=linear` or the legacy `HIKET_PREINIT_LINEAR=1` gives the linear ramp.
> ⚠ On Roihu it needs the `SINGULARITYENV_` prefix to reach R inside the container.

---

## ⭐⭐ STANDING ITEM — forecast robustness: SINK vs SOURCE (raised 2026-08-20, Lorenzo)

**Politically decisive and not yet defensible. Re-check after every change to the input
specification.** Memory: `forecast-sink-robustness`.

Measured on run `20260819_1025*`, 2024→2084, balanced set, tC/ha/yr:

| TP2 | TP3 | Yasso07 | Yasso15 | Yasso20 | SP1 |
|---|---|---|---|---|---|
| **+0.218** | +0.206 | +0.101 | +0.081 | +0.060 (→ SOURCE by the 2070s) | **−0.001** |

**The projection freezes litter input at its 2024 value**, so the entire sink is disequilibrium
relaxation — the soil catching up to *past* input increases. Consequently:

- **r(inferred historical litter rise, projected sink) = +0.870**; r(σ_init, sink) = −0.800,
  Spearman −0.943.
- ⚠⚠ The inferred historical rise (+66–97%) is precisely what Liski (**+14.6%**) and
  NFI+elasticity (**+10.7%**) contradict. **If the input history is 4–6× too steep, so is the
  projected sink.** SP1 — the only model with a defensible input history — is the only one
  projecting no sink.

**To do, in order:** (1) re-measure after the Liski-ramp run lands; (2) treat constant-2024 input as
a *scenario* and test alternatives (the post-2006 litter decline is real); (3) **verify whether
CLIMATE is also frozen in the projection — inputs demonstrably are, climate was never checked**;
(4) report the ensemble spread, never a central estimate; (5) pair with the equifinality result —
fit explains only 10–18% of forecast spread, so a good fit says nothing about a right forecast.


## 🚨 0-NOW. STATE AT 2026-08-17 — READ THIS FIRST

> **⏳ IN FLIGHT: jobs 695142–695147** (SP1/TP2/TP3/Yasso07/Yasso15/Yasso20), launched 2026-08-17
> from commit `ffce254`. **Tests ONE factor: the auxiliary-sigma priors** — widths 0.50→0.25 and
> `sigma_input` centre 1.30→1.08. Everything else is byte-identical to the run it will be compared
> against, so the comparison is clean. Results ~2026-08-18.
>
> **READ IT AS A `sigma_init` RUN, NOT AN MRT RUN.** Pre-registered: `sigma_init` rises (prior and
> likelihood are comparably informative, ~45% prior weight), `sigma_input` falls only slightly (the
> likelihood pins it ~3.6× more sharply than the prior, ~7% weight), and **R² should FALL** — that is
> the intended cost of a physically-bounded prior, not a regression. ⚠ **Magnitudes are NOT
> predictable.** A joint reweight of the previous draws returned **ESS = 1 of 1001**, because
> `sigma_init` has to travel ~6.4 prior SDs into territory the old posterior never visited. Direction
> only. The real question is whether the 2006→2024 sink survives once the depleted-1917 crutch is
> removed, and what it costs in log-likelihood (comparable here — σ stays 0.800).
>
> ⚠ **All three Yassos landed on the SAME node (rc5113)** — the configuration associated with the
> 2026-08-13 triple failure. Left to run deliberately: the node reports 1.33 TB free and peaks are the
> usual 9–14 GB, and if it dies we finally get a co-located failure WITH cgroup telemetry. Check
> around hour 8, historically when the window opens.
>
> **✅ THE CHRONIC OOM IS BROKEN — by NODE SEPARATION, not memory.** Jobs 654390–2 completed with
> peaks 14.5/12.0/14.0 GB against 160 GB and `events:max = 0`. The 16→40→80→160 GB ladder was the
> wrong variable. ⚠ **Never lower the memory request to "save" resources** — a smaller ask lets SLURM
> pack more jobs per node, raising the node-level pressure that does the killing. The big request is
> a de-facto node reservation. See memory `roihu-oom-instrumentation`.
>
> **✅ THE PREVIOUS RUN IS SNAPSHOTTED** to `snapshots/20260817_pre_sigma_tightening/` (823 MB, with
> a MANIFEST holding all headline numbers). `runs/` and `diagnostics/` are gitignored and Roihu
> scratch is purged at 180 days, and `manuscript/figures/run_ids.R` auto-selects the NEWEST
> posterior — so without that snapshot the next run silently erases the comparison.
>
> **✅ REPORTING BASIS DECIDED:** balanced set (n=310), unweighted, whole profile, true obs years —
> what `obs_basis.R` implements. ⚠ Worth ~1.8×: 2006→2024 observed is **+0.117** balanced vs
> **+0.209** pairwise. Every older "+0.209" is the pairwise figure.
>
> **NEXT STEP IF THE POSTERIORS ARE STILL UNREASONABLE:** the correlated-error likelihood, written up
> in full in `manuscript/HIKET_correlated_likelihood_proposal.tex` (compound symmetry, τ fixed at the
> measured ICC, closed form, no new free parameter). ⚠ Verify `τ = 0` reproduces the current
> log-likelihood exactly before trusting any implementation. Memory: `correlated-likelihood-proposal`.
>
> **Still open:** whether to keep F3 now that the merged F4 spans the same window; and
> `appendix_sigma_input.tex` still argues the retired `[0.5, 8.7]` window and repeats the
> "physically impossible" claim that our own records list as dead.

## ⭐⭐ 0-quinquies. THE EQUIFINALITY REACHES THE FORECAST — measured 2026-08-14

**Lorenzo's insight, tested and confirmed:** high-input/short-MRT and low-input/long-MRT fit the
SOC data identically, but they do **not** project identically. Script:
`doublechecks/equifinality_forecast.R` (+ `.png`, `.rds`).

**The number:** within each model the 95% forecast spread is **25–34%** of the median 50-year loss,
and **R² of log-likelihood on the forecast is only 0.10–0.18** — so **82–90% of the projection
spread is among draws the data cannot rank**. The equifinality is not a nuisance beside the
result, it **is** the uncertainty in the result — and saturation is a *timing* question, which is
exactly what MRT sets (after 50 yr a step change is 96% realised at MRT 15 but 78% at MRT 33).

⚠ **A prediction of mine was WRONG, informatively.** I expected the *equilibrium* response to be
degenerate (`C = J/(kξ)` ⇒ warming by factor `f` gives `C/f` whatever the split). It isn't, because
**`f` itself depends on β₁**, which is what sets MRT. corr(log σ_input, equilibrium change) is
**−0.58 for Yasso07** (single ξ: MRT and climate sensitivity are the *same* parameter) but only
**−0.10 / −0.16 for Yasso15/20** (pool-specific ξ, near-degenerate as predicted). So the degeneracy
of the equilibrium response is **itself structural**, splitting on the same single-ξ vs three-ξ line
as the MRT finding.

**CLIMATE SENSITIVITY DIFFERS ~2× ACROSS THE FAMILY, and it was in NO figure.**
`doublechecks/climate_sensitivity_sweep.{R,png}` puts all three on shared axes (0–5 °C, posterior
propagated as bands). Equilibrium SOC change at **+2 °C: Yasso07 −21.7%, Yasso20 −13.7%,
Yasso15 −9.6%**; at +5 °C: −44.7 / −29.6 / −21.9. **The between-model range is 3.1× the mean
within-model 95% band at +2 °C (2.9× at +5 °C)** — for climate sensitivity, structure dominates
calibration uncertainty, and the gap **widens** with warming. This is Thread A ("structure dominates
out of sample") in the climate dimension, quantified.

🚩 **THE ACTUAL QUESTION IS NOT YET ANSWERED, and `climate_sensitivity_sweep.png` does NOT answer
it.** Lorenzo's question is **within-model, between-calibration**: *calibrate the same model under
two input assumptions ⇒ two different MRT posteriors ⇒ how much does climate sensitivity differ?*
The sweep figure shows **between-MODEL** differences with the posterior band of a **single**
calibration — a different quantity. It cannot be answered from what is on disk, because the
restrictive calibration does not exist yet.

⚠ **Slicing the existing posterior along the ridge is NOT a substitute** — a genuine re-calibration
lets every other parameter readjust, which a slice cannot capture. The only partial signal is
corr(log σ_input, equilibrium change) = **−0.58 (Yasso07)** vs **−0.10 / −0.16 (Yasso15/20)**,
suggesting the effect is large where a single ξ ties MRT to climate sensitivity and small where
three pool-specific modifiers decouple them. **A hypothesis to test, not a result.** This is now the
stated primary question of the two-calibration experiment in `HIKET_storyline_note.tex` §Outlook.

⚠ **Both scripts are DEMONSTRATIONS, not projections.** Single-exponential transient on the bulk
MRT (a multi-pool system is a sum of exponentials), a step change rather than a trajectory, and a
common J̄ for absolute scale. A proper version needs forward runs — which is precisely the
two-calibration experiment now recorded in `HIKET_storyline_note.tex` §Outlook as the **priority
future direction, ahead of the stratification test**: restrictive vs permissive σ_input, compared on
RMSE distributions, climate reactivity, sink/source sign, and 2100 trajectories. Publishable either
way. It also notes that **radiocarbon would break the degeneracy directly**, since ¹⁴C constrains
turnover independently of input magnitude — the degeneracy is a property of the observation design,
not of nature.

⚠ **Two traps, both cost time today:** the **chains store UNCONSTRAINED parameters** (the posterior
RDS is physical but carries no likelihood, which is why `intrinsic_mrt.R` reads the latter) — skip
`to_original()` and you silently get MRT = 0; and a **named scalar** taken from the sample matrix
renames the data-frame column.

## ⭐ 0-quater. THE σ_input WINDOW IS ANCHORED ON THE WRONG STATISTIC (2026-08-14)

**Decided: narrow the `flux_pair` window to `[0.05, 4.62]` — but NOT YET, and not bundled.**
Full write-up with tables, references and the boxed proposal:
`manuscript/M&M_parameterization_working_document.tex` §"Prior specification: the litter-input
flux window". Appendix figure: `manuscript/figures/SX_input_vs_literature.png` (⚠ diagnostic,
**not for production**, several comparisons on it are deliberately not like-for-like).

**The defect.** `σ_input` is a single global scalar, so `σ_input × J̄` is a **national mean**. The
current ceiling 8.7 is Gower's Class I evergreen **maximum** (9.12 gC-based) — a global boreal
single-stand maximum. **Bounding a mean with a maximum** is why the constraint has never bound
(posterior effective fluxes 4.99–6.56, all inside). The direction of the bound is fine —
`litter ≤ NPP` is a true identity — the *statistic* and the *population* are wrong.

**⚠ THE UNIT TRAP.** Gower 2001 reports **gC** m⁻² yr⁻¹ (Table 4 caption); Zheng 2004 reports
**dry matter**. Proof is internal: Zheng cites Gower's world-boreal TNPP as 109–1827 (mean 892)
where Gower's own Appendix A gives 218–912 (mean 424) gC — ratio exactly 2.00–2.10. So Zheng's
563 is **2.81**, not 5.63, tC/ha/yr. Mixing them inflates the ceiling 2×; I made this error
mid-session and it produced a "54% above NPP" claim that is wrong.

| option | σ_input ≤ | binds? | implied MRT Y07/Y15/Y20 |
|---|---|---|---|
| 8.70 current (Gower global max) | 3.46 | **no** | 15.2 / 21.9 / 17.5 (unchanged) |
| **4.62 Gower Nordic max ← DECIDED** | **1.84** | all six | 21.3 / 25.8 / 24.0 |
| 2.81 Zheng gridded mean | 1.12 | all six | 35.0 / 42.4 / 39.5 |

*(published for reference: 33.5 / 30.4 / 19.0)*

**Why Nordic max and not Zheng.** Zheng's ceiling (2.81) is **below the prior centre**
(σ_input 1.30 ⇒ flux 3.26), so the centre is not representable and re-centring would be required;
and since 2.81 is only 4% above LUKE's own total litter, it amounts to asserting the inventory as
near-exact. Nordic max keeps the project's "widest defensible" principle, fixes the statistic, and
depends on none of the arguments that failed scrutiny.

⚠ **Narrowing also tightens the PRIOR**, not just the wall: `flux_pair` is a scaled logit onto the
window, so at fixed `sigma_ppm = 0.50` the ±1 ppm range moves from σ_input [0.93, 1.72] to
[1.09, 1.47]. Decide whether that is wanted or whether `sigma_ppm` should rise to compensate.

⚠ **THREE CLAIMS THAT DID NOT SURVIVE — do not revive them.** (1) "our flux is physically
impossible" — exceeds the Nordic *mean* but not the Nordic *range*; (2) "LUKE's litter is biased
high" — its 2.70 is 84% of Gower's Nordic mean, but that range is 2.15–4.62, too wide to support
the inference; (3) "Finnish forests are 2× Gower's productivity" — compared Korhonen's *current
annual increment* against Gower's *mean annual increment* (biomass ÷ stand age, stands averaging
99 yr); the factor is largely definitional and Gower's stands are selected for having complete NPP
budgets, not for representing Finnish managed forest.

**Also established (Zenodo deposit + YaYasso):** `J̄` is a flux (biomass × turnover, Liski 2006
rates); it **includes harvest residues and natural mortality** (both previously assumed missing)
and **excludes understorey**. Component legend: `nwl` = foliage + fine roots, `fwl` = branches +
coarse roots + stem/bark, `cwl` = stumps. The Zenodo DOI **is now live** (CLAUDE.md says it isn't).

**➡ OPEN QUESTION FOR B. TUPEK — one line:** which fine-root biomass model and turnover rate were
used? YaYasso variants span foliage×0.18–2.5 for biomass and ×0.5–0.85 for turnover — a ~20× range
in the dominant `nwl` component, and exactly the term L&H exclude ("The uncertainty in these leaf
mass-to-fine root ratios was not included in our analyses"), which makes their ±10% a floor. One
variant (`fineroot.total.tsum`) is annotated as **including understorey**, which would make the
1.30 centre a partial double-count.

**SEQUENCING: wait for 654390–92 → analyse → then decide.** Do not bundle the window change with
the σ_input prior tightening (§0-ter) or anything else.

## ✅ 0. RESOLVED 2026-08-14 — IDENTIFY AND READ THE ROIHU RUN

> **Do this before anything else, and do not assume what it is.**
>
> Lorenzo believes a **corrected-target production run** was launched in an earlier session and is
> running / has landed on Roihu (as of 2026-08-13). That is **unverified** — this file recorded the
> corrected-target run as *pending* and `Data/` as *not synced*, so the possibility that it launched
> against the **OLD** target is real and would make its numbers actively misleading.
>
> **THE DECISIVE CHECK — the observation count in the launch log:**
> - **1205 observations / 456 calib-ready** ⇒ the CORRECTED target (1985 LM added, true sampling
>   years, 82 undated plot-years dropped). Good: read it as the new baseline.
> - **1269** ⇒ the OLD target. The run repeats what we already have; discard it and relaunch after
>   rsyncing `Data/`.
>
> Confirm at the same time: `Cores per chain: 40` (not 383 — the OOM trap), `5 chains x 50000`, two
> `ERROR MODEL` lines, and all six models reaching `Chain 1/5`. Guard output goes to `.err`, not
> `.out`; restrict globs to the actual job numbers or you will read August's earlier logs.
>
> **What it should show if it IS the corrected target** (pre-registered, so it cannot be
> rationalised afterwards): observations move TOWARD the models over 1985–2024
> (+0.248 → ≈ +0.182); **2006–2024 untouched at +0.107 against a source in all six**; MRT unmoved.
> If 2006–2024 unexpectedly flips sign, stop and re-plan — the priority changes.
>
> **Then:** sync `runs/`, `diagnostics/` AND `Data/model_inputs/` back, run stages 2–4 locally, and
> rebuild the figures. Only after that does the next launch (§0-ter) make sense.
>
> ⚠ Regardless of the outcome: **`Data/` on Roihu may not be current.** Before any NEW launch,
> `rsync -av "<mac-repo>/Data/" roihu:/scratch/project_2019134/HIKET/Data/` —
> `Data/model_inputs/site_attributes.csv` is a **NEW** file and a partial rsync would leave the SOC
> builder without it. ✅ Code is pushed through `391fc4a`.

## ⭐ 0-ter. THE LAUNCH AFTER THAT — tighten the `sigma_input` prior

**One factor, on top of the corrected target.** Do not bundle it with anything else; the discipline
that produced every clear result so far is one change per run (the σ=0.80+t run deliberately
excluded campaign-σ so it could not confound the trend).

`sigma_input` log-SD **0.50 → 0.15–0.20**, justified from **Lehtonen & Heikkinen 2015** — ⚠ *not*
chosen because it lands MRT near the published value. That is the circularity trap, and the 2.4×
inflation is not itself a finding.

**Why it should work, and why that is newly established:** `CLAUDE.md` used to claim the posterior
does not ride the MRT×σ_input ridge — which would predict this lever does nothing. That claim was
**simple-models-only** and is false for Yasso (r = −0.79/−0.73/−0.57). Pinning σ_input drags MRT
along that ridge; the ridge is the mechanism.

**PRE-REGISTERED expected fit cost**, from F14's profile lower bounds (run `20260812_0809*`):

| model | published MRT | expected ll cost to reach it |
|---|---|---|
| Yasso20 | 19.0 | **≈1.6** — nearly free, should move easily |
| Yasso15 | 30.4 | **≈14.7** — visible but survivable penalty |
| Yasso07 | 33.5 | **posterior never reaches it** — expect it to strain hardest, or refuse |

Costs far BELOW these ⇒ our max-over-draws bounds were loose. Far ABOVE ⇒ something outside the
ridge resists. **Watch R² and the log-likelihood, not MRT** — MRT moves by construction if the ridge
is real; what is informative is what it costs.

**Neither open item blocks this launch:** the published-MRT comparator (point vs posterior median)
is a reporting decision, and the Yasso07 profile likelihood is separate Puhti work.

## ⭐ 0-bis. SESSION OF 2026-08-13 — the MRT/fitness work (read with §0)

Nothing here touches the data target or the runs; it is analysis of the **existing** run
`20260812_0809*`. All of it is committed as scripts + figures.

### What landed

1. **F14 is IN the manuscript** (Lorenzo's decision). `manuscript/figures/build_F14_mrt_ridge.R`
   → `F14_mrt_ridge.png`, written into `HIKET_storyline_note.tex` §"What a longer residence time
   would cost", with **two boxes**: a `methodnote` explaining the two layers and a `provisional`
   box recording the planned rework. **F13 was DELETED** as not useful — do not resurrect it; its
   variance panel survives as F14's lower row.
2. **F14 now uses EVERY draw** (~75k, was 6000; one MRT eval is 1e-4 s). Not cosmetic: the
   max-per-cell bias shrinks with draws *in that cell*, so subsampling made the tails look worse
   than they are — biased the convenient way. `NS <- Inf`, `NB <- 60`, `MINN <- 5`.
3. **The ridge claim was mis-scoped and is now corrected** — `ridge_test.R` §1 covers SP1/TP2/TP3
   only. Yasso rides a real ridge (r −0.79/−0.73/−0.57) while the simple models do not
   (−0.26/−0.31); not a run effect. Fixed in CLAUDE.md, memory, and `ridge_test.R`'s own header.
   🚩 **AND THAT CORRECTION IS ITSELF WRONG** — S13 overturns it; see memory `s13-benchmark-ridge`.
4. **NEW STRUCTURAL FINDING** — `doublechecks/xi_published_vs_ours.R`. Yasso07 has **one** ξ for
   all pools ⇒ `MRT = MRT_ref/ξ` exactly; ξ 0.855→1.822 predicts **15.7 yr vs actual 15.2**, the
   whole gap. Yasso15/20 have **three** pool-specific ξ and their humus modifier moved only
   ×1.04/×1.06. So the family's MRT gaps are **not one phenomenon**. Calibration also **inverts
   the ordering**: published 33.5/30.4/19.0 → ours 15.2/22.1/17.5.

### ⚠ Two things to resolve before quoting any of it

- **TWO "published MRT" bases exist and older records mix them.** `intrinsic_mrt.R` reports the
  published **POINT** (`to_original(best_x)`) *and* the published **POSTERIOR** median (FMI `.dat`,
  Yasso15/20 only); MRT is nonlinear many-to-one so they differ. Yasso20: **19.03** vs **25.02**,
  and the fit cost of reaching it goes **1.4 → 9.8**. F14 and the new text use the POINT.
  **Decide which the paper uses.**
- **The costs are max-over-draws lower bounds**, weakest where draws are thin. Corrected values:
  **Yasso07 never reaches its published MRT at all** (0 draws; 10.4–25.4 yr vs published 33.5, so no
  estimate exists), Yasso15 **≤14.7** (741 in band), Yasso20 **≤1.6** (121 141).
  ⚠ **The first version was built on BURN-IN ARTEFACTS** — `getSample()` returns the first retained
  iteration of each internal DEzs chain (15 rows/model, 100–200 ll below the bulk), and those were the
  only draws near the published value for Yasso07/15. Fixed by per-sampler extraction with
  `start = 2`, which also removes the 1-in-3 thinning ⇒ **225 015 draws**, not 75 015.

### The planned rework (Lorenzo: on Puhti, later)

True **profile likelihood** for the colour map: optimise the remaining ~18 parameters at each grid
node. Removes the prior from *coverage* as well as *value*, prices the region past the posterior's
reach, and is smooth enough to interpolate coarsely. Grid over **(β₁, σ_input)** not (MRT, σ_input)
— both are box constraints, whereas fixing MRT is a curved manifold, and β₁ carries the MRT
variation. Warm-start each node from its neighbour. ~7 h/model at 20×20; cost scales as nb².
**Yasso07 first.** Profile the *likelihood*, not the posterior, and show the prior separately —
otherwise a region the data reject cannot be told from one the prior merely disfavours.
Expect the figure to change again once the σ_input prior is tightened (§6).

### Still open from this session

- Is ξ = 1.82 defensible? Needs `beta1` 0.0987→**0.1578** (+60%), ~4.4 prior σ off centre. It is
  essentially Yasso20's *published* β1 (0.1580) — suggestive, but Yasso07's β1 scales all pools and
  Yasso20's scales AWE only, so not strictly comparable. Pair with the FMI warming-rate check.
- Yasso20 is the least ridge-like Yasso (ratio 0.92) yet the most accommodating on MRT. Unexplained.
- Decide whether F7 and S5 are retired (superseded by S11/S12), and whether S12 is promoted to
  main text.

## 0a. The Roihu run in flight (superseded target — diagnostic value only)

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

Added 2026-08-12, because these bear on the data work in §0c:

5. **How much does the Student-t already down-weight 1985?** ν=6 de-weights outliers, and 1985 is
   where the heavy residuals are. Measure the campaign's effective weight under t vs Gaussian — the
   planned `SIGMA_1985_INFL = 2` lands *on top* of this, and we need to know what it adds.
6. **Per-campaign residual spreads under the t.** The 0.97 / 0.62 / 0.54 were measured under a
   Gaussian. If they compress, the empirical case that 1985 is special weakens and the appendix
   must rest entirely on the mechanism argument.
7. **Does the 2006→2024 source survive?** If yes, "it's an error-model artefact" is closed for the
   one window no data fix can reach.

⚠ Log-likelihoods are **not** comparable to any previous run — both family and scale changed.
⚠ Do **not** refresh downstream figures from this run: §0c supersedes its target. Read it as a
diagnostic, then do the data fixes, then one re-run.

---

## 0c. THE IMPLEMENTATION — ✅ DONE 2026-08-12. Only the production re-run remains.

**Commits, each a clean revert point**

| | |
|---|---|
| `3a143cd` | analysis only — target untouched. **Revert here to undo everything.** |
| `8ef5c62` | treatment C + true sampling years + `SIGMA_1985_INFL = 2` |
| `b4bee98` | appendix figures + basis table on the corrected target |
| `247eab4` / `0a25ca4` | the build↔Data_work loop: recorded, then scoped |
| `44215b0` | **loop broken; fixed-point test passes** |

**Converged target (this is what is on disk)**

| | |
|---|---|
| plot-years | **1408** (fixed point; `Data_work → build → Data_work → build` is byte-identical, 0 of 53 numeric columns differing) |
| weighted profile means | VMI8 **66.309**, Biosoil **69.971**, Komeetta **72.444** |
| official validation | 2006 = 59055, 2024 = 61047 — still exact |
| `Data_work` | all checks pass, **456** calib-ready, 1205 observations, median target 64.9 |
| observed rates (balanced, unweighted, whole profile, true intervals 17.1/18/35.1 yr) | 1985→2006 **+0.368**, 2006→2024 **+0.143**, 1985→2024 **+0.254** |

**What was done**
1. **Treatment C** — each plot's own 2006 LM added to its 1985 organic layer (478 of 488 rows; 10
   skipped where OFH is zero/absent so `organic_zero` keeps firing). A/B verified: 1985 moves,
   2006/2024 bit-identical. Switch: `HIKET_ADD_1985_LM=0`.
2. **True sampling years** — `samp_year` → `obs_year` → `soc_obs_year`, indexed by all six
   calibration scripts. `year` stays the campaign key. 82 undated 1985 plot-years dropped
   (71 South / 9 North — stated, not compensated).
3. **`SIGMA_1985_INFL` 1.0 → 2.0**, pre-registered with its rationale in `calib_config.R`.
4. **The build↔Data_work loop broken** — see below.
5. `Data_work.R` and `build_soc_homogenized.R` brought under version control (they were untracked,
   collateral of the blanket `Data/*` rule).

### ⚠ THE LOOP — fixed 2026-08-12, and the invariant that must hold

`build_soc_homogenized.R` read `site_raw.csv`, which is built from `plot_data`, whose row set is
`merge(avg_inputs, avg_SOC)` on `common_plots` — an **inner join on the SOC data the builder
itself produces**. Worse, the GTK soil-class extraction (hence the per-class λ, hence the deep tail
for all three campaigns) ran only over `plot_data`'s plots, so `soil_code` was SOC-gated at source.

It did **not** reach a fixed point in one pass: 1411 vs 1408 plot-years, 30 differing
`soc_profile`, 38 differing λ, means moving ~0.05 Mg/ha. Numerically small; the real defect was
that the target could not be regenerated from a clean checkout.

**Fix:** `Data_work.R` extracts soil types over the SOC-independent universe (`plots_sf`, from
litter-input coordinates) and writes **`Data/model_inputs/site_attributes.csv`** (2719 plots)
before anything SOC-derived. The builder takes every target-affecting lookup from there and keeps
`site_raw.csv` only for descriptive covariates.

> **INVARIANT — do not break it.** Nothing SOC-dependent may be added to `site_attributes.csv`,
> and the **ROW SET** matters as much as the column list: gating those rows on SOC restores the
> loop invisibly. If you ever add a column there, ask first whether it is a raw attribute or a
> modelling result.

### 🔄 RUNNING NOW: local Yasso15 f-sweep

`doublechecks/ablation_logs/Yasso15_SUITE_20260812.log` — 4 configs (`A0`=f 2, `A1`=1, `A2`=1.5,
`A3`=3), 3 chains × 6000, ~70–90 min each, launched on the **converged** target.

**When it lands:**
1. `Rscript doublechecks/quarantine_ablation_runs.R` — ⚠ **mandatory**, or `run_ids.R` picks a
   3×6000 short chain as production and every figure silently rebuilds from it.
2. `Rscript doublechecks/summarise_ablation.R` to compare.
3. **`A1_C5_off` is the attribution arm** — corrected data at f = 1, so A1→A0 isolates what f = 2
   does and the remainder is the data fixes.
4. ⚠ Short chains compare *locations*, not publication posteriors. Do not quote them as if they were.

### ⬜ STILL OUTSTANDING

1. **The production re-run**: six models, Roihu, corrected data, f = 2. Nothing else blocks it.
2. **Then rebuild F2 / F3 / F4 and `appendix_delta_reconciliation`** — deliberately NOT rebuilt
   today, because their posterior bundles are from the old target and would compare old models
   against new observations. S8/S9/S10 *were* rebuilt (observation-only).
3. **Hannu's answer on the LM question.** ⚠ If 1985 already includes LM, treatment C **inverts**;
   revert `8ef5c62` or set `HIKET_ADD_1985_LM=0` and rebuild.
4. Re-run `doublechecks/observed_soc_basis.R` after the production run and refresh the two basis
   annotations in `HIKET_next_session.tex` and the M&M document, which quote current values.
5. ✅ **Denominators are DONE** — `observed_soc_basis.R`, `build_S9_soc_change_by_depth.R` and
   `soc_depth_distribution.R` all compute the true mean interval from `samp_year` (35.1 yr, not 39).
   The observed 1985→2024 rate on the corrected target now reads **+0.18**, against +0.248 at the
   start of 2026-08-12.

### 🔮 WHAT WE EXPECT FROM THIS RUN (written before it lands — 2026-08-13)

**This run is a BASELINE on a corrected, reproducible target. It is not an attempt on the MRT.**

| quantity | expectation | basis |
|---|---|---|
| `sigma_init` | **falls hard**, ~1.0 → ~0.15 | local sweep; A1 says most of it is the DATA fixes, not f |
| pre-run inversion | **resolved** — all sweep values ≪ 0.818 | the defect open since 2026-08-07 |
| `sigma_input` | **unchanged**, ~2.6–2.7 | flat across f ∈ {1, 1.5, 2, 3} in the sweep |
| **MRT** | **unchanged**, ~15 / 22 / 18 | σ 0.72→0.80 + Student-t moved it <2%; f doesn't touch σ_input |
| 1985 level | +3.2 Mg/ha, but 1985 now carries ¼ weight | treatment C + f = 2 partly cancel |
| 2006→2024 | still a source in most models | survived the error-model change already |
| residual sd / kurtosis | ~0.74 / ~7.2 | insensitive to what we assume σ is |

⚠ **If MRT moves materially, something we believe is wrong** — investigate before celebrating.

### ➡ THEN: the σ_input prior, and the question underneath it

If MRT does not move (expected), the remaining identified lever is the `sigma_input` prior width,
log-SD **0.50 → 0.20** (Lehtonen & Heikkinen 2015), predicted to bring Yasso15 to ≈34 yr.
**But that is not a tuning step — it is a discriminating experiment**, and the discussion has to
happen first. The question it answers:

> Is the fast-MRT / high-input solution the only place the model can go given which parameter has
> the loosest prior (**H1**), or is it a real structural inadequacy being hidden in σ_input (**H2**)?

- **H1 — the split is set by the PRIORS, not the data.** The data pin MRT × σ_input (constant to 4%)
  but not the split. Yasso's rates are FIXED, so MRT can only move through the transfer fractions,
  whose priors are tight (logit SD 0.4); σ_input's prior is the loosest thing in the system
  (log-SD 0.50, centre 1), so it absorbs the residual. **Evidence FOR:** the ridge test —
  corr(log MRT, log σ_input) is only −0.30 to −0.37 and sd(log product)/sd(log MRT) ≈ 1. If the
  data constrained only the product, the posterior would be a thin diagonal with corr ≈ −1.
  It is not: each is pinned *separately*, i.e. by its own prior.
  **Prediction under H1:** tightening σ_input forces the fractions to move instead → MRT rises,
  fit roughly holds.
- **H2 — structural inadequacy displaced onto σ_input.** The models cannot reproduce the observed
  pattern, so they manufacture carbon. **Prediction under H2:** MRT rises little and R² visibly
  degrades — which *exposes* the inadequacy instead of hiding it, and is a result in its own right.

**The experiment discriminates them: watch R², not just MRT.**

⚠ **A caution on the prior's CENTRE, not just its width.** σ_input conflates tree-litter model error
with **missing understorey** (the Tupek product excludes it), so part of σ_input > 1 is real input,
not error. Tightening to log-SD 0.20 *around a centre of 1* would assert the understorey is
negligible. Decide the centre explicitly — the effective flux 2.6 × 2.4 ≈ 6.2 tC/ha/yr is still
under the ~9 boreal NPP ceiling, i.e. not physically absurd.

### 📋 TOMORROW, in order

1. **Read the Roihu run 597028–597033** against the §0b checklist. ⚠ Old target — diagnostic only.
2. **Check the local sweep** (`doublechecks/ablation_logs/Yasso15_SUITE_20260812.log`), then
   `quarantine_ablation_runs.R`, then `summarise_ablation.R`. `A1_C5_off` is the attribution arm.
3. **Push + rsync** (see the box at the top).
4. **Launch the production re-run**: six models, corrected data, f = 2.
5. Only then rebuild F2 / F3 / F4 and `appendix_delta_reconciliation`.

### ⚠ OPEN DEFECT FOUND 2026-08-12: build_soc_homogenized.R and Data_work.R form a LOOP

`build_soc_homogenized.R` reads `Data/model_inputs/site_raw.csv`, which `Data_work.R` **writes**;
`Data_work.R` reads the SOC CSVs that `build_soc_homogenized.R` **writes**. `site_raw` supplies the
peat exclusion, `soil_code` (hence the per-GTK-class λ) and region/weights, so the loop closes:
SOC data → `calib_ready`/peat → `site_raw` → plot set and λ → deep tail → SOC data.

**One pass does not reach a fixed point.** Measured: a second build pass (on the `site_raw.csv`
that the first `Data_work` run produced) gives **1408 plot-years instead of 1411**, 30 differing
`soc_profile` values and 38 differing λ, moving the campaign means by ~0.04–0.06 Mg/ha:

| | pass 1 | pass 2 |
|---|---|---|
| VMI8 | 66.2689 | 66.3087 |
| Biosoil | 69.9120 | 69.9712 |
| Komeetta | 72.4195 | 72.4438 |

This is also the origin of the "input drift" seen against the 3 August build: the SOC swap of
4 August re-ran `Data_work` and changed `site_raw`, but `build_soc_homogenized.R` was never re-run,
so the live target had been carrying a stale plot set ever since.

**Current state: pass 1 is on disk**, and it is what the commits, the input bundles and the running
Yasso15 sweep are all built on — deliberately, so everything is mutually consistent. Do not re-run
either script piecemeal.

**THE FIX IS SCOPED (2026-08-12).** The builder needs only these from `site_raw.csv`:
`plot_id`, `peatland`, `soil_code`, `region`, `x_ETRS`, `y_ETRS`, `lon_WGS84`, `lat_WGS84`, plus
the era-tagged stand covariates (`basal_area_85`, `stand_age_85`, `mean_height_85_dm`,
`dev_class_85`, and the 2024 descriptors). **Every one is a raw NFI/GTK attribute — none is a
modelling result.** The contamination is incidental: those attributes arrive in a file that also
carries `calib_ready` / `soc_outlier` / `n_soc_obs` and, critically, a ROW SET shaped by the SOC
data (`site_raw$plot_id = plot_data$plot_id`, and `plot_data` is joined against `obs_count`, which
is aggregated from `SOC_agg`).

**Preferred route:** have `Data_work.R` write a second, minimal `Data/model_inputs/site_attributes.csv`
containing only the raw columns above, emitted from a point in the script that is provably upstream
of anything SOC-derived; point `build_soc_homogenized.R` at that instead. One new write, one changed
path. ⚠ The verification that matters is the ROW SET, not the columns: the new file must be built
before `SOC_agg` is touched, or the loop survives in a subtler form. The weaker alternative — having
the builder re-read the original GTK/NFI sources — duplicates parsing logic and will drift.

**To resolve (before the production re-run):**
1. Establish whether this is a 2-cycle or converges — needs one more `Data_work` → build iteration,
   which must NOT be done while the sweep is running (it would change `input_raw_monthly.csv`
   underneath the later configs).
2. Then either iterate to a fixed point and record the number of passes, or break the loop — e.g.
   have `build_soc_homogenized.R` read the *raw* site table rather than the `Data_work` product,
   so the dependency runs one way only. The second is the real fix.
3. The effect is small (~0.07% on the means) but it makes the target **non-reproducible from a
   clean checkout**, which matters more than the magnitude.

### Original plan, retained as the record of what was decided

Full findings: `manuscript/HIKET_discussion_memo.tex`, memories
[[soc-campaign-comparability]] and [[kramarenko-thesis-vmi8-biosoil]], scripts
`doublechecks/{soc_depth_distribution,litter_in_1985_organic,organic_mineral_boundary,subsoil_offset_mechanism,litter_layer_vs_input}.R`.

### Why (one paragraph)

The 1985 target is not comparable with 2006/2024. Its organic layer is **OFH only** — the LM
(litter+moss) layer, a separately coded layer in the protocol, was measured in 2006/2024 and is
absent from 1985. And "1985" is really **1986–1995**: 19.5% of plots were sampled in 1995, so the
campaign mean represents ~**1989** and every rate denominator is 12% too large. Both are fixable.
What is left after them — LOI-derived mineral C%, relocated subplots, and a physically implausible
depth-inverted subsoil gain — is not, and gets an error-model term.

### The changes, in order

**1. `Data/SOC_homogeneized/build_soc_homogenized.R`** — back up the current CSVs first.
- **Treatment C:** add each plot's own **2006 LM** (BiSo col 247) to its 1985 organic layer, at the
  **layer level** so it propagates to `soc_0_40` and `soc_profile` consistently.
  ⚠ **Compute the organic-quality flags (`organic_missing`/`organic_zero`) BEFORE adding LM** —
  10 VMI8 plots have organic C = 0 and would silently become calibration-ready otherwise.
  ⚠ 113 of 458 VMI8 plots have no 2006 LM, but **77 of those are also the no-year plots dropped in
  step 2**, so only **36** need a fallback. Median LM = **3.187** Mg C/ha.
- **Sampling year:** read BiSo **col 140 (`MAANAYTE`)**, expose as `samp_year`. Two sources agree
  perfectly (BiSo col 140 and sheet `kiv_maat85_95`, 387 overlapping, 0 disagreements; union 391).
  Keep `year` = campaign key (1985/2006/2024); `samp_year` is the true year.
- Re-run the build.

**2. `Data/Data_work.R`**
- Carry the true year through as **`obs_year`** (= `samp_year` for VMI8, = `year` otherwise).
- **Drop the VMI8 observation** for the 84 plots with no recovered year. Only that observation —
  the plots keep 2006/2024. `calib_ready` has **no minimum-obs condition**, so 2-observation plots
  are fine. ⚠ The dropped set is **71 South / 9 North**, so it shifts the first campaign's regional
  balance under the 1:3 weights — quantify and state it.

**3. `Calibration_real_data_transient/calib_config.R`**
- `SIGMA_1985_INFL` **1.0 → 2.0**, with the rationale beside the existing C5-withdrawal note.
  **Pre-registered**: the value is set BEFORE the corrected data are run, and must not be revisited
  on the basis of results — that is the whole defence. Wording agreed:
  > *f = 2 was set a priori as a judgement reflecting documented but unquantified differences in the
  > first campaign's protocol — carbon concentrations largely derived from loss-on-ignition,
  > subplots relocated between campaigns, and a subsoil gain with no physical mechanism. It is not
  > estimated from the data. Its influence was examined with short-chain runs at
  > f ∈ {1, 1.5, 3} for Yasso15.*

  (Wording revised 2026-08-12: the sweep is short-chain and single-model, so it cannot be
  described as if all four were production posteriors.)

**4. The six `run_*_transient_calibration.R`**
- `idx = match(obs_plot$year, clim$year)` → **`obs_plot$obs_year`**.
- Leave `year` as the campaign key, so `sigma_infl = ifelse(obs_plot$year == 1985L, …)` and
  `HIKET_DROP_CAMPAIGN` keep working. That is D3 achieved **without renaming anything**.
  ⚠ If `year` were changed to the true year instead, both would break **silently**.

**5. Reporting only (no re-run needed)** — denominator 39 → **34.7**:
`manuscript/figures/build_S9_soc_change_by_depth.R`, `doublechecks/soc_depth_distribution.R`.
(`build_appendix_delta_reconciliation.R` deliberately keeps 39 — it must match the models' x-axis.)
Corner-cut validated: Δ-of-means/34.7 equals the mean of per-plot rates to **0.1%** (interval length
vs change r = 0.036, p = 0.51).

### Run plan (decided 2026-08-12 — deliberately cheap)

**Production: six models, corrected data, f = 2.** That is all. No full sweep — too expensive.

**Plus one short-chain sweep on Yasso15 only**, which is already coded:
```
Rscript doublechecks/run_ablation.R Yasso15 6000 3 A0_reference,A1_C5_off,A2_C5_1.5,A3_C5_3.0
```
`A0` = f 2 (the new default), `A1` = 1.0, `A2` = 1.5, `A3` = 3.0. Optionally add
`A5_LOCO_no1985`, which withholds 1985 entirely.

⚠ **`A1_C5_off` is not just a sensitivity point — it is the ATTRIBUTION arm.** Production changes
the data *and* the σ at once; A1 has the corrected data at f = 1, so the A1→A0 difference isolates
what f = 2 does and the rest is the data. Without it the two changes are confounded.

⚠ **Run `doublechecks/quarantine_ablation_runs.R` afterwards**, or `run_ids.R` will select an
ablation posterior as production and every figure will silently rebuild from short chains.

⚠ The ablations source the real calibration scripts, so they must run **after** the data fixes.

**RUN THE SWEEP LOCALLY — it does not need Roihu** (checked 2026-08-12). The Yasso20 six-config
suite ran on this Mac at **~70 min/config** (2026-08-05, 00:33→07:41); Yasso07 took ~3 h. Yasso15
shares Yasso20's Fortran, so budget **70–90 min × 4 configs ≈ 5–6 h**, i.e. one overnight, on
12 cores. Prerequisites: the data fixes re-run first (the ablation calls the real calibration
script, which rebuilds its own input bundle from `Data_work.R` output, so it picks up the corrected
data automatically) and `SIGMA_1985_INFL = 2.0` already set, so `A0_reference` *is* the f = 2 arm.

⇒ **Schedule gain:** the local sweep and the Roihu production run are independent and can run
simultaneously from the same corrected data. The attribution arm `A1_C5_off` may therefore be in
hand *before* the production results land — useful, since it is what separates "the data moved it"
from "f = 2 moved it".

### Settled details (agreed 2026-08-12, do not re-open)

- LM fallback for the 36 remaining plots: **global median 3.187** Mg C/ha.
- The 71 South / 9 North imbalance from dropping the no-year plots: **state it, do not compensate.**
- Treatment C applied at the **layer level**, so `soc_obs_tCha_sum` also changes — add an explicit
  flag column marking which 1985 observations carry an imputed LM component.
- Hannu has been asked about the LM question (2026-08-12). ⚠ If he answers that 1985 *includes*
  LM, treatment C **inverts** and this plan must be redone before running.
- f = 2 is the user's decision, to be discussed with coauthors but not blocking the launch.

### ⚠ 0c-bis. THE REPORTING BASIS — settled 2026-08-12. Independent of the data fix.

`doublechecks/observed_soc_basis.R` is now the **canonical table**. Run it before quoting any
observed level or trend. There are **three** axes, not one, and together they span a **4× range**
for the same data.

| 2006→2024 rate | unweighted | weighted |
|---|---|---|
| whole profile, paired | **+0.216** | +0.161 |
| whole profile, balanced | +0.143 | +0.101 |
| measured 0–40, paired | +0.144 | **+0.104** |
| measured 0–40, balanced | +0.076 | +0.049 |

- **weighting** — unweighted = the average *plot in our sample*; weighted = the average *hectare of
  Finland* (North sampled at ⅓ density ⇒ weight 3). Worth 25–30%.
- **plot set** — for 2006→2024, requiring a 1985 observation too drops the rate **+0.216 → +0.143**
  (34%). Plots lacking a 1985 measurement gained more. For 1985→2024 it barely matters (0.316 vs
  0.312). ⚠ Never difference means over *different* subsets — the script prints those rows labelled
  `UNPAIRED (invalid)` so they are recognisable in the wild.
- **depth** — `soc_0_40` (measured; LUKE official and the Hannu–Juha paper) vs `soc_profile`
  (+ modelled deep tail; HIKET's target). Worth ~35% on this interval, because the tail is refitted
  per campaign and so carries its own trend.

**THE RULE.** Model comparison → **paired, unweighted, whole profile**. Anything national →
**weighted, measured 0–40**. Comparing intervals with each other → **balanced** (constant plot set).
State the depth basis whenever a level is quoted.

**Verdict on existing text:** the manuscript's `+0.209` / `+0.312` are internally consistent and are
the *correct* basis for model comparison — nothing to retract. The only hazard is that
`M&M_parameterization_working_document` also quotes LUKE's *weighted* national stocks; those two
must never share a sentence with a model–observation gap.

#### Already fixed (2026-08-12)
- `build_appendix_delta_reconciliation.R`: observed rate was hardcoded from the **retired**
  63/102/105 series (1.077 vs the true 0.312 — **3.5× too high**), and the axis limits were stale
  too (`ylim=c(84,112)` against an actual 63–81, i.e. drawing off-scale). Both now computed.
- **Same figure: model and observed were on different plot sets** (512 vs 316). Aligning them moves
  Yasso15 from +0.422 to **+0.379** and its 1985 level 63.3 → 65.1, i.e. the mismatch was inflating
  the apparent gap by ~40%. Both sides now restricted to the paired plots. Model rates on the
  aligned set: SP1 −0.092, TP2 +0.148, TP3 +0.143, Yasso07 +0.259, Yasso15 **+0.379**, Yasso20
  +0.204, against observed **+0.312** — Yasso15 now *exceeds* the observation.

#### Also done 2026-08-12 (was "still to do", now closed)
1. ✅ **F2 / F3 / F4 aligned.** New shared `manuscript/figures/obs_basis.R` gives all three one
   **balanced** plot set (n = 316) for both the model curves and the observed markers. The F4 cache
   key gained a `_bal_` marker so the pre-alignment cache cannot be silently reused. F4 regenerated
   (spin-up recomputed), then F3 and F2. **F5 was already correct** (paired per-plot differences).
2. ✅ **`HIKET_next_session.tex`** — this work *is* the action H6 called for ("fix one plot set and
   one statistic, applied identically to observed and modelled rates"). H6 is now marked RESOLVED
   there, with the three axes, the rule, and the instruction to use **+0.143** rather than +0.209
   where 2006→2024 sits beside 1985→2024. Existing tables left intact — they are a correct record of
   the paired basis, and overwriting them would have desynchronised their other columns.
3. ✅ **M&M working document** — the trend table now carries a basis footnote giving all four values
   (+0.209 paired-unweighted, +0.143 balanced, +0.161 weighted, +0.111 LUKE official 0–40) and the
   warning not to put a weighted national stock in the same sentence as a model–observation gap.

#### Still to do
4. Re-run `observed_soc_basis.R` after the data fix; every number in it will move. Then re-check
   the two annotations above, which quote current values.

### Numbers to expect

| | now | after |
|---|---|---|
| 1985 profile mean | 61.8 | ~65.1 |
| observed 1985→2024 | +0.248 | ~+0.182 |
| observed 2006→2024 | +0.107 | **unchanged** |

The litter fix pushes the rate down, the year fix pushes it up, litter wins. Observations move
**toward** the models. ⚠ Do not read the net movement as either fix "working".

### What this does NOT fix

**2006→2024 stays at +0.107 observed against a source in all six models**, and the MRT tension
(15.2 / 21.9 / 17.4 vs published 33.4 / 30.3 / 25.0) is untouched. This narrows the problem and
removes the "maybe it's the data" escape route; it does not solve it.

### Decisions already taken, do not re-litigate

- **No subsoil reconstruction.** The 1985 20–40 cm value is *measured*; replacing it would overwrite
  data with a model, and the choice among the three variants is worth 0.058 Mg/ha/yr (~23% of the
  signal). The table (measured 61.75 / shape 63.14 / static 64.10 / regress 65.40) goes in the
  appendix as a stated sensitivity. Litter is different because it is *missing*, not measured.
- **No campaign bias term.** The trend *is* the campaign difference; a free `c_j` eats it
  (δ-offset test).
- **Report the f sweep every time**, never a single number: the C5 ablation moved the 1985–2024
  change from +0.638 to −0.184 as the weighting weakened. It flips the sign of the sink.

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
