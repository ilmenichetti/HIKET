# NEXT SESSION — pending items

**Written 2026-08-04 evening.** Live checklist for the run-up to the Roihu recalibration.
Delete an item only when it is actually done, not when it is merely started.

---

## 0. First thing tomorrow — collect the ablation batch

Running overnight: 4 models (TP2, SP1, Yasso07, Yasso20) x 6 configs, plus a TP2 catch-up
for `A5_LOCO_no1985`. Logs in `doublechecks/ablation_logs/`, progress in the scratchpad
`queue.log`.

```bash
Rscript doublechecks/summarise_ablation.R      TP2   # posterior locations
Rscript doublechecks/ablation_fit_by_campaign.R TP2  # UNWEIGHTED per-campaign fit
Rscript doublechecks/quarantine_ablation_runs.R      # BEFORE any predictive/rsync step
```

**Read it in this order and do not skip the caveat:**

- **Log-likelihoods are NOT comparable across C5 settings.** C5 changes the observation SD and
  therefore the density's normalising constant. Judge only on the unweighted per-campaign
  metrics from `ablation_fit_by_campaign.R`.
- **Noise floor:** 3 chains x 6000 iterations. Differences under ~2% in trusted-campaign
  (2006+2024) RMSE are not differences.

### The C5 decision (rule agreed BEFORE seeing results, 2026-08-04)

User's lean: **drop C5**. Keep it only if either:

1. **`A5_LOCO_no1985`** (fit on 2006+2024, 1985 withheld) predicts 1985 **above** observed by
   ~+8-10 tC/ha or more (>=15% of the 59.4 median) -> VMI8 really does read low; or
2. **`A1_C5_off`** fits 2006+2024 **materially worse** than `A0_reference` -> 1985 is genuinely
   constraining.

Otherwise remove C5 and state in the paper that the 1985 correction was made **in the data**,
not in the likelihood. Rationale: the SOC homogenization already flattened the 1985->2006 rise
from +61% to +12%, so C5 — diagnosed against the OLD data — risks double-correcting.

### Also read off the same batch

- **`sigma_init` vs 0.826.** Above that threshold the reconstructed 1917->1985 pre-run
  *declines*, contradicting the growing-stock history C3 encodes. Check before spending the
  Roihu allocation.
- **`sigma_input`**: expected ~1.0-1.2. If it lands **below 1.0**, the C4b story needs
  rethinking (prior says "missing understorey litter", data would say "too much carbon").

---

## 1. Then: commit (nothing is committed yet)

Deliberate — the tree was held so the data would be right before Roihu.

Tracked and pending: `CLAUDE.md`, `manuscript/HIKET_main_manuscript.{tex,pdf}`,
`calib_config.R`, `calibration_engine_transient.R`, all six `run_*_transient_calibration.R`,
`preinit_input_shape.R`, `Prior_specs/TP3_priors.R`, and five new `doublechecks/` scripts.
Delete the stray `Rplots.pdf`.

**`Data/` is gitignored** — `Data_work.R`, the SOC baseline and the regenerated
`model_inputs/` do **not** travel with git. Backup of the pre-swap script:
`Data/Data_work_pre_SOC_swap_20260804.R`.

---

## 2. Roihu recalibration — pre-flight

- [ ] `git pull` on Roihu **and rsync `Data/`** (use the `roihu:` alias, not the raw hostname).
      The stale-input guard will abort the job if this is forgotten — that is what it is for.
- [ ] Recompile Fortran: `R CMD SHLIB yasso07.f90` and `yasso15.f90` in **separate** calls.
      NB the documented "~69 tC/ha" sanity value predates the new data — re-derive it locally
      first, or the check is ambiguous.
- [ ] Confirm `Cores per chain: 40` in the logs (not 383 — the singularity OOM trap).
- [ ] Watch the two new startup lines: `Input currency check:` and `C5 sigma_1985 inflation:`.
- [ ] Runtime is ~16% longer than the last run (447 -> 520 plots); 36 h should still fit,
      partition `small` allows 72 h if more headroom is wanted.

---

## 3. Documentation still to update (after the recalibration)

- [ ] **`HIKET_calibration.Rmd`** — NOT yet updated for the new SOC target or the litter
      reconstruction. The largest remaining doc gap.
- [ ] **`manuscript/REVISION_PLAN.md`** — C3's status (second-order, quantified) and C5's
      status (in question) are not reflected there.
- [ ] **`manuscript/HIKET_storyline_note.tex`** — currently carries a "levels are provisional"
      reading note. Drop it once figures are rebuilt.
- [ ] **All figures F1-F14 + S1-S7, tables T1/T2, the NextGenC report bundle** — every one
      predates the rebaseline.
- [ ] **Restate every number quoted against a campaign mean**, including the "+26 tC/ha above
      observed 1985" claim. The trajectory changed from +61% to +12%; the paper's core effect
      is real but ~a fifth of what the old figures implied. State it honestly at its true size.

---

## 3b. Input QC sweep (`doublechecks/input_data_qc.R`, run 2026-08-04)

Written after the 1985 artefact, on the principle that it survived every existing check
because none was looking for it. Deliberately probes series *ends*, year-on-year steps,
per-plot constancy and component structure rather than summary statistics.

**Clean:** no negative/NA/non-finite litter, no duplicate plot-year-month rows, complete
coverage (555 plots x 480 rows), all 520 calib plots present, climate ranges sensible
(T -27.6 to 22.7, precip >= 0), SOC obs all in June, no component zero everywhere, and
**no anomalous year-steps** — the 1985/1986 ratio is now 1.05, confirming the backcast.

**Finding 1 — 37 plots with near-zero litter 1985-1995: REAL, not an artefact.**
36 of 37 recover >10x (median 52x), running 0.02-0.18 -> 2.8-5.7 tC/ha/yr. This is the
signature of stands clearcut or very young at the start of the record. *I initially flagged
this as a serious problem on the grounds that `J_t0_mean` (the pre-run endpoint) collapses for
them — that reasoning was WRONG*: it applied a steady-state argument to a transient pre-run,
where humus retains carbon from the earlier, higher-litter years. Tested directly at TP2 prior
centre: these plots are fitted slightly BETTER at 1985 than the rest
(bias **-1.7** vs **+16.9** tC/ha).
**But a genuine structural point survives:** the 1985 residual correlates with `J_t0`
(**r = 0.47**) — the modelled initial state is anchored on *recent* litter, so for
recently-disturbed plots the model carries no memory of the previous rotation. Here that
happens to offset the general over-prediction; it is a coincidence, not correctness, and is
worth a sentence in the limitations.

**Finding 2 — 10 plots have LITERALLY CONSTANT litter for 39 years [CHECK, unresolved].**
Identical to 4 decimal places for 1986-2024 (e.g. plot 19571 = 2.336 every year), constant in
the SOURCE too, and **all 10 are calib_ready**. No stand produces identical litter for 39
consecutive years; the litter model most likely returned a fixed value where the underlying
biomass series was static. Same *class* of problem as the 1985 artefact — plausible-looking,
positive, finite, and wrong. They contribute a spurious "no temporal change" forcing to a
paper whose subject IS the temporal change. Not yet acted on.
Plots: 19571, 35391, 37591, 39251, 41351, 43651, 43771, 45331, 61591, 67631.
**Decide:** flag-and-keep (documented) vs exclude. Worth raising with Boris alongside 1985.

**Minor:** 2 isolated zero plot-years (79651/1985, 23571/1995); 1 plot marginally over the
~9 tC/ha/yr NPP ceiling (33452 at 9.35).

---

## 4. Open questions / people

- [ ] **Boris Tupek** — confirm 1985 is a first-year differencing artefact. Currently an
      inference from its shape, flagged as unconfirmed in three places, but now load-bearing
      in the M&M.
- [ ] **Aleksi Lehtonen** (back ~mid-Aug) — corroborate the stock QC and own the understorey
      question. NB the QC itself is **ours**; his role is review, not execution.
- [ ] **Understorey fraction** still rests on a literature range (~15-35%); pinning it down
      needs the LUKE MUSTIKKA ground-vegetation data, absent locally.

---

## 5. Tests scoped but not run

- **SOC-target ablation** (old vs new target) — the only *destructive* test: it regenerates
  the old bundle and overwrites `Data/model_inputs/`. Do it deliberately, with a backup, not
  inside a batch.
- **Measured-only (0-40 cm) target** — retires the depth-extrapolation question entirely for
  the robustness appendix. Cheap: a column swap.
- **C1 ablation** — the largest round-1 change and still completely unattributed; its headline
  claim (anchoring rates makes sigma_input fall *as a consequence*) is untested. Most invasive.
- **Free post-processing** on posteriors we already have: pre-run direction diagnostic, KL
  information gain per parameter, posterior-predictive coverage *with* observation error
  (the long-standing "Param cov ~0.05" artefact).

---

## 6. Session-state snapshot (written 2026-08-04 evening)

**Running detached** (both drivers have PPID 1 — they survive terminal/session close):
`run_all_ablations.sh` -> TP2, SP1, Yasso07, Yasso20 x 6 configs, then a `queue_tp2_loco.sh`
catch-up adding `A5_LOCO_no1985` to TP2. Progress in the scratchpad `queue.log`; per-config
logs land in `doublechecks/ablation_logs/`. If a driver died, look for the last config log
written and restart from there with:
`Rscript doublechecks/run_ablation.R <MODEL> 6000 3 <comma-separated-configs>`

**NOTHING IS COMMITTED.** Deliberate (get the data right before Roihu), not an oversight.
The working tree is the only copy of today's code changes apart from OneDrive version history;
`Data/` is gitignored throughout, and `Data/Data_work_pre_SOC_swap_20260804.R` is the pre-swap
backup. Committing is safe at any point and does not touch Roihu.

**Two upstream data errors fixed today**, both previously invisible and both present in every
run since April 2026 (including the `20260710_*` production posteriors, which are therefore
stale for two independent reasons):
1. SOC target — missing stoniness correction + drifted per-campaign processing.
2. Litter — artefactual first year (= t0).
Net: the observed 1985->2006 rise fell from **+61% to +12%**. The paper's core accumulation
signal is real but roughly a fifth of what the old figures implied. Every number quoted
against a campaign mean must be restated.

**Guards now in place** so neither class of error can recur silently: `assert_inputs_current()`
(aborts on a stale bundle), the C5 factor logged at startup, ablation posteriors quarantined
away from the predictive stage's newest-file auto-detect, and `doublechecks/input_data_qc.R`
as a repeatable input sweep.

---

## 7. ABLATION RESULTS (2026-08-05) — C5 verdict

### Usable evidence
**TP2 only fully clean** (R-hat 1.009-1.070). **SP1 usable** except A0/A4 (R-hat 1.34/1.38).
**Yasso07 + Yasso20 EXCLUDED** — R-hat 1.40-2.31 at 6000x3; ~19 free parameters need far more
than the setting that worked for TP2's 7. Re-running A0/A1/A5 at **20000x3** (~18 h).
**Lost and re-run:** TP2 `A1_C5_off` — I edited the run scripts WHILE the suite was launching
from them; it read a truncated file (`mcmc_c_chains` = garbled `run_mcmc_chains`) and never
saved a posterior. `A2` crashed at the end but saved a valid posterior.
**LESSON: never edit run scripts while a suite is launching processes from them.**

### Rule 2 — ANSWERED: C5 does no useful work
Trusted-campaign (2006+2024) RMSE vs C5 strength:
| C5 | TP2 | SP1 |
|---|---|---|
| 1.5 | 45.03 | 47.35 |
| 2.0 | 44.63 | 47.43 |
| 3.0 | 45.26 | 47.39 |
TP2 spread **1.4%**, SP1 spread **0.25%**, noise floor ~2%. Halving or tripling C5 changes the
fit to the trusted campaigns by nothing, in two structurally different models independently.

### Rule 1 — models disagree; SP1's answer is unphysical
Bias-controlled 1985 excess (1985 bias MINUS the model's bias on campaigns it DID see — the
control matters: without it any positively-biased model would "prove" VMI8 reads low):
- **TP2 LOCO: +6.4** (below the +8-10 threshold) -> does not support C5
- **SP1 LOCO: +13.0** (above threshold) -> would support C5
BUT SP1's LOCO sits at **sigma_init = 1.126**, i.e. the 1917 flux EXCEEDS the whole-record
mean — past the 0.826 inversion threshold and contradicting the growing-stock history. Its
answer comes from an unphysical corner; likely SP1's single pool straining without the 1985
anchor rather than evidence about VMI8. TP2's LOCO is at 0.780, just under the threshold.

### Verdict (pending the long Yasso runs)
**DROP C5** — rule 2 fails decisively in both models; rule 1 fails in the trustworthy one.
State in the paper that the 1985 correction was made in the DATA, not the likelihood.

### Two findings NOT looked for
- **C3 costs fit, consistently.** `A4_C3_off` gives the BEST trusted RMSE in both models
  (SP1 44.72 vs ~47.4 = 5.7%; TP2 43.26 vs 44.63 = 3%), both above the noise floor, same
  direction. C3 is better motivated physically but slightly worse fitting. State the tension.
- **sigma_input diverges hard by model.** SP1 **0.23-0.61** (tree litter ~2x too MUCH) vs TP2
  **1.2-2.2** (too little). SP1 is in genuine conflict with the C4b prior centred at 1.30 on
  the argument that understorey litter is MISSING. Sharpens the existing sigma_input result.

### Metric caveat (do not repeat this confusion)
`run_*_predictive.R` reports **`cor(obs,hat)^2`**, which ignores bias — that is the source of
the documented "R2 0.05-0.11". The variance-explained form (1 - SSres/SStot) is far harsher and
goes NEGATIVE here because of the ~+10 tC/ha over-prediction. Both are now printed by
`ablation_fit_by_campaign.R`. On the production metric TP2_A0 scores 0.023 — low, but the same
order as before, NOT a collapse. The componentwise posterior median was checked and IS
representative (0.37 sd from the nearest draw vs 3.48 typical), so none of this is an artefact
of evaluating at the median.

---

## 8. STATUS 2026-08-05 afternoon — ROIHU IS RUNNING

**Jobs 474800–474805 launched 13:24, all six models, all RUNNING.** Verified live in the
`.err` files: `Cores per chain: 40` (the OOM trap avoided), `5 chains x 50000 iterations`
(no ablation env vars leaked), and all six reached `Chain 1 / 5` — which means
`assert_inputs_current()` passed, forward sanity passed, and the Yasso `.so` files loaded.

**⚠ Checklist correction:** the guard output goes to **`.err`**, not `.out` — R's `message()`
writes to stderr. `.out` only carries the prior tables. Grep the `.err` files, and restrict
the glob to the job IDs (`*_4748*.err`) or you will be reading July's 98 KB logs.

### Decided and implemented today
- **C5 WITHDRAWN.** `SIGMA_1985_INFL` default 2.0 → 1.0. Rationale in `calib_config.R` at the
  point of use, in `REVISION_PLAN.md` §C5, and in full in the working document.
- **δ-offset tested and also rejected** — identifiable (so the revision plan's
  "non-identifiable" claim was wrong) but a misfit sink.
- **8 constant-litter plots excluded**; calib-ready 520 → 512. Plots 39251/67631 borderline
  and deliberately kept.
- Pre-flight re-cleared on 512 plots with C5 off; Yasso07 blow-ups 2.5% (historical band).
- Committed as `7fc0c82`.

### Local runs still going
Yasso07 ablation `A1` (~15:50), then a watcher stops the driver before `A5` — `A0`/`A1` used
the 520-plot bundle and anything launching now would use 512. Yasso20 leg cancelled. Fold the
`A0`-vs-`A1` result into the working document when it lands.

### Docs status (audited 2026-08-05)
- ✅ current: `CLAUDE.md`, `Data/Data_work.R` (M&M), `manuscript/HIKET_main_manuscript.tex`,
  `manuscript/HIKET_data_and_ablation_tests.tex`, `manuscript/REVISION_PLAN.md`, this file.
- ⬜ **`Calibration_real_data_transient/documentation/HIKET_calibration.Rmd`** — NOT updated for
  the SOC target, the litter reconstruction, the constant-litter exclusion or the C5 withdrawal.
  **The largest remaining doc gap.** Best done after the run, when the numbers are final.
- ⬜ `manuscript/HIKET_storyline_note.tex` — carries the SOC "levels are provisional" note but
  not the litter fix. Its figures are all stale anyway; rewrite once figures are rebuilt.
