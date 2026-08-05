# NEXT SESSION — start here

**Rewritten 2026-08-05 evening**, superseding the running checklist kept during 4–5 August.
Current as of commit `0277ad7`.

---

## 0. State in one paragraph

Two upstream data errors were found and fixed (SOC stocks missing a coarse-fragment
correction; the litter product's first year unusable), eight constant-litter plots were
excluded, and C5 was withdrawn after ablation testing. The six-model Roihu recalibration
launched **2026-08-05 13:24** (jobs **474800–474805**) with all of that in, and is the first
run to include C1–C4. Everything downstream — every figure, table, metric and the NextGenC
bundle — predates it and is stale. Nothing is running locally.

---

## 1. FIRST THING: check the Roihu jobs

Launched 13:24 on 2026-08-05, ~36 h, so expect completion around **01:00 on 2026-08-07**.

```bash
ssh roihu                      # needs a freshly signed cert (24 h validity)
squeue -u menichet
sacct -u menichet --starttime 2026-08-05 --format=JobID,JobName%16,State,Elapsed,ExitCode
```

If finished, check health before trusting anything:

```bash
cd /scratch/project_2019134/HIKET/Calibration_real_data_transient/progress_logs
grep -H -E "R-hat|inf_rate|Wallclock" *_4748*.err
```

> **The guard output goes to `.err`, not `.out`** — R's `message()` writes to stderr; `.out`
> only has the prior tables. And restrict globs to `*_4748*` or you will be reading July's
> 98 KB logs, which is what happened the first time.

**Already verified at launch:** `Cores per chain: 40` (the OOM trap avoided),
`5 chains x 50000 iterations` (no ablation env vars leaked), and all six reached `Chain 1 / 5`
— which means `assert_inputs_current()` passed, forward sanity passed, and the Yasso `.so`
files loaded correctly.

### Then sync back

```bash
# from the Mac. Use the `roihu:` ALIAS -- the bare hostname has no certificate.
rsync -av roihu:/scratch/project_2019134/HIKET/Calibration_real_data_transient/runs/ \
  ./Calibration_real_data_transient/runs/
rsync -av roihu:/scratch/project_2019134/HIKET/Calibration_real_data_transient/diagnostics/ \
  ./Calibration_real_data_transient/diagnostics/
rsync -av roihu:/scratch/project_2019134/HIKET/Data/model_inputs/ ./Data/model_inputs/
```

The third is not optional: the predictive stage hard-loads
`Data/model_inputs/<MODEL>_inputs_<RUN_ID>.rds` with no fallback.

**Before running anything downstream:**

```bash
Rscript doublechecks/quarantine_ablation_runs.R
```

`run_*_predictive.R` picks its posterior by sorting filenames and taking the newest, with no
notion of what kind of run produced it. The local ablation posteriors sit in the same
directory and would be picked up silently.

---

## 2. What to look at in the results

Four numbers carry most of the meaning.

- **σ_input.** Expect ~1.0–1.2. The old pathology (TP2/TP3 at 13–20×, manufacturing carbon)
  should be gone twice over — once from the flux bound, once from an honest target.
  ⚠ **Watch the other direction too**: in the ablations SP1 landed at 0.23–0.61 and Yasso07 at
  0.945 once C5 was removed — i.e. *too much* litter, against a C4b prior centred at 1.30
  precisely because understorey litter is **missing**. Two of three models on the wrong side of
  that prior is a live tension, not a rounding error.
- **σ_init vs 0.826.** Above that threshold the reconstructed 1917→1985 pre-run *declines*,
  contradicting the growing-stock history C3 encodes. Ablations put it at 0.29–0.78 depending
  on configuration; the upper end is close.
- **The 1985 gap.** The old "+26 tC/ha above observed 1985" figure is tied to the superseded
  target and **must be re-derived**. Do not quote it.
- **Skill.** Metric trap: `run_*_predictive.R` reports `cor(obs,hat)^2`, which ignores bias —
  that is the source of the documented "R² 0.05–0.11". Variance-explained R² is much harsher
  and goes negative here because of the ~+10 tC/ha over-prediction. Both are printed by
  `doublechecks/ablation_fit_by_campaign.R`. Do not quote one against the other.

---

## 3. Then: refresh everything downstream

- [ ] Stages 2–4 locally: `Rscript Calibration_real_data_transient/run_hiket_pipeline.R --skip-calibration`
- [ ] Figures **F1–F14** and **S1–S7**, tables **T1/T2** — all predate the rebaseline
- [ ] NextGenC bundle (`Reporting/NextgenC_report/`: `build_soc_matrices.R` → RUN_IDs →
      `plot_soc_trajectories.R` → `build_soc_maps.R`)
- [ ] **Restate every number quoted against a campaign mean.** The observed 1985→2006 rise fell
      from **+61% to +12%**: the effect the paper rests on is real but roughly a fifth of what
      the old figures implied. State it at its true size rather than softening it.

---

## 4. Documentation status (audited 2026-08-05)

**Current:** `CLAUDE.md`; `Data/Data_work.R` (M&M block — the canonical methods text);
`manuscript/HIKET_main_manuscript.tex`; `manuscript/HIKET_data_and_ablation_tests.tex` (the
working record, 9 pp); `manuscript/REVISION_PLAN.md`; this file; memory.

**Still to do:**

- [ ] **`Calibration_real_data_transient/documentation/HIKET_calibration.Rmd`** — not updated
      for the SOC target, the litter reconstruction, the constant-litter exclusion or the C5
      withdrawal. **The largest remaining doc gap.** Deliberately deferred: it is the methods
      document and should carry final numbers rather than be rewritten twice.
- [ ] **`manuscript/HIKET_storyline_note.tex`** — carries the SOC "levels are provisional"
      reading note but not the litter fix. Its figures are all stale; rewrite once figures are
      rebuilt, and drop the provisional caveat then.

---

## 5. Decisions taken 4–5 August (all committed)

| | decision | why |
|---|---|---|
| SOC target | homogenized three-campaign baseline | mineral C over-counted ~1.6× (missing stoniness); campaigns on drifted paths |
| Litter 1985 | reconstructed by 1986–1990 backcast | first year unusable (0.060 vs 1.398); it is t₀, so it corrupted four quantities |
| 8 plots | excluded (`const_litter`) | litter identical to ~15 s.f. for 39 years — not a measured series |
| **C5** | **withdrawn** | inert for prediction (0.25–0.45% across three models), moves σ_init ~4× |
| δ-offset | tested, **not adopted** | identifiable, but a misfit sink |
| C3 | **kept** | justified by the growing-stock record; ~3% fit cost reported honestly |

Full reasoning, including where the analysis went wrong before it went right, is in
`manuscript/HIKET_data_and_ablation_tests.pdf`.

---

## 6. Open questions / people

- [ ] **B. Tupek** — confirm (a) that the 1985 litter collapse is a first-year differencing
      artefact, and (b) the constant-litter plots. Both treatments are defensible without him
      — the values cannot be real either way — but the *stated mechanism* is an inference, and
      it is now load-bearing in the M&M.
- [ ] **A. Lehtonen** (back ~mid-Aug) — corroborate the stock QC, own the understorey question.
      NB the QC itself is **ours** and was performed here; his role is review, not execution.
- [ ] **Understorey fraction** still rests on a literature range (~15–35%); pinning it down
      needs the LUKE MUSTIKKA ground-vegetation data, absent locally. This matters more now
      that two models put σ_input *below* 1.

---

## 7. Tests scoped but not run

- **C1 ablation** — the largest round-1 change and still completely unattributed; its headline
  claim (anchoring rates makes σ_input fall *as a consequence*) is untested. Most invasive.
- **SOC-target ablation** (old vs new) — the only *destructive* test: it regenerates the old
  bundle and overwrites `Data/model_inputs/`. Do it deliberately, with a backup.
- **Measured-only (0–40 cm) target** — retires the depth-extrapolation question for the
  robustness appendix. Cheap: a column swap.
- **Yasso `A5` (LOCO) and the whole Yasso20 leg** — deliberately not run: the bundle changed to
  512 plots mid-batch, so anything launching later would not have been comparable. Re-run on
  the final data if wanted for the record.
- **Free post-processing** on the new posteriors: pre-run direction diagnostic, KL information
  gain per parameter, posterior-predictive coverage *with* observation error (the long-standing
  "Param cov ≈ 0.05" artefact).

---

## 8. Two process lessons worth keeping

1. **Never edit run scripts while a suite is launching processes from them.** TP2's `A1` read a
   half-written file (`mcmc_c_chains` = garbled `run_mcmc_chains`) and was lost.
2. **The run scripts source `calib_config.R` without `local=`**, so its `N_ITER` / `N_CHAINS`
   land in the global environment and will silently overwrite same-named variables in any
   harness that sources the setup. Prefix harness variables (`TEST_*`).
