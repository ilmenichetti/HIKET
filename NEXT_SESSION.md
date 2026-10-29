# NEXT SESSION — start here

## SESSION 2026-10-27/29 (closed) — PULL FIRST (Lorenzo wrote Intro prose on Overleaf in parallel)
- ⚠ **VPN on ⇒ GitHub fails with "SSL certificate problem"** (curl too). Not a git config issue: turn the VPN off.
- ⭐ **Input gap staged as R2's UNEXPECTED result** (`634680c`, `% ROUND 9`): new subsubsection
  `sec:r2inputgap` "Unexpected: the input corner gives way" after the ridge (gap 1.28–1.75×; data-set not
  prior-set; persists in Yasso20; smaller but present under the balanced start, Y15 1.54→1.31, TP2 1.75→1.60;
  independent estimates moved here). Intro plants the product's SCOPE (understorey, mortality, below-ground)
  without announcing. R4: "Two gaps at two scales" (national level = input gap; plot pattern = start state);
  subparagraph retitled "At plot scale, a wrong starting state rather than a missing input".
  ⚠⚠ **The measured candidates do NOT bracket the gap**: understorey 1.13 + mortality 4–6% ≈ 1.18–1.19, NID full
  estimate ≤ 1.26, all below the lowest posterior 1.28 ⇒ the heading "a set that brackets the input gap" is wrong
  unless the below-ground terms carry the remainder. ⬜ Lorenzo decides.
  ⬜ Check on the current run: arm-A figure (1.28–1.50× above centre, earlier setup) and Yasso20 still at its
  published transit time.
- **Intro texts given in chat (Lorenzo pasting on Overleaf; check what landed):** (1) "a small relative change is a
  large flux" — Minasny et al. 2017 *Geoderma* 292:59–86 (8.9 Gt C/yr ÷ 2400 Gt to 2 m = 0.4%, VERIFIED in the OA
  PDF); our 0.26/70 = 0.37%/yr ≈ 4‰. Use as arithmetic only, never endorse the initiative. ⬜ **`Minasny2017` NOT yet
  in `library.bib`** (entry in chat; add if Lorenzo has not). (2) one phrase on why most inventories model ("rely",
  not "relies"). (3) measured vs modelled: `SwedenNID2026`, `Gruneberg2014`, `StatFin2025NID`, `NorwayNID2025`,
  `Didion2016` (CH/AT are secondary via Norway NID fn 50). (4) Two timescales: transit time follows the INPUT, most of
  which leaves within years; adjustment follows the STOCK, held by the slow pool. Worked example (ours, uncited):
  100 in, 90 → 2 yr pool, 10 → 300 yr pool ⇒ stocks 180/3000 (94% slow), MTT 31.8 yr, adjustment centuries.
  Cite `Sierra2018` (stock ages vs transit times) + `Wutzler2007` (slow pool adjusts over centuries).
- NextGenC (not HIKET): Yu-Feng's CTree questions answered; pending items recorded in that repo's CLAUDE.md.

## SESSION 2026-10-26 (closed)
- **S10 saturation figure** (`build_S_saturation.R`): stock panels + sink per 20-yr repetition as grouped bars; benchmarks
  at alpha 0.3 with solid borders. Still "Promote?".
- **Introduction opening restructured** (`% ROUND 9`): inventories matter → report changes, not stocks → mostly modelled →
  start state → soil memory. ⬜ bib: Paris Agreement, EU 2018/841 (+2023/839), IPCC 2006 vol. 4; ⬜ source "soil term =
  largest uncertainty"; ✅ Finnish 2021 example: CUT from the opening, stays in the case.
- **Intro opening, later the same day (Lorenzo writing prose on Overleaf in parallel — PULL FIRST):** 1.1 retitled "a
  modelled soil sink depends on where the soil starts"; abstract plan (1) follows; stage-setting cut to 2–3 lines; soil
  subparagraph = small arc: forests usually sink (Pan2024, ~3.5 PgC/yr ≈ ½ fossil) → not always (deforestation cancels ⅔;
  boreal sink −36%) → boreal special case (Pan2011: 239 vs tropical 242 tC/ha, 272 PgC ≈ ⅓ of forest C, ~60% in soil)
  → high latitudes warm fastest (IPCC2021SPM, Rantanen2022 — both ARCTIC; ⬜ boreal-specific source, Gauthier 2015 unread).
  Bib added + Crossref-checked: Pan2024, Rantanen2022, IPCC2021SPM. ⬜ read the Pan2024 correction (Nature 632:E7).
  ⚠ "a quarter of fossil emissions" was wrong (Pan2011 gives none; ≈⅓ gross 1990–2007) — do not reintroduce.
- **Aleksi's question** (email 2026-10-26: is a steady-state-1970s inventory run outside the measured uncertainty?):
  `doublechecks/ghgi_steady_state_vs_obs.R` uses the inventory's OWN Yasso07 run (LUKE package). Nationally (area-wt.)
  it matches (+0.11) only by CANCELLATION: South measured +0.48/+0.16 vs inventory +0.19/+0.05 (outside), North measured
  −0.15/−0.16 vs inventory +0.13/+0.07. ⬜ reply to Aleksi (Samuli: keep for the paper/final report, not the interim).
- ⭐ **Far North LOSES SOC** (39 plots >66°N, −0.16 [−0.27, −0.04], smooth latitude gradient) → **Fig S2**
  (`manuscript/figures/build_S_soc_by_region.R`) + one paragraph under Fig 2. Lorenzo: paper stays NATIONAL, no figure is
  split. ⬜ `\gap` there: state the caveats (39 plots; first-campaign comparability) or not. Supplement figures after S2
  shifted by one (saturation is now S11). Memory `north-loses-soc`.
- **All draft-v2 figures redrawn at print size** (memory `figures-print-size-convention`). ⬜ offered: benchmark alpha
  treatment in F17, F19, attribution figures.

## ⭐⭐⭐⭐⭐ RESUME HERE (closed 2026-10-22 evening) — Lorenzo reviews the draft v2 skeleton

**State:** everything committed and pushed (`f518a92` + Overleaf merges). `manuscript/HIKET_draft_v2.tex` is now the
MASTER (59 pp, builds clean): M&M in full (round 8), R1–R4 / Intro / Conclusions as nested OCAR-tagged headings with
red `\gap` gists and red `\refs` lines, all figures anchored with captions, Supplement S1–S6 in the same file.
`HIKET_methods_v2.tex` is only the record of the M&M review — do not edit it further. On Overleaf switch the main
document to `manuscript/HIKET_draft_v2.tex`.

**First:** `git pull`, re-read draft_v2 from disk (Lorenzo edits in parallel), collect his notes.

**Decisions waiting for Lorenzo (skeleton stage):**
- F2 (rebuilt: Yasso07 published, balanced start, product litter, 310 plots; +22/+17/+12 above obs, flat) → move to
  R1 as the first picture of the problem? And the wrong Methods pointer to S_baseline (red flag in the SOC section).
- Keep the KL paragraph (no result uses KL)? n_eff 220 (Methods) vs 221 (storyline)?
- Promote F19, F9, the saturation figure (still to build)? Keep harvest / Rh-vs-age figures? Merge F10b with F3alt?
- Name the balanced-start comparison in the objectives? Third-corner run (balanced start + widened rate priors) or a
  limitation?

**Remaining Methods `\gap`s:** exclusion rules, 2024 stone volumes, data statement, region codes, prior-criteria
table `tab:priors` (from `Prior_specs/PRIOR_HOMOGENIZATION_PLAN.md` §2, with the CORRECTED input prior: N(·,0.25) on
the logistic scale of the flux within [0.05, 8.7] ⇒ σ_input median 1.07, 95% 0.75–1.45; σ_init 0.93, 0.58–1.49),
five terms (after the new inputs).

**Other open:** cite Lauri Mehtätalo (which work?); stone effect vs stone correction; ✅ saturation figure BUILT 2026-10-26 (`manuscript/figures/build_S_saturation.R` → `S_saturation.png`, in draft v2 Supplement + numbers in the R3 side-arc gists; promote?);
F9 file name still says "auxiliary"; "log-SD 0.25" statements in CLAUDE.md / memory / v1 / storyline should state the
effective prior (offered, not done). New inputs (Boris understorey litter, Ossi 2006 biomass) → `submit_arms.sh all`
+ `refresh_all_local.sh`, then restate.

**R3b result (done 2026-10-22, F15 third arm):** the start moves the trajectory, hardly the warming sensitivity
(transient − balanced: −0.2 / −0.9 / −1.0 tC/ha at 2084, SSP2-4.5, CIs span 0); only Yasso07 separates from its
published parameters (−2.66 [−4.21, −1.40]).

---

## (closed 2026-10-22) RESUME HERE (2026-10-21 evening) — Lorenzo reviews the new M&M

**TODO tomorrow (2026-10-22), in this order:**
1. ✅ **DONE 2026-10-22 — R3b.** F15 now has a third arm (`HIKET_F15_ARMS=eqinit` merges it into the existing
   .rds; backup `doublechecks/f15_forward_experiment_pre_eqinit_20261022.rds`). SSP2-4.5 at 2084, tC/ha (ours /
   balanced / published): Yasso07 −7.3 / −7.1 / −4.7; Yasso15 −5.5 / −4.6 / −5.0; Yasso20 −6.0 / −4.9 / −6.3.
   Ours − balanced: −0.2 [−2.3, 2.2], −0.9 [−2.1, 0.3], −1.0 [−2.4, 0.4] ⇒ the start moves the trajectory, hardly
   the sensitivity; only Yasso07's published parameters separate. Methods `\gap` replaced; storyline v4 updated.
   OLD TEXT: **While Lorenzo reads the M&M: run the warming response under the balanced start (R3b).** Add a third arm
   `eqinit` to `doublechecks/f15_forward_experiment.R` (draws from `<M>_eqinit_posterior_*`, `sigma_init` = 1,
   the run's own input bundle / pre-run climate as the other arms), Yasso07/15/20 only, 40 draws, ~20–30 min.
   Then show it in F15 (third, dashed line) or as a small table (response at 2084, SSP2-4.5, per model and
   start), and turn the R3b `\gap` in `HIKET_methods_v2.tex` §Simulations into one sentence.
2. **M&M revision cycle — notes from Lorenzo's reading (2026-10-22). He is editing the .tex himself; when he
   hands it over, `git pull`, re-read from disk, then:**
   - ✅ DONE 2026-10-22 (bib `Yasso15code`/`Yasso20code`; files are in YASSOmodel/YASSO15 and /Yasso20, NOT Ryassofortran). **Yasso15/20 climate source — DECIDED: cite the REPOSITORY.** Yasso20 → `Viskari2022` (+ Zenodo
     10.5281/zenodo.5059909). Yasso15 → the posterior sample distributed with the model code (FMI
     Ryassofortran, `data/Yasso15.dat`; state version/commit or access date; use a Zenodo DOI if one covers
     it). Viskari 2022 never mentions Yasso15 — no peer-reviewed source exists. Not "in preparation".
     Add a bib entry, replace the `\gap`. Optional: ask Toni whether a citable Yasso15 reference exists.
     Name it as a remaining asymmetry: Yasso15 widths come from an unpublished calibration.
   - Check whether he applied (answered in chat, NOT edited by me): natural mortality — "that product" =
     Tupek's, it has NO mortality term (wording implied partial); upland AREA held 13.90 Mha → linear rise
     1965–1980 → 15.22 Mha; the pre-run equation $J(t)=J_{1917}+(J_{1985}-J_{1917})\,s(t)$ applied per AWEN ×
     size (composition fixed at 1985–89) and naming $J_{1917}=\sinit J_{1985}$; "what remains is stated with
     each compromise" is honoured only for the rates — add one leftover sentence per subsection (fractions:
     number 0/1/2/12; climate: borrowed response, 1 vs 3 modifiers; auxiliary: none); the ICBM rescaling
     = divide by Yasso07 ξ at Ultuna (5.4 °C, amplitude 10, 520 mm) = 0.94.
3. Lorenzo's review of `HIKET_methods_v2.tex` → settle its `\gap`s → merge into `HIKET_draft_v2.tex`.
3. ✅ DONE 2026-10-22: `HIKET_draft_v2.tex` is the full skeleton (M&M copied whole after round 8; R1–R4 + Intro +
   Conclusions as numbered nested headings to \subparagraph with \gap gists; 12 main figures + input table; Supplement S1–S6
   in the same file). From now on edit the M&M in draft_v2, not in HIKET_methods_v2.tex.

**Headroom: DECIDED 2026-10-21 — the paper quotes the 2024 headroom (21–40%, `forward_scenarios.R`, F19).**
The 1985 figure (`doublechecks/equilibrium_headroom.R`) is NOT important (Liski comparison only). Its script
was FIXED the same day: it built the 1985 state under the plot's 1985–2024 climate instead of the run's
1961–1984 pre-run climate (same bug class as F15), with 5 draws. Now pre-run climate + bundle litter + 100
draws: C_1985 62–67 tC/ha (was 52–58), headroom SP1 12 / TP2 63 / TP3 65 / Yasso07 57 / Yasso15 48 /
Yasso20 41 %. Open item (b) of §0-NOW is closed.

**Decided (Lorenzo, 2026-10-21):** `HIKET_draft_v1.tex` will be DROPPED. A new draft, `manuscript/HIKET_draft_v2.tex`,
is built on storyline v4: for now an EMPTY skeleton (title, authors, coarse sections: Intro, M&M, Results and
discussion R1–R4, Conclusions, Acknowledgements, Data). Next it gets nested headings (many levels), a little text,
figures and captions, for the coauthors. **The M&M is the one section written in full**:
`manuscript/HIKET_methods_v2.tex` (standalone, 19 pp), for Lorenzo's review → then merged into v2.
- M&M v2 = v1's Methods (rounds 2–7) + NEW: pre-run climate (1961–1984 mean, plot mean 2.2 °C), the comparison
  with a balanced start (with the four pre-registered expectations — `\gap` asks whether to keep them in Methods),
  the division of the sink (recalibrated, three parts), and the forward runs (projection under both starts, input
  that holds the sink, warming scenarios with the published arm, Ruosteenoja & Jylhä 2021 added to the bib).
  CHANGED: order setup → calibration → experiments → evaluation; R 4.6.1 / BayesianTools 0.1.9 (v1 said
  4.3.1 / 0.1.8); convergence of all 18 calibrations stated (R-hat ≤ 1.009, ESS ≥ 1716).
- Open `\gap`s in it: exclusion rules not written out; 2024 stone volumes; data statement; region codes;
  Yasso15/20 climate source; the prior-criteria table; the arm-A sensitivity (earlier setup only); five terms;
  warming response under the balanced start (R3b).

## ⭐⭐⭐ RESUME HERE (closed 2026-10-13 evening)

**✅ 1 DONE 2026-10-21:** `eqinit_climheld` landed (RUN_IDs `20261013_13*`, R-hat ≤ 1.007, ESS ≥ 1845), synced,
decomposed, figures `S15/S16/S17_*_recal.png`. Recalibrated three parts (MEANS, six-model): inherited **+11.6**,
litter **+7.3**, climate **−11.0**, sink +8.0. The balanced-start calibration loses carbon in every draw ⇒
inherited > sink in all six. ⭐ **DECIDED (Lorenzo): the fixed-parameter division is DROPPED** — the appendix
(`appendix_attribution.tex`) was rewritten for the recalibrated version only (three calibrations, one assumption
added at a time; joint effects inside the parts by order) and the storyline Methods now present the division as
one more calibration on top of the balanced-start comparison. The recalibrated path now ENFORCES two checks
(balance drift ≤ 1.7e-13, adding-up ≤ 3e-16). The parameter-source test (`HIKET_ATTR_PARAMS=eq`) and the
fixed-parameter joint terms are no longer in any text; their code paths and old S15–S17 figures remain on disk
(still run by `refresh_all_local.sh`). Five terms wait for the new inputs (`submit_arms.sh all`).

~~**1. First: the `eqinit_climheld` arm** (6 jobs, launched by Lorenzo the evening of 2026-10-13 with
`submit_arms.sh eqinit_climheld`): jobs **2265953 SP1, 2265954 TP2, 2265955 TP3, 2265956 Yasso07, 2265957 Yasso15, 2265958 Yasso20** (submission order). ✅ All six started; switches verified in the logs (EQUILIBRIUM INIT + CLIMATE HELD + PRERUN CLIMATE 1961–1984 ON, 40 cores). `squeue -u menichet`; logs with the ONE-line check
below (§ recalibrated decomposition). When landed: rsync `runs/`, `diagnostics/`, `Data/model_inputs/`, then
`HIKET_ATTR_MODE=recalibrated Rscript doublechecks/attribution_decomposition.R` and
`HIKET_ATTR_RECAL=1 Rscript manuscript/figures/build_S15S17_attribution.R` → a THREE-PART recalibrated S16
(temporary main figure). Compare with the expectations written before launch (§ recalibrated decomposition).
Then restate S16 numbers in the storyline (R2 end), the attribution appendix, and decide the appendix's
fixed-parameter vs recalibrated presentation.~~

**2. Where the story stands** — `manuscript/HIKET_storyline_v4.tex`, recomposed with Lorenzo's email to Samuli and
Aleksi (`~/Downloads/The storyline I propose.eml`): general → Finnish case (saturation?) → trilemma in Finland as a
hypothesis; R1 balanced start fails; R2 relaxed works, ridge shifts, and the SCENARIO DECOMPOSITION closes it
(S16 five terms in the email; temporarily three parts once recalibrated); R3 forecasts, "did we get it wrong?",
two side arcs (temperature; inputs & saturation); R4 future development, independent side arm with suggestions.
Load-bearing figures = those in the email; all others kept in the appendix (promote: F19, saturation figure, F4, F9).
**Coauthor feedback:** Samuli likes it ("can become complex if not extremely well written"; offers time to make it
digestible); his F15 question exposed a bug, now fixed (below). Aleksi likes it; Boris is building understorey
litter for 1985/1995/2006/2023 and Ossi could estimate 2006 biomass for BioSoil plots → **new inputs expected in
~2 weeks; Lorenzo: we keep the current inputs, rerun everything once they come** (`submit_arms.sh all` +
`refresh_all_local.sh`). Aleksi also pointed to Peltoniemi et al. 2006 Fig. 4 (below).

**3. Fixed today (all pushed):** (a) **F15 started both arms from the wrong climate** (1985–2024 mean instead of the
run's 1961–1984 pre-run): our arm 6–11 tC/ha low. Fixed in `published_arm_dat.R` + `f15_forward_experiment.R`;
Yasso07 arm difference now −2.66 [−4.21, −1.40]; caption says the published arm is level-fitted to the campaign
means. (b) Setup by-products (input bundles keyed to RUN_IDs without posteriors) cleaned; `refresh_all_local.sh`
now cleans them automatically.

**4. Analyses done today, not yet in the storyline/draft text** (offered, Lorenzo to decide):
- **Parameter-source test of the decomposition** (`HIKET_ATTR_PARAMS=eq`): climate −17…+9%, litter −18…−1% when
  the parameters come from the balanced-start calibration; sum negative in all. In the attribution appendix.
- **Recalibrated inherited part** (production − eqinit sinks) = fixed-parameter one within ±5% (mean 11.5 vs 11.3).
- **Persistence of the inherited sink** (answers Aleksi / Peltoniemi 2006 Fig. 4): half-time Yasso20 12, Yasso15 16,
  SP1 25, Yasso07 30, TP2/TP3 62 yr; share of the 1985 inherited sink left in 2024: 17 / 24 / 28 / 41 / 73–76%.
  Peltoniemi 2006 uses the ORIGINAL Yasso only (no Yasso07 — ask Aleksi for his Yasso07 source). Proposed for
  R2 (decomposition) + the appendix; ⬜ not written yet.
- **Within-model spread of the parts** (inherited driven by σ_init, r −0.85…−0.99; litter and climate tight and
  anticorrelated; "inherited > sink" in every posterior draw). In the appendix.

**5. Still open from the morning** (below): "within 3.5 tC/ha" in the draft does not reproduce; the
equilibrium-headroom reference climate; NextGenC far-south loss patch; arm A/B paragraph on old runs.

## ⭐⭐ 0-NOW (2026-10-13). THE PRE-RUN-CLIMATE PAIR HAS LANDED AND EVERYTHING IS RESTATED

**Runs.** Production `20261012_104221` (SP1, TP2, TP3, Yasso07, Yasso20) / `_104220` (Yasso15); balanced start
`20261012_122348` (SP1, TP2) · `_125950` (TP3) · `_130503` (Yasso07) · `_132602` (Yasso15) · `_132728` (Yasso20).
All 12 jobs COMPLETED in 13–20 h, `PRERUN CLIMATE ON` (1961–1984, 2.18 °C) in all, `EQUILIBRIUM INIT ON` in
the six eq jobs, 40 cores, R-hat ≤ 1.009, 0 warnings, 100% finite, no OOM. Synced; stages 2–4, the eq arm,
every `manuscript/figures/build_*.R`, the residual chain (F11 stand-and-soil and F11b re-copied), the core
doublechecks and the NextGenC bundle rebuilt. Old outputs kept for comparison in `/tmp/claude-501/before/`
(scratch, not durable — the old posteriors themselves stay in `runs/`).

**What moved (new vs 20260831 / 20261008):**

| | before | now |
|---|---|---|
| calib R² / holdout R² | 0.016–0.021 / ≤0.0012 | 0.015–0.020 / ≤0.001 |
| σ_input (flux) | 3.19–4.34 | 1.28–1.75 (3.22–4.40) |
| σ_init (5 models; SP1) | 0.49–0.60; 0.83 | 0.45–0.56; 0.78 |
| C_1917 | 41.9–52.3 | **46–56** (all clear the 45 placeholder) |
| MTT Yasso07/15/20 | 24.6/27.1/22.6 | 23.7/26.8/22.3 (ratios 0.71/0.88/1.02) |
| published-MTT fit cost | 9.7/3.0/0 | 8.4/1.3/0 |
| ridge r (Yasso; SP1/TP2/TP3) | −0.59/−0.53/−0.56; −0.76/−0.47/−0.45 | −0.61/−0.54/−0.54; −0.77/−0.44/−0.42 |
| rate 1985–2024 production | +0.10…+0.35 | **+0.05…+0.30** (SP1 +0.05, Yasso20 +0.10 below obs) |
| rate 1985–2024 balanced start | −0.05…+0.04 | **−0.11…−0.02 (all lose carbon)**; 1985–2006 none reaches obs |
| balanced start 1985 stock | ~72 | 73–76 (obs 64.5) |
| ΔLL eq − prod | within 6 | up to −8 (production better in 5 of 6) |
| attribution (inherited / litter / climate / total) | +9.3 / +8.6 / −7.5 / +10.5 | **+11.3 / +8.9 / −12.2 / +8.0** |
| headroom 2024 (forward) | +26…+42% | +21…+40% (SP1 ≈ 0) |
| equilibrium headroom 1985 (doublecheck) | +28…+49% | +39…+67% ⚠ reference = cold 1961–84 climate |
| sink 2025–44 production / eq | +0.07…+0.22 / −0.06…+0.02 | +0.07…+0.18 / −0.12…+0.02 |
| warming SSP2-4.5, Yasso07 arm diff | −2.08 | **−2.66 [−4.21, −1.40]** (after the F15 fix below) |

⭐ **Story consequences.** (1) The balanced start now LOSES carbon in all six — R1 is sharper. (2) The
inherited part EXCEEDS the whole modelled sink in all six; the warming (vs the colder reference) takes back
more than the litter adds. "Most of the sink is the history" → "the sink is the history, and more". (3) The
size of the climate part rests on the reference climate — state it. (4) The six-model mean modelled sink is
≈0 from ~2006 (S16), against the observed +0.12 (fragile). (5) Residual decay by 1985 basal area unchanged.

⭐ **Attribution vs the calibration's initial-state assumption (Samuli's worry, quantified 2026-10-13).** No run of
the division is reparameterised (one draw = one parameter set; runs change only start and forcing). The assumption
enters only via the calibration: rerunning with the balanced-start posteriors (`HIKET_ATTR_PARAMS=eq` →
`doublechecks/attribution_eqparams.rds`) moves climate −17…+9% and litter −18…−1% (6-model mean −6 / −10%); their
sum −3.2 vs −3.5, negative in every model under both. Exceeds within-model uncertainty only for Yasso15/20. Written
up in `appendix_attribution.tex` §"How much the parts depend on the calibration"; ⬜ manuscript section later.

🚀 **THE RECALIBRATED DECOMPOSITION — MACHINERY READY 2026-10-13, launch tonight (Lorenzo):**
On Roihu: `cd /scratch/project_2019134/HIKET && git pull`, then in a shell with `module load r-env`:
`bash Calibration_real_data_transient/submit_arms.sh eqinit_climheld` (6 jobs; no .f90 change, no recompile).
Check (after `squeue -u menichet` shows R; ONE line -- the `/dev/null` stops grep waiting on stdin):
`grep -H -E 'EQUILIBRIUM INIT\] ON|CLIMATE HELD\] ON|PRERUN CLIMATE\] ON|Cores per chain' Calibration_real_data_transient/progress_logs/eqinit_climheld_*.err /dev/null`.
Landed → rsync, then `HIKET_ATTR_MODE=recalibrated Rscript doublechecks/attribution_decomposition.R` and
`HIKET_ATTR_RECAL=1 Rscript manuscript/figures/build_S15S17_attribution.R` (or the whole
`refresh_all_local.sh`). **After the new litter data: `submit_arms.sh all` + `refresh_all_local.sh`.**
Tests passed (`doublechecks/scenario_switch_tests.R`): switches off = production ll identical (5 models);
held hooks = an independent construction of the decomposition's scenarios, exactly; a stand-in S3 reproduced
the fixed-parameter parts; a 600-iteration SP1 smoke calibration ran end to end, tagged `_eqinit_climheld`.
**Written before launch:** (1) inherited part unchanged (S1, S2 already exist: ~11.5); (2) S3 (balanced,
reference climate) within ±20% of the fixed-parameter S3 under balanced-start parameters (+8.0 six-model mean),
so the climate part ≈ −11.5 ± 2.5 — a much smaller climate part would mean the calibration compensates the
missing warming through σ_input/kinetics, i.e. the climate part rests on the priors; (3) S3 fits at least as
well as S2 (at production's best point, holding the climate RAISED the log-likelihood in 4 of 5 models tested
— SP1 −169.6 → −158.6, Yasso20 −260.1 → −249.6: worth knowing whatever the run shows).

📌 **PLAN (Lorenzo, 2026-10-13):** when `eqinit_climheld` lands → the recalibrated decomposition in THREE PARTS
as a TEMPORARY main S16 (interactions inside the parts, assigned by order, as now). The interactions (five
terms) wait until the input datasets are FROZEN (Boris's understorey litter, Ossi's 2006 biomass): then
`submit_arms.sh all` (adds `eqinit_litheld`, `climheld`) + `refresh_all_local.sh`. The joint-term code path in
the recalibrated mode and the `_recal_joint` builder path are written but UNTESTED end to end (no arms yet) --
expect small fixes then, not a rebuild. Meanwhile the joint sizes (+1.5 / −1.3 tC/ha) come from the
fixed-parameter version (appendix).

⬜ **EVENTUAL TO-DO (Lorenzo, 2026-10-13): a RECALIBRATED decomposition** — every run its own calibration, so
the division becomes "a comparison between assumptions ONLY" (conceptually more robust, easier to explain,
consistent with the rest of the methods). Runs 1 and 2 already exist (production, eqinit); needs run 3 =
balanced start + reference climate held through 1985-2024, recalibrated (new switch + 6 Roihu jobs). Already
known for free: the recalibrated inherited part (production − eqinit sinks) = fixed-parameter one within ±5%
in all six (mean 11.5 vs 11.3). Pre-check first: prior→posterior shift of the climate parameters (S11/S12) —
if prior-pinned, the recalibrated climate part will mostly test the priors. Not now: story first.

⚠ **F15 FIX (2026-10-13, Samuli's catch: "our calibrated models perform worse than published" in panel A).**
`f15_forward_experiment.R` and `published_arm_dat.R` built the 1985 start under the mean of the whole observed
window (1985–2024) instead of the run's own pre-run climate (1961–1984, ~1.2 °C colder): our arm started 6–11
tC/ha low (Yasso07 54.9 vs 66.0 in the predictive stage), while the published arm hid it by refitting σ_input
to the campaign means. Both now load the run's input bundle and use `hiket_ss_clim` (shared `setup()`); the
published arm was refitted with σ_init pinned to the NEW median (its .rds dated from 2 Sep). After the fix: ours
starts at 65.6, sits 1.5–5 tC/ha above the campaign means (the known small positive bias), published ON them by
construction (now stated in the F15 caption and panel label). SSP2-4.5: Yasso loss 5.0–7.3; Yasso07 arm diff
−2.66 [−4.21, −1.40]. ⚠ Lesson: any script that builds a starting state must use the run's pre-run climate.

⬜ **Open from this update:** (a) the draft's "within 3.5 tC/ha of every campaign mean" does not reproduce
(F3 cache: max gap 4.9 old, 5.0 new) — find its basis or restate; (b) the equilibrium-headroom doublecheck
uses the pre-run climate for C_eq → 39–67%, inconsistent with the forward headroom (21–40%) — decide which
reference the paper quotes (the draft quotes the forward one); (c) look at the NextGenC maps (far-south loss
region) before any deposit; (d) arm A vs arm B paragraph in the draft still quotes the 2026-08/09 runs (a
sensitivity on the old basis) — keep as is or rerun.

## ✅ 0-NOW (closed 2026-10-12; DONE 2026-10-13, see above)

**1. Roihu runs — check first.** The pre-run-climate pair (1961–1984 mean) is in flight: production
**2227219–24**, balanced start **2234048–53** (`progress_logs/pclim_*`). Production verified at launch
(`PRERUN CLIMATE ON`, 456 plots × 24 yr, 2.18 °C, 40 cores, 5×50000). ⬜ The balanced-start jobs were pending
(Priority, est. start 14:11–15:00 on 10-12) — grep `pclim_*_eqinit_22340*.err` for `EQUILIBRIUM INIT ON` +
`PRERUN CLIMATE`. When landed: rsync `runs/`, `diagnostics/`, `Data/model_inputs/`; run stages 2–4; then FIRST
run the three scripts edited but only parse-checked (`attribution_decomposition.R`, `build_F4_initialization.R`,
`init_state_plausibility.R` — the latter now uses production's 20-yr pre-run window, its old 1917 stocks are
stale); then restate every number in the storyline, the attribution appendix and the input-to-hold table.

**2. The storyline was rebuilt today with Lorenzo** — `manuscript/HIKET_storyline_v4.tex` (v3 deleted; NEVER
refer to earlier versions). Memory `storyline-trilemma-20261012` has the decisions. In short: general problem
first (soil C carries history; a balanced start assumes it away; that decides the forecast), Finland as the
test field; **the trilemma** (soil series / kinetics / inputs) is the Introduction's connecting object,
**visualised in R1 by F16** and echoed at the ridge (F18); **R1 = opening** (the problem tested the hard way →
first conclusions → justification); **R2 = the main arc, the resolution** (relaxing steady state; σ_init as
the implementation; F3, S16, F12, F18); **the σ's = a nested side arm at the end of R1 (σ_input only) and R2
(both)** — σ_init's need was hidden in R1, not absent; R3 forecast & GHG-balance side arm (headroom,
input-to-hold table, F15); R4 fit (R²) then plots; Discussion "how to do better". The **applied aim is now
declared in the objectives** (future of the Finnish GHG balance), with a scope box. Open questions are listed
at the end of the file. ⬜ Lorenzo to review; ⬜ maybe rename the file to `HIKET_storyline.tex`.

**3. New today (all committed with this note):**
- National GHG-inventory litter received (`Data/GHG_datasets/`, gitignored, extracted). Reconciled with the
  Tupek product: same scope (turnover + harvest residues) agrees in the South and is flat after ~2006 in both;
  the North differs (inventory higher, rising; 50 of 456 plots); the inventory's growth = natural mortality
  (0.03→0.16 tC/ha/yr, South), absent from the product. **S18** (`manuscript/figures/build_S18_…R`) cited in
  the draft's Litter input subsection (new paragraph + `\gap` on region codes). Diagnostic bars:
  `doublechecks/litter_reconciliation_ghgi.R`. ⬜ Confirm region codes with LUKE.
- **Input to hold the sink** (`doublechecks/input_to_hold_sink.R`, exact by linearity from the +20% runs):
  ~+10% litter input for 20 yr (Yasso) holds the 1985–2024 observed sink, half for 10 yr; recent rate needs ~0.
  Lower bounds (no warming). ⬜ Which target; warming-corrected column; national units (area + Aleksi).
- Attribution appendix: joint effects described as measured; three parts vs five terms. `_joint` figure
  variant via `HIKET_ATTR_JOINT=1` (PNGs committed for review; rebuild on the new runs).

## ⭐⭐ 0-NOW (2026-10-09). NEXT = A HEAVY STORY RE-DISCUSSION (Lorenzo), prepared by two new analyses

**Read first: memory `transient-init-rationale-20261009`** — five findings that justify the transient
initialisation. Do NOT rewrite storyline/manuscript before the discussion; some figures will go to the supplement.

1. **✅ Equilibrium-start counterfactual LANDED (2026-10-09, RUN_ID `20261008_134019`).** All six models, recalibrated
   with an equilibrium start, go FLAT over 1985–2024 (−0.05…+0.04 vs observed +0.26 tC/ha/yr); they overshoot the 1985
   stock, compensate through faster and more temperature-sensitive kinetics (σ_input relaxes), and project a sink ≈ 0.
   Results + pre-registration scorecard: `NEXT_RUN_equilibrium_init.md` §Results; figures F16–F19, `T_eqinit.csv`.
   ⬜ Appendix posterior comparison (parameter strain) still to design.
2. **✅ ATTRIBUTION ANALYSIS — BUILT (2026-10-08/09), final method = 3-run SEQUENTIAL split (see below).** Determine which
   components of the models — **inputs, previous C (the initial/inherited state), climate** —
   drive the results **over time**. To discuss before any code. Will feed the story update together
   with item 1.
   Approach (Samuli's note, `Calibration_real_data_transient/documentation/Launiainen_20261008_National_soil_C_sink_dynamics_attribution.pdf`): single pool `F = dC/dt = L(Vol) + I(H,U) − k(T,P)·C`; the total
   differential splits the sink change into growing stock (∂L/∂Vol·dVol), harvest intensity (αU·dH), residue
   management (αH·dU), temperature (−C·∂k/∂T·dT), precipitation (−C·∂k/∂P·dP) and stock feedback (−k·dC).
   Climate term ∝ existing stock. Extends to N pools when ∂k_n/∂T, ∂k_n/∂P are equal across pools; add
   natural mortality as its own term. ⚠ For discussion, not objections: equal climate response is needed only to
   COLLAPSE to the single-pool form; in vector form (dF/dt = b·du/dt + (∂A/∂ξ)(dξ/dt)·C + A·dC/dt) the climate
   term is exact for ANY linear model — Yasso15/20 just get THREE climate terms (AWE, N, H), each weighted by
   its own pools' carbon. Integrated along the actual path it gives an exact additive split of the stock
   change (no interaction term; small 2nd-order remainder with annual steps). Counterfactual runs WOULD
   leave an interaction term — prefer the path integral; our litter product folds
   harvest residues into J (no separate I) and has no natural-mortality term.
   **Intent (Lorenzo, 2026-10-08):** decompose each MODEL to isolate, over time, the effect of the THREE
   components — inputs, historical C (inherited state), climate. Samuli's note is the template; the
   decomposition will need adapting per model.
   ✅ **REVISED 2026-10-08 (Lorenzo: 'historical C' must tell the INITIALISATION story; extend to the pre-run;
   S naming).** Now splits the SINK ITSELF 1917–2084: historical C = pre-run sink (1917–84) + from 1985 the sink
   owed to the deficit (production start − equilibrium start, same forcing; exact by linearity); inputs/climate =
   2-factor Shapley on the equilibrium-start sink vs references 1985–89 litter / 1985–2004 climate. 5 runs per
   plot×draw, 50 draws, checks exact. Figures S15 (6 panels) / S16 (average, annual) / S17 (cumulative, HELD BACK: its projection climate band reads as future warming), each
   `_smooth10` (pre-run and post-1985 smoothed separately). ⭐ Six-model average: pre-run +11.8 tC/ha; 1985–2024
   history +9.3, inputs +8.2, climate −7.0 (total +10.5) ⇒ the 1985–2024 sink is mostly INHERITED; inputs and
   warming nearly cancel. Builder `manuscript/figures/build_S15S17_attribution.R`.
   Methodology for review: `manuscript/appendices/appendix_attribution.tex` (subfile, NOT in the draft).

3. **⬜ DECIDED DIRECTION (Lorenzo, 2026-10-09): pre-run climate = a CONSTANT LONG-TERM PRE-1985 MEAN (not a
   transient series), and attribution figures restricted to 1985–2024.**
   - Why: the pre-run now applies the 1985–2004 mean to 1917–1984 (inherited from `STEADY_STATE_YEARS`: the first
     20 years of the model climate). Plot means: 1961–84 2.13 °C, 1985–2004 2.67, 2005–24 3.83. The attribution's
     climate band is then measured against years the pre-run never saw → it sums to ~0 until ~2004 (artefact).
   - Fix: ONE long pre-1985 mean climate used BOTH for the pre-run and as the attribution's reference/balanced
     state. Then the 1985 soil differs from the balanced one only through the litter history (3 bands, as now),
     and the climate band shows the full warming since 1985. Cost: no weather variability / early-century path in
     the pre-run (modest: the pre-run mostly shapes slow pools, which respond to the mean). A TRANSIENT pre-run
     climate was considered and set aside: it adds a climate legacy to the initialisation (4 components, two
     anchors) and is harder to read.
   - Window options: **1961–1984** = FMI, ALREADY extracted per plot (`Data/model_inputs/climate_monthly.csv`,
     1961–2025, 555 plots; dropped by the merge onto the litter rows in `Data/Data_work.R`) → only wiring, no new
     data. **1917–1984** = FMI 1961–84 + CRU TS 4.09 anomalies on the FMI 1961–1990 climatology for 1917–1960
     (CRU TS 4.09: CEDA, 0.5°, monthly T+P, 1901–2024; Harris et al. 2020 *Sci. Data* 7:109). FMI ClimGrid starts
     1961 (Aalto et al. 2016). ⬜ Tietäväinen, Tuomenvirta & Venäläinen 2010 *Int. J. Climatol.* unverified.
   - Runs: production recalibrates (initial states change) + the equilibrium arm re-run with the SAME mean climate
     (Lorenzo) → ONE PAIR. Current eqinit jobs 2075931–36 stay valid for the current production.
   - Attribution figures S15–S17 now show 1985–2024 only (pre-run gain only as a CSV column).
   - ✅ **CODE READY 2026-10-12 — window 1961–1984 (FMI, already extracted; no new data).** Switch
     `HIKET_PRERUN_CLIM_YEARS=1961-1984` (`Calibration_real_data_transient/prerun_climate.R`; unset =
     production unchanged). One helper `hiket_ss_clim()` replaces the 1985–2004 slice everywhere (engine, six
     calibration + six predictive scripts, `attribution_decomposition.R`); Yasso15/20 `precip_mean` follows.
     The pre-run climate travels in the input bundle (`climate_ss_by_plot`, `PRERUN_CLIM`), so predictive and
     attribution replay the RUN's own climate regardless of the shell.
     Tests (`doublechecks/prerun_climate_test.R`): the new path on the 1985–2004 window reproduces production ll
     EXACTLY in all six; at the production posterior median the 1961–84 window lowers the pre-run ξ ~8% and
     raises the 1985 stock **+9–10%** (SP1/TP2/TP3/Yasso07) / **+5%** (Yasso15/20 — own humus modifier).
     Also with the window ON: prior pre-flight 0% blow-ups in all six (K=300, full N); equilibrium nesting
     test PASS (≤2e-13); end-to-end smoke (tiny TP2 calibration, both arms → predictive with the switch UNSET)
     replays `1961-1984` from the bundle. Updated for post-run use: `build_F4_initialization.R`,
     `init_state_plausibility.R` (⚠ it used the WHOLE 1985–2024 climate for the pre-run ξ, not production's 20 yr —
     its recorded 1917 stocks predate the fix). Attribution figures: `HIKET_ATTR_JOINT=1` → `_joint` variant with
     the two joint effects of climate as hatched bands (main effects at the reference corner; no new runs).
     ⚠ `attribution_decomposition.R`, `build_F4_initialization.R`, `init_state_plausibility.R` were edited but only
     PARSE-checked (Lorenzo: test them on the new production runs) — run them first when the runs land and check
     each prints the pre-run climate `1961-1984`.
     🚀 **PRODUCTION ARM LAUNCHED 2026-10-12** (`main` @ `7d09d8b`): jobs **2227219 SP1, 2227220 TP2, 2227221 TP3,
     2227222 Yasso07, 2227223 Yasso15, 2227224 Yasso20**; logs `progress_logs/pclim_<m>_<jobid>.{out,err}`. Verified:
     `PRERUN CLIMATE ON 1961-1984`, `456 plots x 24 yr | mean T 2.18 C`, 40 cores/chain, 5×50000 in all six.
     Equilibrium arm (`--eq-only`) SUBMITTED ~51 min later: jobs **2234048 SP1, 2234049 TP2, 2234050 TP3, 2234051 Yasso07,
     2234052 Yasso15, 2234053 Yasso20** (logs `pclim_<m>_eqinit_<jobid>.*`); ⬜ check `EQUILIBRIUM INIT ON` once running.
     **LAUNCH (Roihu):** `git pull` (no `.f90` changed → no recompile), then
     `bash Calibration_real_data_transient/submit_prerun_clim.sh --pair` → 12 jobs (`pclim_<m>` +
     `pclim_<m>_eq`). ⚠ Consider staggering the two arms (node co-location). Check:
     `grep -H "PRERUN CLIMATE" progress_logs/pclim_*.err` → `ON -- … 1961-1984` and `456 plots x 24 yr`;
     eq jobs also `EQUILIBRIUM INIT ON`. The production arm becomes the newest untagged posterior ⇒ the figure
     layer auto-picks it (intended). ⬜ After landing: `doublechecks/init_state_plausibility.R` and other
     doublechecks that slice `clim[seq_len(n_ss),]` themselves still use 1985–2004 — switch them to
     `hiket_ss_clim()` before quoting their numbers.
   - ✅ Documentation REWRITTEN IN WORDS (2026-10-09): `manuscript/appendices/appendix_attribution.tex`, 5 pp, for
     Lorenzo's review. FINAL METHOD: three runs (1 model; 2 balanced start; 3 balanced start + reference climate):
     inherited = 1−2, climate = 2−3, litter = 3. Joint effects assigned by order: litter×climate −0.9 tC/ha in the
     climate part, inheritance×climate +0.7 in the inherited part (~7–10 %). 1985–2024 average: inherited +9.3,
     litter +8.6, climate −7.5, total +10.5 tC/ha. Figures S15/S16 (+_smooth10), S17 held back; window 1985–2024.

## ⭐⭐ (2026-10-07) LORENZO REVIEWS THE MANUSCRIPT TEXT EXTENSIVELY — then we answer it

> **Lorenzo is doing a full text review of `manuscript/HIKET_draft_v1.tex`.** Next session starts
> from his notes. Before anything: `git pull` (he may push from Overleaf), **re-read the `.tex` from
> disk, never revert his prose** ([[lorenzo-edits-in-parallel]]), then
> `grep -n "%LORENZO" manuscript/HIKET_draft_v1.tex` and group the notes by section.
> ⚠ **The review had already started on 2026-10-07** (Overleaf merge `893c6da`): Introduction revised
> up to §"Three campaigns on the same plots", marked `%LORENZO: revised until here`. He applied the NID
> wording fixes of 2026-09-30 himself ("one of the reasons", foliage biomass, Luke adds global warming).
> Notes that say **"Discuss" are questions to answer in the reply (propose wording, let him choose),
> not edit instructions. Each change gets a `% ROUND 7` comment. Rebuild (`latexmk -pdf` inside
> `manuscript/`), check 0 undefined refs/citations, commit + push when he says so.
>
> ### What changed in the draft on 2026-10-07 (round 6, `% ROUND 6` comments) — he will review it
> - **sec:residuals rebuilt** — *"What the models leave behind follows the history of each stand"*:
>   **F11 RF heatmap (unchanged) opens it** → basal area emerges → `F11_stand_and_soil` (basal area +
>   fertility/texture/stoniness) → `F11b_time_since_cut` (the development-class lead made explicit) →
>   **reconciliation: one starting-state error** (small 1985 stands median 19 yr old; residual
>   +0.91 → +0.53 → +0.28 in the same plots) → **pointer to pioneer vegetation** + minor site ecology.
> - Methods §residual analysis rewritten; sec:plotscale "+40% cut plots / understorey" trimmed;
>   `\label{sec:init}`; 10 new refs (doi.org-verified). Storyline v3 Strand 3 rewritten to match.
>
> ### Open items from 2026-10-07 (not for the review itself, but keep in view)
> 1. ⬜ **Lauri Mehtätalo (LUKE) must be cited** for pioneer-vegetation modelling — no dedicated
>    published paper found (Rana et al. 2023 + Mehtätalo et al. 2025 as placeholders; `\gap` asks
>    Lorenzo which work/functions). Swedish work to extend (Peichl 2022, Grelle 2023, Dahlgren Lidman
>    2024 in; Fahlvik & Nyström 2006 / Heureka young-stand models verified, not cited).
> 2. ⬜ **How far to develop the pioneer point** — to discuss with Lorenzo. Honest caveat in the text:
>    post-1985 dated cuts show no comparable excess, so inheritance vs pioneer input is not separated.
> 3. ⬜ **Plot-specific initialisation** (cut dated from the 1985 stand age) = natural next model step.
> 4. ⬜ The two residual PNGs in `manuscript/figures/` are COPIES of `doublechecks/figures/` outputs —
>    re-copy after re-running `residuals_stand_and_soil.R` / `residuals_time_since_cut.R`.
> 5. ⬜ Stoniness effect may partly be the stone correction of the observed stocks (CHECK comment in .tex).
> 6. ⬜ Ask **Aleksi**: is NID ground-vegetation litter 50.6/66.6 gC/m²/yr above-ground only? (Our plots,
>    2023 cover → 33 above-ground, (J+u)/J = 1.13.) Ask **Boris**: does the Tupek product see sub-tally
>    trees and cleaning residues?
> 7. ⏸ Samuli's flux-tower Rh comparison parked (`doublechecks/rh_vs_stand_age.R`): proper version needs
>    the tower Rh normalised to a reference temperature first.

---

## ✅ (DONE 2026-09-23) LORENZO'S 72 `%LORENZO` NOTES OF 2026-09-21 — kept for reference

> **Tomorrow's work is the revision: answer the 72 notes.** He reviewed in parallel while the
> 2026-09-16 notes were being closed, and his pass reaches the END of the paper — M&M from
> §"The six models" onward, all of Results and Discussion, and the Conclusions. He also rewrote
> prose himself in §"The six models", §"Transient initialisation" and §Climate. **Re-read the
> `.tex` from disk before touching it; never revert his text** ([[lorenzo-edits-in-parallel]]).
>
> Find them with `grep -n "^%LORENZO" manuscript/HIKET_draft_v1.tex`. Distribution:
>
> | notes | section |
> |---|---|
> | **14** | The error model ⭐ the heaviest, and the most conceptual |
> | 7 | Equal external information across models |
> | 6 | Numerical integration · The models bracket the observed rate · **Conclusions** |
> | 5 | Transient initialisation · The models reach the observed stocks |
> | 4 | Model structure is nearly invisible in the fit |
> | 3 | The six models · The shape of the trajectory separates the six models |
> | 2 | The data fix a stock and little else · Mean transit time · Litter input |
> | 1 each | Study plots · SOC basis · Sampling · Reproducing the trajectory · Between individual plots · Other national inventories · Data availability |
>
> ### ⚠ FOURTEEN notes say "Discuss" — they are NOT edit instructions
> Lines 1368, 1373, 1377, 1447, 1457, 1461, 1480, 1484, 1501, 1516, 1530, 1588, 1647, 1696. He is
> asking to be talked through the maths before anything is rewritten: the matrix exponential
> ("is the algorithm really called a twenty-term Taylor series?"), what "integrated out",
> "combine in quadrature", "a triangular solve" and "factorise" mean in practice, and whether the
> offsets are estimated or prescribed. **Answer in the reply, propose wording, let him choose.**
>
> ### Three notes that challenge a CLAIM, not the prose — settle these first
> 1. **σ_total = 0.800 is circular** (line ~1457): *"Those residuals are estimated after a
>    calibration? If so, it is circular."* He is right that the plug-in route looks circular; the
>    defence is that 0.800 is a DELIBERATE step above the measured 0.71–0.74 and that
>    self-consistency was checked afterwards, not assumed. Needs a straight answer and probably a
>    rewritten justification.
> 2. **"Equal external information" is self-serving** (§Equal external information): *"we cannot
>    give it the same amount of information, and we are not doing so"*, and the "same weakly
>    informative prior" phrasing *"seems a self righteous statement"*. The homogenisation argument
>    has to be restated as what it is — a design choice, not a guarantee of fairness.
> 3. **"The spread between models is a diagnostic rather than an estimate of uncertainty"** —
>    *"a bit self reliant … this is not how science should work"*. ⚠ This is the SAME claim that
>    was walked back on 2026-09-03; his note says the current wording still overreaches.
>
> ### What he wants structurally
> - **Stocks and rates should stop being two subsections** — *"two sides of the same coin"*, the
>   rate discussion mirrors the stock discussion and repeats it.
> - **Reaching the data is not a result** — *"It says that they match the rate because they were
>   calibrated on measurements that match the rate. Uninformative."* Compress, do not defend.
> - **The extrapolation IS the interesting part** — *"when we extrapolate we might have surprises.
>   This is leading to important thing for us."* Give it the weight the fit paragraphs lose.
> - **Conclusions**: say openly that this is a **baseline for further model development**; move the
>   "what it means for other inventories" material out of the Conclusions to just before Limitations.
> - ⚠ **"the models describe a gradient the soil does not have"** reads as "the models do not work";
>   retitle. And **Yasso's chemically-defined pools** (§The six models) is *"controversial"* — he
>   wants a short hook to Bosatta & Ågren's continuous-quality theory instead.
> - ⚠ **"Are the rates exact values or tight priors?"** — must be stated precisely; the ICBM anchor
>   pins the slow rate with a very-informative prior and FIXES the fast one, and the text blurs it.

## ✅ CLOSED 2026-09-21 — the four notes of 2026-09-16, plus the climate source

> All four answered and committed (`27016c4`, `fe098ce`, `f405a71`, `e19b5f8`); dispositions are
> dated source comments. Full record: [[lorenzo-notes-round4-20260921]].
> - **War reparations** → written on the NFI's own record (**Korhonen et al. 2021** §3.5.2: stock
>   fell in the 1950s–60s; logging surge before NFI4). ⚠ Korhonen **2024** does NOT carry those
>   sentences. Mitrunen 2024 (4% of GDP, metal goods) ⇒ **"financed by the forest export economy",
>   never "paid in timber"**. Kochetkova 2020 for ceded Karelia.
> - **"How much time"** → pool-specific timescale on Wutzler & Reichstein 2007; **no own e-folding
>   computation** (decided). ⭐ **Sierra, Hoyt, He & Trumbore 2018** now cited in the Introduction to
>   separate the soil's ADJUSTMENT timescale (slow pool, centuries) from the MEAN TRANSIT TIME
>   (outflow, 20–30 yr) — line 623 had wrongly made the transit time the soil's memory.
> - ⭐ **Peltoniemi et al. 2007** (*Silva Fennica* 41(3):575–602, seven-model review, Mikko first
>   author with Liski, Palosuo, Mäkipää, Wutzler) now cited in §"Steady state as a starting
>   condition": the group itself wrote in 2007 that steady state *"is likely to be violated in most
>   applications"* and that matching an equilibrium to a measured stock needs *"a further correction
>   of the parameters and the pools"*. Frame as lineage.
> - **M&M**: "plot-count-weighted mean sampling year" → the plain average with per-year counts;
>   "balanced plot set" → "the plots measured in all three campaigns (310 plots)".
> - ✅ **CLIMATE SOURCE CONFIRMED** (LUKE Data Support, 2026-09-14): FMI daily gridded observations,
>   **10 km grid** (verified on the NetCDF — YKJ cell centres, uniform 10 000 m spacing), mirrored
>   into the LUKE weather database; the only "modification" is the nearest-cell extraction. Cite
>   **Aalto, Pirinen & Jylhä 2016** (FMI's own guidance). That `\gap` is closed; **16 remain**.
> - **Acknowledgements**: LUKE Data Support goes there; **O.-P. Tikkasalo is a coauthor**, not
>   acknowledged. FMI gridded observations are open data.
>
> ⚠ **Bib housekeeping**: duplicate `Sierra2017` and `BayesianTools` removed; `Menichetti2024`
> (ICBM/3, *Carbon Management* 15) was lost in that edit and **rebuilt from Crossref — worth a
> glance**. `library.bib` now 68 entries, bibtex warning-free.

---

## ✅ 0-NOW. ARM B HAS LANDED — the one-factor test PASSED (2026-09-01)

> All six of jobs **972096–972101** COMPLETED. RUN_IDs **`20260831_162458`** (SP1, TP2) ·
> **`_162456`** (TP3) · **`_162457`** (Yasso07, Yasso15, Yasso20). R-hat 1.000–1.008 (0 warnings),
> ESS 1747–10821, all chains 100% finite, memory peaks 12.3–14.2 GB, `events:max 0`.
> Synced (`runs/`, `diagnostics/`, `Data/model_inputs/`), stages 2–4 run, MRT extracted with
> provenance verified (`intrinsic_mrt.rds` stamped `20260831_162457`).
>
> ### Against the pre-registration
>
> | pre-registered | actual | |
> |---|---|---|
> | σ_input −6…−9% | **−5.0…−8.3%** | ✅ |
> | MRT +7…+9% (24.0/27.6/22.5) | **+5.9…+11.0% (24.57 / 27.11 / 22.56)** | ✅ |
> | flux ≈3.14–4.38, all < 4.62 | **3.19–4.34**, all under | ✅ |
> | RMSE ±0.3 | −0.79…+0.43 (SP1 outside) | ~✅ |
> | σ_init unchanged | **−3.9…−6.7% in all six** | ❌ |
> | C_1917 unchanged | **−2…−4 tC/ha** (43.4–55.9 → 41.9–52.3) | ❌ |
>
> ⭐ **Arm A's inference is now DEMONSTRATED**: the MRT drop rode the σ_input rise along the
> near-conserved ridge, not the Liski ramp. Yasso15 regained overlap with its published posterior;
> Yasso20 sits entirely above its published POINT; Yasso07's shortfall stays structural (single ξ).
>
> ### ⚠⚠ Arm B also FITS BETTER — not pre-registered either way
> Calibration R² up in **all six** (0.0166–0.0209); RMSE better in **five of six**; **bias down in
> five** (SP1 3.16→1.52, Yasso20 2.46→1.17, Yasso15 3.87→2.99). Holdout R² unchanged at ~zero.
> **The lower input anchor costs nothing and is marginally preferred by the data.**
> ⚠ **This is exactly when NOT to decide from the table.** The anchor question is Liski's 27% vs
> Lehtonen & Heikkinen's 8% understorey share, and **Aleksi authored both papers**. Decide on his
> answer, before looking at which arm reads better.
>
> ### ⬜ BINDING TASK — source the 1917 stock floor
> It is the **only** axis favouring arm A: C_1917 41.9–52.3 puts **TP2 42.9 · TP3 41.9 ·
> Yasso07 43.2** below it. But the floor is an **unsourced placeholder reading 45 in
> `doublechecks/init_state_plausibility.R` and 40 in the memory** — so "three models fail" is not yet
> a finding, and it flips on a 5 tC/ha invented threshold. Take the **MINIMUM** from the boreal SOC
> literature, whole profile to `z_cap` (Peltoniemi 2004 the natural anchor), then reconcile both.
> ⚠ The pre-run RATE criterion does NOT discriminate — it passes in all six, both arms.
>
> ### Tomorrow: does the story hold?
> Read [[manuscript-readiness-gaps]] first. The writing-side items are unchanged and none of them
> depends on the anchor decision: the `metrics_calib` fix, the T1 caption, the dead over-prediction
> subsection, and the §"multiplier crosses a physical line" rewrite.

---

## 🚀 0-NOW. ARM B IS IN FLIGHT — jobs 972096–972101, launched 2026-08-31 16:25 EEST

> **Launched** from `78f8ac4`; all six RUNNING within 2 s on separate nodes
> (rc4139/rc5140/rc4106/rc4124/rc5125/rc5126). 36 h wall ⇒ hard deadline **04:25 on 2026-09-02**.
> **Expected landing 09:00–15:30 on 2026-09-01** (from iteration rates at 3 h 41 m, not from
> previous wallclocks): Yasso07 ~16.5 h · TP3 ~17.5 h · SP1/Yasso15 ~20.5 h · TP2 ~22 h ·
> Yasso20 ~23 h.
>
> **Launch verified in all six logs:** `Cores per chain: 40` · `5 chains x 50000` ·
> `[PREINIT SHAPE] liski | mean 0.450 | 1950 0.760 | 1970 0.645` (identical to arm A ⇒ ramp truly
> fixed) · `SOC observations: 1205 | obs CV 0.432` · log-normal ·
> `flux_pair J_bar = 2.511 | sigma_input window [0.020, 3.465]` (the WIDE window). The prior centre
> is not echoed to the log — verified by grepping the six prior files on the cluster after `git pull`.
>
> **Health at 3 h 41 m:** memory peaks 9.4–13.4 GB against 156 GB, `events:max 0`, `oom_kill 0`,
> node free 304 GB–1.46 TB. Nothing to watch.
>
> ⚠ **Three traps when checking the logs, all hit in-session:** `grep -i error` matches the
> `[ERROR MODEL]` banner (false positive); the glob `*_9720*.err` **misses jobs 972100/972101**;
> `memlog_*.csv` is in **MB** with a `#` comment line before the header (`tail -n +3`, /1024).
> ⚠ **The `.out` iteration counter RESETS per chain** — read it with the count of `complete (`
> lines in the `.err`, or a job that is ahead looks stalled.
>
> ### When it lands
> 1. rsync `runs/`, `diagnostics/` AND `Data/model_inputs/` (the predictive stage hard-loads the
>    bundle keyed to each RUN_ID)
> 2. stages 2–4 locally: `run_hiket_pipeline.R --skip-calibration`
> 3. `doublechecks/intrinsic_mrt.R` — **check the RUN_IDs it echoes**
> 4. read against the pre-registration below BEFORE forming a view

---

## 🧪 THE PRE-REGISTRATION (written before launch, 2026-08-31)

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
stated primary question of the forward experiment. ✅ **ANSWERED 2026-09-02** — executed as
published-Yasso vs ours, forward only (F15 + A1); the hypothesis held for Yasso07 and NOT for
Yasso15/20. See `manuscript/HIKET_storyline_v2.tex` §"The consequence, tested".

⚠ **Both scripts are DEMONSTRATIONS, not projections.** Single-exponential transient on the bulk
MRT (a multi-pool system is a sum of exponentials), a step change rather than a trajectory, and a
common J̄ for absolute scale. A proper version needs forward runs — which is precisely the
forward experiment — ✅ **DONE 2026-09-02**, though not in the form sketched here. The
restrictive-vs-permissive σ_input design was **abandoned** ([[sigma-input-physical-bounds]]); the
experiment that ran instead contrasts the *published* parameterisation with ours under ramped IPCC
scenarios, forward only. Compared on RMSE distributions, climate reactivity, and trajectories to
2084. Publishable either
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
   → `F14_mrt_ridge.png`, and its caption now lives in
   `manuscript/figure_captions.tex` (rescued 2026-09-03 when the storyline note was deleted), with
   **two boxes**: a `methodnote` explaining the two layers and a `provisional`
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
