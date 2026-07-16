# HIKET — Revision plan (round 1)

**Status:** planning, pre-recalibration. Created 2026-07-15.
**Trigger:** first coauthor comments on `HIKET_storyline_note` (annotated PDF in
`manuscript/revisions/HIKET_storyline_note_sl.pdf`) + two structural issues raised in
discussion (σ_input height; the "simplicity" story is confounded by prior information).

**The gate.** Several changes alter model structure or priors, so a full six-model
recalibration on Roihu is required. *Everything under §2 must land and be verified
before that job is launched.* §3 (framing/discussion) can follow the run. Nothing here
touches the trusted decomposition rates by *fitting* — where kinetics move, they move to
**independent external calibrations** (Ultuna/ICBM, litterbag), not to the SOC data.

---

## 1. Coauthor comments — inventory & disposition
*(Round-1 comments from **Samuli Launiainen**, coauthor. His candid margin notes — e.g.
"ridiculously flawed" — are internal; keep manuscript prose measured. Aleksi Lehtonen away
until ~mid-Aug 2026; Mikko Peltoniemi approachable meanwhile for litter/stock questions.)*

| # | Pg | Comment (paraphrased) | Disposition |
|---|----|------------------------|-------------|
| 1 | 1 | "nice tension" (opening) | — keep |
| 2 | 2 | Steady-state w/ constant pre-1990 litter can't recover post-1990 accumulation; pre-assumptions "ridiculously flawed" — societally relevant for NGHGI | §3 discussion: frame as the headline non-equilibrium result |
| 3 | 2 | How does this compare with **Peltoniemi et al.** (basis of the steady-state assumption)? Likely they had no repeated stocks → couldn't see the dynamics | §3 discussion: add explicit comparison |
| 4 | 2 | Add phrasing: "recognizing the dynamical nature of forest systems, and the interplay between stand dynamics and soil C stocks" | §3 wording |
| 5 | 2 | **Julius' inversion** of litter inputs from stocks (agric) — do similar for forest at regional/national scale to back the non-steady-state claim | C4c: σ_input external anchor; **possible future work (optional), not now** |
| 6 | 4 | Are T & moisture sensitivities **fixed across models**, or part of "structure"? | **ANSWERED:** C2 (TP3 climate on all pools) + the xi-normalization (C1) — sensitivities are calibrated on a common xi *form* & common climate reference (Ultuna); document |
| 7 | 5 | Is the measured stock dynamic realistic? 1985–2006 sink ≈ 200 gC/m²/yr is very large. How would a **non-linear pre-1985 litter shape** change dynamics? | §2 C3 (growing-stock spin-up) + **D3 (stock QC — CORROBORATED: 184 gC/m²/yr is a mineral artefact)** |
| 8 | 10 | RF residuals: which **direction** does each feature push residuals? Do they rise with basal area? | **§7 ACTION** (post-recalibration): signed residual-driver analysis — sign/slope of residual vs each top driver, esp. basal area |
| 9 | 10 | Transfers already constrained by litterbag data; complexity doesn't help. Reframe deliverable (here or paper 2): (i) re-initialized Yasso07 for NGHGI, (ii) steady-state vs transient over reporting period, (iii) future under MELA scenario | C1 (kinetics anchored) + **D4 (one paper; policy downplayed; operational/MELA study = possible future, not now)** |
| 10| 10 | Sticky note: since transfers can't be calibrated, **keep them fixed and calibrate only input (init) + T & P response** | §2 C1 — adopted: kinetics anchored externally (derived), only σ_input/σ_init/climate calibrated |
| 11| 11 | SOC-change map ±0.6 tC/ha/yr = 60 gC/m²/yr. Compute **regional/national means**; document map construction | **§7 ACTION:** compute regional + national mean SOC & 10-yr change; caption: construction = simple kriging |
| 12| 16 | Baseline litter input? Multiplier ~2.3 to reach current stocks. **Double-check Biosoil & Komeetta stocks — 10 kgC/m² seems high** | §2 C4 stock QC (cheap, do first) |

---

## 2. Changes to implement BEFORE recalibration

### C1 — ICBM-anchored kinetics for the simple models (SP1/TP2/TP3)
**Why.** The intercomparison is "equal degrees of freedom" but *not* equal prior
*information*: Yasso kinetics carry litterbag calibration (tight priors); the simple
models' rates are loosely pinned (log-SD 0.5) around borrowed centres. So "complexity
buys no skill" is confounded with "the simple models are under-informed." Fix: give the
simple models external kinetic constraint of comparable weight, from the **Ultuna bare
fallow ICBM calibration** (Andrén & Kätterer 1997): `k1 = 0.8`, `k2 = 0.00605`,
`h = 0.13` yr⁻¹, reference `r = 1`. Reframe as a *result*: external constraint is
available in proportion to how well a model's **timescale structure** matches an existing
calibration.

**DESIGN — finalized 2026-07-16 (supersedes the "derived per-draw rate" scheme below).**
The guiding principle became *maximum homogeneity with how Yasso is already
parameterized*. Yasso does **not** calibrate its decomposition rates — the whole `a`-vector
is a **fixed constant** injected in `assemble_model_params` (`run_Yasso07_*`:
`FIXED_RATE_NAMES <- c("alpha_A","alpha_W","alpha_E","alpha_N","p_H","alpha_H")`); what
Yasso leaves *free* is the 12 **lateral** inter-AWEN transfer fractions, climate (β),
woody size, and the two σ's. So the homogeneous simple-model parameterization splits into
two classes:
- **Intrinsic decomposition rates → externally anchored, NOT fit to SOC data** (the analog
  of Yasso's fixed `a`-vector). The *fast* rate (`k1` → SP1 fast component, TP2/TP3
  `alpha_A`) is **fixed** (litterbag-grounded, well transferable). The *slow* rate
  (`k2` → `alpha_H`, TP3 `alpha_S`; and SP1's slow-dominated single rate) carries a **very
  informative prior** centred at the ICBM value — free but tightly pinned, so the SOC data
  *can* nudge it and the posterior width becomes a **transferability diagnostic** for the
  one weak link (`k2` was optimized on *arable* Ultuna, the least transferable to boreal
  forest). This hybrid = "fix what's solid, let the data speak a little on what's uncertain."
- **Humification fractions → FREE** (the analog of Yasso's free lateral fractions, *not* of
  its fixed into-humus `p_H`; the simple models have only sequential humification-into-a-
  slower-pool flows, which is the partitioning Yasso calibrates). `p_H` (TP2) and `p_S`,
  `p_H` (TP3) keep the **common Tier-2 logit prior (SD 0.4)**, **re-centred at ICBM
  `h ≈ 0.13`**. Fixing them would leave the simple models with *zero* partitioning freedom
  while Yasso keeps twelve — i.e. *more* constrained than Yasso, the opposite of homogeneity.

**No per-draw derivation, no TP3 single-knob reparam.** Rates are fixed/tight **constants**
computed **once** at the β prior centre (Yasso does not re-normalize its fixed rates each
draw; neither do we). The ICBM `r=1` reference sits at our `xi_Ultuna ≈ 0.94` (Yasso07 xi
form, β-centre, half-range `T_amp≈10`), so the fixed constant is `alpha = k_ICBM /
xi_Ultuna` applied **once** — a ~6% offset, not a per-draw computation. Effective rate at a
plot is then `alpha · xi_FI(β)`, exactly as Yasso scales its fixed `a`-vector by the free xi.

**VERIFIED 2026-07-16 — `doublechecks/icbm_anchor_sanity.R`** (447 real plots, `σ_input=1`,
observed SOC median 71, IQR 55–90). The ICBM anchor reproduces boreal stocks at a *physical*
input multiplier — no escape hatch needed:
- **SP1** bulk MRT 23.8 yr → stock 59.5, `σ_input`≈**1.18**.
- **TP2** MRT 21.5–23.8 yr → 55–60, `σ_input`≈**1.18–1.36**.
- **TP3** *only if `alpha_S` is pinned SLOW (k2-scale), not intermediate.* At
  `alpha_S = k2/(1−p_H) ≈ 0.0070` (S+H = ICBM "old" subsystem) TP3 gives the **identical**
  23.8-yr bulk MRT and `σ_input`≈**1.18** as SP1/TP2 → the "all three share one bulk MRT"
  property emerges. (Intermediate `alpha_S`≈0.05–0.1 starves the cascade — two sequential
  0.13 humification steps → MRT ~6 yr, `σ_input` 4–5×; this was a design-error catch.)
- **C2 (xi on H)** moves TP3 stock by ~0.6 tC/ha — confirmed a fairness/consistency fix, not
  a stock lever.
All models land at `σ_input` 1.05–1.36 = right inside the physical understorey correction
(D2, ~1.3): the kinetics close the σ_input escape hatch, and the D2 re-centre is corroborated.

**Grounded ICBM values (verified against the source paper, 2026-07-15).** Andrén, O.
& Kätterer, T. (1997), *ICBM: the introductory carbon balance model for exploration of
soil carbon balances*, Ecological Applications 7(4):1226–1236. From Table 1 + text:
`k1 = 0.8 yr⁻¹` (young/active; half-life 0.87 yr; assumed for all litter, from
litterbag/¹⁴C — single-exponential k ≈ 0.65 at h=0.13, r=1); `k2 = 0.00605 yr⁻¹`
(old/humus; **optimized from Ultuna bare-fallow topsoil C**, with h=0.13, Y₀=0.3 kgC/m²,
r=1; the paper recommends `k2 ≈ 0.006` as a global constant); humification `h = 0.13`
(bare fallow) / `0.125` (vegetated arable straw treatments) / 0.25 manure / 0.34 sludge;
reference `r = 1` = central Sweden (MAT +5.4 °C, precip 520 mm), bare-fallow r = 1.32.
**Caveats to state in the paper:** (i) these are *arable clay-soil* values (Ultuna) used
as the best-constrained available anchor for boreal *forest* soil — an approximation;
(ii) `r = 1` is the Swedish climate reference, so when we apply our own `xi` the rate
anchor's climate reference must be kept consistent (the coauthor's "are T/moisture
sensitivities fixed across models?" point — C2); (iii) `h` is arable-topsoil
humification, forest may differ. Bulk MRT at reference = `1/k1 + h/k2 = 1.25 + 21.5 ≈
22.7 yr` (→ SP1 rate ≈ 0.044); slow-subsystem MRT `1/k2 ≈ 165 yr` (→ TP3 anchor). Full
citation in memory `icbm-parameters`.

**xi-normalization — ⚠ SUPERSEDED 2026-07-16 by the DESIGN block above.** The per-draw
`alpha = k/xi_Ultuna(β)` derivation below was replaced by a **one-time** fixed-constant
offset (rates no longer leave the parameter vector as "derived"; fast = fixed constant,
slow = very-informative-prior free param, both set once at the β centre). Retained here for
the rationale on why the ICBM `r=1` reference must be placed on our `xi`. Original note:
ICBM's rates are defined at `r=1` = *central-Sweden* climate; our `xi` (Yasso07 form:
`mean(exp(β1 T + β2 T²)) · (1 − exp(γ P/1000))`) has **no normalization constant** — the
absolute scale is absorbed by the rates. So the ICBM rate cannot be dropped in as
`alpha = k` directly; it must be placed at Ultuna's climate on *our* `xi`:
**`alpha = k_ICBM / xi_Ultuna`**, with `xi_Ultuna = compute_xi_yasso07(T_mean=5.4,
T_amp≈11, P=520, β)`. Then Finnish climate scales it: effective rate `= k · xi_FI/xi_Ultuna`.
- **Per-draw, not a fixed constant:** β is calibrated (the "T&P response"), so `xi_Ultuna`
  moves with each draw ⇒ compute `alpha = k/xi_Ultuna(β)` inside the likelihood. This makes
  the simple-model **rates DERIVED, not free** (α leaves the parameter vector; only `p_H`
  stays free in TP3), and β now means "climate sensitivity *relative to Ultuna*." Matches
  the coauthor's "keep kinetics fixed, calibrate T&P response."
- **Numbers at the current β centres (0.095, −0.00014, −1.21):** `xi_Ultuna ≈ 0.976`
  ⇒ `alpha_A ≈ 0.82`, `alpha_H ≈ 0.0062`, `SP1 α ≈ 0.045`. Correction is small (~2.4%) now
  but stays exact as β moves. Self-check: `alpha_A · xi_FI(south) ≈ 0.85` (≈ ICBM 0.8,
  climate-scaled). Reproduce: `doublechecks/` xi-normalization script (to write).
- **Simple-model-specific:** Yasso needs no such step (its rates + `xi` were co-calibrated
  on the same data — internally consistent); only the ICBM graft requires it.
- **Open (minor):** confirm `T_amp` convention (monthly-mean amplitude vs half-range) so
  Ultuna `xi` matches the plot `xi` definition; ~1% on `xi_Ultuna`, not structural. This is
  also the concrete answer to the coauthor's pg-4 "are T/moisture sensitivities fixed across
  models?" — they are calibrated, but on a common `xi` form with a common climate reference.

Per-model parameter status under the finalized design (centres shown at reference; the
one-time `/xi_Ultuna` offset is applied to the rate constants):

| Model | FIXED | Very-informative prior (slow, free-but-pinned) | FREE (logit SD 0.4, centre 0.13) | Old (loose) centres |
|-------|-------------|-----------------------------------------------|----------------------------------|---------------------|
| **SP1** (1 timescale) | — | `alpha ≈ 1/22.7 ≈ 0.044` (bulk MRT; slow-dominated → pinned, not fixed) | — (no split) | 0.09 (≈2× too fast) |
| **TP2** (2 timescales) | `alpha_A = 0.8` | `alpha_H = 0.00605` | `p_H` | 0.73, 0.0015, p_H 0.028 |
| **TP3** (3 timescales) | `alpha_A = 0.8` | `alpha_H = 0.00605`, `alpha_S ≈ k2/(1−p_H) ≈ 0.0070` (both k2-scale = ICBM "old" subsystem) | `p_S`, `p_H` | 0.73, alpha_S 0.10, alpha_H 0.0015, p_S 0.028, p_H 0.50 |

Free set is then **identical across SP1/TP2/TP3** for everything non-kinetic: `{β1, β2, γ,
σ_input, σ_init}`. Complexity adds only free *partitioning* splits (SP1 0 → TP2 1 → TP3 2),
mirroring Yasso's 12 free lateral fractions — the honest "complexity = more partitioning
DoF" ladder. SP1's single rate is slow-dominated (`h/k2 = 21.5` of the 22.7-yr bulk MRT),
so it inherits the `k2` uncertainty → very-informative prior rather than a fixed constant.

**TP3 constraint mechanism → RESOLVED 2026-07-16 = option (b), simplified.** No
reparameterization. `alpha_A` fixed; `alpha_S` and `alpha_H` are **very-informative-
prior** free params **centred at the k2-scale ICBM "old" subsystem** (`alpha_S ≈
k2/(1−p_H_centre) ≈ 0.0070`, `alpha_H ≈ k2`); `p_S`, `p_H` **free** (logit SD 0.4, centre
0.13). The ~165-yr subsystem MRT is thus set by the rate *centres*, not pinned by a hard
identity — the free splits let it drift (partitioning freedom, like Yasso's fractions). The
sanity check confirms the centres give the shared ~24-yr bulk MRT. This keeps the "fix
intrinsic rates, free the fractions" rule uniform with SP1/TP2 (no TP3-only machinery).

**Consequences / cross-checks:**
- All three simple models then share **one common ICBM bulk MRT** (≈22.7 yr at reference)
  and differ *only* in internally non-identified structure → the ideal "does complexity
  buy skill" design. Verify: SP1/TP2/TP3 reference steady-state stock identical for a
  given J.
- **Side benefit on σ_input (not a fitting target):** current TP2/TP3 humification is
  ~5× below Ultuna (`p_S=0.028` vs `h=0.13`), starving the slow pools → bulk MRT ≈ 11 yr
  → forces σ_input up. ICBM humification roughly **doubles** bulk MRT (≈23 yr) →
  ~**halves** required σ_input. This falls out of using the right external constraint.

### C2 — TP3 climate response on ALL pools  *(DECIDED 2026-07-15)*
**Why.** Currently `k_H = alpha_H` (humus pool has **no** xi) while A and S do — see
`tp3_wrapper_transient.R` header lines 9–13. This makes TP3's climate treatment
structurally different from the other models → confounds the comparison (coauthor #6).
**Change:** `k_H = alpha_H * xi` everywhere H turns over.
**Sites (`tp3_wrapper_transient.R`):** `.tp3_step` (`kH` argument / within-year eq.
`Hss`), `tp3_transient_init` (line ~141 `kH <- alpha_H` and the C_init `H`),
`.tp3_steady_state` (`H_ss` at line ~106, add `* xi_mean`). Pure R — no Fortran recompile.
**Note:** with xi on H, the C1 165-yr subsystem anchor is at reference climate; boreal
xi<1 now lengthens the *whole* subsystem consistently (cleaner than before). Add a
verification vs the exact matrix-exp (cf. `doublechecks/test_tp3_exact.R`) — the
lower-triangular eigenvalues change but structure is preserved; watch the
divided-difference "coincident eigenvalue" guard.

### C3 — Historical (non-linear) input interpolation for transient init
**Why.** The 1917→1985 pre-run currently ramps litter **linearly** (`frac <- (i-1)/(n_pre-1)`;
`J <- J_1917 + (J_1985-J_1917)*frac`) — see `sp1/tp2/tp3_wrapper_transient.R` and the
Yasso transient inits. Coauthor #7 + the storyboard "growing-stock-shape" note: replace
the linear ramp with a shape derived from the **historical growing-stock reconstruction**
(Korhonen 2024 / F10b: depleted, near-stationary into the 1970s, then a sustained rise).
This is the physically-honest pre-observation forcing and directly tests whether the
+26 tC/ha 1985 over-prediction shrinks.
**Scope:** all six transient inits. Wire the growing-stock time series → per-year litter
scalar → replace the linear `frac`. Keep `J_1917`/`J_1985` endpoints as calibrated
anchors; only the *shape* between them changes.

### C4 — σ_input robustness & reframe
Rates are trusted, so σ_input height is a **discussion result**, not a thing to fit away.
Make the estimate robust and interpretable:
- **C4a Stock QC (do first, cheap):** verify Biosoil/Komeetta depth & layers; 10 kgC/m²
  is high (coauthor #12). If the target is biased high, σ_input shrinks for free.
- **C4b Re-centre the σ_input prior on a physical expectation >1**, not on 1: the Tupek
  product is a *tree*-litter model. **What it OMITS (corrected 2026-07-16 against the source
  + author confirmation, B. Tupek):** (i) **understorey/ground-vegetation (incl. moss)
  litter** — the dominant missing term, ~15–35% of total boreal litter (larger N); (ii)
  **mycorrhizal mycelial turnover / root exudates** — large but very uncertain; (iii)
  **aboveground tree mortality (deadwood/CWD from whole-tree death)** — a genuine omission,
  but a *small* soil-input term in managed Finnish forests (harvest removes stems, natural
  mortality suppressed) that also enters measured soil C only partly and slowly. **What it
  INCLUDES — roots.** Fine roots (in `nwl`) AND coarse roots (in `fwl`) are covered, along
  with foliage, branches, stem bark, stumps: J is the full living-tree biomass-component
  turnover, not aboveground-only. (The earlier "omits fine-root turnover" here was WRONG —
  Boris confirms roots are in his estimates; only aboveground mortality and understorey are
  out. See `HIKET_data_preparation.Rmd` size-class table: `nwl`=foliage+fine roots,
  `fwl`=branches+coarse roots+stem bark, `cwl`=stumps.) Build an independent estimate of the
  missing fraction (understorey-dominated) and centre σ_input there (~1.3) so the posterior
  is read against *expected total litter*. Converts σ_input from a fudge into an estimate of
  unmodelled input. **NB the Zenodo DOI 10.5281/zenodo.19736499 is not yet registered
  (embargoed 2026 deposit) — this comparison is against our local derivation + the author's
  statement; re-verify against the public M&M when it lands.**
- **C4c (paper-2 / optional):** external input anchor via stock→input inversion or
  NPP-allometric total litter (coauthor #5, "Julius inversion").

### C5 — Down-weight the suspect 1985 (VMI8) stocks via inflated observation variance
*(DECIDED 2026-07-15 — this is how we handle the D3 mineral-soil artifact; supersedes
dropping/ignoring VMI8.)* The 1985 mineral measurement is systematically suspect (D3), so
rather than drop it we **inflate its observation error** and let the trusted 2006/2024
anchors + the transient dynamics carry the initial state.
- **Mechanism (trivial — hook exists):** likelihood is `dnorm(soc_obs, SOC_hat, sd_vec)`
  with `sd_vec = SOC_hat * sigma_obs_fixed` (`calibration_engine_transient.R` ~line 132).
  Add a per-obs factor: `sd_vec = SOC_hat * sigma_obs_fixed * meta$sigma_infl`, populate
  `sigma_infl` when building plot meta. **Key on campaign year == 1985**, NOT `meta$is_first`
  (not every plot has VMI8).
- **Which points → the whole 1985 campaign** (systematic method bias, not random plots;
  per-plot flagging would wrongly penalize large real changes).
- **How much → ~2× (`sigma_1985 ≈ 0.9–1.0`)**, grounded in the ~30% (≈0.36 log) 1985-low
  discrepancy vs base `sigma_obs ≈ 0.47`. Puts the discrepancy inside 1σ, makes 1985 a weak
  anchor. **Fixed factor, not calibrated** (a calibrated 1985 offset is non-identifiable
  against `sigma_init`). Run a **sensitivity (1.5×/2×/3×)**.
- **Clean:** no double-count with `sigma_init` (that's model-state spread; this is
  measurement distrust). Complementary with C3. Limitation: total-SOC likelihood, so the
  good 1985 *organic* info is also down-weighted (minor; a layer-resolved likelihood could
  isolate mineral later).

---

## 3. Manuscript / discussion changes (no recalibration)
- **Thread A reframe:** informational-asymmetry / timescale-matching as a *result*
  ("complexity buys no skill *given equal external constraint*"), not an apology.
- **σ_input:** present as input-completeness + stock QC, tension resolved by external
  anchoring; report effective flux inside the NPP envelope.
- **Peltoniemi et al.** comparison (why steady state was assumable pre-repeated-stocks).
  **Framing note: Aleksi Lehtonen (lead of the steady-state inventory paper, Lehtonen et al.
  2016) is a HIKET coauthor** — present the steady-state lineage collegially (the group
  evolving its own then-reasonable method with new repeated-stock data), NOT as a critique;
  keep prose measured (no "flawed"). Lehtonen = natural owner for D3 stock QC & D2 understorey.
- Add the "dynamical nature / stand–soil interplay" phrasing (coauthor #4).
- **RF residuals:** report *direction*/slope vs basal area, not just importance (#8).
- **Maps:** regional + national means; methods note on kriging construction (#11).
- **Policy tone → DOWNPLAY to a side-thread in ONE paper (see D4; NOT a split):** keep policy
  as motivation, remove recommendations; render it as a minor thread (like `[Bayes]`/`[History]`)
  riding on the existing wall-to-wall product + forecast-divergence spread; rename Landing
  "policy deliverable" → "a national estimate and its uncertainty"; strip prescriptive language
  from abstract + discussion; add explicit "no sink-magnitude / inventory-method
  recommendations" line. A full operational/MELA study = possible future successor, not now (#9).

---

## 4. Implementation order, verification, launch
1. **C4a stock QC** — DONE (local sanity check; see D3). Authoritative QC = someone else.
1b. **ICBM anchor validation** — DONE 2026-07-16 (`doublechecks/icbm_anchor_sanity.R`;
   see C1 VERIFIED block). Confirms σ_input lands 1.05–1.36 (physical) and `alpha_S` must
   be k2-scale.
2. **C1** ICBM kinetics (SP1/TP2/TP3 `Prior_specs/*_priors.R` + each `assemble_model_params`):
   **fast rate fixed** (`alpha_A`/SP1 fast) as a constant like Yasso's `fixed_rates`;
   **slow rate very-informative-prior** free (`alpha_H`, TP3 `alpha_S`; SP1 single rate);
   **humification fractions free**, logit SD 0.4, re-centred at `h≈0.13` (`p_H`; TP3 `p_S`,
   `p_H`). One-time `/xi_Ultuna` offset on the rate constants; NO per-draw derivation, NO
   TP3 reparam.
3. **C2** TP3 climate-on-all-pools (`tp3_wrapper_transient.R`).
4. **C3** historical input interpolation (all six transient inits).
5. **C4b** σ_input prior re-centre (≈1.3) — corroborated by 1b.
6. **C5** 1985 observation-variance inflation (~2×) in `calibration_engine_transient.R`
   + plot-meta `sigma_infl`.
7. **Verify (local `doublechecks/`):** ICBM steady-state stock ≈ shared across SP1/TP2/TP3
   (DONE, 1b); TP3 exact integrator still matches matrix-exp with xi on H; interpolation
   shape sanity; slow-rate very-informative prior pushforward (no forward blow-ups); flux
   stays in NPP envelope.
8. **Sync + recompile** on Roihu (recompile only if any `.f90` changed — C1–C5 are R-only,
   so no Fortran rebuild expected; still confirm `.so` currency per CLAUDE.md).
9. **Launch** six-model calibration on Roihu; watch `Cores per chain: 40`.
10. Post-run: stages 2–4, refresh figures/tables/NextGenC report, then §3 edits.

## 5. Decisions — sharpened for tomorrow's session

**D1. TP3 parameterization → REVISED 2026-07-16 (supersedes the 2026-07-15 single-knob
reparam).** Under the finalized "fix intrinsic rates, free the fractions" rule (homogeneous
with Yasso, which fixes rates + frees its lateral fractions), TP3 needs **no
reparameterization**: `alpha_A = k1 = 0.8` **fixed**; `alpha_H ≈ k2` and `alpha_S ≈
k2/(1−p_H_centre) ≈ 0.0070` are **very-informative-prior** free params (both k2-scale = ICBM
"old" subsystem); `p_S`, `p_H` **free** (logit SD 0.4, centre `h≈0.13`). The ~165-yr
subsystem MRT is set by the rate *centres*, not pinned by a hard identity — the free splits
let it drift (the partitioning DoF, exactly like Yasso's fractions). Rationale for dropping
the single-knob reparam: it *fixed* `p_S` and derived `alpha_S`, leaving TP3 with only one
free kinetic DoF and the simple models *more* constrained than Yasso; freeing both splits
restores parity. The k2-scale `alpha_S` is essential (verified 1b: an intermediate `alpha_S`
starves the cascade). (Interacts with C2: with xi on H, the k2-scale centres are the
reference-climate rates; boreal xi<1 lengthens the whole subsystem consistently.)

**D2. σ_input prior centre → SET ≈ 1.3 (log), literature-grounded 2026-07-15;
CORROBORATED 2026-07-16.** The ICBM-anchor sanity check (1b) independently lands σ_input at
**1.05–1.36** across all three simple models to reach observed SOC — the same window the
missing-understorey argument predicts. Two independent lines (litter-completeness physics +
the kinetic anchor) agree on ~1.3, so the re-centre is well-founded, not a tuning.
The litter product is **tree-litter** (full living-tree biomass-component turnover:
foliage, **fine roots**, branches, **coarse roots**, stem bark, stumps — so **roots ARE
included**, confirmed by B. Tupek and the size-class table). **Not** represented, in order
of magnitude: **understorey/ground-vegetation (incl. moss)** — the dominant term;
**mycorrhizal mycelial turnover / exudates** — large but very uncertain; and **aboveground
tree mortality (deadwood/CWD)** — genuinely omitted but a *small* soil-input term in managed
Finnish forest (harvest removes stems; natural mortality suppressed; slow/partial entry to
measured soil C). Boreal-Finland litterfall studies: understorey ≈ **15% (south) to 33%
(north)** of aboveground litter (up to ~50% of total in Lapland), i.e. understorey/tree
ratio ≈ 0.2–0.35 → a physically-expected multiplier ≈ **1.2–1.35** from ground vegetation
alone; mycorrhiza and the small mortality term add a further (mostly unquantified) upward
push. **Recommend centring σ_input at ≈ 1.3** (understorey-dominated; keep the flux-pair NPP
bound on top). Caveat: the fraction rises strongly S→N, so a single national
multiplier cannot capture it — a limitation that motivates stratifying inputs by site
class (future work). *To pin exact citation:* the 15/33% figures are from the boreal
site-type-gradient litterfall literature — confirm the primary ref before the paper.

**D3. Stock QC → LOCAL SANITY CHECK DONE 2026-07-15 (authoritative QC = someone else).**
Campaign means (calibration target RDS): VMI8 63.3, Biosoil 102.0, Komeetta 105.3 tC/ha.
Implied sinks: **VMI8→Biosoil 184 gC/m²/yr** (implausibly large; typical boreal soil sink
10–50), Biosoil→Komeetta **18 gC/m²/yr** (plausible). Layer breakdown of the **paired**
1985–2006 data (`Data/SOC/soilC1985_2006.csv`, same 414 plots both years, both covering
organic + 0–40 cm mineral):
- **Organic layer CONSISTENT:** 20.1 (1985) vs 19.7 (2006) tC/ha — comparable, reassuring.
- **Mineral soil ~DOUBLES:** 0–40 cm = 34.9 (1985) vs 63.1 (2006); incl. **deep 20–40 cm
  11.6 → 23.3** — deep mineral C *cannot* double in 21 yr → this is a between-campaign
  **mineral-soil method/measurement inconsistency**, not real accumulation.
**Implications (big):** (a) the 1985→2006 leg of the accumulation narrative is largely
artefactual; the 2006→2024 near-plateau is real-ish. (b) The **+26 tC/ha "1985
over-prediction" likely means VMI8 mineral C is biased LOW, not the model high** — flips
that interpretation ([[transient-phase-open-question]]). (c) σ_input is partly chasing the
(possibly inflated) 2006/2024 level, entangling D2 with this. **For the authoritative QC:**
focus on 1985-vs-2006 *mineral* sampling (bulk density, coarse-fragment correction, the
0–5/5–20 vs 0–10/10–20 increment split, 1 m extrapolation — cf. `Data/SOC/extrap_*` figs);
the organic layer looks fine. We proceed "with what we have" but the paper must state this
caveat prominently.
**External support (Peltoniemi et al. 2004, GChB — now in `literature/`):** their careful
chronosequence found **mineral soil C does not change with stand age** (organic layer does),
and a careful total of **6.8 kgC/m² to 1 m** (≈ our VMI8 ~63 tC/ha; makes Biosoil/Komeetta
~100+ look high). Both independently corroborate that the VMI8→Biosoil *mineral* jump is a
measurement artifact, not accumulation — cite them for the D3 caveat.

**D4. Paper scope → DECIDED: ONE PAPER, policy downplayed to a side-story/thread
(2026-07-15; reversed the earlier "split" — that was premised on policy being a FULL
deliverable, which the downplay removes).**
- **Single paper.** Spine = the methodological contribution: initialization problem as sole
  protagonist; transient init under a fair frame across all six models; proof of concept
  (accumulation→saturation); structural intercomparison; non-identifiability + σ_input reveal;
  wall-to-wall national product (current stocks + between-model uncertainty). Round-1 revisions
  C1–C5 serve this.
- **Policy = a light THREAD, not a section or protagonist** — woven through like the existing
  `[Bayes]` and `[History]` threads. Keep policy as *motivation* (incipient saturation, the
  NGHGI must represent it, initialization changes the answer); cut *recommendations* (no
  sink-magnitude claims, no "the inventory should adopt X"). Justified by the data: the sink
  magnitude leans on the artefactual 1985→2006 rise (D3), the future is an un-adjudicable
  forecast divergence, σ_input/σ_init carry real uncertainty.
- **Why one paper now works cleanly (no new machinery):** the cautious policy thread rides on
  material ALREADY in the paper — the wall-to-wall product, the forecast-divergence panel
  (ensemble spread to 2084), and the "honest output is the spread" framing. No MELA runs, no
  operational re-init, no second paper needed.
- **Landing reframed as methodological:** (1) initialization matters, (2) a calibrated national
  estimate + its irreducible structural uncertainty, (3) narrow it via better inputs/history,
  not more structure. Wall-to-wall product = research output WITH caveats (incl. stock-QC),
  not an inventory prescription.
- **Kept SEPARATE (decide later, independently):** the locally-stratified-calibration
  extension is a *methods* question (Yasso-only, RF-guided), NOT policy — it may stay an
  optional `\ifstrat` strand or become a future paper on its own merits; the policy downplay
  says nothing about it. A full operational/MELA study also remains a *possible* future
  successor, not a commitment, and would need the clean authoritative stock data first.
- **Authorship/scope stance (user, 2026-07-15):** the user is **open to authoring** a future
  policy paper — **mildly interested ("mildly yes"), just not a priority** — so it's a genuine
  possible successor, not something to offload. But it stays OUT of this paper: policy here is
  downplayed to a thread, and the figures already downplay it enough. The wall-to-wall product
  is "just a simple kriging" (also the answer to coauthor #11 "how were maps built"). So:
  record the tone for the write-up, keep the policy-paper door open, don't over-invest now.
- **Prose to-do (§3):** rename Landing "the policy deliverable" → e.g. "a national estimate
  and its uncertainty"; strip prescriptive language from abstract + discussion; add an explicit
  line that the paper makes no sink-magnitude or inventory-method recommendations; render
  policy as a minor thread, not the payoff.

---

## 6. Deferred figure items (carried from the retired storyboard)
The `FIGURE_STORYBOARD_DRAFT` served its purpose (F1–F14 built, OCAR schema now the
manuscript skeleton). Its substantive content is already in the storyline note + main
manuscript; only these minor figure to-dos remain and are preserved here so the
storyboard can be retired:
- **F13(b) labile-fraction map:** add the between-Yasso *spread* of the labile share
  (the field is model-structural; a spread panel quantifies that).
- **S7 parameter marginals:** currently the two identified σ's; extend to the full
  parameter set when wanted.
- **F8 density estimator:** ridge panels use a Scott full-bandwidth-matrix KDE; swap to a
  plug-in selector (`ks::Hpi`) if a reviewer wants the rigorous bandwidth (one-line change).
- **Growing-stock spin-up shape:** = C3 above (now promoted from "idea" to a planned change).
- **OCAR schema** is a working map; collapse to standard IMRaD at cleanup.

(Storyboard `.tex`/`.pdf` left in place; safe to delete once you're satisfied nothing
else is needed — all unique content is captured here.)

---

## 7. Comment-driven figure/analysis actions (post-recalibration)
These need actual analysis/figure work (not just prose) but depend on the new posteriors,
so they belong in the post-run figure pass, not the pre-recalibration batch.

- **[F11 — coauthor #8, priority] Signed residual-driver analysis.** The RF heatmap shows
  *importance* but not *direction* — the coauthor asks which way each driver pushes the
  residual, and specifically whether residuals rise with basal area. **Action:** for each
  model, compute the signed relationship of the plot residual vs each top driver (sign +
  slope of a univariate fit, or Spearman ρ; and/or an RF partial-dependence curve for the
  top 2–3, esp. `basal_area_85`). **Represent** either as a signed companion to F11 (e.g.
  diverging colour for sign alongside importance, or a small "residual vs basal area" panel
  with fitted slopes per model). **Expected & why it matters:** residuals likely *increase*
  with basal area (under-prediction at productive, high-input plots) → basal area is a
  productivity/input proxy, so this *reinforces* the thesis that the missing signal is in the
  **inputs**, not the kinetics ([[residual-rf-basal-area]]). Scripted in the residual stage
  (`run_residual_analysis.R` / `run_multimodel_comparison.R`).
- **[F13 — coauthor #11] Regional + national means.** Add regional (e.g. latitude terciles /
  regions) and national mean SOC and mean 10-yr change to accompany the maps; state the map
  construction in the caption = **simple kriging** (per user). NextGenC report scripts
  (`Reporting/NextgenC_report/build_soc_maps.R`).
- **[Peltoniemi — coauthor #2/#3] Lit action → DONE 2026-07-16.** All three papers in
  `literature/` (catalogue `literature/README.md`; memory `steady-state-convention-lineage`):
  **Peltoniemi et al. 2004 (GChB)** — read in full = the DIRECT ANTECEDENT (it *identified*
  the non-equilibrium init problem in 2004: "long-term trend cannot be distinguished … unless
  the measurements cover at least two rotations"; found mineral soil C unchanging → supports
  D3; total 6.8 kgC/m² to 1 m → Biosoil/Komeetta look high; Lehtonen acknowledged in it).
  **Lehtonen et al. 2016 (GMD)** = the steady-state convention verbatim + understorey
  underestimation (D2). **Palosuo 2008 thesis** = "no repeated stocks → too expensive" basis.
  **Discussion framing (§3):** position HIKET as the successor Peltoniemi 2004 *called for* —
  repeated stocks + calibrated init finally separate the landscape non-equilibrium trend from
  within-rotation recovery. Collegial (Lehtonen coauthor). Peltoniemi 2004 GChB not paywall-
  blocked anymore.
