# HIKET — literature (crucial references)

Local library of the load-bearing references for the HIKET manuscript. **PDFs are
git-ignored** (copyright + size); this catalogue is tracked. Drop new PDFs here and add a
row + anchors below. Started 2026-07-15.

---

## 1. Andrén & Kätterer (1997) — ICBM
`Andren_Katterer_1997_ICBM_EcolAppl.pdf`
Andrén, O. & Kätterer, T. (1997). *ICBM: the introductory carbon balance model for
exploration of soil carbon balances.* Ecological Applications 7(4):1226–1236.

**Role:** the external anchor for the simple-model (SP1/TP2/TP3) kinetics (revision C1).
**Anchors:** two-pool Y/O; `k1=0.8`, `k2=0.00605 yr⁻¹` (Ultuna bare fallow), `h=0.13`,
`r=1` = central Sweden. Full parameter memory: `icbm-parameters`.

---

## 2. Lehtonen et al. (2016) — Finnish national soil-C inventory  ★ KEYSTONE (the foil)
`Ortiz_etal_2016_Finnish_forest_soilC_inventory_GMD.pdf`
Lehtonen, A., Linkosalo, T., Peltoniemi, M., Sievänen, R., Mäkipää, R., Tamminen, P.,
Salemaa, M., Nieminen, T., Ťupek, B., Heikkinen, J., Komarov, A. (2016). *Forest soil
carbon stock estimates in a nationwide inventory: evaluating performance of the ROMULv and
Yasso07 models in Finland.* Geosci. Model Dev. 9:4169–4183.

**Role:** documents the **steady-state initialization convention that HIKET challenges**,
*and* the **understorey-litter underestimation** that grounds our σ_input centre (D2).
**Anchors:**
- **Steady-state convention (verbatim):** "Carbon stocks were estimated by running Yasso07
  and ROMULv models into a steady state (i.e. a state where carbon input for the model
  equals carbon flux due to decomposition). … **If we assume that this average level of
  inputs and climate has remained steady over centuries, then our soils should approach
  steady-state conditions.** … For Yasso07, steady state was simulated by running the model
  10 000 years, after which relative change of carbon stock was less than 1:10 000." → this
  is exactly the assumption HIKET argues is false (managed forests out of equilibrium; a
  ~70% growing-stock rise → rising litter).
- **Understorey underestimated (D2):** "the role of understorey litter input was
  underestimated when the Yasso07 model was parameterised, especially in northern Finland."
  Bryophytes = 20–60% of NDVI at high northern latitudes (Yuan et al. 2014). Table 2 =
  understorey turnover rates (dwarf shrubs, grasses, herbs, bryophytes, lichen). Supports
  centring σ_input > 1 and the S→N gradient argument.
- **Depth:** Yasso07 "calibrated with litter input estimates and soil carbon measurements
  down to 1 m soil depth" — relevant to the D3 depth question.
- Inventory favours **simple** models ("transparent and verifiable" for GHG reporting) —
  supports HIKET's "Yasso is the right inventory tool" framing.

---

## 3. Palosuo (2008) — YASSO thesis  ★ KEYSTONE (lineage + the "no repeated stocks" point)
`Palosuo_2008_thesis_YASSO_DissForestales61.pdf`
Palosuo, T. (2008). *Soil carbon modelling as a tool for carbon balance studies in
forestry.* Dissertationes Forestales 61, University of Helsinki. (Supervisors Liski,
Sievänen.)

**Role:** develops/evaluates YASSO; states *why* dynamics were modelled rather than
measured — the exact basis for the coauthor's "no repeated stocks" hypothesis.
**Anchors:**
- **Confirms the coauthor (verbatim, §1.3):** "Repeated sampling would be the most
  straightforward way to assess such changes, but the method is often considered to be too
  expensive and to require excessive effort. **Such is the case with forest soils.**" → soil
  C dynamics were inferred indirectly (models, chronosequences) *because* repeated stocks
  were unavailable. HIKET's VMI8/Biosoil/Komeetta repeats are what finally allow a *direct*
  test — subject to the D3 between-campaign QC caveat.
- "A typical feature of the dynamic model is the **memory**: the state … affects its state
  in the following moments." — the humus-memory idea behind transient init.
- **Study V** applied YASSO to Finland **1922–2004** (a historical/transient national run) —
  precedent for transient *forward* runs; HIKET's novelty is calibrating the *initial state*.
- Older YASSO structure (extractives/celluloses/lignin/Humus1/Humus2), pre-AWEN.

---

## 4. Peltoniemi et al. (2004) — stand-age chronosequence  ★★ KEYSTONE (the direct antecedent)
`Peltoniemi_etal_2004_soilC_standage_GCB.pdf` (read in full 2026-07-16; DOI 10.1111/j.1365-2486.2004.00881.x)
Peltoniemi, M., Mäkipää, R., Liski, J., Tamminen, P. (2004). *Changes in soil carbon with
stand age — an evaluation of a modelling method with empirical data.* Global Change Biology
10(12):2078–2091. **Aleksi Lehtonen is acknowledged in it** — direct lineage to the HIKET coauthor.

**Role:** far more than "the steady-state basis" — it is the paper that **identified the
non-equilibrium initialization problem HIKET solves**, and it independently corroborates
our D2 and D3 findings. Method: Motti stand simulator → biomass/allometry + turnover →
litter → Yasso; evaluated against a **stand-age chronosequence** of 64 southern-Finland NFI
plots (space-for-time, NOT repeated same-plot measurements).
**Anchors:**
- **Non-equilibrium, called out in 2004 (p2087, "Effect of site history"):** sites may sit
  *below* equilibrium from historical **slash-and-burn** ending ~early 20th c.; and — key
  quote — "**the long-term trend in carbon accumulation cannot be distinguished from the
  measured data unless the measurements cover at least two rotations.**" → HIKET is the
  direct successor this called for: repeated campaigns (VMI8/Biosoil/Komeetta, ~40 yr) + a
  *calibrated* initial state finally separate the long-term (landscape non-equilibrium) trend
  from within-rotation recovery.
- **Mineral soil C does NOT change with stand age** (p2085–86: "The amount of carbon in
  mineral soil did not change"; "Changes in the mineral soil carbon were not found"). →
  **strong external support for D3**: our VMI8→Biosoil near-doubling of *mineral* C (incl.
  deep 20–40 cm) is not physically plausible = a between-campaign method artifact.
- **Level (D3/D2 "stocks seem high"):** careful measured total (F/H + mineral to **1 m**) =
  **6.8 ± 2.5 kgC/m²** (sim 7.0). Matches our VMI8 (~63 tC/ha); makes Biosoil/Komeetta
  (~102–105 tC/ha = 10+ kgC/m²) look high — and ours are only to 40 cm, so even more so.
- **Initialization sensitivity (prefigures HIKET):** "Estimation of the accumulation rate …
  was highly dependent on the initial stock value" (tested ±20%). HIKET makes the initial
  state a *calibrated* parameter with propagated uncertainty — the next step.
- **Understorey litter (D2):** includes ground-vegetation litter (Tables 2–4); understorey
  turnover rates (Table 3: dwarf shrubs 0.25/0.33, bryophytes 0.33, herbs/grasses 1.0/0.33,
  lichen 0.1); southern-Finland total litter ~2.6–3.35 tC/ha/yr incl. understorey.
- **Layer lumping:** notes Yasso does not separate organic vs mineral, and that the changes
  are in the organic layer not mineral — the same limitation behind our layer-resolved-
  likelihood idea for the 1985 down-weighting (C5).

---

## 5. Korhonen, Räty et al. (2024) — Forests of Finland / NFI growing stock
`Korhonen_etal_2024_ForestsOfFinland_SilvaFennica_24045.pdf`
Korhonen K.T., Räty M., Haakana H., Heikkinen J., Hotanen J.-P., Kuronen M., Pitkänen J.
(2024). *Forests of Finland 2019–2023 and their development 1921–2023.* Silva Fennica
58(5) art. 24045. https://doi.org/10.14214/sf.24045 (open access).

**Role:** the [History] thread's quantitative backbone — the **rising, non-stationary
growing-stock base** that physically justifies a below-equilibrium 1917 start
(`sigma_init < 1`) and rising litter inputs; the source for the **F10b** figure and the
**C3** growing-stock-shaped pre-run interpolation.
**Anchors:** total growing stock 1.4 G m³ (NFI1, recalculated) → 2.6 G m³ (NFI13),
**+84% over 100 yr**; "most of the increase … after the end of the 1960s" (§3.5.2, Fig 10a).
Series lives *only* as Fig 10a + endpoints (no per-NFI table) → digitized into
`Data/forest_history/nfi_growing_stock.csv` (endpoints exact, NFI2–NFI10 read ±~50; see
that folder's README). Supplement S2 Table 55 = NFI13 detailed breakdown (not the series).

---

## The one-line answer to the coauthor (comments #2/#3)
Yes: the steady-state convention rests on a lineage (Peltoniemi 2004 chronosequence → YASSO
development, Palosuo 2008 → operational inventory, Lehtonen 2016) that **lacked repeated soil
C stocks** and so inferred the dynamics indirectly (space-for-time, or an assumed
centuries-steady input/climate). HIKET's repeated campaigns let us test the dynamics
*directly* and drop the equilibrium-start assumption — the contribution. (Caveat: the
between-campaign mineral-soil QC issue, D3, partly compromises that repeated-stock advantage
and must be stated.)

---

## `input_estimates/` — independent NPP, for bounding the litter-input flux (added 2026-08-14)

PDFs are untracked (as elsewhere here); this records what they are and what they settle.

- **Gower et al. 2001**, *Ecol. Appl.* — NPP and carbon allocation of boreal forest ecosystems.
  **Reports in g C m⁻², explicitly (Table 4 caption), and its NPP includes an "Understory NPP" row.**
  Nordic Class I: total NPP mean **321**, range **215–462** g C m⁻² yr⁻¹ (3.21; 2.15–4.62 tC/ha/yr).
  Class I evergreen (global): mean 387, range 214–912 — **the 912 is where the current 8.7 flux
  window ceiling comes from**, i.e. a global single-stand maximum.
  ⚠ Gower's own caveat: mycorrhizal NPP is excluded, so "all the total NPP estimates are likely to
  be underestimates" (Vogt et al. 1982 put mycorrhizae at ~15% of NPP).

- **Zheng et al. 2004**, *J. Veg. Sci.* 15:161–170 — gridded NPP for Finland and Sweden from field
  inventory + AVHRR. Above-ground NPP mean **408** (172–1091); **total NPP mean 563** (252–1426).
  ⚠⚠ **THESE ARE DRY MATTER, NOT CARBON.** Proof is internal: Zheng cites Gower's world-boreal TNPP
  as 109–1827 (mean 892) where Gower's own Appendix A gives 218–912 (mean 424) g C — a ratio of
  exactly 2.00–2.10. So Zheng's 563 is **≈2.81**, not 5.63, tC/ha/yr. Mixing the two inflates any
  NPP-based ceiling twofold.
  ⚠ Not fully independent of the inventory lineage: it runs on NFI field data, so it is a different
  *estimator*, not a different *data source*.

**What they were used for:** M&M §"Prior specification: the litter-input flux window" — the decision
to narrow the `flux_pair` window to `[0.05, 4.62]` (Gower's *Nordic* maximum) because the present
ceiling bounds a **national mean** with a **global single-stand maximum**.
