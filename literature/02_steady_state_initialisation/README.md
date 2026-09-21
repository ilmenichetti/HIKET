# 02 — Steady-state initialisation across model lineages

**Serves:** `HIKET_storyline_v3.tex` §"The convention: soils initialised at steady state" (`:193`),
flagged there as a **blocking** dependency — *"the claim of near-universality cannot be made without
it"* (Lorenzo, annotation p3, twice).

**The claim this review must support (or weaken):**
> The dominant practice across inventory applications is to start the soil at the equilibrium implied
> by contemporary litter input — and the assumption is rarely tested.

## ⭐⭐ START HERE: Liski et al. 2006 is the reference this section is really about

**Lorenzo, 2026-09-04:** *"Liski did something similar to us, pretty much same approach even if
implemented differently — that's the most important reference."* Confirmed from the PDF
(`../Liski_2006_…pdf`), and it is a **closer antecedent than Wutzler & Reichstein**: same country,
same model family, same inventory-derived litter, same idea — run the national soil forward over the
management history instead of reporting an equilibrium.

Verbatim, §2:

> *"The soil and litter carbon pools at the beginning of the study period were calculated by
> **assuming a steady state** with mean litter input between 1922 and 1936 and mean temperature
> between 1901 and 1930. **Starting from this steady state in 1922, the model was run using annually
> varying values of litter input and temperature.**"*

**That sentence positions the whole paper.** Liski et al. did the transient *run*; the initial state
was **assumed at equilibrium**. HIKET keeps the run and **calibrates the state**. The difference is
one assumption — the one the repeated campaigns can now inform. Write the lineage that way: evolving
the group's own method (Aleksi and Mikko are coauthors on it), not correcting it.

⭐ **They already say the soil is far from equilibrium, in their own discussion:**

> *"soil carbon stock responded slowly to the increased litter production. On the other hand, for
> this same reason, carbon would still accumulate in the litter and soil with no further expansion of
> the forested area if the production of litter is only maintained at the level of 2004 and,
> **centuries later, these carbon stocks would stabilize at a 38% higher level than in 1922**."*

⇒ **This is the lag citation folder 04 §C was missing**, from the most authoritative possible
source: same country, same model, coauthors in common. Centuries to equilibrate, ending 38% higher.

⚠⚠ **Their result also quantifies what an equilibrium start costs.** Per-area soil carbon goes
**6.1 → 6.3 kg C m⁻²** over 1922–2004 — about **+0.024 tC ha⁻¹ yr⁻¹**, an order of magnitude below
the **+0.259** our repeated campaigns show for 1985–2024. A model started full has little room left
to accumulate. ⬜ **Verify before use**: theirs is a national mean over a *changing forest area* on a
"soil and litter" basis, and the periods differ — but if it survives checking it is the most direct
demonstration of this paper's thesis, in the antecedent's own numbers.

⬜ Also extract: **6.1 kg C m⁻² = 61 tC ha⁻¹ in 1922** is a candidate anchor for the currently
**unsourced 1917 stock floor** (`init_state_plausibility.R`'s `STOCK_FLOOR = 40` is a placeholder).
See `../04_finnish_management_history/README.md` §F.

---

**Verdict so far: the first half is supportable and the second half needs softening.** The convention
is real and documentable, but it has been *criticised repeatedly since 2007*, by people inside the
Yasso lineage itself. HIKET is therefore **not the first to say the assumption is wrong** — it is the
first (in this data setting) able to *replace* it with a calibrated transient start and test the
result against repeated national inventories. **Write it that way**; claiming novelty for the
critique would be checkable and wrong.

⚠ Provenance discipline in the tables below: **✅ = verbatim quote or number obtained from the source**
(quotable now); **⬜ = from secondary/summary sources, primary PDF still needed** before it goes in the
manuscript.

---

## A. How each lineage initialises

| Lineage | Initialisation convention | Evidence | ✅/⬜ |
|---|---|---|---|
| **Yasso07 / Yasso15** (Finnish NGHGI) | Run to steady state under mean contemporary litter + climate. Lehtonen et al. 2016: *"Carbon stocks were estimated by running Yasso07 and ROMULv models into a steady state … **If we assume that this average level of inputs and climate has remained steady over centuries, then our soils should approach steady-state conditions.** … For Yasso07, steady state was simulated by running the model 10 000 years, after which relative change of carbon stock was less than 1:10 000."* | local PDF, `../Ortiz_etal_2016_…pdf` | ✅ |
| **Yasso (Liski et al. 2006)** ★★ | **The direct antecedent — see the section above.** Transient national run 1922–2004, but *"assuming a steady state"* for the 1922 initial pools. | local PDF `../Liski_2006_…pdf` | ✅ |
| **Yasso (earlier)** | Same convention in the national/regional applications; Palosuo 2008 documents the lineage and states *why* dynamics were modelled rather than measured (repeated stocks too expensive). | local PDF, `../Palosuo_2008_…pdf` | ✅ |
| **RothC** | Standard practice is an equilibrium ("spin-up") run to stabilise pools, with the annual C input back-solved so equilibrium SOC matches the measured stock. Contreras et al. 2026: *"In the standard RothC, the model is run under long-term equilibrium conditions in order to estimate both the steady-state pool sizes and the annual carbon inputs from plant residues."* | SOIL 12:773–790 (2026) | ✅ |
| **AMG** (French inventory model) | Pools initialised by a fixed default stable fraction (C_S/C_0 = 0.65) rather than by spin-up — an equilibrium-free but equally untested convention. Kanari et al. 2022 measure what it costs. | BG 19:375–387 | ✅ |
| **Century / DayCent** | Long spin-up (typically millennia) to steady state under assumed historic management, then a land-use-change sequence. | ⬜ need Parton et al. 1987 + a DayCent inventory application | ⬜ |
| **ICBM** | Analytical steady state is available in closed form (two-pool linear); the Ultuna parameterisation is anchored on bare-fallow decay, not on an equilibrium fit. Relevant because HIKET's simple models are ICBM-anchored. | local PDF `../Andren_Katterer_1997_…pdf` | ⬜ (quote the initialisation sentence) |
| **CBM-CFS3** (Canadian NGHGI) | **The instructive exception.** Spin-up simulates *repeated rotations of growth and disturbance* until humified pools stabilise, then applies one final stand-replacing "last-pass disturbance" before the reporting period — so the reported initial state is explicitly **not** an equilibrium with current litter. Reported as yielding *"non-equilibrium soil conditions that reflect changes in disturbance regime, management, or species relative to historical conditions."* | Kurz et al. 2009, Ecol. Model. 220:480–504 | ⬜ (secondary; get the primary) |
| **ESMs / global land models** | Spin-up to steady state is standard; Luo et al. 2016 note that *because external forcing is never at steady state, steady-state carbon storage always deviates from realistic storage.* | GBC 30:40–56 | ⬜ |

**The CBM-CFS3 row is the one that decides the tone of the paragraph.** Canada's inventory already
does something transient. So the honest claim is **not** "everyone starts at equilibrium" but:
> *equilibrium initialisation is the dominant convention in the Yasso/RothC/Century inventory
> lineages, including Finland's own; where it has been relaxed (CBM-CFS3's rotation spin-up,
> historical spin-ups) the initial state is imposed from an assumed disturbance history rather than
> **calibrated against observed stock change** — which is what HIKET adds.*

---

## B. The critiques: what has already been shown, and by whom

| Source | What it establishes | Number to quote |
|---|---|---|
| **Wutzler & Reichstein 2007**, BG 4:125–136 ★★ | The clearest *theoretical* statement of the problem, and it uses Yasso — the closest *methodological* antecedent is Liski et al. 2006, above. *"Parameters of these models are often determined in a way that the steady state of the model matches observed carbon stocks. The underlying simplifying assumption is that observed carbon stocks are near equilibrium. This assumption is challenged by observations of very old soils that do still accumulate carbon."* Two consequences named: equilibrium calibration **overestimates the decay rate of the slowest pool**, and spin-up **overestimates stocks of recently disturbed sites**. | small current accumulation ⇒ theoretical equilibrium stocks *"virtually approach infinity"*; transient correction at a beech site = **+5.7 ± 1.5 tC/ha over 100 yr** |
| **Carvalhais et al. 2008**, GBC 22:GB2007 | Relaxing the carbon-cycle steady-state assumption in a model–data fusion improves fit **and removes parameter bias** — i.e. the assumption contaminates the *calibration*, not only the initial state. This is HIKET's mechanism stated for a different model family. | model efficiency **+21%**, normalised average error **−92%** |
| **Luo et al. 2016**, GBC 30:40–56 | Steady-state spin-up is structurally wrong under non-stationary forcing; proposes traceability analysis of C input × residence time — the same decomposition as HIKET's `MRT × σ_input` ridge. | — |
| **Lee & Viscarra Rossel 2020**, NCA 116:245–255 | Initialisation choice alone biases RothC predictions. | **−2.0 to −4.5 Mg C ha⁻¹** over the simulation |
| **Kanari et al. 2022**, BG 19:375–387 | Measures the cost of a default (non-calibrated) pool initialisation and halves it with a data-driven one; explicitly names both spin-up shortcomings: unreliable historical forcing, and *"the unrealistic assumption of SOC equilibrium at simulation start given actual site disturbance histories."* | RMSE **5.95 → 3.60** (best attainable 2.12) Mg C ha⁻¹ |
| **Peltoniemi et al. 2006**, FEM 232:75–85 ★★ **READ 2026-09-11** | **The one citation the Finnish inventory rests its equilibrium initialisation on.** It is a **precision** decomposition, not an accuracy test: the initial state's share of the SD of the annual soil sink is **77% (1989) → 30% (1998) → 21% (2004)**, so "cancel out" in the NID's paraphrase means *no longer dominant*. ⭐ **The initial state they perturbed WAS ITSELF A STEADY STATE** (§2.2: *"At the start of the calculation period (1988), the soil model was in steady state with the input and climate of the first year"*; Appendix A derives it *"based on inputs and model parameters"*) — so they propagated the uncertainty **OF** an equilibrium value, never **OF** the assumption. §2.4: *"our results represent the **precision, not accuracy**, of national level calculation."* ⚠ **Frame as scope, never contradiction — M. Peltoniemi is a coauthor.** Local PDF in this folder; his 2007 thesis (open access, same folder) is Study III's own synthesis and states the caveat in his words. | 77 / 30 / 21% |
| **Le Noë et al. 2023**, Commun. Earth Environ. 4:158 | Review of ~250 SOC models: a *critical lack of independent validation against observed time series*. HIKET's three repeated campaigns are exactly that. Useful for the Introduction's "what the inventory asks a model to do". | ~**60%** of models are not built for prediction at all |
| **Rantakari et al. 2012**, FEM 286:137–147 ★ | Yasso07 tested against a **repeated** Finnish soil inventory (organic layer): stock and accumulation within the error limits of measurement, but accumulation **slightly underestimated**. Nearest precedent to HIKET's target, and its "slight underestimation" is the same direction as our input displacement. | — (get the numbers from the PDF) |
| **Peltoniemi et al. 2004** ★★ | Already in `../README.md`: named the non-equilibrium problem and said *"the long-term trend in carbon accumulation cannot be distinguished from the measured data unless the measurements cover at least two rotations."* | — |
| **Contreras et al. 2026**, SOIL 12:773–790 | Recent, and adopts a historical spin-up *instead of* equilibrium — evidence the convention is actively being abandoned in 2026. | — |

---

## C. Draft text — §"The convention: soils initialised at steady state"

> Inventory-scale soil carbon models are almost always started from an equilibrium. The state at the
> beginning of the simulated period is taken to be the one at which decomposition balances a
> contemporary estimate of litter input, obtained either analytically or by running the model forward
> for thousands of years under recycled climate and constant inputs. Finland's own national inventory
> is explicit about both the procedure and its premise: Yasso07 and ROMULv were *"run into a steady
> state"* over 10 000 years, on the stated assumption that *"this average level of inputs and climate
> has remained steady over centuries"* (Lehtonen et al., 2016). The same convention is standard for
> RothC, where the equilibrium run is additionally used to back-solve the annual carbon input from the
> measured stock (Contreras et al., 2026), and for the spin-up procedures of global land models
> (Luo et al., 2016).
>
> The convention is not unexamined. Wutzler and Reichstein (2007), working with Yasso, showed that if
> a soil is in fact still accumulating, its implied equilibrium stock diverges — small accumulation
> rates push the theoretical equilibrium toward infinity — so that calibrating to current stocks under
> an equilibrium assumption systematically overestimates the decay rate of the slowest pool, while
> spin-up overestimates the stocks of recently disturbed sites. Carvalhais et al. (2008) showed the
> same assumption biases *parameters* retrieved by model–data fusion, not merely the initial state.
> Later work has measured what initialisation choices cost in prediction error (Lee and Viscarra
> Rossel, 2020; Kanari et al., 2022), and one national inventory model, CBM-CFS3, already initialises
> through repeated rotations of growth and disturbance rather than an equilibrium (Kurz et al., 2009).
>
> What has been missing is not the criticism but the test. Where the assumption has been relaxed, the
> replacement initial state has been *imposed* — from an assumed disturbance history, a thermal
> proxy, or a prescribed correction — rather than **inferred from observed change in the stock
> itself**. Repeated national soil inventories are rare, and their absence is the reason the dynamics
> were modelled rather than measured in the first place (Palosuo, 2008); Peltoniemi et al. (2004)
> stated the requirement precisely, that the long-term accumulation trend cannot be separated from the
> data unless the measurements span at least two rotations. Finland now has three campaigns over
> roughly four decades. That is what makes the initial state estimable here, and it is the position
> this paper takes: not that equilibrium initialisation is a novel criticism, but that it can now be
> replaced by a calibrated transient state and the difference measured.

⚠ **Tone note for the writing pass.** The paragraph above deliberately *concedes* the critique's
history. The v3 text currently reads as if the convention were unchallenged; that is not defensible
and a Yasso-lineage reviewer would know Wutzler & Reichstein immediately.

---

## D. What is still missing (blocking → nice-to-have)

1. **Primary sources for the ⬜ rows** — Kurz et al. 2009 (CBM-CFS3), Parton et al. 1987 (Century),
   Coleman & Jenkinson (RothC guide), so each lineage is quoted in its own words rather than through
   a secondary summary. **Blocking** for the near-universality sentence.
2. **Rantakari et al. 2012 PDF** — the repeated-inventory Yasso07 test. Paywalled (Elsevier); Aleksi
   Lehtonen is a coauthor and can supply it. Its "slightly underestimated accumulation" needs its
   actual numbers, because it may be the *same* finding as ours from the other direction.
3. **A count, if the claim stays "dominant".** If the manuscript says *dominant* rather than
   *common*, the defensible form is a tally over a defined set (e.g. every LULUCF Tier-3 soil method
   in the EU NIR submissions) — otherwise soften to *"the convention in the model lineages used for
   forest soil reporting, including Finland's."*
4. Le Noë et al. 2023 supplementary — it may already contain the tally in point 3.

## E. Local PDFs

- `Wutzler_Reichstein_2007_soils_apart_from_equilibrium_BG.pdf` (OA, Copernicus) ★★
- `Kanari_etal_2022_robust_initialization_BG.pdf` (OA, Copernicus)
- Related, already in the parent folder: Lehtonen 2016, Palosuo 2008, Peltoniemi 2004, Andrén & Kätterer 1997.
