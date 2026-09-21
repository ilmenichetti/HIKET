# 01 — The GHG inventory chain: how the forest soil term is produced, reported and used

**Serves:** `HIKET_storyline_v3.tex` §"Why it matters for the national greenhouse gas balance"
(`:183`) — *"the GHG inventory chain: how the soil term is produced, reported and used; where our
models sit in it"*.

**Primary source, and it settles most of this topic on its own:** Finland's National Inventory
Document (NID) submitted to the UNFCCC/EU, *Greenhouse gas emissions in Finland 1990 to 2022*,
Statistics Finland, March 2024. Local copy `FI_NID_2024.pdf` (546 pp; gitignored). Extract text with
`pdftotext -layout FI_NID_2024.pdf FI_NID_2024.txt`. Everything quoted below is verbatim from §6.4
and §6.1, with the text-dump line numbers given so it can be re-found.

⚠ **This is the document our models are ultimately arguing with.** Not a paper — the reporting
instrument itself. Cite the NID, not a secondary summary, for every statement about the method.

---

## A. Where our models sit in the chain (verbatim)

| Link | What the NID says | Line |
|---|---|---|
| **Method class** | *"The methodology for estimating carbon stock changes in soil, litter and dead wood in mineral soils is a **Tier 3** approach and builds on the research by **Liski et al. (2006)**. This method combines forest inventory data, biomass models, litter turnover rates and the dynamic soil carbon model."* | 15974 |
| **The model** | *"For Forest Land Remaining Forest Land, the **Yasso07** model (Tuomi et al. 2011b) was applied."* | 15977 |
| **What is reported** | DOM+SOM **aggregated**, because *"the Yasso07 soil carbon model estimates carbon stock change for the total"*; splitting it *"would be artificial"*. | 15982 |
| **Depth** | *"The Yasso07 model has been defined to estimate carbon stock change to a depth of **one metre**."* | 15984 |
| **Litter input** | Living trees (biomass × turnover coefficients, Table 6.4-2), ground vegetation (Muukkonen et al. 2006: **50.6** S / **66.6** N gC m⁻² a⁻¹ aboveground, Table 6.4-3), harvesting residues and unrecovered natural losses. Fine roots via Helmisaari et al. (2007). | 16005–16050 |
| **Resolution** | Two regions only — *"simulations were made separately for the mineral soils of Southern and Northern Finland"*; weather as **30-year moving averages**. | 15990, 16060 |
| **⚠ INITIALISATION** | *"iii) Estimating the initial values of the model state variables based on NFI6 data (1971–1976) (**so-called spin-up runs to obtain a steady state for the model**)"* … *"The model used the given litter and mean weather data for 1960 to 1990 **as the steady state**. Earlier research has shown that **approximately ten years of simulation since spin-up is enough to cancel out the effect of the spin-up level** (Peltoniemi et al. 2006)."* | 15999, 16069 |
| **Uncertainty** | Mineral soils **31.5%**, *"described in Appendix 6h and by **Lehtonen and Heikkinen (2015)**"*; organic soils 76.0%; **combined Forest Land Remaining Forest Land 83.9%** against a total change of 1.9 Mt C. | 16320, 16336 |

### ⚠⚠ The single most important sentence in this folder

> *"approximately ten years of simulation since spin-up is enough to cancel out the effect of the
> spin-up level (Peltoniemi et al. 2006)"*

This is the inventory's **own written justification** for equilibrium initialisation, and it is
exactly the proposition HIKET tests — over a **35-year** observed window, with the initial state
calibrated rather than assumed. It belongs in the Introduction verbatim. Two cautions before it is
used as a foil:

1. **Read Peltoniemi et al. 2006 first** (*For. Ecol. Manage.* 232:75–85, "Factors affecting the
   uncertainty of sinks and stocks of carbon in Finnish forests soils and vegetation"). The claim is
   about the *spin-up level's effect on the simulated **change*** in a specific setting; quoting it as
   a general claim that initialisation does not matter would be unfair, and Mikko Peltoniemi is a
   close collaborator. **⬜ Blocking: get and read that paper.**
2. Note the NID initialises on **NFI6 (1971–76)** litter, i.e. mid-history, not on today's — a
   detail worth stating accurately rather than caricaturing.

---

## B. Why the soil term matters: the reported numbers (NID 2024, year 2022)

| Quantity | Value | Line |
|---|---|---|
| LULUCF sector | net sink 1990–2017; **net source of 4.4 Mt CO₂ eq in 2022**, 28% more than 2021 | 322, 417 |
| Forest Land net removals | **32.0 Mt CO₂ (1990) → 7.4 (2021) → 7.2 (2022)**; the sink is **down 84%** since 1990 | 15275, 15817 |
| Living biomass | net sink 12.6 Mt CO₂ | 15834 |
| **Mineral soils (DOM+SOM)** | **net sink 4.8 Mt CO₂** — i.e. the soil term is roughly **a third of the size of the biomass sink** | 15835 |
| Organic soils (drained peat) | source 10.1 Mt CO₂ | 15835 |
| Stated causes of the decline | fellings up, NFI increment down, organic-soil emissions up, *"carbon sink of mineral soils has decreased"* | 15260 |

**The argument this supports:** the mineral-soil term is (i) large enough to matter against the
biomass sink, (ii) entirely model-produced — there is no measured soil flux in the chain — and
(iii) carries **31.5%** uncertainty, the second-largest contributor to Forest Land's **83.9%**. When
the national sink is near zero, a term of that size and that uncertainty decides the sign.

⬜ **Policy layer, still unverified.** EU LULUCF Regulation (2018/841 as amended by 2023/839), the
"no-debit" rule, Finland's 2021–2025 reference level (~21 Mt CO₂ eq per a secondary NGO source) and
the 2026–2030 national removal target. **Do not quote numbers from news or NGO pages** — take them
from the Regulation's Annex and the Commission's implementing decision.

---

## C. Draft text — §"Why it matters for the national greenhouse gas balance"

> Finland reports the carbon balance of its forest mineral soils to the UNFCCC with a Tier 3 method:
> forest inventory data and biomass models supply a litter input, and the Yasso07 decomposition model
> converts it into an annual stock change for dead organic matter and soil organic matter combined,
> to one metre depth, resolved into a southern and a northern region (Statistics Finland, 2024,
> §6.4; the approach follows Liski et al., 2006). No soil flux is measured anywhere in this chain.
> The reported number is a model output, and its uncertainty is stated as 31.5%, against 16.3% for
> the tree biomass change it accompanies.
>
> The term is not small. In 2022 mineral forest soils were reported as a sink of 4.8 Mt CO₂, beside a
> living-biomass sink of 12.6 Mt CO₂, while the forest land sink as a whole had fallen 84% since 1990
> and the LULUCF sector had become a net source of 4.4 Mt CO₂ eq. When the national balance sits this
> close to zero, a model-derived term of that magnitude, carrying that uncertainty, is capable of
> deciding its sign.
>
> The same reporting chain also fixes the initial state of the soil by a spin-up to steady state — in
> the Finnish inventory, on the average litter input of the sixth National Forest Inventory
> (1971–1976) and 1960–1990 mean weather — on the argument that "approximately ten years of
> simulation since spin-up is enough to cancel out the effect of the spin-up level" (Statistics
> Finland, 2024, §6.4, citing Peltoniemi et al., 2006). That proposition is testable, and testing it
> over four decades of repeated soil measurements is what this paper does.
>
> The quantity the inventory needs from a soil model follows from this position in the chain. It is
> not the ability to rank plots against one another, but the ability to reproduce a national
> trajectory over decades and project it forward. That is the quantity evaluated here, and stating it
> now prevents a later result — near-zero plot-level skill in every model — from being read as
> failure at a task the inventory never asks of these models.

---

## D. What is still missing

1. **⬜ Peltoniemi et al. 2006** (FEM 232:75–85) — blocking, see §A.
2. **⬜ Liski et al. 2006** — already in the parent folder (`../Liski_2006_…pdf`) and already reviewed
   in memory `liski-2006-input-benchmark`; cross-link rather than re-read. It is the *method's own
   founding paper* and our input benchmark at the same time, which is worth saying once in the text.
3. **⬜ The EU policy layer** from primary legal sources (§B).
4. **⬜ Lehtonen & Heikkinen 2015** — already the σ_input anchor (parent folder); note it is *also*
   the source of the inventory's 31.5% soil uncertainty. Same paper doing double duty; say so.
5. ⬜ The 2026 NID (`FI_NID_UNFCCC_BTR2_2024_2026-04-15.pdf`, data to 2024) — too large for automated
   fetch here; download manually if the paper needs the newest figures. The 2024 edition is enough
   for the method, which has not changed.
