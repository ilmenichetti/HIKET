# 10 — How national inventories obtain the forest mineral-soil term

**Question (Lorenzo, 2026-09-15):** the Introduction says *"few countries have a survey of that kind,
so the soil term is more often estimated with a model"* — uncited. Source it from the UNFCCC
submissions, go broader where an aggregated source exists, and aggregate what is online about the
Finnish inventory's own soil term over time. Built 2026-09-15/16 from the sources listed below;
✅ = read at the primary source, ⬜ = secondary only, needs the primary.

## 1. Aggregated sources

| source | what it gives | status |
|---|---|---|
| **Blujdea, Abad Viñas, Federici & Grassi 2015**, *Carbon Management* 6:247–259, doi 10.1080/17583004.2016.1151504 (CC-BY) | Comparative analysis of LULUCF methods across the 28 EU member states, from the NIRs. Search-engine summary of the forest-soil part: *five member states use models, covering ~18% of EU forest area, driven by NFI biomass at plot level (Austria, Finland, Sweden among them)*. Abstract: "the obvious trend is the move toward statistical sampling covering all land categories … and slight increasing use of models." | ⬜ **PDF blocked to scripts (Cloudflare) — download by browser**; the forest-soil table must be read before the "18%" is quoted. NB the Sweden row is ~2013 vintage: Sweden reported with the Q model then and measures now. |
| **Didion, Blujdea, Grassi, Hernández, Jandl, Kriiska, Lehtonen & Saint-André 2016**, *Models for reporting forest litter and soil C pools in national greenhouse gas inventories: methodological considerations and requirements*, *Carbon Management* 7, doi 10.1080/17583004.2016.1166457 | The methods review for exactly our question — models for the DOM/soil pools in NGHGIs. **Aleksi is a coauthor**; ask him for the PDF and, more usefully, for the current picture. | ⬜ abstract elided by the publisher; not read |
| **Hernández et al. 2017** (Jandl, Blujdea, **Lehtonen**, Kriiska, … Didion), *Sci. Total Environ.* 599–600:1171–1180, doi 10.1016/j.scitotenv.2017.03.298 | Yasso07 across European forest conditions, country cases with differing data. Abstract (verbatim): *"changes in soil organic carbon (SOC) occur slowly and these changes may not be captured through repeated soil inventories. Simulation models may be used as alternatives to SOC measurement"*; *"The obstacles encountered when applying the Yasso07 model reflect a lack of available input data. Future research should focus on improving our knowledge of C inputs from compartments such as shrubs, herbs, coarse woody debris and fine roots."* | ✅ abstract (PubMed E-utilities). ⚠ that last sentence is our missing-flux argument, stated by the model's own user community in 2017. |
| **Jandl et al. 2014**, *Sci. Total Environ.* 468–469:376–383, doi 10.1016/j.scitotenv.2013.08.026 | Status of SOC monitoring; BioSoil as the harmonised European attempt and its difficulties. Not about inventory methods per se. | ⬜ abstract |
| **Ortiz et al. 2013**, *Ecol. Model.* 251:221–231 (already in library.bib) | Swedish measured estimate vs Yasso07 and Q for the same soils: 1994–2000 change 6.6 (±7) Tg C/yr measured, 1.7 (±8.8) Yasso07, −3.2 (+10.5/−16.9) Q. | ✅ in library |

## 2. Country by country — forest land remaining forest land, mineral soils

| country | approach | model / data | initialisation | soil measurements in the chain? | source |
|---|---|---|---|---|---|
| **Finland** | model, Tier 3 | Yasso07; litter from NFI biomass × turnover; DOM+SOM combined to 1 m | steady state at NFI6 litter (1971–76) and 1960–90 climate; defended by Peltoniemi et al. 2006 | none in the chain; evaluation only (Lehtonen et al. 2016) | ✅ `FI_NID_2025.txt` 16413–16480; `../01_ghg_inventory_chain/` |
| **Sweden** | **measured**, Tier 3 | Swedish Forest Soil Inventory on the NFI permanent plots, since 1993, ten-year cycle; stock-change with pedotransfer functions for bulk density and stoniness | n/a | yes — it *is* the chain | ✅ `SE_NID_2026.txt` 21425–21560 (§6.4.2.1, §6.4.2.4.1) |
| **Norway** | model, Tier 3 | Yasso07 per NFI plot, "total SOC" = dead wood + litter + mineral soil, then split | **two-step**: 5000-yr spin-up to equilibrium at the **mean 1986–2016 input** by species × site index, then a constructed 1951–1990 backcast series "to reduce the effect of the equilibrium assumption". Verbatim: *"In a managed forest landscape neither forest biomass, soil nor DOM pools can be assumed to be in an equilibrium state"* | none yet: soil sampling on NFI plots began 2022 (388 plots); the 1988–92 survey "does not qualify as a baseline"; first national change estimate "earliest in 2035"; six-site comparison project GJENFERD instead | ✅ `NO_NID_2025.txt` 21929–22660 |
| **Germany** | **measured**, national forest soil inventory | NFSI/BZE: ~1 800 plots on an 8×8 km grid, 1987–92 and 2006–08; organic layer + mineral soil to 30 cm | n/a | yes | ✅ abstract of Grüneberg, Ziche & Wellbrock 2014, *GCB* 20:2644–2662, doi 10.1111/gcb.12558: *"provides the Greenhouse Gas Reporting in Germany with a quantitative assessment"*; **mineral soil sequestered 0.41 Mg C/ha/yr**, organic layer stable |
| **Switzerland** | model | Yasso07 for dead wood, litter and mineral-soil change; NFI-driven | (not read) | Swiss forest soil inventory used for evaluation (Yasso20, Biogeosciences 2025) | ⬜ secondary (Didion et al.; Norway NID fn 50) |
| **Austria** | model | Yasso07 | (not read) | Austrian Forest Soil Survey for evaluation (Hernández 2017) | ⬜ Norway NID fn 50: "Yasso07 is used for UNFCCC reporting for forests in Finland, Switzerland, and Austria" |
| **Canada** | model, Tier 3 | CBM-CFS3 (Kurz et al. 2009); ~2.7 million inventory records; DOM and mineral-soil pools | spin-up of repeated growth–disturbance cycles to a steady state, then a dated "last pass disturbance" before 1990 | none | ⬜ secondary; primary = Canada NIR Part 2 (publications.gc.ca) and Kurz 2009 (in library) |
| **USA** | empirical stock model | FIA soil measurements (since 2001, 20 cm) scaled to 100 cm by an empirical model (Domke et al. 2017); stock-change from the FIADB; **no decomposition model, no spin-up**. Alaska / territories: Tier 1, *"no net carbon stock change reported"* | n/a | yes (FIA soil subplots) | ✅ `US_GHGI_2024_ch6_LULUCF.pdf` pp. 6-30, "Carbon in Forest Soil" |
| **Japan** | model | CENTURY-jfos, with a national forest soil carbon inventory (NFSCI) for stocks | (not read) | NFSCI for evaluation | ⬜ NIES presentation (Shirato) |
| **Denmark, Estonia, Latvia, Lithuania, France, UK** | — | not established | | | ⬜ not read |

**What this supports.** The sentence in the draft is defensible as *"most countries that report the term
use a model; a few measure it"*, with Sweden and Germany as the measuring cases and Finland, Norway,
Switzerland, Austria, Canada, Japan as modelling cases. It is **not** yet defensible as a count or a share
until Blujdea 2015 / Didion 2016 are read.

**Two things found on the way that matter more than the count:**
1. **Norway states our premise in its NID** and acts on it with a backcast pre-run — the same structure
   as HIKET's transient initialisation (equilibrium + prescribed ramp), except that Norway equilibrates at
   the *modern* mean input and has **no soil measurement** to estimate the start from (first change
   estimate 2035). HIKET is the case where the measurement exists.
2. **Germany's measured rate, +0.41 Mg C/ha/yr (mineral soil to 30 cm, 1987–92 → 2006–08)**, is the
   same order as our +0.399 for 1989 → 2006 on the whole profile. Different depth, different basis — but
   worth one clause where the Finnish rate is first given.

## 3. The Finnish inventory's own soil term over time

Reported mineral-soil (DOM + SOM) term on forest land remaining forest land, as the submissions state it.
Mineral forest-land area for the per-hectare conversion: **15.9–16.0 Mha** (NID 2025 Table 6.4-1:
15 998 kha in 1990, 15 846 in 2014). 1 Mt CO₂ = 0.2727 Mt C.

| period | reported term | ≈ tC ha⁻¹ yr⁻¹ | source |
|---|---|---|---|
| early 2000s | sink of ~**10 Mt CO₂ eq** | **+0.17** | Luke news 15.1.2025 (2025 submission, recalculated) |
| 2022 | sink of **4.8 Mt CO₂**; "net sink during the whole time series" | +0.08 | NID 2024 (`../01_ghg_inventory_chain/FI_NID_2024.txt` 331, 15835) |
| **2021 onwards** | **source**; 0.4 Mt CO₂ emitted in 2023 | −0.007 | NID 2025 (`FI_NID_2025.txt` 16206–16208, 17037) |
| 2024 (preliminary) | forest land net sink 0.1 Mt CO₂ eq, "decreased accumulation of carbon into the growing stock, mineral soil and forest litter" | | Statistics Finland, 2024 preliminary |

⚠ **The 2024 and 2025 submissions disagree about the sign in 2021–2022** because the 2025 submission
recalculated the litter input (new BCEFs, NFI13, drain data → *"decreased total litter input to the mineral
soil"*, Fig. 6.4-2). The draft's §"What the Finnish chain has reported" was written from the 2024 text and
**must be restated on the 2025 submission** (done 2026-09-16, see the draft).

**The inventory's own explanation of the decline** (verbatim, NID 2025 and Luke 15.1.2025):
- *"The decrease in the litter input is one of the reasons which turned mineral soils from a net sink to a
  net source from 2021 onwards."* Growing stock still rose NFI12 → NFI13, so the cause is **foliage biomass**
  falling between the two inventories (smaller BCEFs) and slash recovered for energy.
- *"The long-term change in the carbon sink in mineral soils from −10 Mt CO₂ eq at the beginning of the
  2000s into an emission source can be explained by the decrease in the litter input generated by the growing
  stock … as well as global warming."*
- Mechanism of the sink itself, NID 2024: *"Increase in the volume of roundwood removals increases the
  litter input to the soil, thus increasing the carbon storage of combined pools of DOM and SOM."*

**Against the measurements** (balanced set, n = 310, whole profile, true years — `obs_basis.R`):
1989 → 2006 **+0.399**, 2006 → 2024 **+0.117 ± 0.066**, 1989 → 2024 **+0.259 ± 0.040** tC ha⁻¹ yr⁻¹.
So the modelled term went from about half the measured early rate to a source, while the measured
soil went from +0.40 to +0.12 (interval including zero). Same direction, larger swing in the model.
⚠ Bases differ: the inventory is area-weighted over all mineral forest land and includes dead wood; ours is
an unweighted balanced plot set of the soil. Quote as an order-of-magnitude comparison, not a test.

**The Finnish NFI and its soil campaigns** (for the Introduction's "three campaigns" paragraph):
NFI since 1921 (NFI1 1921–24; NFI6 1971–76 is the spin-up basis; NFI8 1986–94 carried the first soil
campaign on its permanent plots, fieldwork 1986–1995; BioSoil 2006 on a subset (EU/ICP Forests Level I);
third campaign 2024 (Komeetta); NFI13 2019–2023 is the current inventory used for reporting).
Growing stock 1.4 → 2.6 billion m³ over the century (Korhonen et al. 2024, in library).

## 4. Draft text

**§ The soil term is often modelled rather than measured** — replace the uncited clause with:

> Countries obtain the soil term of their inventories in more than one way. Sweden and Germany measure it:
> Sweden from a soil inventory carried on the permanent plots of its national forest inventory since 1993 on a
> ten-year cycle \citep{SwedenNID2026}, Germany from two national forest soil inventories two decades apart
> \citep{Gruneberg2014}. Repeating a national soil survey often enough to yield a stock change is expensive,
> and the more common route is a model driven by the forest inventory: Finland, Norway, Switzerland and
> Austria run Yasso07 on their inventory plots \citep{StatFin2025NID,NorwayNID2025,Didion2016}, Canada
> runs CBM-CFS3 \citep{Kurz2009}, and Japan a version of CENTURY. Norway's inventory states the
> difficulty the modelling route inherits: ``in a managed forest landscape neither forest biomass, soil nor DOM
> pools can be assumed to be in an equilibrium state'' \citep{NorwayNID2025} --- and, having no repeated
> soil measurement to start from, expects its first measured national change estimate in 2035.

**§ What the Finnish chain has reported, and why** — restate on the 2025 submission (applied in the draft).

## 5. Gaps
- ⬜ Blujdea 2015 PDF (browser download; CC-BY) — read the forest-soil table, replace "five member states … 18%".
- ⬜ Didion 2016 — ask Aleksi; also the natural person to confirm the country table above.
- ⬜ Switzerland, Austria, Canada, Japan: read the primary NIDs before citing them as such (secondary now).
- ⬜ The CRT series behind NID 2025 Fig. 6.4-x for mineral soils, 1990–2023, to draw the reported term
  against the measured rate (one panel, Introduction or Discussion).
- ⬜ Confirm with Luke colleagues the −10 Mt CO₂ eq "beginning of the 2000s" figure (news release wording).

## Files
- `SE_NID_2026.pdf/.txt` — Sweden, NID 2026 (497 pp)
- `NO_NID_2025.pdf/.txt` — Norway, NID 2025 (626 pp)
- `FI_NID_2025.pdf/.txt` — Finland, NID 2025 (1990–2023; 558 pp). The 2024 submission is in `../01_ghg_inventory_chain/`.
- `US_GHGI_2024_ch6_LULUCF.pdf` — USA, inventory 2024, chapter 6
- `refs.bib` — entries for the sources above, ready to merge into `manuscript/library.bib`
