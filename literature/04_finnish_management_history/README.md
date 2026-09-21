# 04 — Finnish forest management history and the soil-carbon lag

**Serves:** `HIKET_storyline_v3.tex` §"Finnish forests are a managed landscape, not an equilibrium
one" (`:217`). Lorenzo's annotation (p5): *"This is load bearing and will also need proper literature
revision."* — and the storyline itself calls it *"the load-bearing claim of the whole paper."*

**The claim to support, in three parts.** (i) Finnish forests were heavily depleted by the early
20th century and have been recovering under management since; (ii) that recovery raised litter input
over the period our models integrate; (iii) soil carbon lags a change in input by long enough that a
soil integrating this history is **not** at equilibrium with today's litter.

**Status: (i) and (ii) are well sourced. (iii) is the weak leg and needs the most work.**

---

## A. The depleted starting point (~1920)

| Source | What it gives | ✅/⬜ |
|---|---|---|
| **Aakala, Kulha & Kuuluvainen 2023**, *Landscape Ecology* 38:2417–2431 ★ | The best single citation for the state of the forest at NFI1. *"In eastern Finland almost all productive forest land was used in slash-and-burn cultivation."* Woodland grazing *"efficiently prevented forest regeneration"*; high-grading of sawn timber from the mid-1800s; shipbuilding, tar and potash. Southern Finland *"virtually absent of old trees"* outside inaccessible areas, against a north that retained them — a **striking north–south contrast**, with **human population density the strongest predictor** of depletion. *"The wood volume in Finnish forests reached a low-point in the beginning of the twentieth century."* | ✅ |
| **Korhonen et al. 2024**, *Silva Fennica* 58:24045 | The quantitative backbone, already local (`../Korhonen_etal_2024_…pdf`): growing stock **1.4 → 2.6 G m³** (NFI1 → NFI13), **+84% in 100 yr**, most of it after the late 1960s. Source of F10b and of the pre-run ramp. | ✅ |
| Slash-and-burn extent | ~4 Mha of forest land affected *cumulatively*; widest at the turn of the 18th/19th centuries; abandoned through the 19th century, ended in the 1940s. Tar: Finland the leading European exporter, 17th to mid-19th century. | ⬜ **secondary web sources only — do not cite as is.** Kaila (1932, Silva Fennica) on tar burning is a primary lead; a forest-history monograph is needed for the area figures. |

⚠ **The north–south contrast in Aakala et al. is not decoration — it is a testable prediction for
us.** If depletion was population-driven and southern, then the below-equilibrium initial state
should be *deeper in the south*. Our σ_init is a **single global scalar**, so we cannot express that;
worth stating as a limitation, and it is a natural next-study hook alongside the residual strand.

---

## B. The managed landscape is below its own potential (contemporary evidence)

| Source | What it gives | ✅/⬜ |
|---|---|---|
| **Kumpu, Peltoniemi, Forsius & Mäkelä 2025**, *Silva Fennica* 59, art. 24072 ★★ | 27 old **unmanaged** stands vs 213 managed Biosoil sites: mineral soil **0–10 cm 2.44 vs 1.78 kg C m⁻²** and **10–20 cm 1.62 vs 1.14** (both p<0.001) — i.e. **+37% and +42%** in the mineral layers; organic layer 2.58 vs 2.06, n.s. And at **resampled** sites (2001–03 → 2021) the unmanaged stocks are **static**: mean change **−692 kg C ha⁻¹ (p = 0.75)**. | ✅ |
| **Peltoniemi et al. 2004** (parent folder) | Names the mechanism for our period: sites may sit *below* equilibrium from slash-and-burn ending in the early 20th century; and *"the long-term trend in carbon accumulation cannot be distinguished from the measured data unless the measurements cover at least two rotations."* | ✅ |
| **Liski et al. 2006** (parent folder) | Finland's forests 1922–2004 as an accumulating system — the inventory method's founding paper *and* our input benchmark. | ✅ |

**This pair is the strongest single argument in the folder.** Managed soils hold ~40% less mineral
soil carbon than comparable unmanaged ones, and the unmanaged ones are not changing. A depleted,
still-recovering managed landscape is exactly the premise a transient initialisation encodes and an
equilibrium initialisation denies.

⚠ Two cautions before it is used. Unmanaged sites are **not a chronosequence of our plots** — they
differ in site type selection and possibly in soil sampling depth basis; check comparability against
our whole-profile target before quoting the percentages side by side with our stocks. And "unmanaged
is stable" is evidence about *old* stands, not proof that managed stands are climbing toward that
level.

---

## C. The lag: how long soil carbon takes to catch up  ⬜ **the weak leg**

What we can say now, and what it rests on:

- Post-fire boreal chronosequences: accumulation fastest in the first ~50 yr, ecosystem carbon at
  ~90% of maximum by ~60 yr, an asymptotic form over the first ~80 yr; rates **−0.6 to +1.7**, median
  **0.15 Mg C ha⁻¹ yr⁻¹**; but a 5000-yr chronosequence gives only **~0.05 Mg C ha⁻¹ yr⁻¹**, and a
  314-yr sequence shows all pools still rising linearly. ⬜ **All of this is from search summaries —
  every number needs its primary source before use.** Leads: Wiley GCB 2024 *"The biological controls
  of soil carbon accumulation following wildfire"*; *Global Biogeochemical Cycles* 2020GB006612
  (paywalled, 403 here); Wardle's Swedish island chronosequence.
- **Wutzler & Reichstein 2007** (folder 02) is the theoretical statement of the same thing and is
  already ✅: soils disturbed *several centuries* ago remain in a transient state because the slowest
  pool is still accumulating; small current accumulation implies an enormous implied equilibrium.
- ⚠ **Note the tension to handle honestly.** A ~60-yr, 90%-of-maximum recovery after fire would argue
  that a century of regrowth is nearly *enough* — that soils should be close to equilibrium by now.
  The reconciliation is pool-specific: the fast pools equilibrate in decades, the slow humus pool
  does not, and it is the slow pool that carries the stock. **Our own result speaks to this** (the
  models' intrinsic MRT is 22–27 yr, yet the fitted transient still matters over 35 yr), but the
  literature paragraph must not paper over it.

---

## D. Draft text — §"Finnish forests are a managed landscape, not an equilibrium one"

> The Finnish forest of the early twentieth century was a heavily used one. Slash-and-burn
> cultivation had passed over most of the productive forest land of the east, woodland grazing
> suppressed regeneration, and tar burning and the high-grading of sawn timber had drawn down the
> accessible stands; by the first National Forest Inventory, old trees were effectively absent from
> southern Finland outside the least accessible areas, and forest volume was at its historical low
> (Aakala et al., 2023). What followed was a century of managed regrowth: growing stock rose from
> roughly 1.4 to 2.6 billion m³ between the first and the thirteenth inventories, an increase of 84%,
> concentrated after the late 1960s (Korhonen et al., 2024).
>
> Litter input to the soil is tied to that standing biomass, so the input our models integrate has
> been rising for the whole period they are asked to reproduce. Whether the soil has kept pace is a
> separate question, and the contemporary evidence says it has not: old unmanaged stands in southern
> and eastern Finland hold about 40% more carbon in the mineral soil than managed stands of
> comparable site types, while their own stocks are statistically unchanged over a twenty-year
> resampling (Kumpu et al., 2025). A managed landscape below the level of its unmanaged counterparts,
> under a still-rising input, is the definition of a soil away from equilibrium.
>
> The consequence for initialisation was stated twenty years ago in this literature. Peltoniemi et
> al. (2004) attributed part of the soil carbon of Finnish forest plots to a history of slash-and-burn
> ending in the early twentieth century, and concluded that the long-term accumulation trend cannot be
> separated from measured data unless those data span at least two rotations; Wutzler and Reichstein
> (2007) showed that soils disturbed centuries earlier remain in a transient state, because it is the
> slowest pool that is still filling. A soil integrating a century of managed recovery is therefore
> not the soil an equilibrium spin-up produces, and the difference is precisely the quantity the
> national inventory needs.

---

## E. What is still missing

1. **⬜ The lag evidence (§C), with primary sources.** Blocking — it is the leg that turns "inputs
   rose" into "the soil is still catching up", and it currently rests on search summaries.
2. **⬜ A forest-history source for the slash-and-burn and tar figures** that is citable in an
   ecosystem journal (Aakala et al. 2023 covers the qualitative claim already; the area numbers may
   simply be dropped rather than sourced).
3. **⬜ Kumpu et al. 2025 comparability check** (§B caution) before the percentages go in.
4. Consider whether the **north–south** structure of the depletion becomes a stated limitation of the
   single global σ_init (§A).

---

## F. ⭐ The lag leg is now sourced — from Liski et al. 2006 (added 2026-09-04)

§C above was the weak leg. The best citation for it turns out to be the antecedent paper itself
(`../Liski_2006_…pdf`, discussion §4):

> *"soil carbon stock responded slowly to the increased litter production. On the other hand, for
> this same reason, carbon would still accumulate in the litter and soil with no further expansion of
> the forested area if the production of litter is only maintained at the level of 2004 and,
> **centuries later, these carbon stocks would stabilize at a 38% higher level than in 1922**."*

Same country, same model family, coauthors in common: **centuries** to equilibrate, ending **38%**
above the 1922 level, even with litter frozen at 2004. That is the lag claim, stated by the people
whose method this paper extends — and it is far better than a generic chronosequence citation.

⚠ Note what it implies about our own forward arm: a projection with inputs frozen at 2024 that shows
continued accumulation is doing the *same thing they described*, not something new. Keep the two
consistent (memory `forecast-sink-robustness`).

### ⬜ And a lead on the unsourced 1917 stock floor

Liski et al. report a mean soil-and-litter carbon content of **6.1 kg C m⁻² (61 tC ha⁻¹) in 1922**,
rising to 6.3 by 2004. Our calibrated models start 1917 at **30–40 tC ha⁻¹** (runs before the
correlated likelihood) or **44–56** (after it) — i.e. **below** the antecedent's estimate for
essentially the same epoch, on what memory `yasso-depth-basis-matches` says is a comparable
compartment.

**Before this is used as the floor, check three things:** (i) the area basis — 848 Tg / 6.1 kg m⁻²
implies ~14 Mha, so it is not obviously today's forest-land area; (ii) whether "soil and litter"
includes the same deep tail as our `soc_profile`; (iii) that it is an **equilibrium** estimate, hence
an estimate of where the soil *would* sit, not a measurement — which makes it a soft floor at best.
Still the best lead we have; `STOCK_FLOOR = 40` in `init_state_plausibility.R` remains a placeholder.

## Added 2026-09-16 — the post-war pressure (Lorenzo: "no mention of the war reparation pressure")

What the sources support, and one nuance that matters for how it is written:

- **Reparations 1944–1952** (Moscow armistice → last delivery Sept 1952): on average **4% of GDP per year** in goods to the USSR. ⚠ They were paid **mostly in metal and engineering products** — about one third in paper and timber products; the USSR had little interest in wood, which it had itself. So the reparations did not draw down the forest *directly*; they forced the build-up of an export-earning industry, and the forest sector was the earner. Source: **Mitrunen 2024**, *Q. J. Econ.* 140(1):521–584, doi 10.1093/qje/qjae036 (abstract read).
- **Ceded Karelia 1944**: over **12% of Finland's forests**, ~70 sawmills, ~20% of cellulose capacity (secondary: encyclopaedia sources; primary needed).
- **Resettlement / land reform 1940s–50s**: land distributed to ~**half a million** evacuees; selective logging **banned 1947**, clearcutting became the norm; "After the war, Finland needed to rebuild its economy and it had to pay the Soviet Union compensation" — **Kröger 2025**, ch. 8 "Finland's clearcutting forestry", in *Clearcut: Political Economies of Deforestation*, Cambridge UP (fetched summary).
- **The peer-reviewed forest-history source for the whole century: Henttonen, Nöjd, Suvanto, Heikkinen & Mäkinen 2020**, *Size-class structure of the forests of Finland during 1921–2013: a recovery from centuries of exploitation, guided by forest policies*, *Eur. J. For. Res.* 139:279–293, doi 10.1007/s10342-019-01241-y, **CC-BY**. Springer blocks scripted download — ⬜ **download by browser**; it should carry the 1940s–50s fellings-vs-growth and the HKLN/MERA programmes of the 1960s that the "most of it after the late 1960s" clause rests on.

Suggested clause for §"A managed forest landscape": after "…had reached its historical low point (Aakala 2023)":
> The pressure did not end with the war: Finland ceded more than a tenth of its forests in 1944, resettled half a million evacuees through land reform, and paid reparations to the Soviet Union until 1952 that had to be earned largely through the forest-based export industries \citep{Mitrunen2024,Kroger2025}.

### Update 2026-09-21 — the search re-run, and what it changed

- **The NFI record itself says it (Korhonen et al. 2021, *Silva Fennica* 55(5) art. 10662, §3.5.2 and §3.5.3, open access):**
  *"In the 1930's and in the 1950's and 1960's, the use of wood and construction reduced the volume of growing stock. The increase in growing stock volume in the 1940s is explained by the fact that logging was low during the war years."* and *"The slowdown in increment between NFI3 and NFI4 is partly explained by the strong increase in logging in the years before NFI4."* ⭐ This is the citation to use for the post-war drawdown — peer-reviewed, the inventory's own authors, and it needs no political-history claim. ⚠ Korhonen et al. **2024** (the one we cite for the 84%) does NOT carry these sentences; only the 2021 report does.
- **Mitrunen 2024** (*QJE* 140(1):521–584, CC-BY-NC-ND): abstract confirms *"From 1944 to 1952 … on average, 4% of its yearly GDP in industrial products"*, and that *"the more advanced heavy industry … received the majority of state assistance."* Composition (one third forest/paper, two thirds metal) is from secondary sources (Wikipedia; Jensen-Eriksen 2007 is the Finnish primary). Prime Minister's Office press release 343/2012 (2012-11-07): 500 M USD (1953 dollars), 141 490 wagons, last train 19 Sept 1952.
- **Ceded Karelia:** Kochetkova, E. (2020) *Finnish Forestry: from the periphery to the centre of the European forestry*, Encyclopédie d'histoire numérique de l'Europe (EHNE), ISSN 2677-6588, published 22/06/2020, https://ehne.fr/en/node/12493 — *"over 12 percent of all Finnish forests … more than 70 sawmills and 20 percent of Finland's capacity in the cellulose production."* Secondary but signed and dated; acceptable for an Introduction clause.
- Henttonen et al. 2020 fetched (via the Springer HTML): it **jumps from the 1920s to the 1970s** and says nothing about the war decades — drop it as a lead for this point.
- Myllyntaus & Mattila 2002, *Ecological Economics* 41:271–288, *"Decline or increase? The standing timber stock in Finland, 1800–1997"* — the environmental-history reconstruction of the whole two centuries; abstract not retrievable by script. ⬜ optional.

**Decision on wording:** lead with the NFI record (stock fell in the 1950s–60s, logging surged), then name the three drivers (reparations financed through the forest-based export industries, ceded Karelia, resettlement) in one clause with Mitrunen + Kochetkova. Reparations were NOT paid in wood — write "financed", never "paid in timber".
