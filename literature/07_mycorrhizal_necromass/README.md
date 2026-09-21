# 07 — Mycorrhizal necromass as a precursor of boreal soil carbon

**Serves:** `HIKET_storyline_v3.tex` §"What might be missing from the input flux" → *Mycorrhizal
necromass* (`:407`): *"increasingly implicated as a dominant precursor of stable soil carbon in
boreal systems, and entirely outside the litter product"* — with the Clemmensen line of evidence
named in the annotation.

**⭐ This is the strongest of the four candidates**, because unlike the others it has (i) a
high-profile boreal result asserting *dominance*, (ii) production estimates in the right units, and
(iii) an unambiguous absence from the inventory chain: mycelium is not a tree biomass compartment, so
no turnover coefficient can produce it.

---

## A. The evidence

| Source | What it establishes |
|---|---|
| **Clemmensen et al. 2013**, *Science* 339:1615–1618 ★★ | The headline. Bomb-¹⁴C modelling across a chronosequence of boreal forested islands (N Sweden): **50–70% of stored soil carbon derives from roots and root-associated microorganisms**, not from aboveground litter. *"Mycorrhizal mycelium plays a fundamental role in carbon sequestration by boreal forests."* Fungal biomarkers indicate preservation of fungal residues in late-successional forests. |
| **Clemmensen et al. 2015**, *New Phytologist* | Carbon sequestration related to **mycorrhizal fungal community shifts** — i.e. the flux is not a constant per unit biomass, which matters for a model that would want one. |
| **Ekblad et al. 2013**, *Plant and Soil* | ECM mycelial production in coniferous forests **"several hundred kg ha⁻¹ yr⁻¹"** in the upper soil. ⬜ get the stated range. |
| **Hagenbo et al. 2017**, *New Phytologist* | Up to **1.2 kg mycelium ha⁻¹ day⁻¹** in young Scots pine (≈ **438 kg ha⁻¹ yr⁻¹** mycelium mass, i.e. roughly **0.2 tC ha⁻¹ yr⁻¹** at ~50% C — ⬜ verify the carbon conversion). Across a *P. sylvestris* chronosequence, **turnover rather than production** regulates standing mycelial biomass — so a young, growing forest is not simply a bigger version of an old one. |
| **Ectomycorrhizal necromass turnover** (recent, *Plants People Planet* / DiVA) | Necromass turnover reported as **one-third of** biomass turnover — the step from mycelium production to *necromass input* is not one-to-one. ⬜ get the primary. |

---

## B. Why it fits our result better than the alternatives

1. **It is absent by construction**, like exudates — but unlike exudates it is *particulate necromass*,
   chemically closer to what a litter model represents, so adding it to `J` is less of a category
   error.
2. **Its magnitude is in the right band.** Several hundred kg ha⁻¹ yr⁻¹ of production is 0.2–0.5 tC
   ha⁻¹ yr⁻¹, a meaningful fraction of the **0.68–1.83 tC ha⁻¹ yr⁻¹** gap
   (`../MISSING_FLUX_BUDGET.md`) though not, by itself, all of it.
3. **It has a stabilisation argument attached**, which no other candidate has: Clemmensen's evidence
   is specifically about carbon that *persisted*, from bomb-¹⁴C, not about a gross flux.

⚠ **The trap to avoid.** Clemmensen's 50–70% is a statement about the **origin of stored carbon**,
not about a missing input to a litter model. Some of that root-derived carbon *is* already in our
input as fine-root litter. Quoting the 50–70% as though it were the size of a missing flux would be a
straightforward misreading, and an expert reader will catch it. The defensible sentence is: *the
dominant precursor of stable boreal soil carbon is root-derived and partly fungal, and the fungal
component is entirely outside inventory litter models.*

⚠ Second trap: the mycelium/necromass distinction (§A last row) — production is not input.

---

## C. Draft text — *Mycorrhizal necromass*

> Where boreal soil carbon comes from has been substantially revised over the last decade. Using
> bomb-radiocarbon across a chronosequence of forested islands, Clemmensen et al. (2013) attributed
> 50 to 70 percent of the stored carbon to roots and root-associated organisms rather than to
> aboveground litter, and identified mycorrhizal mycelium as a fundamental agent of the sequestration.
> Production of extramatrical mycelium in coniferous stands is measured in hundreds of kilograms per
> hectare per year (Ekblad et al., 2013; Hagenbo et al., 2017), of which only part becomes necromass
> input.
>
> None of this appears in an inventory litter model. Such models are assembled from tree biomass
> compartments and their turnover coefficients, and fungal mycelium is not one of them; the flux is
> therefore missing by construction rather than by underestimation. Its magnitude — of the order of a
> few tenths of a tonne of carbon per hectare per year — is a substantial fraction of the input
> displacement inferred here without accounting for all of it, and the pathway carries the additional
> property, unusual among the candidates, that the evidence for it is evidence about carbon that
> persisted rather than about carbon that entered.

---

## D. What is still missing

1. ⬜ **Clemmensen et al. 2013** — read in full before quoting; paywalled (*Science*), but Lorenzo's
   institution will have it. **Blocking** for this subsection.
2. ⬜ **Ekblad et al. 2013** and **Hagenbo et al. 2017** — the actual production ranges, in carbon.
3. ⬜ The **necromass fraction** of production (§A last row), which is what a litter model would need.
4. ⬜ Check whether Yasso's own parameterisation datasets contain any fungal material — if the model
   was calibrated on litterbag data only, the fungal pathway is missing from *both* the input and the
   decomposition parameters, which is a sharper statement than the input one alone.
