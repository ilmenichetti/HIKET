# 08 — Early-succession and pioneer vegetation litter

**Serves:** `HIKET_storyline_v3.tex` §"What might be missing from the input flux" →
*Early-succession vegetation* (`:416`): *"Pioneer species after disturbance contribute litter that
inventory-based models, keyed to standing merchantable biomass, do not represent."*

**⚠ Read the claim carefully before reviewing it: it is narrower than the others, and partly wrong as
stated.** The Finnish inventory **does** carry a ground-vegetation litter term (NID 2024, Table
6.4-3: **50.6** gC m⁻² yr⁻¹ aboveground in the south, **66.6** in the north, from Muukkonen et al.
2006). What our *own input product* lacks is understorey altogether — which is the **understorey**
paragraph, not this one. This subsection is specifically about the **post-disturbance transient**:
the years after a clear-cut, when pioneer and early-successional vegetation is at its most productive
and the stand-keyed biomass is near zero.

---

## A. What the literature gives (⬜ largely unverified)

| Source | What it gives |
|---|---|
| **Palviainen et al.** (several, Finnish, 2004–2010) | Ground vegetation after clear-cutting: total biomass **falls 46–65%** immediately, then recovers; *"a relatively large proportion of the annual ground vegetation production returns to the soil as litterfall"*, and after clear-cutting the vegetation **recycles nutrients faster** than in old-growth because a higher share of annual production is shed each year. Also: *Logging residues and ground vegetation in nutrient dynamics of a clear-cut boreal forest*; *Development of ground vegetation biomass and nutrient pools in a clear-cut, disc-plowed boreal forest*. |
| Ground vegetation responses to clear-cutting, first 7 years | Species-level aboveground biomass and nutrient contents through the recovery window. |

⚠ **The direction is not obvious, and this is the subsection most at risk of wishful reading.**
Clear-cutting *reduces* ground vegetation biomass at first; the pioneer flush comes later. Whether
the integral over a rotation is larger or smaller than a stand-keyed model assumes is an open
question — and the *sign* is what the paragraph needs.

⬜ **The decisive check is ours, not the literature's**, and it is cheap: our plots have management
history and stand age. **Do the residuals of the six models depend on time since disturbance?** If
the models under-predict on young stands specifically, that is direct evidence for this mechanism
in our own data — far stronger than a citation. It also connects this subsection to the residual
strand (`F11_rf_heatmap`, where `basal_area_85` is already the top driver in five of six models).

---

## B. Draft text — *Early-succession vegetation*

> Inventory litter models are keyed to standing tree biomass, and the years immediately after a
> stand-replacing harvest are precisely the years in which that key is least informative: tree biomass
> is near zero while the site is neither bare nor unproductive. Ground vegetation biomass falls
> sharply at clear-cutting and then recovers through a pioneer phase in which a large share of annual
> production is returned to the soil each year, faster than in the old-growth stands the coefficients
> are calibrated on (Palviainen et al., 2005). Whether the integral over a rotation exceeds what a
> stand-keyed model assumes is not established, and we raise it as a structural gap in the timing of
> the input rather than as a claimed addition to its total.

---

## C. What is still missing

1. ⬜ **The Palviainen papers themselves** — currently only search summaries. They are Finnish, in
   *Plant and Soil* / *Ecological Research*, and Mikko Peltoniemi will know them.
2. ⬜ **The sign of the effect over a rotation** (§A). Without it this paragraph cannot say anything
   quantitative and should stay short.
3. ⬜ **The residual-vs-stand-age check on our own data** — the strongest available evidence, and it
   is a few hours of work, not a literature search.
4. ⚠ Keep this paragraph **separate from the understorey one**: the inventory has an understorey term
   and our input product does not, which is a different (and better evidenced) argument.
