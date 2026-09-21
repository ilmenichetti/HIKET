# 05 — Fine root mean residence time in boreal conifers

**Serves:** `HIKET_storyline_v3.tex` §"What might be missing from the input flux" → *Belowground
turnover* (`:400`): *"fine root MRT in boreal conifers: isotopic vs minirhizotron estimates, and the
factor-of-several disagreement between methods."*

**The quantity at stake.** The Finnish inventory (and, through the same lineage, the litter product we
use) builds fine-root litter as `biomass × turnover rate`, with the **turnover rate fixed at 0.85
yr⁻¹ for every species and region** (NID 2024, Table 6.4-2; see `../01_ghg_inventory_chain/`) and
fine-root biomass from root:leaf coefficients (Helmisaari et al. 2007). Both factors are uncertain,
and fine roots are the largest belowground litter term. See `../MISSING_FLUX_BUDGET.md` for the size
of gap this would have to close.

---

## A. The methodological disagreement

| Method | What it gives for boreal conifers | Source |
|---|---|---|
| **Minirhizotron** (direct observation of individual roots) | Median lifespan **< 2.5 yr** for Norway spruce | Strand et al. 2008, *Science* — the "irreconcilable differences" paper |
| **Radiocarbon** (age of C in root tissue) | Apparent ages **from recently grown to ~14 yr**, same species | Strand et al. 2008; Gaudinski et al. |
| **Reconciliation attempt** | Two populations: **60 ± 2%** turning over within **0.75 ± 0.10 yr**, the remainder within **8.4 ± 0.2 yr**; overall MRT **3.80 ± 0.16 yr** | Hansson et al., "Reconcilable differences: a joint calibration of fine-root turnover", *New Phytologist* |
| **The likely artefact** | Root cellulose contains carbon older than the year of formation — roots are built partly from **stored reserves**, so ¹⁴C age overestimates root *age* | "Old carbon in young fine roots in boreal forests"; Sah et al., *Plant and Soil* |

⚠⚠ **THIS EVIDENCE POINTS THE WRONG WAY FOR US, AND THAT MUST BE SAID.** A reconciled MRT of ~3.8 yr
implies a turnover of **~0.26 yr⁻¹**, roughly **a third** of the 0.85 yr⁻¹ the inventory applies.
Taken at face value it would make fine-root litter *smaller*, **widening** the gap our calibration is
trying to explain rather than closing it.

**Two honest readings; the section must carry both.**
1. A mean lifespan is the wrong summary for a *flux*. The flux is dominated by the short-lived 60%
   (0.75 yr ⇒ ~1.3 yr⁻¹), so 0.85 yr⁻¹ may be a defensible flux-weighted value even where the mean
   MRT is 3.8 yr. **⬜ Check this arithmetic against Hansson's own numbers before writing it down.**
2. If the long-MRT reading is right, fine roots are *not* where the missing flux is, and the section
   should say so and move to the other candidates. **That would be a better paper than one in which
   every mechanism conveniently helps.**

⚠ The *biomass* factor is independently uncertain (root:leaf allometry; whether deep fine roots are
captured at all) and can move the flux in either direction.

---

## B. Additional Finnish/Nordic anchors (⬜ all need the primary PDF)

- **Leppälammi-Kujansuu et al.** — Norway spruce fine-root turnover and litter production under
  long-term temperature and nutrient manipulation: warming and fertility **cut median longevity and
  raised belowground litter production three- to fourfold**, particularly into the *mineral* soil.
  ⭐ Directly relevant to a warming-driven, non-stationary input — i.e. to the forward arm — and it is
  Helmisaari's group, so it is close to the source of the coefficients we use.
- *Fine-root turnover rates of European forests revisited* (*Plant and Soil*, 2013) — may supply a
  defensible **range** to replace the single 0.85.
- Fine-root longevity and below/aboveground litter production in a boreal *Betula pendula* forest.

---

## C. Draft text — *Belowground turnover*

> Fine roots are the largest belowground litter term in the inventory chain, and they enter it as a
> product of two uncertain factors: a biomass estimated from root-to-leaf allometry, and a turnover
> rate applied as a single constant, 0.85 yr⁻¹, across species and regions (Statistics Finland, 2024).
> The lifespan behind that constant is among the least settled quantities in boreal ecology. Direct
> observation with minirhizotrons puts the median lifespan of Norway spruce fine roots below two and a
> half years, while radiocarbon dating of the same species returns apparent ages up to fourteen
> (Strand et al., 2008) — a disagreement too large for both methods to be describing the same
> population. Part of it resolves once roots are recognised to be built partly from stored carbon, so
> that their carbon is older than they are; part once the population is split into a short-lived
> majority turning over within a year and a long-lived remainder persisting for most of a decade.
>
> The consequence for a litter model is not only that the flux is uncertain, but that a single mean
> lifespan is the wrong summary of a population whose flux is carried by its fastest members, and that
> the mass and rate factors are uncertain in different ways. We flag it as the largest identified
> uncertainty in the belowground input, not as a resolution of the displacement inferred here.

---

## D. What is still missing

1. ⬜ **Strand et al. 2008** (*Science* 319:456–458) — the anchor citation, not yet read.
2. ⬜ **Hansson et al.** — needed for the arithmetic in §A, which decides whether this mechanism is a
   candidate at all. **Blocking for this subsection.**
3. ⬜ The provenance of **0.85 yr⁻¹** in the Finnish chain: whose measurement, on what, with what
   stated uncertainty. Our input product depends on it. **Boris Tupek or Aleksi Lehtonen can answer
   this in one email, faster than the literature can.**
