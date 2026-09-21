# 06 — Root exudates and rhizodeposition

**Serves:** `HIKET_storyline_v3.tex` §"What might be missing from the input flux" → *Rhizodeposition*
(`:404`). **Rewritten 2026-09-04** after Lorenzo pointed to the ICBM/Bolinder convention — the first
version of this file called exudates "the weakest of the four candidates". **That was wrong, and the
correction is recorded in §E.**

**Why it is a candidate.** The inventory's litter input is built from **biomass compartments ×
turnover coefficients** (NID 2024, Table 6.4-2: needles, branches, bark, stumps, coarse roots, fine
roots). Rhizodeposition is not a biomass compartment, so no coefficient can express it. It is missing
by construction.

---

## A. ⭐⭐ There is an established quantitative convention, and forest inventories do not apply it

**Bolinder et al. 2007**, *Agric. Ecosyst. Environ.* 118:29–42 — the standard approach for cropland C
input accounting, and the one used in ICBM applications — assigns coefficients to **four** C pools:
crop product, straw, **root biomass**, and **extra-root material**. The extra-root term is set by
*"a coefficient to estimate total belowground C input when assuming that extra-root material
represents **65% of standing root biomass at harvest**"* — i.e. the **× 1.65** multiplier on
root-derived C that Lorenzo used with Kätterer and Bolinder.

⚠ Its own uncertainty is stated in the same literature: rhizodeposition coefficients range from
**below 0.5 to 2.0**, and 0.65 is offered as an approximation, not a measurement.

**And the stabilisation argument is not speculative — it is measured, in the same lineage.**
**Kätterer, Bolinder, Andrén, Kirchmann & Menichetti 2011**, *AEE* 141:184–192 (Ultuna long-term
experiment): the optimised humification coefficient for **root-derived carbon is ~2.3× that of
aboveground residues**. Root-derived carbon is not a low-efficiency input; it is the *high*-efficiency
one.

⭐ Note the position this puts us in: **the convention exists, is decades old, is applied routinely in
cropland carbon accounting — and the forest inventory chain applies nothing analogous.** That is a
much stronger statement than "exudates are a poorly quantified flux".

---

## B. What it would be worth here

Our litter product is resolved by component (`Data/LitterData/tree_litter_per_site_year_by_component_29.04.26.csv`).
Summed over the record, the component shares of tree litter C are:

| component code | share of total litter C |
|---|---|
| 10 | **44.9%** |
| 4 | 33.5% |
| 7 | 8.8% |
| 3 | 7.4% |
| 1 | 4.1% |
| 6 | 1.3% |

✅ **The component key is now resolved without asking anyone (2026-09-04).** The Yasso input bundle
carries an independent woody share, `(fwl+cwl)/J = 0.218`. Components **1, 3, 6, 7 sum to 0.217** —
a match to three decimals — so those four are the **woody** classes (stem, branches, stump, coarse
roots) and the two large non-woody ones are **4 = foliage (33.5%)** and **10 = fine roots (44.9%)**.
⬜ Worth one line of confirmation from Boris, but the arithmetic is not in doubt.

If it holds, a Bolinder-style extra-root factor applied to our `J̄ = 2.511` gives:

| applied to | added flux (tC ha⁻¹ yr⁻¹) |
|---|---|
| fine roots only (0.449 × 0.65) | **0.73** |
| fine + coarse roots (0.449 + 0.088 = 0.537, × 0.65) | **0.88** |

against a gap of **0.68–1.83** (`../MISSING_FLUX_BUDGET.md`). **That is the single largest candidate
in the review** — larger than understorey, larger than mycorrhizal necromass.

### ⚠⚠ But the cropland factor cannot be carried over at face value, and this is the crux

In **cropland**, root biomass is a **snapshot at harvest**. Everything that grew and died during the
season — turnover, root hairs, sloughed cells — *plus* exudates is missing from that snapshot, and
the 65% is compensating for **all of it**. In **forest**, the chain already applies an explicit
**fine-root turnover of 0.85 yr⁻¹** to a standing biomass, so within-year root necromass **is already
counted**. The forest analogue of the extra-root term therefore covers only what remains — exudates
proper, root hairs, mycorrhizal transfer — and **must be smaller than 0.65**.

⇒ **The honest form of the argument is not "apply 1.65 to forest roots" but: cropland accounting
recognises a belowground input beyond root biomass and prices it at 65%; forest inventory accounting
recognises none beyond root turnover; the truth is somewhere between, and the gap our calibration
infers is consistent with it being non-zero.**

### ✅ ANSWERED 2026-09-04 — read from the paper (`Bolinder_2007_AEE.pdf`, §3.3.3)

**How the 0.65 is derived, verbatim:** *"Recent reviews of tracer studies indicates that roughly
**33% of C allocated below-ground in wheat and barley is released by living roots and remains in
soil** (including in the soil microbial biomass), and that **50% of below-ground C remains in roots**
(Kuzyakov and Domanski, 2000; Kuzyakov and Schneckenberger, 2004). This implies that
C_E ≈ 0.65 (i.e., 33/50) × C_R."*

**What C_E contains, verbatim:** *"This material includes **exudates, as well as root hairs and fine
roots sloughed off during the growing season which, because of sampling difficulties, are not
included in the 'root' fraction**."*

⇒ **The cropland-vs-forest objection is weaker than the first version of this file stated.** The 65%
is not mostly the within-season root turnover a forest chain already counts. It is (i) exudation
proper and (ii) fine roots and root hairs that **soil sampling misses** — and forest fine-root biomass
is estimated from soil cores and root:leaf allometry, which has *the same* sampling problem. The
forest analogue is therefore smaller than 0.65 by an unknown and possibly modest amount, not smaller
by construction.

⚠ Two things that keep it honest. The 33/50 comes from **tracer studies on wheat and barley**, and
transfer to boreal trees is not established. And Bolinder et al. say so themselves: *"Our estimates of
extra-root C have high uncertainty reflecting the variability of measured values … a first
approximation."* ⭐ Note also that **0.65 is the conservative end** — earlier work assumed
C_E = C_R (a factor of 2.0), and Kuzyakov & Domanski put pastures at the same 0.65.

⬜ **What remains:** a boreal-forest tracer estimate of the same ratio, and whether Helmisaari's
root:leaf coefficients (our chain's fine-root biomass) come from cores that miss sloughed material.
The second is a question for Boris or Aleksi, not for the literature.

⚠ Also still live: **priming.** Exudates can accelerate decomposition of existing SOM, which would
move the stock the other way. Kätterer's 2.3× is a *net* humification coefficient from a long-term
experiment, so it already contains whatever priming occurred there — which is a good reason to lean
on that number rather than on gross exudation fluxes.

---

## C. Field flux estimates (⬜ unverified, and now secondary)

| Quantity | Value |
|---|---|
| Exudation as a fraction of NPP, forests | **1–17%**; syntheses quote "up to 10%" |
| Boreal Scots pine, upper 15 cm | **~9 g C m⁻² yr⁻¹ = 0.09 tC ha⁻¹ yr⁻¹** (1–2% of NPP) |
| ¹³C labelling, non-boreal | **118.5 g C m⁻² yr⁻¹ = 1.19 tC ha⁻¹ yr⁻¹** (16.7%) |

These span more than an order of magnitude and are worse evidence than §A. **Lead with the
accounting convention and the humification coefficient; use the flux measurements to bound them.**

---

## D. Draft text — *Rhizodeposition*

> Litter models built from inventories are assembled compartment by compartment: a biomass is
> estimated for needles, branches, bark, coarse and fine roots, and each is multiplied by a turnover
> coefficient. Carbon released from living roots as exudates belongs to no compartment and cannot
> appear in such a model at all. This is not an obscure omission. Cropland carbon accounting has
> treated it as a standard term for two decades: the widely used approach of Bolinder et al. (2007)
> carries an explicit "extra-root material" pool, set at 65% of standing root biomass, precisely
> because root biomass alone understates belowground input. Nothing analogous appears in the forest
> inventory chain.
>
> The omission matters more than its size suggests, because root-derived carbon is preferentially
> stabilised. In the Ultuna long-term experiment the humification coefficient optimised for
> root-derived carbon was about 2.3 times that for the same mass of aboveground residues (Kätterer et
> al., 2011). A belowground input that is both unaccounted and disproportionately efficient at forming
> stable organic matter is a candidate of exactly the right kind for a displacement inferred from
> stocks rather than from fluxes.
>
> The cropland coefficient cannot be transferred directly. It compensates for everything missing from
> a single harvest-time snapshot of root mass, including the within-season turnover that a forest
> model already represents explicitly through a fine-root turnover rate. The forest analogue is
> therefore smaller than 65%, and how much smaller is not established. We raise rhizodeposition as a
> structurally missing pathway with a quantitative precedent, not as a calibrated correction.

---

## E. ⚠ Correction to the first version of this file

The first draft ranked exudates last, on the grounds that *"most exudate carbon is respired within
days"*. **That reasoning does not survive contact with Kätterer et al. 2011:** what matters is the
humification coefficient, and for root-derived carbon it is measured to be **higher**, not lower,
than for aboveground litter. The "mostly respired" objection applies to the gross flux and says
nothing about the stabilised fraction — which is the quantity a stock-level displacement depends on.

The residual objection that *does* survive is narrower and is now in §B: exudates enter as
low-molecular-weight, water-soluble material, so representing them as a **uniform multiplier on the
existing AWEN-partitioned litter vector** — which is what `σ_input` is — misstates their quality.
`σ_input` is the right *diagnostic* of a missing input and the wrong *representation* of this
particular one.

---

## F. What is still missing

1. ⬜ **The turnover/exudation split inside Bolinder's 65%** (§B) — the highest-value item.
2. ⬜ **Confirm the component key** with Boris Tupek (§B) — everything quantitative rests on it.
3. ⬜ Read Bolinder et al. 2007 and Kätterer et al. 2011 in full (the second is Lorenzo's own paper;
   the first has an OA copy at `verdeterreprod.fr`).
4. ⬜ A boreal-forest rhizodeposition estimate to check the cropland-derived range against.
