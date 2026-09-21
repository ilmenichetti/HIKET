# The σ_input understorey anchor — resolved 2026-09-04

**Question (Open decisions, v3):** the inventory's ground-vegetation litter is $0.51$–$0.67$
tC ha⁻¹ yr⁻¹, but Lehtonen & Heikkinen appeared to imply only $0.19$, and our `sigma_input` prior
centre (**1.08**, arm B, the *decided* anchor) is built on the smaller number. Which is right?

## ⚠ The premise was wrong. The two sources do not disagree.

**Lehtonen & Heikkinen 2015 carry the same understorey number as the inventory — it is in their own
Table A5.** Caption: *"Litter turnover rates for tree biomass by component and species … **For
understorey litter, biomass estimates and their CVs are provided at the bottom of the table.**"*
Bottom rows, under the header `Biomass (kg C·ha⁻¹)`:

| | | |
|---|---|---|
| Understorey | Southern | **506** |
| Understorey | Northern | **666** |

i.e. understorey litter = biomass, turnover 1 yr⁻¹ ⇒ **0.506 / 0.666 tC ha⁻¹ yr⁻¹**, sourced to
Muukkonen & Mäkipää 2006 — **the same numbers as NID Table 6.4-3** (50.6 / 66.6 gC m⁻² a⁻¹). Liski's
implied understorey (2.88 − 2.27 = 0.61) sits inside that range too. **All three agree.**

## Where 1.08 actually came from: a mixed basis

`Prior_specs/*_priors.R` derives arm B as **L&H's absolute total (2.70) ÷ our J̄ (2.511)**. Those are
two *different litter products*. If L&H's total of 2.70 includes their own understorey of ~0.59, then
their tree-only subtotal is ~**2.11** — about 16% below our 2.511 — so the ratio 1.08 is the
understorey correction **minus** the offset between the two products' tree litter. It measures both
at once.

⚠ **The prior file already warns against exactly this**, for Liski: *"Matching their absolute total
instead would give 2.88/2.511 = 1.15 … The RATIO is the right invariant here because sigma_input
exists to correct OUR J̄ for the understorey gap."* The rule was applied to arm A and **not** to arm B.

## The correction, computed directly

`sigma_input` multiplies **our** J̄, so the anchor should add **an understorey estimate** to **our
tree-only flux** — no second product involved:

| south share of plots | understorey (tC ha⁻¹ yr⁻¹) | centre = (2.511 + u)/2.511 |
|---|---|---|
| 90% (our balanced set: 330 / 35) | 0.522 | **1.208** |
| 60% | 0.570 | 1.227 |
| 50% | 0.586 | 1.233 |
| 40% | 0.602 | 1.240 |

**⇒ Anchor ≈ 1.21–1.24 from Lehtonen & Heikkinen, against 1.27 from Liski.** The two anchors that the
prior file records as *"DISAGREE … 27% vs 8% … Unresolved"* **agree to within 5%** once both are used
as ratios. **The 18% disagreement was ours, not theirs, and that note can be closed.**

## What follows

1. **Arm B's 1.08 is an artefact of the mixed basis.** The decision of 2026-09-02 (memory
   `decisions-20260902`: *"the input anchor is Lehtonen & Heikkinen … i.e. arm B"*) was a decision
   about *which source*, and that choice survives — but the *number* it implies is ~1.21, not 1.08.
2. **Arm A and arm B stop being a sensitivity pair.** They were run as competing anchors
   (`20260820_1554*` vs `20260831_1624*`); on the corrected basis they are the same anchor twice.
3. ⚠ **It shrinks the input finding.** The posterior effective correction is 1.27–1.73. Against a
   prior centre of ~1.21 rather than 1.08, most of the ensemble sits **at or just above its prior**,
   and the honest statement becomes *"the product omits a term we can name, and the calibration
   recovers it"* rather than *"the calibration demands more than the product supplies."*
4. ⚠⚠ **Do not adopt 1.21 because it moves the posterior anywhere.** It is adopted because it is the
   arithmetic the prior file already prescribes. The same discipline as *"never re-tune to land MRT on
   a published value."*

## ⬜ Two checks before this is acted on

1. **The provenance of `2.70`.** It is not findable in the L&H text dump; it was probably read off
   their Fig. 1 (total litter input). ⚠ If 2.70 *excludes* understorey, their internal ratio is
   (2.70+0.59)/2.70 = **1.22** — the conclusion is unchanged either way, which is why it is safe to
   state now, but the number should be sourced before it appears in a methods section.
2. **Region coding.** `site_attributes.csv` gives 330 / 35 for `region` 1 / 2; the table above assumes
   **1 = South**. Consistent with the NFI's 1/3 northern sampling density, but confirm.

**Owner for both, and for the 2.70 itself: Aleksi Lehtonen — author of the paper *and* coauthor here.**
