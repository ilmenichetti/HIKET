# The missing-flux budget — what folders 05–08 have to add up to

Cross-cutting note for `HIKET_storyline_v3.tex` §"What might be missing from the input flux"
(`sec:missingflux`). Each of folders **05** (fine roots), **06** (exudates), **07** (mycorrhizal
necromass) and **08** (early succession) proposes a source for the same displacement. This file states
the size of that displacement once, so no subsection has to guess it, and so the section can be
checked for the obvious failure mode: **candidates that sum to far more than the gap**.

---

## The gap

| Quantity | Value | Source |
|---|---|---|
| Litter input supplied to the models, `J̄` | **2.511** tC ha⁻¹ yr⁻¹ | Tupek product; tree-only, **includes** harvest residues and natural mortality, **excludes** understorey (`build_T1_model_summary.R`) |
| Effective flux after calibration, `σ_input × J̄` | **3.19 – 4.34** | posterior medians, six models, arm B run `20260831_1624*` |
| Implied correction `σ_input` | **1.27 – 1.73** | ditto |
| **THE GAP** | **≈ 0.68 – 1.83 tC ha⁻¹ yr⁻¹** | difference of the two rows above |

⚠ **The gap is a national mean, not a maximum**, and it is inferred **jointly with turnover** — it is
one end of the `MRT × σ_input` ridge (`F14_mrt_ridge`). A mechanism that supplies 0.5 tC ha⁻¹ yr⁻¹
does not "prove" anything; it makes one point on the ridge physically reachable. Say it that way.

⚠ Independent corroboration of the *level*: **Liski et al. 2006** implies `σ_input ≈ 1.27` for the
same country, period and model family — the **floor** of our range (memory `liski-2006-input-benchmark`).

---

## The candidates, against that gap

| # | Candidate | Best current estimate | Covers | Confidence |
|---|---|---|---|---|
| **—** | **Understorey** (in the *inventory*, absent from *our* product) | **0.51 (S) – 0.67 (N)** tC ha⁻¹ yr⁻¹ — and this is now ✅ **agreed across all three sources**: NID Table 6.4-3, Lehtonen & Heikkinen 2015 Table A5, and Liski's implied 0.61, all from Muukkonen & Mäkipää 2006. On our plot distribution, **0.52** | ~75% of the low end | ⭐⭐ **settled** — see `UNDERSTOREY_ANCHOR.md`; the apparent 3× disagreement was **our** mixed basis, not the sources |
| **07** | Mycorrhizal necromass | production "several hundred kg ha⁻¹ yr⁻¹" ⇒ **~0.2–0.5** tC ha⁻¹ yr⁻¹, of which only part is necromass input | up to ~50% of the low end | ⭐⭐ strongest mechanism; magnitude needs primary sources |
| **06** | **Rhizodeposition** ⭐⭐ **(re-rated 2026-09-04)** | A **cropland accounting convention already exists**: Bolinder et al. 2007 prices extra-root material at **65% of root biomass** (the ×1.65 used in ICBM work). Applied to our belowground litter share (≈45% fine roots, ≈54% all roots ⬜) that is **0.73–0.88** — ⚠ the forest analogue is smaller, but by **less than assumed** — ✅ read from the paper 2026-09-04, the 65% is exudates *plus* fine roots and root hairs that **soil sampling misses**, and forest fine-root biomass has the same sampling problem | plausibly the largest single term | ⭐⭐ has a convention AND a measured humification coefficient (Kätterer et al. 2011: root C ~**2.3×** aboveground into refractory SOM) |
| **05** | Fine-root turnover | ⚠ **sign uncertain, and the literature reading points the wrong way**: reconciled MRT ~3.8 yr ⇒ turnover ~0.26 yr⁻¹ against the 0.85 yr⁻¹ in use | possibly **negative** | ⚠ may *widen* the gap — see folder 05 §A |
| **08** | Early-succession vegetation | timing, not total; sign over a rotation unknown | unquantified | ⚠ weakest as stated |

---

## What this table is for

1. **It stops the section from over-claiming.** Understorey + mycorrhizal alone plausibly reach the
   low end of the gap (~0.7). Adding exudates at their high end would overshoot the *high* end (1.83)
   — so the section must not simply list four mechanisms as if they were additive.
2. **It makes the fine-root subsection honest.** One candidate points the other way. A section in
   which every mechanism helps is not credible, and the storyline's own warning box
   (*"This section must not become an argument that we are right"*) is aimed exactly here.
3. **✅ RESOLVED 2026-09-04 — see `UNDERSTOREY_ANCHOR.md`.** The sources agree (0.51–0.67 in all
   three); the 1.08 anchor came from dividing L&H's **absolute total** by **our** J̄, mixing two
   litter products. Corrected on our own plot distribution the anchor is **≈1.21**, against Liski's
   1.27 — so arm A and arm B were never a sensitivity pair. ⬜ Two checks remain (provenance of the
   2.70; region coding), both for Aleksi Lehtonen.


---

## Does the literature justify our sigmas?  (asked 2026-09-04)

**Short answer: it justifies the LOW half of `σ_input`, it does NOT yet justify `σ_init`, and the two
must not be answered together.**

### `σ_input` — yes, at the bottom of the range, and that is a real result

Name the mechanisms and add up only the ones with numbers on a Finnish basis:

| | tC ha⁻¹ yr⁻¹ |
|---|---|
| Understorey, from the **inventory's own** ground-vegetation litter (NID Table 6.4-3) | 0.51 – 0.67 |
| Mycorrhizal necromass, from ECM mycelium production | ~0.2 – 0.5 |
| Rhizodeposition, at the cropland convention ⚠ upper bound for forest | 0.73 – 0.88 |
| **Nameable total** | **≈ 0.7 – 1.2 without rhizodeposition; 1.4 – 2.1 with it at full cropland weight** |
| **The gap** | **0.68 – 1.83** |

So the **low end of our inferred correction is now nameable biology rather than a fudge factor** —
and `σ_input ≈ 1.27` (the bottom of our 1.27–1.73) is exactly where **Liski et al. 2006** lands
independently. Three sources meeting at the same floor is the strongest thing in this review.

⚠⚠ **REVISED 2026-09-04 — the top of the range is now reachable too, and that creates the opposite
problem.** Rhizodeposition, priced at the cropland convention, adds **0.73–0.88** on its own (folder
06 §B). Understorey + mycorrhizal + rhizodeposition would then be **1.4–2.1** against a gap of
**0.68–1.83** — i.e. the named mechanisms would **overshoot**. Two readings, and the section must not
duck between them:

- The forest rhizodeposition term is genuinely smaller than the cropland 65%, because forest chains
  already count fine-root turnover explicitly (folder 06 §B). Then the candidates land inside the
  gap and the account is coherent.
- Or the terms are not additive as stated — some of what the understorey estimate covers is already
  in the mycorrhizal or root terms.

**Either way, the correct presentation is a set of candidates that BRACKETS the gap, not a sum that
matches it.** A budget that adds up exactly would be the more suspicious result.

### ⚠⚠ Three conditions, and the first one bites

1. **This makes the claim SMALLER, not bigger.** If the understorey term is really 0.51–0.67 and not
   0.19, then the `σ_input` **prior centre** should be near **1.25–1.30**, not the 1.08 we ran (arm B).
   The posterior would then sit close to its prior, and the finding stops being *"the calibration
   demands more than the product supplies"* and becomes *"the product omits a term we can name, and
   the calibration recovers it."* That is a **better** paper — anchored in biology, which is the
   weight budget's own instruction — but it is a **weaker** claim than the current framing. Decide
   which one is being made before writing §missingflux.
2. **Do not re-derive the centre to match the posterior.** The centre must come from litter evidence
   alone; whatever distance remains is then a result. This is the same discipline already recorded
   for MRT (*"never re-tune it to land MRT on a published value"*) and for the 1985 weighting (the C5
   circularity). ⬜ Resolve the 0.19-vs-0.51 discrepancy from the sources, not from the fit.
3. **Input and turnover are one package.** On the ridge, justifying more input is *the same act* as
   justifying faster turnover — `MRT × σ_input` is what the level pins. So "the literature supports
   our `σ_input`" is also a statement in support of our short MRT, which is the more contested half.
   Never present the input justification without saying so.

### `σ_init` — no, and it is still the open one

Nothing here supports the low 1917 states. The opposite, if anything: **Liski et al. 2006 report
6.1 kg C m⁻² ≈ 61 tC ha⁻¹ of soil and litter carbon in 1922**, where our models start 1917 at
**30–40** tC ha⁻¹ (pre-correlated-likelihood) or **44–56** (after it) — see
`../04_finnish_management_history/README.md` §F. Their number is an *equilibrium* estimate, so it is
not a floor in the strict sense, and its area and depth basis need checking. But it is the first
quantitative external figure we have for that epoch, it points the same way as the existing concern,
and it comes from the antecedent paper. **The 1917 stock floor stays open; this is its best lead.**
