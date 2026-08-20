# Litter input to Finnish forest soil, 1922–2004 — extracted from Liski et al. (2006)

## Source

> Liski J., Lehtonen A., Palosuo T., Peltoniemi M., Eggers T., Muukkonen P., Mäkipää R. (2006).
> **Carbon accumulation in Finland's forests 1922–2004 – an estimate obtained by combination of
> forest inventory data with modelling of biomass, litter and soil.**
> *Annals of Forest Science* **63**: 687–697. DOI [10.1051/forest:2006049](https://doi.org/10.1051/forest:2006049)

PDF archived at `literature/Liski_2006_C_accumulation_Finland_f6070 2.pdf`.
**A. Lehtonen and M. Peltoniemi are coauthors of that paper and collaborators on HIKET** — the
numbers here can and should be checked with them rather than relied on blindly.

## Why we need it

HIKET's transient initialisation ramps litter input from 1917 to 1985 along a fixed **shape**. That
shape was derived from *national total growing stock*, which is the wrong driver twice over: it is a
standing **stock** rather than an **input**, and it is a **national total** rather than **per
hectare**. Liski et al. reconstructed the input to soil directly, including harvest residues — the
term that buffers the trajectory. In their own words (§4.4):

> *"Large harvests … showed a decreasing effect on tree carbon stock but a temporary **increasing**
> effect on litter and soil carbon stock because the residues of the harvests were an important
> source of litter and soil carbon. As a result of these contrasting effects, the compounded carbon
> balance … was **less variable than that of any of the components alone**."*

See `manuscript/M&M_parameterization_working_document.tex` §`sec:preinitshape` and
§`sec:liskicompare`.

## What is in here

| file | contents |
|---|---|
| `liski2006_fig5_input_to_soil.csv` | annual series 1922–2004, Tg C yr⁻¹ |
| `extract_liski_fig5.py` | the extraction (re-runnable; needs `pdftocairo`) |

Columns: `tree_litter`, `harvest_residues`, `ground_vegetation`, `natural_mortality`,
`luc_transfer`, plus two totals:

- ⭐ **`total_input_tree_basis_TgC_yr`** = tree litter + harvest residues + natural mortality.
  **THIS IS THE ONE TO USE** — see Decision 1 below.
- `total_input_all_TgC_yr` = the above plus ground vegetation. Reference only.

⚠ `luc_transfer` is excluded from both — it is a transfer between land uses, not a flux from
vegetation, and our per-plot model has no land-use change.

## How it was obtained — NOT a digitisation

Figure 5 is **vector** in the publisher PDF (confirmed: `pdfimages -list` reports no raster on that
page). The curve vertices are therefore read **exactly** out of the drawing. The only approximation
is the polyline simplification the original plotting software applied (tree litter is drawn with 49
vertices over 83 years, harvest residues with 77).

The five series are separated by **stroke grey level and width**, which are distinct — no visual
judgement is involved:

| series | grey | width |
|---|---|---|
| tree litter | 65.5% | 1.189 |
| harvest residues | 22.4% | 1.189 |
| ground vegetation (dotted) | 13.7% | 1.586 |
| natural mortality (dashed) | 13.7% | 0.793 |
| land-use transfer | 31.0% | 0.396 |

**Axis calibration is also exact**, taken from the tick coordinates in the vector path that draws the
frame: y-ticks at 247.74/227.07/206.06/185.44/164.42/143.80/122.79/102.16 → −5…30 Tg C yr⁻¹;
x-ticks at 347.51…569.18 → 1922…2002. (Page 5 of the PDF = journal p. 691; `pdftocairo` emits
`matrix(1,0,0,-1,0,799)`, a pure y-flip, which the script undoes.)

## Validation — four independent checks

Figure 7 of the same paper gives the 1990s means **per unit area** (kg C m⁻² yr⁻¹), so ratios
between components test the extraction **without needing any area**:

| | extracted | Fig. 7 |
|---|---|---|
| residues / tree litter | 0.391 | 0.400 |
| ground vegetation / tree litter | 0.370 | 0.389 |
| natural mortality / tree litter | 0.035 | 0.038 |

And a fourth, which closes the loop on the absolute scale:

- extracted litter production (trees + ground vegetation) 33.32 Tg ÷ Fig. 7's 0.219 kg C m⁻²
  ⇒ **15.22 M ha**
- their litter+soil stock 959 Tg ÷ their stated density 6.3 kg C m⁻² (2004)
  ⇒ **15.22 M ha**

Two entirely separate routes give the same area to four significant figures.

## The area series, derived the same way

Liski et al. report litter+soil carbon as both a national total and a per-area mean, so the ratio
gives the area their soil results run over:

| year | total | density | implied area |
|---|---|---|---|
| 1922 | 848 Tg | 6.1 kg m⁻² | **13.90 M ha** |
| 2004 | 959 Tg | 6.3 kg m⁻² | **15.22 M ha** |

**+9.5%** — smaller than the "+16% forested area" they quote for *all* forest land, as it should be:
this is the **upland** soil basis, which is our basis. Peatland drainage inflated total forest area
far more than the upland soil area.

## Results

Normalised pre-run shape, 1917–1985 (0 = the 1917 level, 1 = the 1985 level):

| year | ours (growing stock) | Liski total | **Liski per hectare** | linear |
|---|---|---|---|---|
| 1930 | 0.000 | 0.000 | 0.000 | 0.194 |
| 1950 | 0.208 | 0.352 | **0.651** | 0.493 |
| 1960 | 0.182 | 0.689 | **1.000** | 0.642 |
| 1970 | 0.213 | 0.622 | **0.870** | 0.791 |
| 1980 | 0.751 | 1.000 | 1.000 | 0.940 |
| **mean** | **0.194** | 0.336 | **0.461** | 0.500 |

⭐ **The true shape is close to LINEAR** (mean 0.461 vs linear 0.500), not the strongly convex 0.194
we currently impose. It even peaks around 1960 and dips — the 1950s–60s heavy-harvesting decades put
*more* carbon into the soil as residues, exactly when our growing-stock shape claims input was at its
1917 minimum.

⭐ **The implied per-hectare input rise 1922→1985 is +14.6%** on the recommended tree basis
(+11.3% if ground vegetation is included), against **+10.7%** from the wholly independent
NFI-volume-plus-elasticity route. Our posteriors infer **+66 to +97%** (all but SP1).

## Decisions (2026-08-20, validated by Lorenzo)

### Decision 1 — use the TREE BASIS; ground vegetation is excluded

**Not because it is small, but for CONSISTENCY ACROSS 1985.** The post-1985 driver is the Tupek
product: tree litter covering foliage, branches, stem+bark, stumps, coarse and fine roots,
**including** harvest residues and natural mortality, **excluding** understorey. The tree basis here
has exactly that composition. So `sigma_input` applies the *same* understorey correction on both
sides of the join.

Including ground vegetation would make `sigma_input` mean one thing before 1985 and another after —
a discontinuity in a fitted parameter, at the join, that **no diagnostic we run would reveal**.
It also removes the mixed-basis problem, since ground vegetation is the only upland-restricted
series in Fig. 5 while the rest cover all forest land.

(For the record, the shape effect is minor either way: mean 0.461 → 0.487. The consistency argument
is the reason; the insensitivity is incidental.)

⚠ **Residual bias, stated not fixed.** The tree terms still cover all forest land including
peatlands, while we divide by the *upland* soil area. Peatland forests went from largely undrained
in 1922 to ~26% of national volume today, so peatland tree litter grew from near-nothing to a
substantial share, inflating the apparent rise. **Therefore +14.6% is an UPPER BOUND on the true
upland per-hectare rise.** The bias runs against our own argument, which makes the claim
conservative — the real trajectory is flatter, and the gap to the models' +66–97% wider.

### Decision 2 — area timing: 1965–1980 as the working choice, "no correction" as the bound

`none` is **not** the neutral option: constant area is a claim their own numbers falsify
(13.90 → 15.22 M ha). The real choice is between a *sourced* timing and *invented* ones, and only
1965–1980 has a source (Korhonen et al. 2024: "since the mid 1960s … most of this change was before
the 1980s"). Linear and step are conveniences.

Sensitivity on the tree basis:

| timing | mean shape | 1950 | 1970 | rise |
|---|---|---|---|---|
| **1965–1980 (sourced)** | **0.487** | 0.737 | 0.845 | **+14.6%** |
| linear 1922–2004 | 0.330 | 0.429 | 0.566 | +15.7% |
| step at 1970 | 0.462 | 0.737 | 0.439 | +14.6% |
| none (bound) | 0.367 | 0.440 | 0.636 | +24.5% |

Report the envelope **0.367–0.487**: `none` is the state we would retreat to if the area proved
wrong, so it brackets the real uncertainty rather than an arbitrary one.

⭐ **The amplitude is robust to all of this: +14.6 to +24.5% across every timing.** The path detail
depends on the choice; the conclusion that the models' +66–97% is unsupportable does not.

## Caveats

⚠ **Mixed basis, unresolved.** Their §2.2: *"Our results for trees cover all forest land including
upland forests and peatlands, whereas our results for soil, litter and ground vegetation are for
upland forests only."* So `tree_litter` is over a larger area than `ground_vegetation`. This affects
the *composition* of the total, not the validation above (which is ratio-based within Fig. 7's own
basis). **Ask A. Lehtonen.**

⚠ **Polyline simplification.** Tree litter carries 49 vertices for 83 years, so sub-decadal detail is
the plotting software's, not the model's. Fine for a pre-run shape; do not read annual variability
off it.

⚠ **The area growth timing is assumed**, not extracted: 13.90 → 15.10 M ha concentrated 1965–1980,
following Korhonen et al. (2024) — *"since the mid 1960s the area of productive forest has increased
by 1.6 million hectares due to drainage … most of this change was before the 1980s."* Only the
timing matters, since σ_init absorbs the amplitude.

⚠ **Their 1922 state is a steady state** (mean litter 1922–36, climate 1901–30) — the equilibrium
convention this manuscript argues against. It does not contaminate the *input* series, only their
soil results.

⚠ These are **national aggregates**. HIKET is per plot; this series can inform the pre-run *shape*
only, never a per-plot input.
