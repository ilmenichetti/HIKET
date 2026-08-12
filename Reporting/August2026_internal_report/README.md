# HIKET internal report — August 2026

Working folder for the SOC **data-quality** assessment (the campaign-comparability
thread), kept separate from the manuscript figure set in `manuscript/figures/`.

## Contents

| File | What it is |
|---|---|
| `build_bland_altman_campaigns.R` | builds the figure below; run from the repo root |
| `bland_altman_campaigns.png` | Bland–Altman agreement between the three SOC campaigns, 4 layers × 2 campaign pairs |

## `bland_altman_campaigns.png`

Asks whether two campaigns measuring the same plot agree well enough for their
difference to be read as carbon change. Neither campaign is a gold standard, so
each panel plots the **difference** (later − earlier) against the **average** of
the two, not against one of them — regressing on the earlier value alone would
manufacture a slope out of measurement error, because the difference contains
`−earlier`.

Per panel: solid line = constant bias; dashed = 95% limits of agreement
(mean ± 1.96 SD); tilted grey line = **proportional** bias, the fingerprint of a
scale error (bulk density, sampled volume, a calibration factor) rather than an
added amount of carbon.

Litter (LM) is excluded from the humus row — it exists for 2006/2024 but not for
1985, so including it would compare different quantities. The whole-profile row
is the live calibration target.

### Results (balanced panel, n = 334 plots)

| layer | pair | bias | slope | limits of agreement | band as % of stock |
|---|---|---|---|---|---|
| humus OFH | 1985→2006 | +0.51 | +0.061 **n.s.** | −13.8 … +14.8 | 150% |
| humus OFH | 2006→2024 | −2.49 | −0.089 **n.s.** | −17.4 … +12.4 | 165% |
| mineral 0–20 | 1985→2006 | +1.92 | +0.145 * | −14.4 … +18.2 | 141% |
| mineral 0–20 | 2006→2024 | +3.66 | +0.333 * | −17.6 … +24.9 | 164% |
| mineral 20–40 | 1985→2006 | +2.64 | +0.370 * | −10.8 … +16.0 | 219% |
| mineral 20–40 | 2006→2024 | −0.47 | +0.112 * | −15.5 … +14.5 | 225% |
| whole profile | 1985→2006 | +6.51 | +0.225 * | −30.5 … +43.6 | 113% |
| whole profile | 2006→2024 | +1.93 | +0.146 * | −39.0 … +42.9 | 118% |

**Two readings.**

1. **The humus layer is the only one free of proportional bias**, in either pair
   (n.s. both times), while every mineral layer shows it. An independent
   statistic pointing the same way as Kramarenko (2012 §4.4): organic-layer
   samples are volumetric, mineral samples are mostly spade-cut from a pit wall.
   ⚠ Note this means proportional bias is **not** specific to 1985 — it is a
   property of the mineral measurement, so it should not be used as the argument
   for down-weighting 1985 in particular.
2. **The limits of agreement are the operative number.** Even for the whole
   profile the 95% band is ~115% of the mean stock, and for the subsoil ~220%.
   Plot-level change is not resolvable in any layer or any interval — which is
   an argument about the error model, not about a correction to the mean.

Companion analyses: `doublechecks/soc_depth_distribution.R`,
`doublechecks/subsoil_offset_mechanism.R`,
`doublechecks/organic_mineral_boundary.R`,
`doublechecks/litter_in_1985_organic.R`. Manuscript versions: figures S8–S10.
