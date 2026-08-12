# Homogenized SOC baseline for HIKET

Procedural note for `build_soc_homogenized.R`. This folder builds **the** soil-organic-carbon
baseline for HIKET: one internally consistent source for all three Finnish campaigns
(VMI8 1985, Biosoil 2006, Komeetta 2024).

## Why this exists

The previous SOC inputs (`Data/SOC/soilC1985_2006.csv` + `Data/Komeetta/Komeetta_mitatut hiilet.xlsx`)
**over-counted mineral carbon by ~1.6×** because they omitted the coarse-fragment (stoniness)
correction — Finnish forest mineral soils are ~40 % stones by volume. The calibration target was
therefore ~1.7× the official LUKE stock once depth extrapolation was added on top. This rebuild
switches to the LUKE stocks, which are already stoniness/bulk-density corrected and homogeneous
across campaigns.

## Source

`Data/Komeetta/Hannu/Komeetta 150526hi--.xlsx` — Hannu Ilvesniemi's LUKE workbook.

| Sheet | Campaign(s) | What is read |
|---|---|---|
| `BiSo` | 2006 + 2024 | per plot-layer; C stock in **kg/m²** (cols FT/FU), litter in kg/ha (IM/IN), `ORGANIC_LAYER_WEIGHT` (V), plus Biosoil-2006 covariates (texture, BD, coarse-frag, pH, basal area, height, age). Header on row 3. |
| `Data_1985` | 1985 | per plot; organic + 0–5 / 5–20 / 20–40 cm stocks in kg/ha. |

The 2006/2024 cleaning follows J. Heikkinen's `Data/Komeetta/Juha/BiosoilKomeetta15052026.R` and
**reproduces his official weighted means exactly** (2006 = 59 055, 2024 = 61 047 kg C ha⁻¹; see the
validation block printed at build time).

## Processing rules

1. **Peat excluded** — plots with `CODE_LAYER` H01/H12 (Heikkinen's rule), plus `site_raw$peatland`
   for 1985-only plots not in BiSo. (Peat gets its own extension, not this mineral-soil baseline.)
2. **Layers kept in their original intervals** — organic + 0–5/5–20/20–40 (1985) and
   organic + 0–10/10–20/20–40 (2006/24). The **measured 40–80 cm layer (Krs 204)** is kept *out of the
   homogeneous stack* (it exists only for 2006) but retained as `soc_40_80_meas` to validate the tail.
3. **Depth: measured + extrapolated, capped at each plot's soil depth.**
   - `soc_0_40` — measured organic + mineral to the deepest measured layer (≤40 cm).
   - `soc_profile` (**primary**) — `soc_0_40` + an exponential extrapolation of the deep tail, integrated
     **only down to `z_cap = min(100 cm, plot soil depth)`**, so it never extends past actual soil.
   - `soc_1m_uncapped` — the old fixed-100 cm version, kept for comparison.
   - **Per-plot soil depth** = `soil_depth_reached` (`LAYER_LIMIT_INF` of the deepest augered layer):
     80 cm = reached the standard target (soil ≥ 80 → extrapolate to 100 cm); 10/20/40 cm = augering
     **refusal** (bedrock/stones → soil ≈ that depth → **no** deep extrapolation). This removes the
     ~14 Mg ha⁻¹ of carbon that a fixed-depth extrapolation was inventing below 38 thin-soil plots.
   - Extrapolation ported from `Data_work.R`: one λ per GTK soil class (`site_raw$soil_code`), global
     fallback (λ ≈ 0.032), C₀ per profile. **Validated** against the measured 40–80 cm layer (Krs 204,
     n=460): median predicted/observed **0.96**, mean bias −2.5 Mg ha⁻¹ — slightly conservative, and
     *under*-predicts (not over) the deep layer in carbon-rich soils. Adds ~10 Mg ha⁻¹ on deep plots.
4. **Incomplete-2024 profiles dropped** (abandoned / all-C%-missing / missing organic or topmost
   mineral) rather than zero-filled — zero-filling would deflate the 2024 mean. Their valid 2006/1985
   rows are **kept**, with `abandoned_2024` / `unmeasured_2024` / `missing_C2024` flags explaining why.
5. **Litter** is folded into the organic layer for 2006/24 (Heikkinen). 1985 has no litter component
   (minor, sub-Mg ha⁻¹) — a documented, small cross-campaign inconsistency.

## Covariates — year-resolved (stand properties vary over time!)

Stand properties (basal area, age, height, development class) change between campaigns, so
**every covariate column carries its measurement era as a suffix** (`_1985` / `_2006` / `_2024`) —
the time is explicit in the header. Time-varying stand properties get one column per era (wide);
they are plot-level, so they repeat across a plot's campaign rows (a 1985 SOC row still shows the
plot's 2006 and 2024 stand columns, and vice-versa). Column → source → coverage:

| Column | Source | Coverage |
|---|---|---|
| `basal_area_1985`, `dev_class_1985` | site_raw NFI8 (`basal_area_85`, `dev_class_85`) | 83% |
| `stand_age_1985`, `mean_height_1985_m` | site_raw NFI8 (`stand_age_85`, `mean_height_85_dm`/10) | 15% (sparse) |
| `basal_area_2006`, `stand_age_2006`, `mean_height_2006_m` | BiSo `Ppa` / `IKA` / `Keskipituus` | 95% |
| `mgmt_op_2006`, `mgmt_yrs_since_2006`, `mgmt_residue_2006` | BiSo management | 95% |
| `dev_class_2024` | Kohdekuvaukset `KEHLK` | 99% |
| `site_main_class_2024`, `site_fertility_2024`, `dominant_species_2024` | Kohdekuvaukset | 95% |
| `clay_2006`, `silt_2006`, `sand_2006`, `bd_est_2006`, `coarse_frag_2006`, `pH_CaCl2_2006`, `texture_class_2006` | BiSo soil physics (≈ static) | 73–100% |
| coordinates, `soil_code` | site_raw | 95% |

No 2006 dev-class (`Puuston_tila` unreliable) and no 2024 basal-area/age/height (Komeetta measured
no stand inventory), hence those gaps. Caveats: (1) 1985 `stand_age`/`mean_height` are **sparse**
(~98 plots) and 1985 height looks low — verify the dm unit with LUKE. (2) `dev_class` uses different
scales across eras (1985 NFI vs 2024 KEHLK) — compare within an era, not across.

## Region — the cross-sampling consistency issue

Region is **not** in the SOC workbook, and the classification differs by data source, so it is
resolved to **one value per plot** with an explicit provenance column (`region_source`):

1. `biosoil_design` — `Data/Komeetta/Juha/region.csv` (LUKE Biosoil design; 1 = South w = 1,
   2 = North w = 3). Covers all 2006/24 plots (~508).
2. `site_raw_region` — `site_raw$region` where the design file is silent.
3. `latitude_ETRS` / `latitude_sitekey` — northing ≥ 7 300 000 m ⇒ North, else South, using
   `site_raw$y_ETRS` then the `soil_litter_site_key.csv` northing (keyed by `koealatunnus_BIOSOIL`; ~59).
4. `default_South` — last resort (≈10 plots, all 1985-only with no geographic info anywhere).

`region_conflict = TRUE` marks the 2 plots where the design and `site_raw` disagree (design wins,
because the design weights reproduce the official means). **Weights** (`weight` = 3 North / 1 South)
are the Biosoil design; the same weights are assumed for 1985 (documented assumption — the 1985 NFI
design weights were not available separately).

## Outlier flags (thresholds are named constants at the top of the script)

Two documented flags; nothing is silently dropped.

- **`soc_outlier`** (per plot-year) — `soc_profile > SOC_OUTLIER_MAX` (**250 Mg ha⁻¹**): physically
  implausible for boreal mineral soil to 1 m (paludified/gley sites or data errors). **6 plot-years,
  4 plots.** Intended to be **excluded from calibration** (wired in `Data_work` when the baseline is
  adopted). E.g. plot 33631 sits at 290–417 Mg ha⁻¹ in all three campaigns.
- **`high_change`** (per plot) — any consecutive-campaign change `> HIGH_CHANGE_RATE` (**3 tC ha⁻¹ yr⁻¹**):
  **23 plots.** **Flagged but KEPT** — the campaigns don't re-core the identical soil volume, so large
  apparent jumps are largely resampling noise, not errors; hard-dropping them would bias the set toward
  artificial stability. Left for the likelihood to absorb; revisit only if they dominate the misfit.

Change either threshold in one line at the top of `build_soc_homogenized.R` and rebuild.

## Outputs

| File | Contents |
|---|---|
| `soc_homogenized_layers.csv` | long, original layer intervals — drop-in replacement for `soilC1985_2006.csv` (now with 2024). `plot_id, year, campaign, layer, depth_lower_cm, depth_upper_cm, C_kgha, C_Mgha`. |
| `soc_homogenized_plot.csv` | per plot × campaign. Stocks (kg/ha **and** `*_Mgha`): `organic`, `mineral_0_40`, `soc_0_40` (measured), **`soc_profile`** (primary, depth-capped whole profile), `soc_deep` (its extrapolated addition), `soc_1m_uncapped` (old fixed-100 cm), `soc_40_80_meas` (measured Krs 204, 2006), `soc_40_80_pred` (extrapolated, for validation). Plus `soil_depth_reached`, `z_cap`, `lambda`, `fit_ok`; `soc_outlier` / `high_change` (see Outlier flags); region/weight/provenance; year-tagged covariates; static soil + site descriptors; and all data flags. |
| `soc_homogenized.rds` | bundle `{layers, plot, meta}` (meta carries λ values + validation). |
| `plots/*.png` | 01 overall trajectory (balanced), 02 by region, 03 layers, 04 depth contribution, 05 coverage/provenance, 06 stock distributions, 07 extrapolation validation vs measured 40–80 cm, 08 new vs old (inflated) target, 09 per-plot ΔSOC, 10 maps of Finland, 11 SOC vs drivers (T/age/stoniness/fertility), 12 organic-layer share, 13 outliers (highest stocks + largest changes). |

## Headline numbers (built baseline)

- Validation vs Heikkinen: **exact** (2006 = 59 055, 2024 = 61 047 kg ha⁻¹, n = 446).
- Extrapolation validated vs measured 40–80 cm (Krs 204, n = 460): median pred/obs **0.96**, bias −2.5 Mg ha⁻¹.
- Balanced trajectory (present in all 3 campaigns, n = 344), weighted Mg ha⁻¹:
  - 0–40 cm (measured): **51.8 → 58.9 → 59.9**
  - `soc_profile` (depth-capped): **62.1 → 70.5 → 72.5**  (accumulation → incipient saturation)
- Per-plot depth cap zeroes the deep extrapolation for 38 thin-soil plots (reached ≤ 40 cm), removing
  ~14 Mg ha⁻¹ each of previously-invented sub-soil carbon; aggregate effect small (deep plots unchanged).
- Coverage: 577 plots (VMI8 460 / Biosoil 503 / Komeetta 445; 344 with all three).

## Plugging into the pipeline — ✅ DONE (2026-08-04)

This baseline is now **wired into `Data/Data_work.R`** and is the live calibration target.
What changed there (backup of the previous version: `Data/Data_work_pre_SOC_swap_20260804.R`):

- **§1** reads `soc_homogenized_layers.csv` instead of `./Data/SOC/soilC1985_2006.csv`
  (same `plot_id / year / layer / C_kgha` schema, all three campaigns, so the organic-quality
  flags in §1.1 work unchanged). `soilC1985_2006.csv` is no longer read.
- **§1.0b removed** — the separate Komeetta 2024 ingest (and its mislabelled layer map) is gone;
  the layers file already carries 2024 on the corrected basis.
- **§1.3b** no longer fits its own depth model. It merges `soc_profile_Mgha` from
  `soc_homogenized_plot.csv` as `soc_obs_tCha`. `soc_obs_tCha_sum` (measured organic + 0–40)
  is retained and reproduces this file's `soc_0_40_Mgha` to **0.000 tC/ha** — the check that
  both tables describe the same data.
- **`soc_outlier` is applied as a calibration exclusion** in `calib_ready`, documented at the
  point of use. `high_change` is carried into `site_raw` but those plots are **retained**.

Verified after the swap: Data_work runs clean, **66/66 sanity checks pass**, calibration-ready
**520 plots** (416 calibration / 104 holdout, up from 447 — the baseline covers more plots),
no plot-year falls back to the raw sum, and `soc_obs_tCha` matches `soc_profile_Mgha` exactly.
Campaign medians of the live target: **59.4 / 66.7 / 67.4 tC/ha** (previously 63 / 102 / 105).
The observation CV used for the fixed likelihood σ drops from 0.472 to **0.440**.

Prior-pushforward pre-flight re-run on the new target (see `preflight_prior_pushforward.R`,
K=300 draws at full N): forward sanity **PASS in all six models**, blow-up rate **0 %** for
SP1/TP2/TP3/Yasso15/Yasso20 and **0.6 %** for Yasso07 (1 of 156 constraint-passing draws).
**TP3 improves from 10–12 % to 0 %.** So the ICBM rate anchors and the `flux_pair` σ_input
bounds — both tuned on the old, inflated target — remain safe on the gentler one, and the
feared σ_input lower-bound clamp did not appear (window [0.021, 3.647], J̄ = 2.386 tC/ha/yr).
This is the SOC change that triggers the six-model Roihu re-run.

## Rebuild

```bash
Rscript Data/SOC_homogeneized/build_soc_homogenized.R   # run from repo root
```
Deterministic (`set.seed(2025)`). Requires `readxl`, `dplyr`, `tidyr`.
