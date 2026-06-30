<!--
Zenodo deposit description for HIKET_SOC_maps_Finland_v1.zip.
Paste the body below into the Zenodo "Description" field. CC-BY-4.0 is set
separately in Zenodo's License field. Keep in sync with SOC_maps/README.md and
SOC_maps_README.tex. Last updated 2026-07-01 (EPSG:3035 / EEA-snapped re-krige).
-->

# Gridded soil organic carbon (SOC) maps of Finland — HIKET / NextGenC

Kriged maps of soil organic carbon stock for Finnish forest land, 1985–2024, from
six soil decomposition models (SP1, TP2, TP3, Yasso07, Yasso15, Yasso20) and their
equal-weight ensemble. Produced by the HIKET Bayesian model-intercomparison project
for the NextGenC reporting effort.

**This is a demonstration / illustrative product** — see *Interpretation &
limitations* below before using the maps quantitatively.

## Contents

- 14 multiband GeoTIFFs: for each model `M` and the `ENSEMBLE`, `M_SOC_mean.tif`
  (mean SOC, tC ha⁻¹) and `M_SOC_sd.tif` (total uncertainty, 1 SD). 40 bands each
  = years 1985–2024 (band *n* = year 1984 + *n*; year also stored as the layer
  name, e.g. `y1985`).
- `ENSEMBLE_SOC_delta10yr.tif` — single-band annual SOC change (tC ha⁻¹ yr⁻¹) over
  2015–2024, the per-cell OLS slope of the ensemble mean (~76 % of cells gaining,
  mean +0.13).
- Per-plot value tables (the numbers behind the maps): `<MODEL>_SOC_mean.csv` /
  `_sd.csv` (447 plots × 40 years) and `Finland_SOC_matrices.ods` (all 12 matrices).
- `thumbnails/` — one PNG preview per raster (1985–2024 time-average; the delta is
  the 2015–2024 change). Peat = black, water/sea = white.
- `README.md` (full file documentation) and `LICENSE`. A methods note is provided
  alongside the archive as `SOC_maps_README.pdf`.

## Raster specification

- CRS: **ETRS89-LAEA Europe (EPSG:3035)**, equal-area. The 2 km grid is snapped to
  the **EEA reference grid** (cell edges on multiples of the resolution from the
  LAEA origin).
- Resolution 2000 m; grid 330 × 589, masked to Finnish forest land; Float32;
  DEFLATE-compressed; sea / outside-Finland = NoData; units tonnes C per hectare.

## Method

Posterior-predictive SOC (mean + SD per plot per year) from each model's Bayesian
calibration is interpolated by ordinary kriging in log space (one pooled variogram
per model), back-transformed to tC ha⁻¹. Per-model SD fuses the kriging
(interpolation) variance with the interpolated model posterior uncertainty; the
ensemble SD applies the law of total variance (within-model + between-model
structural spread). Kriging is performed natively in EPSG:3035. Mineral-soil forest
SOC is targeted: peat/water plots are dropped from the input and cells >50 %
peatland (Luke MS-NFI `paatyyppi`) or >50 % water (SYKE Ranta10 lakes) are masked.

## Interpretation & limitations

Spatial detail is limited by the plot network: plots are ~27 km apart on average and
SOC is spatially weakly autocorrelated at that spacing (variogram ~85–90 % nugget).
The maps therefore reproduce the national mean field with plot-anchored local
deviations; smooth sub-regional gradients are not resolvable from these points, and
detail between plots is interpolation, not measurement. Read the `_sd` layer as the
honest statement of confidence — it inflates between plots and in data-sparse
regions. Point SOC predictability is intrinsically low in this system (calibration
R² ≈ 0.05–0.11). Note: the per-plot CSV/ODS `_sd` is the *parameter* spread on the
predictive mean (observation error not added). The TP3 layers use the forward-Euler
posterior, superseded by an exact matrix-exponential re-calibration; SOC levels
change negligibly.

## Credits & attribution

Project: HIKET — Bayesian calibration and structural intercomparison of SOC models
for the Finnish greenhouse-gas inventory, produced for NextGenC. Litter inputs:
Tupek et al. (Zenodo DOI 10.5281/zenodo.19736499). Mask source data (not
redistributed; please credit if reused): peat — Luke Multi-source NFI (`paatyyppi`),
© Natural Resources Institute Finland, CC-BY-4.0; water — SYKE Ranta10 / Shoreline
1:10 000 lakes, © Finnish Environment Institute (SYKE), CC-BY-4.0. License: Creative
Commons Attribution 4.0 International (CC-BY-4.0). Contact: ilmenichetti@gmail.com ·
GitHub: ilmenichetti/HIKET.
