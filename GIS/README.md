# HIKET GIS — national layers and point extraction

Raw GIS data live on the external drive, in the **GIS store** `$HIKET_GIS_DIR`
(default `/Volumes/NextGenC_SS/HIKET_GIS`). Only the per-plot extractions come
into the project, in `Data/GIS_points/` (gitignored, like the rest of `Data/`).

Plot universe: the 2719 plots of `Data/model_inputs/site_attributes.csv` — the
SOC-independent set, so the extractions do not depend on the calibration layer.
Coordinates are the 2021-23 GPS plot centres (median accuracy 3.2 m, 99th
percentile 9 m).

| step | script | source | store folder | output |
|---|---|---|---|---|
| 01 | `01_download_national.sh` | Metsäkeskus forest use declarations + KEMERA (by region); Luke TWI 16 m | `raw/metsakeskus_*`, `raw/twi_16m` | — |
| 02 | `02_plots_and_dtw_tiles.R` | plot layer; DTW tiles touching a plot (+20 m) | `plots/` | — |
| 03 | `03_download_dtw_tiles.sh` | Luke DTW 2 m, 2023, 5 thresholds, plot tiles only (~999 tiles each) | `raw/dtw_2m` | — |
| 04 | `04_download_emep.sh` | EMEP MSC-W rv5.3 trend runs (Zenodo 12580842) | `raw/emep_ndep` | — |
| 05 | `05_download_history.sh` | Aakala et al. 2023 figshare maps (slash-and-burn, population 1925, ...) | `raw/aakala2023_history` | — |
| 10 | `10_extract_history.R` | | | `history_1925.csv` |
| 11 | `11_extract_wetness.R` | | | `wetness.csv` |
| 12 | `12_extract_emep.R` | | | `n_deposition.csv` |
| 13 | `13_extract_declarations.R` | | | `harvest_declarations.csv` |
| 14 | `14_extract_kemera.R` | | | `kemera_works.csv` |
| 19 | `19_combine_points.R` | | | `gis_point_predictors.csv` |

All downloads resume and skip complete files, so any step can be re-run.

**Caveats found while building (2026-10-05):**
- `site_attributes.csv` repeats 6 plot_ids at two locations each (83471-83473, 83511,
  83552, 83553; none in the calibration set). `19_combine_points.R` drops them →
  **2695 plots × 62 variables**.
- EMEP open archive covers **2010-2022 only**; the 1990-2019 trend runs are on request
  from EMEP MSC-W (david.simpson@met.no).
- Forest use declarations are harvest *intentions*; completion year is empty in the
  open data. Dense from 2004.
- TWI file stores TWI × 1000 as integer (rescaled in `11_`).
- macOS `xargs -I` caps the command at 255 bytes — keep per-item logic in a function.
- Store size: ~73 GB (DTW 48 GB; Metsäkeskus 6.3 GB zipped + ~12 GB unpacked GeoPackages; EMEP 4.5 GB; TWI 2.2 GB).
Licences: Metsäkeskus, Luke, EMEP and figshare data are CC BY 4.0 — cite the
source when used (see each script header).
