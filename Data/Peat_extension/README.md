# Peat extension — peatland layer for masking the SOC maps

Drives the peatland exclusion in the NextGenC SOC maps
(`Reporting/NextgenC_report/build_soc_maps.R`). Only the scripts/readme in this
folder are tracked in git; the downloaded raster and derived cache are ignored
(large + reproducible).

## Source layer

**Luke Multi-source NFI (MVMI / MS-NFI), variable `paatyyppi` — site main class.**

| | |
|---|---|
| File | `paatyyppi_vmi1x_1519.tif` (2019; `1519` = VMI year code) |
| Resolution / CRS | 16 m, ETRS-TM35FIN (EPSG:3067) — same as the SOC maps |
| Coverage | whole country, single GeoTIFF (~249 MB) |
| Classes | **1 = mineral soil, 2 = spruce mire, 3 = pine mire, 4 = treeless peatland**; 32766 = non-forest/nodata |
| Peat definition | organic layer is peat **or** >75% peatland vegetation → captures drained peatland forest (*turvekangas*) |
| Accuracy | ~84% overall (mineral 95% user acc.) |
| License | **CC-BY-4.0** — "© Natural Resources Institute Finland, 2019" |

Why Luke and not CORINE/GTK: CORINE is land *cover* and misses drained forested
peat; GTK soil map catches it via a depth threshold. Luke `paatyyppi` catches it
by *site type*, the natural target for excluding peat from a forest SOC map, and
is already in EPSG:3067 + CC-BY-4.0.

## Files

| File | Tracked | What |
|---|---|---|
| `fetch_peatland_mask.R` | ✓ | Reproducible download via the Luke INSPIRE Atom service (service feed → year feed → `paatyyppi` URL). Arg = year (default 2019). |
| `README.md` | ✓ | this file |
| `paatyyppi_vmi1x_1519.tif` | ✗ (gitignored) | the 16 m source raster |
| `peat_fraction_2000m.tif` | ✗ (gitignored) | cached 2 km peat fraction (built by `build_soc_maps.R`) |

## How it's used

1. **Download:** `Rscript Data/Peat_extension/fetch_peatland_mask.R 2019`
2. `build_soc_maps.R` then:
   - **drops** peat/water calibration plots from the kriging input
     (`site_raw$soil_code ∈ {Tu, Ve}`), and
   - builds a 2 km **peat fraction** = peat pixels (classes 2/3/4) ÷ **all** 16 m
     pixels in the cell (water included in the denominator, so lake-dominated cells
     don't read as peat), cached as `peat_fraction_2000m.tif`, and **masks cells
     >50% peat** (`PEAT_CUTOFF`) to NA in every output map.

To change the year, peat-fraction cutoff, or excluded soil codes, edit the
`PEAT_*` constants at the top of `build_soc_maps.R` (and re-run the fetch for a
different year). Deleting `peat_fraction_2000m.tif` forces a rebuild.
