# Water extension — lake/water layer for masking the SOC maps

Drives the water exclusion in the NextGenC SOC maps
(`Reporting/NextgenC_report/build_soc_maps.R`). Only the scripts/readme here are
tracked in git; the downloaded shapefile and derived cache are ignored (large +
reproducible).

## Source layer

**SYKE Ranta10 / Shoreline 1:10 000 — lake polygons (`Jarvi10`).**

| | |
|---|---|
| File | `Jarvi10.shp` (from `ranta10jarvet.zip`) |
| CRS | ETRS-TM35FIN (EPSG:3067) — same as the SOC maps |
| Content | ~214,000 lake polygons (>200 m²), total ~37,600 km² (Finland's lake area) |
| License | **CC-BY-4.0** — "© SYKE (Ranta10)" |

The MS-NFI peat layer cannot supply water — it sets water to nodata, lumped with
fields/built — so a dedicated waterbody layer is needed. Natural Earth is far too
coarse for Finland (captures <10% of lake area). Ranta10 is the authoritative
national source and is already in EPSG:3067.

## Files

| File | Tracked | What |
|---|---|---|
| `fetch_water_mask.R` | ✓ | Downloads + unzips Ranta10 lakes from SYKE (direct URL). |
| `README.md` | ✓ | this file |
| `Jarvi10.*` | ✗ (gitignored) | the lake shapefile (~480 MB) + `ranta10jarvet.zip` |
| `water_fraction_2000m.tif` | ✗ (gitignored) | cached 2 km water fraction (built by `build_soc_maps.R`) |

## How it's used

1. **Download:** `Rscript Data/Water_extension/fetch_water_mask.R`
2. `build_soc_maps.R` rasterizes the lakes to a 2 km **water fraction** (cached as
   `water_fraction_2000m.tif`) and masks cells **>50% water** (`WATER_CUTOFF`) to
   NA. The sea is already outside the land polygon. In the rendered thumbnails,
   **water shows white** and **peat shows black** (see `make_thumbnails.R` /
   `make_delta_map.R`).

To change the cutoff edit `WATER_CUTOFF` in `build_soc_maps.R`. Deleting
`water_fraction_2000m.tif` forces a rebuild.
