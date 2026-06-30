#!/usr/bin/env Rscript
# === Fetch SYKE Ranta10 lakes (Jarvi10) — water mask source for the SOC maps ===
#
# Source : SYKE (Finnish Environment Institute) Ranta10 / Shoreline 1:10 000.
# Layer  : Jarvi10 = lake polygons (>200 m2 incl. man-made), ETRS-TM35FIN (EPSG:3067).
#          ~200,000 lakes; the comprehensive national waterbody dataset.
# License: CC-BY-4.0. Attribution: "© SYKE (Ranta10)".
# Used by build_soc_maps.R to mask water cells (lakes; the sea is already outside
# the land polygon) so bodies of water render white, peat renders black.
#
# Usage: Rscript Data/Water_extension/fetch_water_mask.R

PROJ   <- "/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling"
OUTDIR <- file.path(PROJ, "Data", "Water_extension")
URL    <- "https://wwwd3.ymparisto.fi/d3/gis_data/spesific/ranta10jarvet.zip"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

zipf <- file.path(OUTDIR, "ranta10jarvet.zip")
shp  <- file.path(OUTDIR, "Jarvi10.shp")
if (file.exists(shp)) {
  cat("already present:", shp, "\n")
} else {
  options(timeout = 3600)
  cat("downloading Ranta10 lakes (~340 MB) ...\n")
  download.file(URL, zipf, mode = "wb", quiet = FALSE)
  unzip(zipf, exdir = OUTDIR)
  cat("unzipped to:", OUTDIR, "\n")
}
cat("lake polygons:", shp, "\n")
cat("Mask use: rasterize to the SOC grid as water fraction, mask cells > cutoff.\n")
