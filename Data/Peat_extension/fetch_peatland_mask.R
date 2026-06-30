#!/usr/bin/env Rscript
# === Fetch Luke MS-NFI 'paatyyppi' (site main class) raster — peatland mask source ===
#
# Reproducible download via the Luke INSPIRE Atom service (machine-readable):
#   service feed  -> per-year dataset feed -> 'paatyyppi' GeoTIFF download URL.
#
# Source : Luke Multi-source National Forest Inventory (MVMI / MS-NFI).
# Layer  : 'paatyyppi' = Kasvupaikan päätyyppi / site main class.
#          16 m, ETRS-TM35FIN (EPSG:3067), single country-wide GeoTIFF (~249 MB).
# Values : 1 = mineral soil, 2 = spruce mire, 3 = pine mire, 4 = treeless peatland.
# Peat definition: organic layer is peat OR >75% peatland vegetation
#          -> captures drained peatland forest (turvekangas). ~84% overall accuracy.
# License: CC-BY-4.0. Attribution: "© Natural Resources Institute Finland, <year>".
#
# PEATLAND MASK for the SOC maps:
#   p  <- terra::rast(dest)                 # 16 m site-main-class
#   peat <- p >= 2                          # 1 = peat (classes 2/3/4), 0 = mineral
#   # aggregate to the 2 km SOC grid as peat fraction, then threshold (e.g. >0.5):
#   pf <- terra::resample(peat, soc_grid, method = "average")
#   soc_masked <- terra::mask(soc_grid, pf > 0.5, maskvalues = TRUE)
#
# Usage: Rscript Data/Peat_extension/fetch_peatland_mask.R [YEAR]
#   YEAR in {2006,2009,2011,2013,2015,2017,2019,2021,2023}; default 2019.

suppressPackageStartupMessages(library(xml2))

args    <- commandArgs(trailingOnly = TRUE)
YEAR    <- if (length(args) >= 1) as.integer(args[1]) else 2019L
SERVICE <- "https://kartta.luke.fi/inspireatom/mvmi.xml"
PROJ    <- "/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling"
OUTDIR  <- file.path(PROJ, "Data", "Peat_extension")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# namespace-agnostic helpers (Atom feeds use the atom: default namespace)
lf <- function(node, name) xml_find_all(node, sprintf(".//*[local-name()='%s']", name))

# 1. service feed -> the dataset feed for YEAR
svc     <- read_xml(SERVICE)
entries <- lf(svc, "entry")
titles  <- vapply(entries, function(e) xml_text(lf(e, "title")[[1]]), character(1))
i <- grep(as.character(YEAR), titles)
if (!length(i))
  stop("No MVMI feed entry for ", YEAR, " (have: ", paste(titles, collapse = ", "), ")")
links  <- lf(entries[[i[1]]], "link")
ds_url <- xml_attr(links, "href")[grepl("\\.xml$", xml_attr(links, "href"))][1]
cat("MVMI", YEAR, "dataset feed:\n  ", ds_url, "\n", sep = "")

# 2. dataset feed -> 'paatyyppi' .tif URL
ds     <- read_xml(ds_url)
dhrefs <- xml_attr(lf(ds, "link"), "href")
tif    <- dhrefs[grepl("paatyyppi", dhrefs, ignore.case = TRUE) &
                 grepl("\\.tif$", dhrefs)][1]
if (is.na(tif))                                   # fallback: match by entry title
  for (e in lf(ds, "entry")) {
    t <- xml_text(lf(e, "title")[[1]])
    if (grepl("p.{0,2}tyyppi|site (main|fertility) class", t, ignore.case = TRUE)) {
      h <- xml_attr(lf(e, "link"), "href"); tif <- h[grepl("\\.tif$", h)][1]; break
    }
  }
if (is.na(tif)) stop("Could not locate 'paatyyppi' .tif in dataset feed: ", ds_url)
cat("paatyyppi download URL:\n  ", tif, "\n", sep = "")

# 3. download (skip if already present)
dest <- file.path(OUTDIR, basename(tif))
if (file.exists(dest) && file.size(dest) > 1e8) {
  cat(sprintf("already present: %s (%.0f MB)\n", dest, file.size(dest) / 1e6))
} else {
  options(timeout = 3600)
  cat("downloading ~249 MB ...\n")
  download.file(tif, dest, mode = "wb", quiet = FALSE)
}
cat(sprintf("saved: %s (%.0f MB)\n", dest, file.size(dest) / 1e6))
cat("\nMask rule: terra::rast(dest) >= 2  (keep 1=mineral; drop 2/3/4=peat)\n")
