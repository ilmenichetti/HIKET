# === NextGenC: thumbnails of the SOC map layers ===
# One PNG per GeoTIFF, each the 1985-2024 TIME-AVERAGE of that band stack.
# Small, print-quality previews (~2000 px tall, 300 dpi). -> SOC_maps/thumbnails/

suppressPackageStartupMessages(library(terra))

PROJ   <- "/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling"
MAPDIR <- file.path(PROJ, "Reporting", "NextgenC_report", "SOC_maps")
OUTDIR <- file.path(MAPDIR, "thumbnails")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# Finland outline (light overlay for context)
fin <- tryCatch(vect(sf::st_union(readRDS(
  file.path(PROJ, "Data", "model_inputs", "maakunta_sf.rds")))),
  error = function(e) NULL)

# peat fraction (rendered BLACK); water/sea are NA in the data -> render white
peat_cache <- file.path(PROJ, "Data", "Peat_extension", "peat_fraction_2000m.tif")
peat_r <- if (file.exists(peat_cache)) rast(peat_cache) else NULL

tifs <- list.files(MAPDIR, pattern = "_SOC_(mean|sd)[.]tif$", full.names = TRUE)

H <- 2000L                                   # target height (px)
for (f in tifs) {
  r    <- rast(f)
  tavg <- mean(r, na.rm = TRUE)              # average across the 40 year-bands

  base  <- sub("[.]tif$", "", basename(f))
  is_sd <- grepl("_sd$", base)
  model <- sub("_SOC_(mean|sd)$", "", base)
  pal   <- if (is_sd) hcl.colors(100, "Purples", rev = TRUE)   # calm sequential
           else        hcl.colors(100, "YlGnBu",  rev = TRUE)
  lab   <- if (is_sd) "uncertainty (SD)" else "mean SOC"
  ttl   <- sprintf("%s — %s, 1985–2024 avg", model, lab)

  # width from raster aspect + room for legend/margins
  asp <- nrow(tavg) / ncol(tavg)
  W   <- round(H / asp + 0.45 * H)

  png(file.path(OUTDIR, paste0(base, ".png")),
      width = W, height = H, res = 300)
  par(mar = c(2, 2, 3, 4))
  plot(tavg, col = pal, axes = FALSE, main = "",
       mar = c(2, 2, 3, 4), plg = list(title = "tC/ha"))
  # left-aligned title over the map panel only (avoids the legend column)
  mtext(ttl, side = 3, line = 0.6, adj = 0, font = 2, cex = 1.05)
  # overlay peat-masked cells in black (water/sea stay white). peat cache + fin are
  # in 3067; project() onto the map raster (now EPSG:3035) — handles either CRS.
  if (!is.null(peat_r)) {
    pk <- project(peat_r, tavg, method = "bilinear") > 0.5
    plot(ifel(pk, 1, NA), col = "black", legend = FALSE, axes = FALSE, add = TRUE)
  }
  if (!is.null(fin)) lines(project(fin, crs(tavg)), col = "grey30", lwd = 0.4)
  dev.off()

  cat(sprintf("  %-26s %d x %d px  %.0f KB\n", paste0(base, ".png"), W, H,
              file.size(file.path(OUTDIR, paste0(base, ".png"))) / 1024))
}

cat("\nWrote", length(tifs), "thumbnails to:\n", OUTDIR, "\n")
