# === NextGenC: ensemble SOC change (delta) map, last 10 years ===
# Per-cell annual SOC change (tC/ha/yr) over 2015-2024 from the ensemble mean,
# as the OLS slope across the 10 yearly layers. Writes a GeoTIFF + a thumbnail.
# Diverging palette centred on zero (loss = red, gain = blue) since the field
# has both gains (~76%) and losses (~24%). -> SOC_maps/ + thumbnails/

suppressPackageStartupMessages(library(terra))

PROJ   <- "/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling"
MAPDIR <- file.path(PROJ, "Reporting", "NextgenC_report", "SOC_maps")
THUMB  <- file.path(MAPDIR, "thumbnails")
NYEARS <- 10L

fin <- tryCatch(vect(sf::st_union(readRDS(
  file.path(PROJ, "Data", "model_inputs", "maakunta_sf.rds")))),
  error = function(e) NULL)
peat_cache <- file.path(PROJ, "Data", "Peat_extension", "peat_fraction_2000m.tif")
peat_r <- if (file.exists(peat_cache)) rast(peat_cache) else NULL

r      <- rast(file.path(MAPDIR, "ENSEMBLE_SOC_mean.tif"))
last   <- r[[(nlyr(r) - NYEARS + 1L):nlyr(r)]]            # final NYEARS layers
yrs    <- as.integer(sub("^y", "", names(last)))
tbar   <- mean(yrs); w <- (yrs - tbar) / sum((yrs - tbar)^2)
delta  <- sum(last * w)                                   # OLS slope, tC/ha/yr
names(delta) <- sprintf("delta_%d_%d", min(yrs), max(yrs))

writeRaster(delta, file.path(MAPDIR, "ENSEMBLE_SOC_delta10yr.tif"),
            overwrite = TRUE, gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3"))

# --- thumbnail: diverging palette, symmetric limits centred on 0 ---
v   <- values(delta); v <- v[is.finite(v)]
lim <- as.numeric(quantile(abs(v), 0.98))                 # robust symmetric cap
# RdYlBu: loss = red, gain = blue, zero = PALE YELLOW (not white) so near-zero
# change stays visibly "data" instead of blending into the white NoData bg.
pal <- colorRampPalette(c("#a50026","#d73027","#f46d43","#fdae61","#fee090",
                          "#ffffbf","#e0f3f8","#abd9e9","#74add1","#4575b4",
                          "#313695"))(100)
ttl <- sprintf("ENSEMBLE — SOC change, %d–%d", min(yrs), max(yrs))

# clamp (squish) extremes into the colour range so out-of-range cells render as
# the end colour, not transparent white (terra plots out-of-range as NoData).
# The saved GeoTIFF keeps the true unclamped values; only the display is clamped.
ddisp <- clamp(delta, -lim, lim, values = TRUE)

H <- 2000L; asp <- nrow(delta) / ncol(delta); W <- round(H / asp + 0.45 * H)
png(file.path(THUMB, "ENSEMBLE_SOC_delta10yr.png"), width = W, height = H, res = 300)
par(mar = c(2, 2, 3, 4))
# range padded 0.2% beyond the clamp so the clamped extremes (= ±lim) fall
# strictly inside it; terra renders cells exactly on the range edge as white.
plot(ddisp, col = pal, range = c(-lim, lim) * 1.002, axes = FALSE, main = "",
     mar = c(2, 2, 3, 4), plg = list(title = "tC/ha/yr"))
mtext(ttl, side = 3, line = 0.6, adj = 0, font = 2, cex = 1.05)
if (!is.null(peat_r)) {                                    # peat black; water white
  # peat cache + fin are in 3067; project onto the delta raster (now EPSG:3035).
  pk <- project(peat_r, delta, method = "bilinear") > 0.5
  plot(ifel(pk, 1, NA), col = "black", legend = FALSE, axes = FALSE, add = TRUE)
}
if (!is.null(fin)) lines(project(fin, crs(delta)), col = "grey30", lwd = 0.4)
dev.off()

cat(sprintf("delta %d-%d: mean %.3f, median %.3f tC/ha/yr; %.0f%% gaining\n",
            min(yrs), max(yrs), mean(v), median(v), 100 * mean(v > 0)))
cat("wrote ENSEMBLE_SOC_delta10yr.tif + thumbnail (palette cap +/-", round(lim, 2),
    "tC/ha/yr)\n")
