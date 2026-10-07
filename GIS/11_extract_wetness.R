# =============================================================================
# 11_extract_wetness.R — topographic wetness at each plot
#   TWI 16 m (Luke 2016)       twi16_point, twi16_mean3x3 (48 m window); the
#                              file stores TWI x 1000 as integer, rescaled here
#   DTW 2 m (Luke 2023), for each stream-initiation threshold v in
#     050/1/2/4/10 ha (very moist ... dry; 2 ha = average conditions):
#     dtw<v>_point     value at the plot centre (cm)
#     dtw<v>_mean10    mean within 10 m (GPS error: median 3.2 m, 99th pct 9 m)
#     dtw<v>_min10     minimum within 10 m
#     dtw<v>_wet10     share of the 10 m disc with DTW < 100 cm ("wet", Luke)
# Reads the GIS store ($HIKET_GIS_DIR); writes Data/GIS_points/wetness.csv.
# DTW tiles are opened only where plots fall (tile index from 02_).
# Run from repo root:  Rscript GIS/11_extract_wetness.R
# =============================================================================
suppressPackageStartupMessages({library(sf); library(terra)})
G   <- Sys.getenv("HIKET_GIS_DIR", "/Volumes/NextGenC_SS/HIKET_GIS")
pts <- st_read(file.path(G, "plots/hiket_plots.gpkg"), quiet = TRUE)
out <- data.frame(plot_id = pts$plot_id)

# --- TWI 16 m -----------------------------------------------------------------
twi <- rast(file.path(G, "raw/twi_16m/TWI_16m_Finland_NA_lakes_int.tif")) / 1000   # stored x1000 as integer
v   <- vect(st_transform(pts, crs(twi)))
out$twi16_point   <- extract(twi, v)[, 2]
out$twi16_mean3x3 <- extract(focal(crop(twi, ext(v) + 100), w = 3, fun = mean, na.rm = TRUE), v)[, 2]

# --- DTW 2 m, tile by tile -----------------------------------------------------
disc <- st_buffer(pts, 10)
for (th in c("050", "1", "2", "4", "10")) {
  dir <- file.path(G, sprintf("raw/dtw_2m/DTW_INT_CMv2_%s", th))
  idx <- st_read(file.path(dir, sprintf("DTW_INT_CMv2_%s.shp", th)), quiet = TRUE)
  hit <- st_intersects(pts, idx)                     # tile containing each centre
  res <- matrix(NA_real_, nrow(pts), 4)
  for (tile in unique(idx$location[unlist(hit)])) {
    f <- file.path(dir, tile); if (!file.exists(f)) next
    r <- rast(f)
    k <- which(vapply(hit, function(h) tile %in% idx$location[h], logical(1)))
    res[k, 1] <- extract(r, vect(pts[k, ]))[, 2]
    e <- extract(r, vect(disc[k, ]))                 # all pixels in each disc
    s <- split(e[, 2], e$ID)
    res[k, 2] <- vapply(s, function(x) mean(x, na.rm = TRUE), numeric(1))
    res[k, 3] <- vapply(s, function(x) suppressWarnings(min(x, na.rm = TRUE)), numeric(1))
    res[k, 4] <- vapply(s, function(x) mean(x < 100, na.rm = TRUE), numeric(1))
  }
  res[!is.finite(res)] <- NA
  out[paste0("dtw", th, c("_point", "_mean10", "_min10", "_wet10"))] <- res
  cat(sprintf("DTW %s ha: %d plots with values\n", th, sum(is.finite(res[, 1]))))
}
write.csv(out, "Data/GIS_points/wetness.csv", row.names = FALSE)
cat("Wrote Data/GIS_points/wetness.csv\n"); print(summary(out[, c("twi16_point", "dtw2_point", "dtw2_wet10")]))
