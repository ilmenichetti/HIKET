# =============================================================================
# 02_plots_and_dtw_tiles.R — plot point layer + the DTW tiles that cover it
#
# Plot universe = Data/model_inputs/site_attributes.csv (2719 plots, the
# SOC-INDEPENDENT set; see CLAUDE.md "The data layer must not depend on the
# calibration layer"). Coordinates x_ETRS/y_ETRS (EPSG:3067) are the 2021-23
# GPS plot centres (median reported accuracy 3.2 m, 99th pct 9 m), so 2 m
# rasters are usable with a small buffer.
#
# Writes to the GIS store ($HIKET_GIS_DIR, default /Volumes/NextGenC_SS/HIKET_GIS):
#   plots/hiket_plots.gpkg            plot points, EPSG:3067
#   plots/dtw_tiles_needed.txt        DTW tile names touching a plot + 20 m buffer
# Run from repo root:  Rscript GIS/02_plots_and_dtw_tiles.R
# =============================================================================
suppressPackageStartupMessages(library(sf))
G  <- Sys.getenv("HIKET_GIS_DIR", "/Volumes/NextGenC_SS/HIKET_GIS")
sa <- read.csv("Data/model_inputs/site_attributes.csv")
pts <- st_as_sf(sa[, c("plot_id", "x_ETRS", "y_ETRS")], coords = c("x_ETRS", "y_ETRS"),
                crs = 3067, remove = FALSE)
st_write(pts, file.path(G, "plots/hiket_plots.gpkg"), delete_dsn = TRUE, quiet = TRUE)

idx   <- st_read(file.path(G, "raw/dtw_2m/DTW_INT_CMv2_1/DTW_INT_CMv2_1.shp"), quiet = TRUE)
hit   <- st_intersects(st_buffer(pts, 20), idx)
tiles <- sort(unique(idx$location[unlist(hit)]))
writeLines(tiles, file.path(G, "plots/dtw_tiles_needed.txt"))
cat(sprintf("%d plots -> %d DTW tiles (of %d); plots without a tile: %d\n",
            nrow(pts), length(tiles), nrow(idx), sum(lengths(hit) == 0)))
