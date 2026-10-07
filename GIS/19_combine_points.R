# =============================================================================
# 19_combine_points.R — merge all per-plot GIS extractions into one table
# Reads every Data/GIS_points/*.csv except the combined file itself, joins on
# plot_id, writes Data/GIS_points/gis_point_predictors.csv (one row per plot).
# ⚠ site_attributes.csv repeats 6 plot_ids at TWO locations each (83471-83473,
# 83511, 83552, 83553 — the deferred ambiguous BIOSOIL->MUSTIKKA mappings, none
# in the calibration set). Their location is unknown, so they are DROPPED here.
# Run from repo root:  Rscript GIS/19_combine_points.R
# =============================================================================
sa  <- read.csv("Data/model_inputs/site_attributes.csv")
loc <- unique(sa[, c("plot_id", "x_ETRS", "y_ETRS")])
amb <- unique(loc$plot_id[duplicated(loc$plot_id)])
cat("Dropping", length(amb), "plot_ids with more than one location:", amb, "\n")
fs  <- setdiff(list.files("Data/GIS_points", "\\.csv$", full.names = TRUE),   # events/ not included
               "Data/GIS_points/gis_point_predictors.csv")
out <- data.frame(plot_id = setdiff(unique(sa$plot_id), amb))
for (f in fs) {
  d <- read.csv(f); d <- d[!d$plot_id %in% amb & !duplicated(d$plot_id), ]
  out <- merge(out, d, by = "plot_id", all.x = TRUE)
  cat(sprintf("%-28s %3d columns\n", basename(f), ncol(d) - 1))
}
stopifnot(!anyDuplicated(out$plot_id))
write.csv(out, "Data/GIS_points/gis_point_predictors.csv", row.names = FALSE)
cat("Wrote Data/GIS_points/gis_point_predictors.csv:", nrow(out), "plots x", ncol(out) - 1, "\n")
