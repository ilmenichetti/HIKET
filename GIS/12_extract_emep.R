# =============================================================================
# 12_extract_emep.R — atmospheric N deposition to FOREST at each plot
# Source: EMEP MSC-W rv5.3 runs with CAMS-REG emissions (Simpson, met.no;
# Zenodo 10.5281/zenodo.12580842, CC BY 4.0), 0.1 x 0.1 deg, "tot" scenario,
# annual files. ⚠ The open archive covers 2010-2022 only; the 1990-2019 trend
# runs exist but are released on request (EMEP MSC-W).
# Forest total N = WDEP_OXN + DDEP_OXN_m2Forest + WDEP_RDN + DDEP_RDN_m2Forest
# (dry deposition per m2 of forest; wet deposition is grid-uniform), mgN/m2 ->
# kg N/ha/yr (/100). Value of the 0.1 deg cell containing the plot.
# Per plot: ndep_mean_2010_22, ndep_oxn_mean, ndep_rdn_mean (kg N/ha/yr),
#   ndep_slope_per_decade (OLS over the 13 years), ndep_2010, ndep_2022.
# Writes Data/GIS_points/n_deposition.csv.
# Run from repo root (after GIS/04_download_emep.sh):  Rscript GIS/12_extract_emep.R
# =============================================================================
suppressPackageStartupMessages({library(sf); library(terra)})
G   <- Sys.getenv("HIKET_GIS_DIR", "/Volumes/NextGenC_SS/HIKET_GIS")
dir <- file.path(G, "raw/emep_ndep/EMEP_Files")
if (!dir.exists(dir)) system2("unzip", c("-o", "-q", file.path(G, "raw/emep_ndep/EMEP_Files.zip"),
                                         shQuote("EMEP_Files/*_tot_year_*"), "-d", file.path(G, "raw/emep_ndep")))
pts <- st_read(file.path(G, "plots/hiket_plots.gpkg"), quiet = TRUE)
ll  <- vect(st_transform(pts, 4326))

YRS <- 2010:2022
get <- function(sp, y) {
  f <- file.path(dir, sprintf("emepv5p3deps_%s_tot_year_%d.nc", sp, y))
  w <- rast(f, subds = paste0("WDEP_", sp)); d <- rast(f, subds = sprintf("DDEP_%s_m2Forest", sp))
  crs(w) <- crs(d) <- "EPSG:4326"
  (extract(w, ll)[, 2] + extract(d, ll)[, 2]) / 100        # kg N/ha/yr
}
oxn <- sapply(YRS, get, sp = "OXN"); rdn <- sapply(YRS, get, sp = "RDN")
tot <- oxn + rdn
slope <- apply(tot, 1, function(v) if (all(is.finite(v))) 10 * coef(lm(v ~ YRS))[2] else NA)
out <- data.frame(plot_id = pts$plot_id,
                  ndep_mean_2010_22 = rowMeans(tot), ndep_oxn_mean = rowMeans(oxn),
                  ndep_rdn_mean = rowMeans(rdn), ndep_slope_per_decade = slope,
                  ndep_2010 = tot[, 1], ndep_2022 = tot[, ncol(tot)])
write.csv(out, "Data/GIS_points/n_deposition.csv", row.names = FALSE)
cat("Wrote Data/GIS_points/n_deposition.csv\n"); print(summary(out[, -1]))
cat(sprintf("Spearman: mean vs latitude %.2f\n",
            cor(out$ndep_mean_2010_22, st_coordinates(pts)[, 2], method = "s", use = "pair")))
