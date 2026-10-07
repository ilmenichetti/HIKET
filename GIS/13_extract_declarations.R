# =============================================================================
# 13_extract_declarations.R — harvest history 1997-2026 from Metsäkeskus forest
# use declarations (metsänkäyttöilmoitukset), point-in-polygon at each plot.
#
# A declaration is a harvest INTENTION valid 3 years (completion year is empty
# in the open data); the realised area differs by ~0-10%. Small fellings (mean
# diameter <= 13 cm, household use) need no declaration. Dense from 2004.
# Cutting codes (cuttingrealizationpractice, Metsäkeskus code list):
#   thinning      2,3,12,13            regeneration  1,4,5,6,7,8,17,19
#   clear-cut     5                    damage-driven 20-25 (storm, insects,
#                                       other) or forestdamagequalifier set
# Per plot: mki_n, mki_n_thin, mki_n_regen, mki_n_clearcut, mki_n_damage,
#   mki_first_year, mki_last_year, mki_last_regen_year, mki_yrs_since_regen
#   (to 2024; NA = none declared), mki_any_regen_since_2006.
# Each region's GeoPackage is unpacked once next to its zip (the zip member
# names carry non-ASCII characters that break unzip).
# Writes Data/GIS_points/harvest_declarations.csv (per-plot summary) and
# Data/GIS_points/events/harvest_declaration_events.csv (one row per declaration).
# Run from repo root:  Rscript GIS/13_extract_declarations.R
# =============================================================================
suppressPackageStartupMessages({library(sf); library(dplyr)})
dir.create("Data/GIS_points/events", showWarnings = FALSE)
G   <- Sys.getenv("HIKET_GIS_DIR", "/Volumes/NextGenC_SS/HIKET_GIS")
dir <- file.path(G, "raw/metsakeskus_declarations")
pts <- st_read(file.path(G, "plots/hiket_plots.gpkg"), quiet = TRUE)

unpack <- function(zip) {                     # -> path of the unpacked .gpkg
  out <- sub("\\.zip$", ".gpkg", zip)
  if (!file.exists(out)) system2("python3", c("-c", shQuote(sprintf(
    "import zipfile,shutil; z=zipfile.ZipFile(%s); m=[i for i in z.infolist() if i.filename.endswith('.gpkg')][0]; shutil.copyfileobj(z.open(m), open(%s,'wb'))",
    deparse(zip), deparse(out)))))
  out
}

THIN <- c(2, 3, 12, 13); REGEN <- c(1, 4, 5, 6, 7, 8, 17, 19); DAMAGE <- 20:25
hits <- list()
for (z in list.files(dir, "\\.zip$", full.names = TRUE)) {
  g <- unpack(z)
  # region extent from the R-tree, then only the plots inside it
  rt <- st_read(g, quiet = TRUE, query = paste(
    "SELECT MIN(minx) xmin, MIN(miny) ymin, MAX(maxx) xmax, MAX(maxy) ymax",
    "FROM rtree_forestusedeclaration_geometry"))
  inreg <- pts[st_coordinates(pts)[, 1] >= rt$xmin & st_coordinates(pts)[, 1] <= rt$xmax &
               st_coordinates(pts)[, 2] >= rt$ymin & st_coordinates(pts)[, 2] <= rt$ymax, ]
  if (!nrow(inreg)) next
  d <- st_read(g, quiet = TRUE, wkt_filter = st_as_text(st_union(st_buffer(inreg, 1))),
               query = paste("SELECT geometry, cuttingrealizationpractice AS crp,",
                             "forestdamagequalifier AS fdq, declarationarrivaldate AS dt",
                             "FROM forestusedeclaration"))
  j <- st_join(inreg[, "plot_id"], d, join = st_within, left = FALSE)
  hits[[basename(z)]] <- st_drop_geometry(j)
  cat(sprintf("%-32s plots in extent %4d, declarations hitting plots %5d\n",
              basename(z), nrow(inreg), nrow(j)))
}
H <- bind_rows(hits) |> distinct() |>
  mutate(yr = as.integer(substr(dt, 1, 4)), crp = as.integer(crp))
# Event list (one row per declaration hitting a plot) for time-resolved analyses
write.csv(transmute(H, plot_id, date = substr(dt, 1, 10), year = yr, cutting_code = crp,
                    damage_code = fdq,
                    type = case_when(crp %in% DAMAGE | !is.na(fdq) ~ "damage",
                                     crp %in% REGEN ~ "regeneration", crp %in% THIN ~ "thinning",
                                     TRUE ~ "other")),
          "Data/GIS_points/events/harvest_declaration_events.csv", row.names = FALSE)
S <- H |> group_by(plot_id) |>
  summarise(mki_n = n(), mki_n_thin = sum(crp %in% THIN), mki_n_regen = sum(crp %in% REGEN),
            mki_n_clearcut = sum(crp == 5, na.rm = TRUE),
            mki_n_damage = sum(crp %in% DAMAGE | !is.na(fdq)),
            mki_first_year = min(yr, na.rm = TRUE), mki_last_year = max(yr, na.rm = TRUE),
            mki_last_regen_year = suppressWarnings(max(yr[crp %in% REGEN], na.rm = TRUE)),
            .groups = "drop") |>
  mutate(mki_last_regen_year = ifelse(is.finite(mki_last_regen_year), mki_last_regen_year, NA),
         mki_yrs_since_regen = 2024 - mki_last_regen_year,
         mki_any_regen_since_2006 = !is.na(mki_last_regen_year) & mki_last_regen_year >= 2006)
out <- data.frame(plot_id = pts$plot_id) |> left_join(S, by = "plot_id")
cnt <- c("mki_n", "mki_n_thin", "mki_n_regen", "mki_n_clearcut", "mki_n_damage")
out[cnt] <- lapply(out[cnt], function(x) replace(x, is.na(x), 0L))
out$mki_any_regen_since_2006[is.na(out$mki_any_regen_since_2006)] <- FALSE
write.csv(out, "Data/GIS_points/harvest_declarations.csv", row.names = FALSE)
cat("Wrote Data/GIS_points/harvest_declarations.csv\n"); print(summary(out[, -1]))
