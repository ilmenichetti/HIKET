# =============================================================================
# 14_extract_kemera.R — COMPLETED subsidised forestry works at each plot
# Source: Metsäkeskus KEMERA / Metka open data (CC BY 4.0), regional
# GeoPackages; only the completiondeclaration_stand_* layers (work reported as
# DONE, with its real end date), point-in-polygon at the plot.
# Work codes are grouped by their name in the official code list
# (docs/kemera-koodisto-ja-tietokantakuvaus.xlsx, field workcode):
#   tending      Taimikon / Nuoren metsän hoito, varhaishoito, Pystykarsinta
#   energywood   Energiapuun korjuu, Haketustuki  (whole trees / residues removed)
#   fertilise    Terveyslannoitus                 (remedial fertilisation)
#   ditching     Kunnostusojitus, Suometsänhoito
#   regenerate   Metsänuudistaminen, Luontainen uudistaminen, Kulotus (burning)
#   afforest     *metsitys (incl. Pellonmetsitys = afforested FIELD)
#   moose        Hirvivahinko*                    (moose-damage compensation)
#   other        everything else (nature management, roots rot control, ...)
# Per plot: kem_n_<group>, kem_last_<group> (year), kem_first_year, kem_n.
# Writes Data/GIS_points/kemera_works.csv (per-plot summary) and
# Data/GIS_points/events/kemera_events.csv (one row per completed work).
# Run from repo root:  Rscript GIS/14_extract_kemera.R
# =============================================================================
suppressPackageStartupMessages({library(sf); library(dplyr); library(readxl); library(tidyr)})
dir.create("Data/GIS_points/events", showWarnings = FALSE)
G   <- Sys.getenv("HIKET_GIS_DIR", "/Volumes/NextGenC_SS/HIKET_GIS")
dir <- file.path(G, "raw/metsakeskus_kemera")
pts <- st_read(file.path(G, "plots/hiket_plots.gpkg"), quiet = TRUE)

cl <- as.data.frame(suppressMessages(read_excel(
  file.path(G, "docs/kemera-koodisto-ja-tietokantakuvaus.xlsx"), "Koodisto", skip = 1)))[, 2:4]
names(cl) <- c("field", "code", "name")
cl <- cl[tolower(cl$field) == "workcode", ]
grp <- function(n) dplyr::case_when(
  grepl("taimikon|nuoren mets|varhaishoito|pystykarsinta", n, ignore.case = TRUE) ~ "tending",
  grepl("energiapuu|haketus", n, ignore.case = TRUE)                              ~ "energywood",
  grepl("lannoitus", n, ignore.case = TRUE)                                       ~ "fertilise",
  grepl("ojitus|suometsänhoito", n, ignore.case = TRUE)                           ~ "ditching",
  grepl("uudistaminen|kulotus", n, ignore.case = TRUE)                            ~ "regenerate",
  grepl("metsitys", n, ignore.case = TRUE)                                        ~ "afforest",
  grepl("hirvi", n, ignore.case = TRUE)                                           ~ "moose",
  TRUE                                                                            ~ "other")
cl$group <- grp(cl$name)

unpack <- function(zip) {
  out <- sub("\\.zip$", ".gpkg", zip)
  if (!file.exists(out)) system2("python3", c("-c", shQuote(sprintf(
    "import zipfile,shutil; z=zipfile.ZipFile(%s); m=[i for i in z.infolist() if i.filename.endswith('.gpkg')][0]; shutil.copyfileobj(z.open(m), open(%s,'wb'))",
    deparse(zip), deparse(out)))))
  out
}

hits <- list()
for (z in list.files(dir, "\\.zip$", full.names = TRUE)) {
  g  <- unpack(z)
  ly <- grep("^completiondeclaration_stand_", st_layers(g)$name, value = TRUE)
  for (l in ly) {
    rt <- st_read(g, quiet = TRUE, query = sprintf(
      "SELECT MIN(minx) xmin, MIN(miny) ymin, MAX(maxx) xmax, MAX(maxy) ymax FROM rtree_%s_geometry", l))
    if (!is.finite(rt$xmin)) next
    xy <- st_coordinates(pts)
    inr <- pts[xy[, 1] >= rt$xmin & xy[, 1] <= rt$xmax & xy[, 2] >= rt$ymin & xy[, 2] <= rt$ymax, ]
    if (!nrow(inr)) next
    d <- st_read(g, layer = l, quiet = TRUE, wkt_filter = st_as_text(st_union(st_buffer(inr, 1))))
    if (!nrow(d)) next
    # first non-empty date: actual end > actual start > completion report > project end
    cand <- intersect(c("realenddate", "realstartdate", "arrivaldate",
                        "completiondeclarationarrivaldate", "projectenddate"), names(d))
    dd <- as.data.frame(lapply(st_drop_geometry(d)[cand], as.character))
    d$date <- apply(dd, 1, function(x) x[!is.na(x) & nzchar(x)][1])
    d <- d[, c("workcode", "date")]
    j <- st_join(inr[, "plot_id"], d, join = st_within, left = FALSE)
    if (nrow(j)) hits[[paste(basename(z), l)]] <- st_drop_geometry(j) |> mutate(layer = l)
  }
  cat(sprintf("%-32s done\n", basename(z)))
}
H <- bind_rows(hits) |> distinct() |>
  mutate(yr = as.integer(substr(as.character(date), 1, 4)),
         group = cl$group[match(workcode, as.numeric(cl$code))],
         group = ifelse(is.na(group), "other", group))
write.csv(transmute(H, plot_id, date = substr(as.character(date), 1, 10), year = yr, workcode, group),
          "Data/GIS_points/events/kemera_events.csv", row.names = FALSE)
cat("\nCompleted works hitting plots, by group:\n"); print(table(H$group))
cat("Years:", range(H$yr, na.rm = TRUE), " records without a date:", sum(is.na(H$yr)), "\n")
W <- H |> group_by(plot_id, group) |>
  summarise(n = n(), last = suppressWarnings(max(yr, na.rm = TRUE)), .groups = "drop") |>
  mutate(last = ifelse(is.finite(last), last, NA)) |>
  pivot_wider(names_from = group, values_from = c(n, last), names_glue = "kem_{.value}_{group}")
A <- H |> group_by(plot_id) |> summarise(kem_n = n(), kem_first_year = suppressWarnings(min(yr, na.rm = TRUE)),
                                          .groups = "drop") |>
  mutate(kem_first_year = ifelse(is.finite(kem_first_year), kem_first_year, NA))
out <- data.frame(plot_id = pts$plot_id) |> left_join(A, by = "plot_id") |> left_join(W, by = "plot_id")
nn <- grep("^kem_n", names(out), value = TRUE)
out[nn] <- lapply(out[nn], function(x) replace(x, is.na(x), 0L))
write.csv(out, "Data/GIS_points/kemera_works.csv", row.names = FALSE)
cat("Wrote Data/GIS_points/kemera_works.csv:", sum(out$kem_n > 0), "plots with completed works\n")
