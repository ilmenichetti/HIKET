# =============================================================================
# 10_extract_history.R — early-20th-century human impact at each plot
# Source: Aakala, Kulha & Kuuluvainen 2023 (figshare 23257562), see
# GIS/05_download_history.sh. Per plot:
#   kaski_share, kaski_1860, kaski_1913  slash-and-burn of the plot's parish
#                                        (Heikinheimo 1915 classes; higher = more)
#   pop1925_5km/10km/20km                persons within the radius (1925 points)
#   dist_rail1925_km                     distance to the 1925 railway network
#   state_forest1925                     plot inside a 1925 state-owned forest
# Writes Data/GIS_points/history_1925.csv (plot_id key).
# Run from repo root:  Rscript GIS/10_extract_history.R
# =============================================================================
suppressPackageStartupMessages(library(sf))
G   <- Sys.getenv("HIKET_GIS_DIR", "/Volumes/NextGenC_SS/HIKET_GIS")
H   <- file.path(G, "raw/aakala2023_history")
pts <- st_read(file.path(G, "plots/hiket_plots.gpkg"), quiet = TRUE)

kb  <- st_read(file.path(H, "slash_and_burn.gpkg"), quiet = TRUE)
kb  <- st_transform(kb, 3067)
j   <- st_join(pts, kb[, c("kaskimaita", "kasket1860", "kasket1913")], join = st_within, left = TRUE)
j   <- j[!duplicated(j$plot_id), ]

pop <- st_transform(st_read(file.path(H, "population_1925.gpkg"), quiet = TRUE), 3067)
cnt <- function(r) vapply(st_is_within_distance(pts, pop, dist = r * 1000),
                          function(i) sum(pop$Asukkaita2[i]), numeric(1))
rail <- st_transform(st_read(file.path(H, "railroads_1925.gpkg"), quiet = TRUE), 3067)
sf25 <- st_transform(st_read(file.path(H, "state_forests_1925.gpkg"), quiet = TRUE), 3067)

out <- data.frame(
  plot_id          = pts$plot_id,
  kaski_share      = j$kaskimaita[match(pts$plot_id, j$plot_id)],
  kaski_1860       = j$kasket1860[match(pts$plot_id, j$plot_id)],
  kaski_1913       = j$kasket1913[match(pts$plot_id, j$plot_id)],
  pop1925_5km      = cnt(5), pop1925_10km = cnt(10), pop1925_20km = cnt(20),
  dist_rail1925_km = as.numeric(apply(st_distance(pts, rail), 1, min)) / 1000,
  state_forest1925 = lengths(st_intersects(pts, st_make_valid(sf25))) > 0)
write.csv(out, "Data/GIS_points/history_1925.csv", row.names = FALSE)
cat("Wrote Data/GIS_points/history_1925.csv:", nrow(out), "plots\n"); print(summary(out[, -1]))
