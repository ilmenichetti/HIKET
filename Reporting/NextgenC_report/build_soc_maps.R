# === NextGenC reporting: kriged SOC maps of Finland ===
# Interpolate per-plot posterior-predictive SOC (mean + SD) to a regular
# 2 km grid over Finland, one layer per year (1985-2024), for each of the 6
# models, plus a six-model ensemble. Outputs multiband GeoTIFFs.
#
# Products (14 files) -> Reporting/NextgenC_report/SOC_maps/:
#   <MODEL>_SOC_mean.tif / _sd.tif   (40 bands each, 6 models)
#   ENSEMBLE_SOC_mean.tif / _sd.tif  (40 bands each)
#
# Method: ordinary kriging in LOG space (SOC is positive, right-skewed,
# multiplicative error model). One pooled variogram per model (fit on the
# per-plot temporal-mean log field), reused for all years. Total per-model SD
# fuses kriging variance + interpolated model uncertainty (option C). Ensemble
# SD fuses within-model + between-model (structural) spread (law of total var).
# See SOC_maps_README.pdf for the full rationale.

suppressPackageStartupMessages({
  library(terra)
  library(gstat)
  library(sf)
  library(sp)
})

set.seed(2025)

# --- config ---
RES_M <- 2000L                 # grid resolution (m); storage ~ 1/res^2
NMAX  <- 50L                   # local kriging neighbourhood (points)
EPSG  <- 3067L                 # ETRS-TM35FIN

# --- peatland exclusion (Luke MS-NFI 'paatyyppi'; see Data/Peat_extension/) ---
PEAT_CUTOFF   <- 0.50          # mask a 2 km cell if peat fraction exceeds this
PEAT_SOILCODE <- c("Tu", "Ve") # site_raw soil_codes dropped from kriging input
                               #   Tu = peat (Turvekerrostuma), Ve = water (Vesistö)

# --- water exclusion (SYKE Ranta10 lakes; see Data/Water_extension/) ---
# Lakes masked to NA (no forest soil under water); the sea is already outside the
# land polygon. Peat renders black, water white (see make_thumbnails/make_delta).
WATER_CUTOFF  <- 0.50          # mask a 2 km cell if lake fraction exceeds this

PROJ   <- "/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling"
RUNS   <- file.path(PROJ, "Calibration_real_data_transient", "runs")
INDIR  <- file.path(PROJ, "Reporting", "NextgenC_report")
OUTDIR <- file.path(INDIR, "SOC_maps")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# --- model -> predictive bundle (keep in sync with build_soc_matrices.R) ---
# NOTE: TP3 is the pre-fix (Euler) bundle; swap when the exact-integrator
# re-calibration is synced from Puhti, then re-run.
bundles <- c(
  SP1     = "SP1_posterior_predictive_20260608_015026.rds",
  TP2     = "TP2_posterior_predictive_20260608_020212.rds",
  TP3     = "TP3_posterior_predictive_20260608_020212.rds",
  Yasso07 = "Yasso07_posterior_predictive_20260611_032825.rds",
  Yasso15 = "Yasso15_posterior_predictive_20260611_032825.rds",
  Yasso20 = "Yasso20_posterior_predictive_20260611_033140.rds"
)
MODELS <- names(bundles)

# --- plot coordinates (EPSG:3067) + peat/water plots to exclude from input ---
site <- read.csv(file.path(PROJ, "Data", "model_inputs", "site_raw.csv"),
                 stringsAsFactors = FALSE)
coords <- site[, c("plot_id", "x_ETRS", "y_ETRS")]
coords <- coords[stats::complete.cases(coords), ]
excl_plots <- site$plot_id[site$soil_code %in% PEAT_SOILCODE]
cat(sprintf("Excluding %d peat/water plots from kriging input (soil_code %s)\n",
            length(excl_plots), paste(PEAT_SOILCODE, collapse = "/")))

# --- Finland land mask + 2 km prediction grid ---
fin     <- readRDS(file.path(PROJ, "Data", "model_inputs", "maakunta_sf.rds"))
fin_v   <- vect(sf::st_union(fin))                       # dissolved boundary
templ   <- rast(fin_v, resolution = RES_M)               # bbox grid, 2 km
crs(templ) <- paste0("EPSG:", EPSG)
mask_r  <- rasterize(fin_v, templ)                       # 1 = land, NA = sea
land    <- which(!is.na(values(mask_r)))                 # land cell indices
gxy     <- xyFromCell(mask_r, land)
pred_grid <- data.frame(x = gxy[, 1], y = gxy[, 2])
coordinates(pred_grid) <- ~ x + y
blank   <- mask_r; values(blank) <- NA_real_             # empty land-shaped layer

# vector over land cells -> one raster layer
to_layer <- function(vals) { r <- blank; r[land] <- vals; r }

# --- peat fraction on the 2 km grid (Luke MS-NFI 16 m -> aggregate; cached) ---
# Mask cells where peat (classes 2/3/4) exceeds PEAT_CUTOFF. The 16 m -> 2 km
# aggregation is heavy, so cache it; resample-on-load guarantees grid alignment.
PEAT_RAST  <- file.path(PROJ, "Data", "Peat_extension", "paatyyppi_vmi1x_1519.tif")
PEAT_CACHE <- file.path(PROJ, "Data", "Peat_extension",
                        sprintf("peat_fraction_%dm.tif", RES_M))
if (file.exists(PEAT_CACHE)) {
  peat_frac <- resample(rast(PEAT_CACHE), mask_r, method = "bilinear")
} else if (file.exists(PEAT_RAST)) {
  cat("Building peat fraction from 16 m MS-NFI (one-off, ~minutes)...\n")
  pr   <- rast(PEAT_RAST)
  fct  <- round(RES_M / res(pr)[1])
  # peat = classes 2/3/4 (mires) -> 1; everything else (mineral, non-forest) -> 0;
  # water is NaN. Fraction = peat pixels / TOTAL block pixels (fct^2), so the
  # denominator includes water -> lake-dominated cells do NOT read as ~100% peat.
  peat1   <- classify(pr, matrix(c(2,1, 3,1, 4,1), ncol = 2, byrow = TRUE), others = 0)
  pcoarse <- aggregate(peat1, fact = fct, fun = "sum", na.rm = TRUE) / (fct * fct)
  writeRaster(pcoarse, PEAT_CACHE, overwrite = TRUE,
              gdal = c("COMPRESS=DEFLATE"))
  peat_frac <- resample(pcoarse, mask_r, method = "bilinear")
} else {
  stop("Peat raster missing: run Data/Peat_extension/fetch_peatland_mask.R first")
}
peat_drop <- peat_frac > PEAT_CUTOFF                      # TRUE = mask this cell

# --- water fraction on the 2 km grid (SYKE Ranta10 lakes; cached) ---
WATER_CACHE <- file.path(PROJ, "Data", "Water_extension",
                         sprintf("water_fraction_%dm.tif", RES_M))
if (file.exists(WATER_CACHE)) {
  water_frac <- resample(rast(WATER_CACHE), mask_r, method = "bilinear")
} else if (file.exists(file.path(PROJ, "Data", "Water_extension", "Jarvi10.shp"))) {
  cat("Building water fraction from Ranta10 lakes (one-off, ~minute)...\n")
  lk   <- vect(file.path(PROJ, "Data", "Water_extension", "Jarvi10.shp"))
  fine <- rast(ext(mask_r), resolution = 200, crs = paste0("EPSG:", EPSG))
  wf   <- aggregate(rasterize(lk, fine, field = 1, background = 0),
                    fact = round(RES_M / 200), fun = "mean")
  writeRaster(wf, WATER_CACHE, overwrite = TRUE, gdal = c("COMPRESS=DEFLATE"))
  water_frac <- resample(wf, mask_r, method = "bilinear")
} else {
  stop("Water raster missing: run Data/Water_extension/fetch_water_mask.R first")
}
water_drop <- water_frac > WATER_CUTOFF                   # lakes -> mask to NA (white)
drop_mask  <- peat_drop | water_drop                      # union: both masked to NA

cat(sprintf("Grid: %d x %d @ %d m | %d land cells | %d models | peat %.1f%% + water %.1f%% masked\n",
            ncol(mask_r), nrow(mask_r), RES_M, length(land), length(MODELS),
            100 * global(peat_drop,  "mean", na.rm = TRUE)[1, 1],
            100 * global(water_drop, "mean", na.rm = TRUE)[1, 1]))

# --- accumulators for the ensemble (per-model linear-space stacks) ---
mean_stacks <- vector("list", length(MODELS)); names(mean_stacks) <- MODELS
sd_stacks   <- vector("list", length(MODELS)); names(sd_stacks)   <- MODELS
years_ref   <- NULL

for (m in MODELS) {
  pp <- readRDS(file.path(RUNS, bundles[[m]]))
  s  <- pp$posterior_summary[, c("plot_id", "year", "soc_mean", "soc_sd")]
  s  <- merge(s, coords, by = "plot_id")
  s  <- s[is.finite(s$soc_mean) & s$soc_mean > 0, ]
  s  <- s[!(s$plot_id %in% excl_plots), ]                 # drop peat/water plots

  # log-space fields: mean and a multiplicative (log) model SD
  s$zlog <- log(s$soc_mean)
  cv     <- pmax(s$soc_sd, 0) / s$soc_mean
  s$lsd  <- sqrt(log1p(cv^2))                            # log-scale model SD

  years <- sort(unique(s$year))
  if (is.null(years_ref)) years_ref <- years

  # --- pooled variogram: per-plot temporal mean of the log field ---
  # NB: temporal-mean log-SOC is ~88% nugget at the ~27.5 km plot spacing (weak
  # spatial autocorrelation). Inits are derived from the empirical variogram
  # (nugget ~ shortest-lag gamma, structured sill ~ field var - nugget) so the
  # fit reflects that reality rather than over-imposing structure. Auto-fit
  # warns on this near-pure-nugget field; we suppress and sanity-clamp the range.
  pooled <- aggregate(zlog ~ plot_id + x_ETRS + y_ETRS, data = s, FUN = mean)
  coordinates(pooled) <- ~ x_ETRS + y_ETRS
  v0    <- stats::var(pooled$zlog)
  vemp  <- variogram(zlog ~ 1, pooled, cutoff = 4e5)    # 400 km
  nug0  <- min(vemp$gamma)
  psill0 <- max(v0 - nug0, 0.05 * v0)
  vfit  <- suppressWarnings(fit.variogram(vemp,
            vgm(psill = psill0, model = "Exp", range = 8e4, nugget = nug0),
            fit.method = 7))
  # sanity guard: singular or implausible range -> honest mostly-nugget model
  if (isTRUE(attr(vfit, "singular")) ||
      vfit$range[2] <= 0 || vfit$range[2] > 4e5) {
    vfit <- vgm(psill = psill0, model = "Exp", range = 1e5, nugget = nug0)
  }

  mean_lyrs <- vector("list", length(years))
  sd_lyrs   <- vector("list", length(years))

  for (i in seq_along(years)) {
    dy <- s[s$year == years[i], ]
    sp_d <- dy; coordinates(sp_d) <- ~ x_ETRS + y_ETRS

    # krige log-mean: prediction + kriging variance (log space)
    km <- krige(zlog ~ 1, sp_d, pred_grid, model = vfit, nmax = NMAX,
                debug.level = 0)
    pred_log <- km$var1.pred
    v_krige  <- pmax(km$var1.var, 0)                     # spatial uncertainty

    # krige model log-SD field (interpolated model uncertainty)
    ks <- krige(lsd ~ 1, sp_d, pred_grid, model = vfit, nmax = NMAX,
                debug.level = 0)
    lsd_grid <- pmax(ks$var1.pred, 0)

    # option C: total log-variance = spatial + model uncertainty
    tot_logvar <- v_krige + lsd_grid^2

    # back-transform to tC/ha (median surface; lognormal SD)
    mean_v <- exp(pred_log)
    sd_v   <- mean_v * sqrt(pmax(expm1(tot_logvar), 0))

    mean_lyrs[[i]] <- to_layer(mean_v)
    sd_lyrs[[i]]   <- to_layer(sd_v)
  }

  mean_stk <- rast(mean_lyrs); sd_stk <- rast(sd_lyrs)
  names(mean_stk) <- paste0("y", years); names(sd_stk) <- paste0("y", years)

  # mask peat- and water-dominated cells to NA in both stacks
  mean_stk <- mask(mean_stk, drop_mask, maskvalues = TRUE)
  sd_stk   <- mask(sd_stk,   drop_mask, maskvalues = TRUE)

  writeRaster(mean_stk, file.path(OUTDIR, sprintf("%s_SOC_mean.tif", m)),
              overwrite = TRUE, gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3"))
  writeRaster(sd_stk,   file.path(OUTDIR, sprintf("%s_SOC_sd.tif", m)),
              overwrite = TRUE, gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3"))

  mean_stacks[[m]] <- mean_stk
  sd_stacks[[m]]   <- sd_stk
  cat(sprintf("  %-8s done (%d years)\n", m, length(years)))
}

# --- six-model ensemble (combine on the grid, equal weight) ---
M <- length(MODELS)
m_sum  <- Reduce(`+`, mean_stacks)
m_ens  <- m_sum / M                                      # ensemble mean
within <- Reduce(`+`, lapply(sd_stacks, function(x) x^2)) / M
between <- Reduce(`+`, lapply(mean_stacks,
                  function(x) (x - m_ens)^2)) / (M - 1)  # structural spread
s_ens  <- sqrt(within + between)                         # law of total variance

names(m_ens) <- paste0("y", years_ref); names(s_ens) <- paste0("y", years_ref)
writeRaster(m_ens, file.path(OUTDIR, "ENSEMBLE_SOC_mean.tif"),
            overwrite = TRUE, gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3"))
writeRaster(s_ens, file.path(OUTDIR, "ENSEMBLE_SOC_sd.tif"),
            overwrite = TRUE, gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3"))

cat("\nWrote 14 GeoTIFFs to:\n", OUTDIR, "\n", sep = "")
