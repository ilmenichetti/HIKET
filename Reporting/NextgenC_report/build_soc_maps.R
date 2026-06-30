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

# --- plot coordinates (EPSG:3067) ---
site <- read.csv(file.path(PROJ, "Data", "model_inputs", "site_raw.csv"),
                 stringsAsFactors = FALSE)
coords <- site[, c("plot_id", "x_ETRS", "y_ETRS")]
coords <- coords[stats::complete.cases(coords), ]

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

cat(sprintf("Grid: %d x %d @ %d m | %d land cells | %d models\n",
            ncol(mask_r), nrow(mask_r), RES_M, length(land), length(MODELS)))

# --- accumulators for the ensemble (per-model linear-space stacks) ---
mean_stacks <- vector("list", length(MODELS)); names(mean_stacks) <- MODELS
sd_stacks   <- vector("list", length(MODELS)); names(sd_stacks)   <- MODELS
years_ref   <- NULL

for (m in MODELS) {
  pp <- readRDS(file.path(RUNS, bundles[[m]]))
  s  <- pp$posterior_summary[, c("plot_id", "year", "soc_mean", "soc_sd")]
  s  <- merge(s, coords, by = "plot_id")
  s  <- s[is.finite(s$soc_mean) & s$soc_mean > 0, ]

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
