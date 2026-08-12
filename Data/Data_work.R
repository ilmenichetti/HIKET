# ========================================================================================
# HIKET Data Preparation
# Lorenzo Menichetti
#
# Produces three files in ./Data/model_inputs/:
#   input_raw_monthly.csv  -- monthly time series (FILTERED to plots with full climate)
#   climate_monthly.csv    -- standalone monthly climate
#   site_raw.csv           -- comprehensive site attributes and stratification metadata
#
# site_raw includes:
#   - GTK soil type (API, section 1.3)
#   - Biosoil mineral + OFH properties (section 1.4)
#   - Cajander fertility class KA and SYKE biogeographic zone kasvyo_syke (section 1.4)
#   - Measured mineral profile depth profile_depth_cm (section 1.4)
#   - Climate normals + variability + Köppen-Geiger class (section 4.2 / 5)
#   - Derived soil chemistry: CN_ratio, CEC, base_saturation (section 5)
#   - Litter quality: AWEN fractions, woody share, conifer share (section 4.2 / 5)
#   - Peatland flag; peatland plots (KA 11-13) excluded via calib_ready = FALSE
#
# SOC observations (soc_obs_tCha) come from the homogenized three-campaign baseline
# in ./Data/SOC_homogeneized/ (section 1, merged in section 1.3b). Full rationale in
# the METHODS & MATERIALS block at the end of this script. In brief:
#   - VMI8 1985 / Biosoil 2006 / Komeetta 2024 from ONE LUKE source, so the campaigns
#     share a basis; stocks are bulk-density AND coarse-fragment (stoniness) corrected
#     and reproduce LUKE's official national figures exactly (2006 = 59.1, 2024 = 61.0
#     Mg/ha, org + 0-40 cm, weighted). The former source lacked the stoniness
#     correction and over-counted mineral C by ~1.6x.
#   - Target = soc_profile: OFH (measured) + mineral 0-40 cm (measured) + deep tail
#     modelled as C(z) = C0*exp(-lambda*z), integrated only to z_cap = min(100 cm,
#     depth augering reached). Per-plot VARIABLE depth, not a fixed 1 m -- shallow
#     soils no longer receive carbon below bedrock.
#   - soc_obs_tCha_sum retained as the measured 0-40 cm + OFH sum for records.
#
# Organic layer exclusions (section 1.1):
#   - organic_missing: OFH layer never recorded -> calib_ready = FALSE
#   - organic_zero:    OFH recorded but C_kgha == 0 in any year -> calib_ready = FALSE
#   Both sets verified non-overlapping.
#
# SOC outlier exclusion (section 1.3b / section 5):
#   - soc_outlier: whole-profile stock > 250 Mg/ha in any campaign -> calib_ready = FALSE
#     (4 plots). Companion flag high_change is recorded but plots are RETAINED.
#
# Constant-litter exclusion (section 2 / section 5):
#   - const_litter: annual litter identical (rel. SD < 1e-6) over 1986-2024 -> calib_ready
#     = FALSE (8 plots). A fixed repeated value, not a measured series.
#
# NOTE for residual analysis covariates:
#   - MeanBulkDensity has 84% NAs -- prefer EstimatedBulkDensity for RF
# ========================================================================================

library(readxl)


# ========================================================================================
# 1. Load and match data sources
# ========================================================================================

# --- SOC source: homogenized three-campaign baseline (see Data/SOC_homogeneized/) ---
# Replaces the former ./Data/SOC/soilC1985_2006.csv + separate Komeetta 2024 append.
# Single homogeneous source for VMI8 1985 / Biosoil 2006 / Komeetta 2024, built from
# Hannu Ilvesniemi's LUKE workbook (Komeetta 150526hi--.xlsx). Key properties:
#   - LUKE pre-computed stocks: bulk-density AND coarse-fragment (stoniness) corrected.
#     The former source lacked the stoniness correction and over-counted mineral C ~1.6x
#     (verified against Juha Heikkinen's official LUKE computation, which this baseline
#     reproduces EXACTLY: org+0-40, weighted, 2006 = 59055 / 2024 = 61047 kg/ha, n = 446).
#   - Peat layers (CODE_LAYER H01/H12) already excluded.
#   - All three campaigns on one corrected basis (1985 verified homogeneous with 2006/2024).
# Layer vocabulary matches the old file exactly (organic, 0-5cm, 5-20cm, 0-10cm,
# 10-20cm, 20-40cm), so the organic-quality flags in 1.1 work unchanged.
SOC    <- read.csv("./Data/SOC_homogeneized/soc_homogenized_layers.csv")
inputs <- read.csv("./Data/LitterData/tree_litter_per_site_year_awen_29.04.26.csv", sep = ";")
inputs_by_component <- read.csv("./Data/LitterData/tree_litter_per_site_year_by_component_29.04.26.csv", sep = ";")
keys   <- read.csv("./Data/soil_litter_site_key.csv", sep = ";")

keys_clean <- keys[!is.na(keys$koealatunnus_BIOSOIL) &
                     !is.na(keys$koealatunnus_MUSTIKKA), ]

inputs$plot_id <- keys_clean$koealatunnus_BIOSOIL[
  match(inputs$site, keys_clean$koealatunnus_MUSTIKKA)]
inputs_by_component$plot_id <- keys_clean$koealatunnus_BIOSOIL[
  match(inputs_by_component$site, keys_clean$koealatunnus_MUSTIKKA)]
SOC$plot <- SOC$plot_id   # baseline is plot_id-native; keep `plot` for back-compatibility

inputs_matched              <- inputs[!is.na(inputs$plot_id), ]
inputs_by_component_matched <- inputs_by_component[!is.na(inputs_by_component$plot_id), ]
SOC_matched                 <- SOC[SOC$plot_id %in% unique(inputs_matched$plot_id), ]

cat("=== COVERAGE ===\n")
cat("Unique plots in inputs (matched):              ", length(unique(inputs_matched$plot_id)), "\n")
cat("Unique plots in SOC (matched):                 ", length(unique(SOC_matched$plot_id)), "\n")
cat("Plots with BOTH inputs & SOC:                  ",
    length(intersect(unique(inputs_matched$plot_id), unique(SOC_matched$plot_id))), "\n")


# SOC_dedup: deduplicated matched observations (1985 & 2006).
# Defined here so the Komeetta 2024 append (section 1.0b) can row-bind
# before the organic-layer quality flags in section 1.1 are computed.
SOC_dedup <- unique(SOC_matched)




# ========================================================================================
# ========================================================================================
# 1.0b Komeetta 2024 SOC data -- REMOVED 2026-08-04 (superseded)
# ========================================================================================
# The former separate ingest of ./Data/Komeetta/Komeetta_mitatut hiilet.xlsx (with its
# own layer map, plausibility filters and rbind into SOC_dedup) is GONE: the homogenized
# baseline read in section 1 already contains the Komeetta 2024 campaign on the same
# stoniness/BD-corrected basis as 1985 and 2006. Keeping a second, differently-processed
# 2024 ingest here is what made the campaigns mutually inconsistent.
#
# Two defects of that block are fixed by construction in the baseline:
#   - KOM_LAYER_MAP mislabelled the Komeetta mineral codes as the 1985 protocol
#     (201 = 0-5 cm, 202 = 5-20 cm). Komeetta in fact uses the 2006 protocol
#     (201 = 0-10 cm, 202 = 10-20 cm), which corrupted the depth extrapolation.
#   - the mineral stocks carried no coarse-fragment (stoniness) correction.
#
# Per-layer plausibility filtering, peat exclusion and outlier flagging now happen once,
# in Data/SOC_homogeneized/build_soc_homogenized.R, for all three campaigns alike.
# ========================================================================================




# ========================================================================================
# 1.1 Aggregating SOC measurements
# ========================================================================================


# --- Organic layer quality flags ---
# Flag 1: plots where organic layer is consistently absent (never recorded in any year)
has_organic_yr <- aggregate(layer ~ plot_id + year,
                            data = SOC_dedup,
                            FUN  = function(x) "organic" %in% x)
names(has_organic_yr)[3] <- "has_organic"
organic_wide <- reshape(has_organic_yr, idvar = "plot_id",
                        timevar = "year", direction = "wide")
organic_missing_plots <- organic_wide$plot_id[
  apply(organic_wide[, -1], 1, function(x) !any(x, na.rm = TRUE))]

# Flag 2: plots with zero organic C in any year (measurement issue, not absence)
zero_organic_plots <- unique(SOC_dedup$plot_id[
  SOC_dedup$layer == "organic" & SOC_dedup$C_kgha == 0])

# Both sets are mutually exclusive (verified) — combined exclusion set
organic_exclude_plots <- union(organic_missing_plots, zero_organic_plots)

cat("Organic layer quality exclusions:\n")
cat("  Consistently absent (no organic layer):    ", length(organic_missing_plots), "\n")
cat("  Zero organic C in at least one year:       ", length(zero_organic_plots), "\n")
cat("  Total excluded (non-overlapping):          ", length(organic_exclude_plots), "\n")

# --- SOC aggregation (C_kgha > 0 filter; organic zeros excluded via plot-level flags) ---
SOC_clean <- aggregate(C_kgha ~ plot_id + year + layer,
                       data = SOC_dedup[SOC_dedup$C_kgha > 0, ], FUN = mean)

has_deep      <- aggregate(layer ~ plot_id, data = SOC_clean,
                           FUN = function(x) any(x == "20-40cm"))
shallow_plots <- has_deep$plot_id[!has_deep$layer]
cat("\nShallow plots (no 20-40cm layer):", length(shallow_plots), "\n")

SOC_agg <- aggregate(C_kgha ~ plot_id + year, data = SOC_clean, FUN = sum)
SOC_agg$shallow         <- SOC_agg$plot_id %in% shallow_plots
SOC_agg$organic_missing <- SOC_agg$plot_id %in% organic_missing_plots
SOC_agg$organic_zero    <- SOC_agg$plot_id %in% zero_organic_plots

# soc_obs_tCha_sum: raw 0-40 cm + OFH sum — kept for records throughout
# soc_obs_tCha:     will be replaced by 1m extrapolation in section 1.3b
SOC_agg$soc_obs_tCha_sum <- SOC_agg$C_kgha / 1000
SOC_agg$soc_obs_tCha     <- SOC_agg$soc_obs_tCha_sum   # placeholder; overwritten in 1.3b

cat("Total plot-year observations:", nrow(SOC_agg), "\n")
cat("Unique plots:                ", length(unique(SOC_agg$plot_id)), "\n")


# ========================================================================================
# 1.2 plot_data construction
# ========================================================================================

inputs_total <- aggregate(Cha ~ plot_id + year, data = inputs_matched, FUN = sum)
common_plots <- intersect(unique(inputs_total$plot_id), unique(SOC_agg$plot_id))

avg_inputs <- aggregate(Cha ~ plot_id,
                        data = inputs_total[inputs_total$plot_id %in% common_plots, ],
                        FUN = mean)
avg_SOC <- aggregate(C_kgha ~ plot_id,
                     data = SOC_agg[SOC_agg$plot_id %in% common_plots, ], FUN = mean)
avg_SOC$C_Mgha <- avg_SOC$C_kgha / 1000

plot_data         <- merge(avg_inputs, avg_SOC, by = "plot_id")
plot_data$shallow <- plot_data$plot_id %in% shallow_plots

plot_region <- aggregate(region ~ plot_id, data = inputs_matched, FUN = function(x) x[1])
plot_data   <- merge(plot_data, plot_region, by = "plot_id")

inputs_by_species <- aggregate(Cha ~ plot_id + species, data = inputs_matched, FUN = sum)
dominant_species  <- do.call(rbind, lapply(split(inputs_by_species, inputs_by_species$plot_id),
                                           function(x) x[which.max(x$Cha), c("plot_id", "species")]))
plot_data <- merge(plot_data, dominant_species, by = "plot_id")

species_labels <- c("1" = "Scots pine", "2" = "Norway spruce", "3" = "Birch",
                    "4" = "Aspen",      "5" = "Grey alder",    "6" = "Black alder",
                    "7" = "Other broadleaved")
plot_data$species_name <- species_labels[as.character(plot_data$species)]

species_totals <- aggregate(Cha ~ plot_id + species, data = inputs_matched, FUN = sum)
plot_totals    <- aggregate(Cha ~ plot_id, data = inputs_matched, FUN = sum)
species_totals <- merge(species_totals, plot_totals, by = "plot_id", suffixes = c("_sp", "_total"))
species_totals$frac <- species_totals$Cha_sp / species_totals$Cha_total

species_wide <- reshape(species_totals[, c("plot_id", "species", "frac")],
                        idvar = "plot_id", timevar = "species", direction = "wide")
names(species_wide) <- sub("^frac\\.", "sp_frac_", names(species_wide))
species_wide[is.na(species_wide)] <- 0
plot_data <- merge(plot_data, species_wide, by = "plot_id", all.x = TRUE)

obs_count <- aggregate(soc_obs_tCha ~ plot_id, data = SOC_agg, FUN = length)
names(obs_count)[2] <- "n_soc_obs"
plot_data <- merge(plot_data, obs_count, by = "plot_id", all.x = TRUE)

cv_litter <- aggregate(Cha ~ plot_id, data = inputs_total,
                       FUN = function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE))
names(cv_litter)[2] <- "litter_cv"
plot_data <- merge(plot_data, cv_litter, by = "plot_id", all.x = TRUE)


# ========================================================================================
# 1.3 Soil type extraction (GTK API)
# ========================================================================================

library(sf)
library(httr)
library(jsonlite)

plot_coords <- unique(inputs_matched[, c("plot_id", "x", "y")])
plot_coords <- plot_coords[!is.na(plot_coords$x) & !is.na(plot_coords$y), ]

plots_sf <- st_as_sf(plot_coords, coords = c("x", "y"), crs = 2393)
plots_sf <- st_transform(plots_sf, crs = 3067)

# ⚠ 2026-08-12 -- LOOP BREAK. This used to be restricted to plot_data$plot_id, and
# plot_data is `merge(avg_inputs, avg_SOC)` keyed on `common_plots`, i.e. an INNER
# join on the SOC data. That made soil_code (and hence the per-GTK-class lambda in
# build_soc_homogenized.R) depend on which plots had SOC observations -- while the
# SOC baseline itself reads soil_code back from site_raw.csv. A closed loop: drop a
# 1985-only plot's observation and the plot vanished from site_raw, changing the peat
# set and lambda on the NEXT build, which changed the SOC data again. Measured: it
# did not reach a fixed point in one pass (1411 vs 1408 plot-years, 38 lambdas).
# plots_sf comes from inputs_matched coordinates, which are SOC-independent, so
# extracting over the whole of it breaks the cycle. The downstream merge into
# plot_data is all.x = TRUE, so a superset is harmless there.
coords_all <- st_coordinates(plots_sf)
plot_ids   <- plots_sf$plot_id

rest_url_1M <- "https://gtkdata.gtk.fi/arcgis/rest/services/Hasu/maapera/MapServer/3/query"

get_soil_type_1M <- function(x, y) {
  resp   <- GET(rest_url_1M, query = list(
    geometry = paste0(x, ",", y), geometryType = "esriGeometryPoint",
    inSR = 3067, spatialRel = "esriSpatialRelIntersects",
    outFields = "CODE,TEKSTI", returnGeometry = "false", f = "json"))
  result <- fromJSON(content(resp, "text"))
  if (length(result$features) == 0) return(data.frame(CODE = NA, TEKSTI = NA))
  result$features$attributes
}

soil_types_1M <- do.call(rbind, lapply(seq_len(nrow(coords_all)), function(i) {
  res <- get_soil_type_1M(coords_all[i, 1], coords_all[i, 2])
  res$plot_id <- plot_ids[i]
  res[, c("plot_id", "CODE", "TEKSTI")]
}))

plot_data <- plot_data[, !names(plot_data) %in% c("CODE", "TEKSTI", "soil_label")]
plot_data <- merge(plot_data, soil_types_1M[, c("plot_id", "CODE", "TEKSTI")],
                   by = "plot_id", all.x = TRUE)

soil_labels <- c(
  "HM" = "HM: Glaciofluvial", "Jk" = "Jk: Residual/lag",
  "Ka" = "Ka: Bedrock",       "KM" = "KM: Rocky moraine",
  "KSa" = "KSa: Coarse sand", "LSHk" = "LSHk: Sorted fine sand",
  "Mr" = "Mr: Till (moraine)", "Sa" = "Sa: Clay",
  "SHM" = "SHM: Sandy moraine", "SiMr" = "SiMr: Silty moraine",
  "SrHk" = "SrHk: Gravel",    "Tu" = "Tu: Peat", "Ve" = "Ve: Water")
plot_data$soil_label <- soil_labels[as.character(plot_data$CODE)]


# ========================================================================================
# ========================================================================================
# 1.3b SOC calibration target: whole-profile stock from the homogenized baseline
# ========================================================================================
# The in-place exponential depth extrapolation that used to live here (pooled per-GTK-class
# lambda fit, per-profile C0 rescaling, integration 40 -> 100 cm) has been REMOVED. The same
# extrapolation is now done once, upstream, in Data/SOC_homogeneized/build_soc_homogenized.R,
# for all three campaigns on the stoniness-corrected basis, and with one substantive
# improvement: the profile is capped at the depth augering actually reached.
#
# Why the cap matters: the old code extrapolated every plot to a flat 100 cm, including
# 38 plots where augering hit refusal at 10/20/40 cm. Those plots were assigned ~14 tC/ha
# of carbon below bedrock. The baseline instead integrates from the deepest measured layer
# down to z_cap = min(100 cm, soil depth reached), so shallow soils stay shallow.
#
# The extrapolation itself is validated, not assumed: Biosoil layer Krs 204 is a MEASURED
# 40-80 cm layer on 501 plots, and the predicted 40-80 stock tracks it with median pred/obs
# 0.96 (bias -2.5 Mg/ha, i.e. mildly conservative).
#
# Target column: soc_profile_Mgha -- organic + measured mineral 0-40 cm + modelled
# (deepest measured -> z_cap). This is a per-plot VARIABLE depth, matching how each plot's
# carbon was actually measured, rather than a fixed standardized depth.
#
# soc_obs_tCha_sum (raw organic + 0-40 cm sum, computed in 1.1) is retained for records.
# ========================================================================================

soc_plot_baseline <- read.csv("./Data/SOC_homogeneized/soc_homogenized_plot.csv")

cat("\n=== SOC target: homogenized baseline (plot level) ===\n")
cat(sprintf("Baseline plot-years: %d across %d plots (years: %s)\n",
            nrow(soc_plot_baseline), length(unique(soc_plot_baseline$plot_id)),
            paste(sort(unique(soc_plot_baseline$year)), collapse = ", ")))

SOC_agg <- merge(
  SOC_agg,
  soc_plot_baseline[, c("plot_id", "year", "soc_profile_Mgha", "soc_0_40_Mgha",
                        "soc_deep_Mgha", "z_cap", "soil_depth_reached",
                        "lambda", "fit_ok", "soc_outlier", "high_change",
                        "samp_year", "lm_added_1985", "lm_imputed_1985")],
  by = c("plot_id", "year"), all.x = TRUE)

# Primary calibration target. Fall back to the raw measured sum where the baseline has no
# profile value (should be none for calib-ready plots; asserted below).
SOC_agg$soc_obs_tCha <- ifelse(
  !is.na(SOC_agg$soc_profile_Mgha),
  SOC_agg$soc_profile_Mgha,
  SOC_agg$soc_obs_tCha_sum)

# ---------------------------------------------------------------------------
# TRUE OBSERVATION YEAR (2026-08-12). The campaign labelled 1985 was actually
# sampled 1986-1995: 19.5% of plots in 1995, and the campaign mean represents
# ~1989. `year` REMAINS the campaign key (1985/2006/2024) because it keys
# SIGMA_1985_INFL and HIKET_DROP_CAMPAIGN downstream; `obs_year` is the year the
# model must be evaluated at. Do not merge the two.
SOC_agg$obs_year <- ifelse(SOC_agg$year == 1985L & is.finite(SOC_agg$samp_year),
                           SOC_agg$samp_year, SOC_agg$year)

# Plots whose first-campaign sampling year could not be recovered are DROPPED for
# that campaign only (they keep 2006/2024). Assigning them 1985 would keep a
# known dating error of up to 10 years; assigning them the median would invent
# one. ⚠ The dropped set is ~71 South / 9 North, so it shifts the regional
# balance of the FIRST CAMPAIGN ONLY -- stated, not compensated.
.drop85 <- SOC_agg$year == 1985L & !is.finite(SOC_agg$samp_year)
cat(sprintf("\nTrue observation year: %d of %d 1985 plot-years dated (%s); %d dropped for want of a year\n",
            sum(SOC_agg$year == 1985L & is.finite(SOC_agg$samp_year)),
            sum(SOC_agg$year == 1985L),
            paste(names(table(SOC_agg$obs_year[SOC_agg$year == 1985L])),
                  table(SOC_agg$obs_year[SOC_agg$year == 1985L]), sep = ":", collapse = " "),
            sum(.drop85)))
SOC_agg <- SOC_agg[!.drop85, ]

SOC_agg$soc_extrap_ok <- !is.na(SOC_agg$soc_profile_Mgha)

n_missing <- sum(!SOC_agg$soc_extrap_ok)
cat(sprintf("SOC_agg matched to baseline: %d plot-years; %d fell back to raw 0-40 sum\n",
            sum(SOC_agg$soc_extrap_ok), n_missing))

# Consistency check: our 1.1 layer sum should reproduce the baseline 0-40 stock.
chk <- SOC_agg[SOC_agg$soc_extrap_ok, ]
d40 <- chk$soc_obs_tCha_sum - chk$soc_0_40_Mgha
cat(sprintf("Layer-sum vs baseline 0-40: median diff %.3f tC/ha (max |diff| %.3f)\n",
            median(d40, na.rm = TRUE), max(abs(d40), na.rm = TRUE)))

cat("\nCalibration target soc_obs_tCha (tC/ha), by year:\n")
print(round(tapply(SOC_agg$soc_obs_tCha, SOC_agg$year, median, na.rm = TRUE), 1))
cat("\nOverall summary:\n"); print(summary(SOC_agg$soc_obs_tCha))

# --- SOC outlier plots (excluded from calibration; see section 5) --------------------
# Rule (set in build_soc_homogenized.R): soc_profile > 250 Mg/ha in ANY campaign year.
# These are implausible for Finnish forest mineral soil and are almost certainly source
# errors rather than real profiles (e.g. plot 33631 reads 290-417 Mg/ha in all three
# campaigns). Flagged at the plot level so a plot is excluded consistently across years.
# NOTE: high_change (|consecutive rate| > 3 tC/ha/yr) is FLAGGED but NOT excluded --
# those are read as resampling noise, not data errors.
soc_outlier_plots <- unique(soc_plot_baseline$plot_id[
  !is.na(soc_plot_baseline$soc_outlier) & soc_plot_baseline$soc_outlier])
cat(sprintf("\nSOC outlier plots (soc_profile > 250 Mg/ha, excluded from calib): %d\n",
            length(soc_outlier_plots)))
if (length(soc_outlier_plots) > 0)
  cat("  plot_ids:", paste(sort(soc_outlier_plots), collapse = ", "), "\n")

high_change_plots <- unique(soc_plot_baseline$plot_id[
  !is.na(soc_plot_baseline$high_change) & soc_plot_baseline$high_change])
cat(sprintf("High-change plots (|rate| > 3 tC/ha/yr, FLAGGED but retained): %d\n",
            length(high_change_plots)))


# --- Implied MRT filter ---
# Implied MRT = mean(1m SOC) / mean(annual litter input).
# Plots with MRT > 100y have litter inputs so low that the model predicts
# near-zero SOC at any parameter value, causing -Inf likelihood at the
# pre-MCMC sanity check. The 100y threshold matches a clear gap in the MRT
# distribution (see diagnostics_inputs.R output 01_mrt_histogram.png).
# This filter also catches plots with zero litter after AWEN mapping (MRT = Inf)
# that escape the raw zero-litter filter, which operates on inputs_matched
# before compound-fraction assignment.
MRT_EXCL_THRESHOLD <- 100L
mean_litter_mrt <- aggregate(Cha ~ plot_id, data = inputs_total, FUN = mean)
mean_soc_mrt    <- aggregate(soc_obs_tCha ~ plot_id, data = SOC_agg,    FUN = mean)
mrt_filter      <- merge(mean_litter_mrt, mean_soc_mrt, by = "plot_id")
mrt_filter$implied_mrt <- mrt_filter$soc_obs_tCha / mrt_filter$Cha
high_mrt_plots  <- mrt_filter$plot_id[
  !is.finite(mrt_filter$implied_mrt) | mrt_filter$implied_mrt > MRT_EXCL_THRESHOLD]
cat(sprintf("Implied MRT > %dy or Inf: %d plots flagged for exclusion\n",
            MRT_EXCL_THRESHOLD, length(high_mrt_plots)))


# ========================================================================================
# 1.4 Biosoil soil properties + Cajander fertility + biogeographic zone + profile depth
# ========================================================================================

biosoil_data  <- read.csv("../../Datasets/BioSoil_maaperäaineisto_2006/soil_data_main.csv")
biosoil_index <- read.csv("../../Datasets/BioSoil_maaperäaineisto_2006/plot_index.csv")
koeala_1985   <- read.csv("./Data/PysyvätKoealat/1985/koeala.csv")

# --- Link PlotIndex -> plot_id via koealatunnus_VANHA ---
biosoil_index$plot_id <- keys$koealatunnus_BIOSOIL[
  match(biosoil_index$KOEALA, keys$koealatunnus_VANHA)]

cat("Biosoil plots matched to plot_id:", sum(!is.na(biosoil_index$plot_id)),
    "of", nrow(biosoil_index), "\n")

# --- Enrich biosoil_index with all zone variants from koeala_1985 ---
# koeala_1985$bio_soil_id == biosoil_index$PlotIndex
# Zone variants (semantics, partial guesses for sk and eliomk):
#   kasvyo_syke  - SYKE (administrative) biogeographic zones, 5 levels S->N
#   kasvyo_ahti  - Ahti (1968) classical Nordic vegetation zonation, alternative
#                  boundaries (more floristic)
#   kasvyo_sk    - third zone variant; meaning to confirm with Finnish colleagues
#   alavyo_syke  - SYKE sub-zones (finer than kasvyo_syke)
#   suovyohyke   - mire vegetation zone (peatland-oriented; many NAs expected
#                  on mineral calib_ready plots)
#   alavyohyke   - sub-zone refinement (finer geographic granularity)
#   eliomk       - biogeographic refinement; meaning to confirm
zone_vars <- c("kasvyo_syke", "kasvyo_ahti", "kasvyo_sk",
               "alavyo_syke", "suovyohyke",  "alavyohyke", "eliomk")
biosoil_index <- merge(
  biosoil_index,
  koeala_1985[, c("bio_soil_id", zone_vars)],
  by.x = "PlotIndex", by.y = "bio_soil_id",
  all.x = TRUE
)

# --- Add plot_id to biosoil_data via PlotIndex ---
biosoil_data$plot_id <- biosoil_index$plot_id[
  match(biosoil_data$PlotIndex, biosoil_index$PlotIndex)]

# --- Mineral topsoil (M01 + M12): average numeric variables across layers ---
mineral_vars <- c("ClayContent", "SiltContent", "SandContent",
                  "MeanBulkDensity", "EstimatedBulkDensity", "CoarseFragments",
                  "pH.CaCl2.", "pH.H2O.", "OrganicCarbon", "TotalNitrogen",
                  "ExchangeableAcidity", "ExchangeableAl", "ExchangeableCa",
                  "ExchangeableFe", "ExchangeableK", "ExchangeableMg",
                  "ExchangeableMn", "ExchangeableNa", "FreeHAcidity", "OrganicMatter")

mineral_top <- biosoil_data[biosoil_data$LayerCode %in% c("M01", "M12") &
                              !is.na(biosoil_data$plot_id), ]

biosoil_mineral <- aggregate(
  mineral_top[, mineral_vars],
  by  = list(plot_id = mineral_top$plot_id),
  FUN = function(x) mean(x, na.rm = TRUE)
)

# TexturalClass: most common value across M01 + M12
textural_mode <- do.call(rbind, lapply(
  split(mineral_top, mineral_top$plot_id), function(d) {
    tc <- d$TexturalClass[d$TexturalClass != "" & !is.na(d$TexturalClass)]
    data.frame(
      plot_id       = d$plot_id[1],
      TexturalClass = if (length(tc) > 0) names(sort(table(tc), decreasing = TRUE))[1]
      else NA_character_
    )
  }))

biosoil_mineral <- merge(biosoil_mineral, textural_mode, by = "plot_id", all.x = TRUE)

# --- Organic layer (OFH only): depth and weight ---
ofh <- biosoil_data[biosoil_data$LayerCode == "OFH" & !is.na(biosoil_data$plot_id), ]

biosoil_ofh <- aggregate(
  ofh[, c("UpperDepthLimit", "LowerDepthLimit", "OrganicLayerWeight")],
  by  = list(plot_id = ofh$plot_id),
  FUN = function(x) mean(x, na.rm = TRUE)
)
names(biosoil_ofh)[names(biosoil_ofh) == "UpperDepthLimit"]    <- "ofh_upper_cm"
names(biosoil_ofh)[names(biosoil_ofh) == "LowerDepthLimit"]    <- "ofh_lower_cm"
names(biosoil_ofh)[names(biosoil_ofh) == "OrganicLayerWeight"] <- "ofh_weight_kgm2"

# --- Measured mineral profile depth (all layers: M01/M12/M24/M48) ---
mineral_biosoil_all <- biosoil_data[grepl("^M", biosoil_data$LayerCode) &
                                      !is.na(biosoil_data$plot_id), ]
profile_depth <- aggregate(LowerDepthLimit ~ plot_id, data = mineral_biosoil_all,
                           FUN = max, na.rm = TRUE)
names(profile_depth)[2] <- "profile_depth_cm"

cat("Profile depth distribution:\n")
print(table(profile_depth$profile_depth_cm))

# --- Cajander fertility class and biogeographic zone ---
ka_lookup <- unique(biosoil_index[!is.na(biosoil_index$plot_id),
                                  c("plot_id", "KA", zone_vars)])
# A small number of plots have multiple rows with conflicting KA or zone values
# (e.g. plot 83552 has KA=12 and KA=2). Keep the first row per plot: the
# conservative choice, since any peatland KA (11-13) in the first row will
# correctly trigger the downstream exclusion filter.
ka_lookup <- ka_lookup[!duplicated(ka_lookup$plot_id), ]

cat("KA distribution in Biosoil plots:\n")
print(table(ka_lookup$KA, useNA = "always"))
cat("kasvyo_syke distribution in Biosoil plots:\n")
print(table(ka_lookup$kasvyo_syke, useNA = "always"))

# =========================================================================================
# site_attributes.csv -- the SOC-INDEPENDENT site table (added 2026-08-12)
# -----------------------------------------------------------------------------------------
# build_soc_homogenized.R must NOT read site_raw.csv: that file is built from plot_data,
# whose row set is an inner join on the SOC data the builder itself produces. Everything the
# builder actually needs is a RAW attribute -- GTK soil class, Cajander peat class, region,
# coordinates -- none of which is a modelling result. They are emitted here, over the plot
# universe defined by the litter-input coordinates, BEFORE anything SOC-derived is used.
# ⚠ Nothing SOC-dependent may be added to this table. The invariant to preserve is the ROW
# SET, not just the columns: gating the rows on SOC would restore the loop invisibly.
# =========================================================================================
site_attributes <- merge(
  data.frame(plot_id = plots_sf$plot_id,
             x_ETRS  = st_coordinates(plots_sf)[, 1],
             y_ETRS  = st_coordinates(plots_sf)[, 2]),
  soil_types_1M[, c("plot_id", "CODE", "TEKSTI")], by = "plot_id", all.x = TRUE)
names(site_attributes)[names(site_attributes) == "CODE"]   <- "soil_code"
names(site_attributes)[names(site_attributes) == "TEKSTI"] <- "soil_type"
site_attributes <- merge(site_attributes, plot_region,               by = "plot_id", all.x = TRUE)
site_attributes <- merge(site_attributes, ka_lookup[, c("plot_id", "KA")],
                         by = "plot_id", all.x = TRUE)
site_attributes$peatland <- !is.na(site_attributes$KA) & site_attributes$KA %in% c(11L, 12L, 13L)
write.csv(site_attributes, "./Data/model_inputs/site_attributes.csv", row.names = FALSE)
cat(sprintf("\nsite_attributes.csv: %d plots (SOC-independent), %d peat, %d with a soil class\n",
            nrow(site_attributes), sum(site_attributes$peatland),
            sum(!is.na(site_attributes$soil_code))))
cat("kasvyo_ahti distribution in Biosoil plots:\n")
print(table(ka_lookup$kasvyo_ahti, useNA = "always"))

# --- Combine all Biosoil-derived tables ---
biosoil_plot <- merge(biosoil_mineral, biosoil_ofh,    by = "plot_id", all = TRUE)
biosoil_plot <- merge(biosoil_plot,    profile_depth,  by = "plot_id", all.x = TRUE)
biosoil_plot <- merge(biosoil_plot,    ka_lookup,      by = "plot_id", all.x = TRUE)

cat("Plots with Biosoil data:        ", nrow(biosoil_plot), "\n")
cat("Plots with SOC and Biosoil:     ",
    sum(unique(SOC_agg$plot_id) %in% biosoil_plot$plot_id), "\n")
cat("SOC plot coverage:              ",
    round(100 * mean(unique(SOC_agg$plot_id) %in% biosoil_plot$plot_id), 1), "%\n")


# ========================================================================================
# 1.5 Coordinates
# ========================================================================================

library(rnaturalearth)

all_plots_sf <- st_as_sf(
  unique(inputs_matched[!is.na(inputs_matched$x), c("plot_id", "x", "y")]),
  coords = c("x", "y"), crs = 2393)
all_plots_sf <- st_transform(all_plots_sf, crs = 3067)

all_plots_wgs84 <- st_transform(all_plots_sf, crs = 4326)
coords_wgs84    <- st_coordinates(all_plots_wgs84)
coords_etrs     <- st_coordinates(all_plots_sf)

coords_lookup <- data.frame(
  plot_id   = all_plots_sf$plot_id,
  x_KKJ    = unique(inputs_matched[!is.na(inputs_matched$x), c("plot_id", "x", "y")])$x,
  y_KKJ    = unique(inputs_matched[!is.na(inputs_matched$x), c("plot_id", "x", "y")])$y,
  x_ETRS   = coords_etrs[, 1],
  y_ETRS   = coords_etrs[, 2],
  lon_WGS84 = coords_wgs84[, 1],
  lat_WGS84 = coords_wgs84[, 2]
)

write.csv(coords_lookup, "NFI_plots_litter_inputs.csv", row.names = FALSE)


# ========================================================================================
# 1.6 NFI permanent plot stand attributes (kuvio.csv across 1985, 1990, 1995)
# ========================================================================================
# Pulls stand-level attributes from PysyvätKoealat (Finnish NFI permanent plots).
# Three measurement waves available, each with different scope:
#   1985: full inventory + biological surveys; cuttings recorded ≤10y back,
#         soil prep ≤30y back (full retrospective)
#   1990: only NEW cuttings/treatments since 1985
#   1995: complete re-measurement; cuttings 1990-1995 only
# Combined: management indicators covering ~1955-1995 (max).
#
# Linkage: kuvio$koeala == keys$koealatunnus_VANHA -> plot_id
#
# Variable codes documented in PysyvätKoealat/{year}/kuvio.md
# Notes on selected variables:
#   mets_ika_1     = mean stand age, layer 1 (years)
#   kuvion_pohja_1 = stand basal area (m²/ha)
#   keskipituus_1  = mean tree height (dm; *not* cm)
#   keskilapim_1   = mean DBH (cm)
#   vall_puul_1    = dominant species (0=treeless, 1=Scots pine, 2=Norway spruce,
#                    3=silver birch, 4=downy birch, 5=aspen, 6=grey alder,
#                    7=black alder, 8=other coniferous, 9=other broadleaved)
#                    -- DIFFERENT coding from species_code (which uses litter convention)
#   keh_luokka     = development class (0=open, 1=seed-tree, 2=small sapling,
#                    3=mature sapling, 4=young thinning, 5=mature thinning,
#                    6=regeneration-ready, 7=shelterwood)
#   kasvup_tyyppi  = site fertility (1-8, finer than KA's 4 mineral classes)
#   alaryhma       = main site type (1=upland, 2=fen, 3=bog, 4=open mire, 5=rich fen)
#   ojitustilanne  = drainage status (0=undrained upland, 1=drained upland tree-growth,
#                    2=undrained peatland, 3=drained peatland tree-growth, etc.;
#                    letter codes for "no tree growth" variants)
#   per_tapa       = origin (1=natural old land, 2=failed planting old, 3=successful
#                    planting old, 4-6 same on new land)
#   tuhon_laatu    = damage type (0=none, 1=standing death, 2=fallen/broken,
#                    3=decay, 4=stem defect, 5=dead crown, 6=other crown,
#                    7=needle loss, 8=discoloration; A-F variants for >5y old)
#   tuhon_syy      = damage cause (1=wind, 2=snow, 3=climate, 4=plant competition,
#                    5=harvesting, 7=vole, 8=moose, 9=pine weevil, B=blister rust, ...)
#   veroluokka     = taxation class (productivity proxy: 0=IA, 1=IB, 2=II, 3=III, 4=IV)
#   lamposumma     = temperature sum (Finnish forestry heat-sum index)
#   maanp_muoto    = landform (0=flat/gentle, 1=slope, 2=hilltop, 3=valley)
#   hakk_laatu     = cutting type (0=none in 10y, 1=sapling thinning, 2=overstory
#                    removal, 3=first thinning, 4=other thinning, 5=selective,
#                    7=artificial regen, 8=natural regen, 9=clearing)
#   mets_h_toim   = treatment without removal (0=none, 1=planting, 2=seeding,
#                    3=natural regen, 5=pruning; A,B = combined)
#   aiempi_toim    = soil preparation (0=none, 1=plowing, 2=plow+drainage,
#                    3=harrowing, 4=mounding, 5=controlled burning, 6=drainage,
#                    7=re-drainage, 8=ditch maintenance) -- 1985 covers up to 30y back
# ========================================================================================

cat("\n=== KUVIO (NFI stand) attribute extraction ===\n")

kuvio_85 <- read.csv("./Data/PysyvätKoealat/1985/kuvio.csv")
kuvio_90 <- read.csv("./Data/PysyvätKoealat/1990/kuvio.csv")
kuvio_95 <- read.csv("./Data/PysyvätKoealat/1995/kuvio.csv")

# --- Helper: take primary compartment per plot (vmikuvio == 0; fallback to first row)
take_main_kuvio <- function(df) {
  if ("vmikuvio" %in% names(df)) {
    main   <- df[!is.na(df$vmikuvio) & df$vmikuvio == 0, ]
    others <- df[!df$koeala %in% main$koeala, ]
    others <- others[!duplicated(others$koeala), ]
    rbind(main, others)
  } else {
    df[!duplicated(df$koeala), ]
  }
}

k85 <- take_main_kuvio(kuvio_85)
k90 <- take_main_kuvio(kuvio_90)
k95 <- take_main_kuvio(kuvio_95)

cat(sprintf("Primary kuvio rows: 1985=%d, 1990=%d, 1995=%d\n",
            nrow(k85), nrow(k90), nrow(k95)))

# --- Helper: extract a variable with NA fallback if column missing
get_var <- function(df, var) {
  if (var %in% names(df)) df[[var]] else rep(NA, nrow(df))
}

# --- Static stand attributes from 1985 (most informative wave) ---
stand_attrs <- data.frame(
  koeala            = k85$koeala,
  stand_age_85      = get_var(k85, "mets_ika_1"),
  basal_area_85     = get_var(k85, "kuvion_pohja_1"),
  mean_height_85_dm = get_var(k85, "keskipituus_1"),
  dev_class_85      = get_var(k85, "keh_luokka"),
  kasvup_tyyppi     = get_var(k85, "kasvup_tyyppi"),
  alaryhma          = get_var(k85, "alaryhma"),
  ojitustilanne     = get_var(k85, "ojitustilanne"),
  temp_sum_NFI      = get_var(k85, "lamposumma"),
  elevation_m       = get_var(k85, "kork_merenp"),
  stringsAsFactors  = FALSE
)

# --- Stand age projected to 2006 ---
age_lkup <- function(df) setNames(get_var(df, "mets_ika_1"), as.character(df$koeala))
a85 <- age_lkup(k85); a90 <- age_lkup(k90); a95 <- age_lkup(k95)

stand_attrs$stand_age_2006_est <- vapply(
  as.character(stand_attrs$koeala),
  function(pid) {
    if (!is.na(a95[pid])) return(a95[pid] + 11)
    if (!is.na(a90[pid])) return(a90[pid] + 16)
    if (!is.na(a85[pid])) return(a85[pid] + 21)
    NA_real_
  },
  numeric(1)
)

# --- Management history flags across 1985, 1990, 1995 ---
is_event <- function(x) !is.na(x) & x != 0 & x != "0" & x != ""

flag_for <- function(df, var, ids) {
  evt_plots <- df$koeala[is_event(get_var(df, var))]
  ids %in% evt_plots
}

all_plots_kuvio <- unique(c(k85$koeala, k90$koeala, k95$koeala))

cut_85 <- flag_for(k85, "hakk_laatu",  all_plots_kuvio)
cut_90 <- flag_for(k90, "hakk_laatu",  all_plots_kuvio)
cut_95 <- flag_for(k95, "hakk_laatu",  all_plots_kuvio)
trt_85 <- flag_for(k85, "mets_h_toim", all_plots_kuvio)
trt_90 <- flag_for(k90, "mets_h_toim", all_plots_kuvio)
trt_95 <- flag_for(k95, "mets_h_toim", all_plots_kuvio)

mgmt <- data.frame(
  koeala          = all_plots_kuvio,
  any_cut_85_95   = cut_85 | cut_90 | cut_95,
  n_cuts_85_95    = as.integer(cut_85) + as.integer(cut_90) + as.integer(cut_95),
  any_trt_85_95   = trt_85 | trt_90 | trt_95,
  soil_prep_pre85 = flag_for(k85, "aiempi_toim", all_plots_kuvio)
)

kuvio_summary <- merge(stand_attrs, mgmt, by = "koeala", all.x = TRUE)
kuvio_summary$plot_id <- keys_clean$koealatunnus_BIOSOIL[
  match(kuvio_summary$koeala, keys_clean$koealatunnus_VANHA)]
kuvio_summary <- kuvio_summary[!is.na(kuvio_summary$plot_id), ]
kuvio_summary$koeala <- NULL

cat(sprintf("\nKuvio summary: %d plots linked to plot_id\n", nrow(kuvio_summary)))
cat("Stand age 1985 summary:\n");  print(summary(kuvio_summary$stand_age_85))
cat("Stand age 2006 (est) summary:\n"); print(summary(kuvio_summary$stand_age_2006_est))
cat("Cutting events 1985-1995 distribution:\n")
print(table(kuvio_summary$n_cuts_85_95, useNA = "always"))
cat("Treatment any 1985-1995 (TRUE/FALSE):\n")
print(table(kuvio_summary$any_trt_85_95, useNA = "always"))
cat("Drainage status (1985):\n")
print(table(kuvio_summary$ojitustilanne, useNA = "always"))
cat("Development class (1985):\n")
print(table(kuvio_summary$dev_class_85, useNA = "always"))




# ========================================================================================
# 2. Build monthly model input table
# ========================================================================================

cat("\n=== BUILDING MONTHLY INPUT TABLE ===\n")

litter_long <- aggregate(Cha ~ plot_id + year + size + compound,
                         data = inputs_matched[!is.na(inputs_matched$size), ], FUN = sum)
litter_long$Cha <- pmax(litter_long$Cha, 0)
litter_long$col <- paste0("C_", litter_long$size, "_", toupper(litter_long$compound))

# NOTE (Apr 2026 audit): the previous x10 multiplier ("kg m^-2 -> tC ha^-1") has been
# removed. Audit showed source values are already in tC/ha/yr (median 2.4 across plots,
# bang in the boreal range 1.5-4.5). The x10 was a guess for a bias of unknown origin
# which has since been resolved at the data source.

input_wide <- reshape(litter_long[, c("plot_id", "year", "col", "Cha")],
                      idvar = c("plot_id", "year"), timevar = "col", direction = "wide")
names(input_wide) <- sub("^Cha\\.", "", names(input_wide))

litter_cols <- grep("^C_", names(input_wide), value = TRUE)
input_wide[, litter_cols][is.na(input_wide[, litter_cols])] <- 0

# NOTE: C_cwl_* (coarse woody litter = stumps) comes directly from the Tupek et al.
# dataset where cwl is a proper size class. A previous placeholder (cwl = 0.30 × fwl)
# has been removed. All three size classes (nwl, fwl, cwl) are present in the source data.

# --- 1985 litter: RECONSTRUCTED (first year of the source series is unusable) -----------
# The Tupek litter product's FIRST year is near-zero and not a real ecological signal:
# source median total litter is 0.060 tC/ha/yr in 1985 against 1.398 in 1986, with an
# identical row count (2805 site-years) in both, so it is not missing data -- the values
# themselves collapse. This is the expected signature of litter derived from between-
# inventory biomass increments, where the first year has no predecessor to difference
# against. (TO CONFIRM with B. Tupek before this wording goes in the manuscript.)
#
# Why it must be fixed rather than tolerated: 1985 is t0. The artefactual year contaminates
# BOTH pre-run litter anchors and the flux bridge --
#   J_t0_mean   = mean of the first 5 years  -> 1/5 of the window   (~22% too low)
#   J_full_mean = mean of the whole series   (1917 pre-run anchor)
#   J_total_mean= mean of the first 20 years -> feeds J_bar, the flux_pair units bridge
# -- and it is also the forward run's own first year, landing exactly on the VMI8
# observation. Left uncorrected the 1917->1985 spin-up terminates ~22% below the true
# 1985 litter level, biasing C_init low in all six models.
#
# METHOD: per plot and per AWEN component, fit a linear trend over 1986-1990 and backcast
# ONE year to 1985. Chosen over the plain 1986-1990 mean (1.949) because litter is RISING
# through that window (1986 = 1.84 -> 1990 = 2.17): a flat mean would place 1985 ABOVE
# 1986, contradicting the growing-stock history used for the pre-run shape (C3). The
# backcast (1.754) keeps 1985 just below 1986, as the forest history implies. Carry-back
# of 1986 (1.838) gives nearly the same answer -- all three options agree within ~3% on
# J_t0_mean, so the result is not sensitive to this choice.
BACKCAST_YEARS <- 1986:1990
stopifnot(all(BACKCAST_YEARS %in% input_wide$year))

.backcast_1985 <- function(df) {
  fit_yrs <- BACKCAST_YEARS
  for (cc in litter_cols) {
    w <- reshape(df[df$year %in% fit_yrs, c("plot_id", "year", cc)],
                 idvar = "plot_id", timevar = "year", direction = "wide")
    M   <- as.matrix(w[, -1, drop = FALSE])
    # slope per plot from the 1986-1990 trend; backcast one year below 1986
    sl  <- apply(M, 1, function(v)
      if (all(is.finite(v))) unname(coef(lm(v ~ fit_yrs))[2]) else 0)
    y86 <- M[, 1]
    pred <- pmax(y86 - sl, 0)                      # non-negative litter
    names(pred) <- as.character(w$plot_id)
    ix  <- df$year == 1985L
    df[ix, cc] <- unname(pred[as.character(df$plot_id[ix])])
  }
  df
}

j85_before <- median(rowSums(input_wide[input_wide$year == 1985L, litter_cols, drop = FALSE]))
input_wide <- .backcast_1985(input_wide)
j85_after  <- median(rowSums(input_wide[input_wide$year == 1985L, litter_cols, drop = FALSE]))
cat(sprintf(paste0("1985 litter reconstructed by 1986-1990 backcast: ",
                   "median %.3f -> %.3f tC/ha/yr (1986 = %.3f)\n"),
            j85_before, j85_after,
            median(rowSums(input_wide[input_wide$year == 1986L, litter_cols, drop = FALSE]))))
if (!is.finite(j85_after) || j85_after < 0.5)
  stop("1985 backcast failed -- median litter still implausibly low.")

# --- Plots with a CONSTANT litter series (detected, not hard-coded) --------------------
# Ten plots carry litter that is IDENTICAL to four decimal places for every year 1986-2024,
# in the source as well as here. No stand produces identical litter for 39 consecutive
# years; this is the litter model returning a fixed value where the underlying biomass
# series was static. It is the same class of defect as the artefactual 1985 first year --
# plausible-looking, positive, finite, and wrong -- and it survived every earlier check for
# the same reason: nothing was testing per-plot constancy.
#
# Why they are EXCLUDED rather than kept and flagged: a constant series injects a spurious
# "no temporal change" forcing into a study whose entire subject is the temporal change.
# The cost is 10 of 520 plots (2%). Detected here rather than listed, so the rule survives
# any change in the litter product; the plot IDs are echoed for the record.
LITTER_CONST_CV_TOL <- 1e-6   # relative SD below which a series counts as constant

.ann_litter <- aggregate(
  list(J = rowSums(input_wide[, litter_cols, drop = FALSE])),
  by = list(plot_id = input_wide$plot_id, year = input_wide$year), FUN = sum)
.post85 <- .ann_litter[.ann_litter$year >= 1986L, ]   # exclude the reconstructed 1985
.cv <- tapply(.post85$J, .post85$plot_id,
              function(v) if (mean(v) > 0) sd(v) / mean(v) else NA_real_)
const_litter_plots <- as.integer(names(.cv)[!is.na(.cv) & .cv < LITTER_CONST_CV_TOL])
cat(sprintf("Constant-litter plots (CV < %.0e over 1986-2024, EXCLUDED from calib): %d\n",
            LITTER_CONST_CV_TOL, length(const_litter_plots)))
if (length(const_litter_plots))
  cat("  plot_ids:", paste(sort(const_litter_plots), collapse = ", "), "\n")

input_monthly <- input_wide[rep(seq_len(nrow(input_wide)), each = 12), ]
input_monthly$month <- rep(1:12, times = nrow(input_wide))
rownames(input_monthly) <- NULL

litter_cols_all <- grep("^C_", names(input_monthly), value = TRUE)
input_monthly[, litter_cols_all] <- input_monthly[, litter_cols_all] / 12

# --- Extend input_monthly to 2024 for Komeetta SOC placement ---
# Litter inputs are not available for 2024 (MUSTIKKA ends 2023). Since models
# are annual and SOC stocks change on decadal timescales, carrying forward the
# 2023 litter row for each plot introduces negligible error. This avoids NA
# litter in the calibration engine and keeps the 2024 Komeetta observation at
# its correct year label rather than displacing it to 2023.
rows_2023 <- input_monthly[input_monthly$year == 2023, ]
rows_2024 <- rows_2023
rows_2024$year <- 2024L
# soc_obs_tCha and soc_obs_tCha_sum NOT assigned here —
# input_monthly doesn't have them yet; they are initialised
# to NA for all rows (including 2024) in the placement block below.
input_monthly <- rbind(input_monthly, rows_2024)
cat(sprintf("input_monthly extended to 2024: %d plots × 12 months added\n",
            nrow(rows_2024) / 12L))

# soc_obs_tCha: 1m-extrapolated values (from SOC_agg, updated in section 1.3b)
# Placed in June of each observation year.
# Komeetta 2024 observations land in year 2024 (rows added above).
input_monthly$soc_obs_tCha     <- NA_real_
input_monthly$soc_obs_tCha_sum <- NA_real_
# The observation still SITS at its campaign-key year row; soc_obs_year records the
# year the model must be evaluated at (1986-1995 for VMI8). See section 1.3b.
input_monthly$soc_obs_year     <- NA_integer_

for (i in seq_len(nrow(SOC_agg))) {
  idx <- input_monthly$plot_id == SOC_agg$plot_id[i] &
    input_monthly$year    == SOC_agg$year[i]    &
    input_monthly$month   == 6L
  if (any(idx)) {
    input_monthly$soc_obs_tCha[idx]     <- SOC_agg$soc_obs_tCha[i]
    input_monthly$soc_obs_tCha_sum[idx] <- SOC_agg$soc_obs_tCha_sum[i]
    input_monthly$soc_obs_year[idx]     <- as.integer(SOC_agg$obs_year[i])
  }
}

cat("SOC observations placed (1m extrapolated):", sum(!is.na(input_monthly$soc_obs_tCha)), "\n")
cat("  1985:", sum(!is.na(input_monthly$soc_obs_tCha) & input_monthly$year == 1985), "\n")
cat("  2006:", sum(!is.na(input_monthly$soc_obs_tCha) & input_monthly$year == 2006), "\n")
cat("  2024:", sum(!is.na(input_monthly$soc_obs_tCha) & input_monthly$year == 2024), "\n")

input_monthly$temp_air <- NA_real_
input_monthly$precip   <- NA_real_
input_monthly$evap     <- NA_real_

input_monthly <- input_monthly[, c(
  "plot_id", "year", "month",
  "C_nwl_A", "C_nwl_W", "C_nwl_E", "C_nwl_N",
  "C_fwl_A", "C_fwl_W", "C_fwl_E", "C_fwl_N",
  "C_cwl_A", "C_cwl_W", "C_cwl_E", "C_cwl_N",
  "temp_air", "precip", "evap", "soc_obs_tCha", "soc_obs_tCha_sum",
  "soc_obs_year")]

input_monthly <- input_monthly[order(input_monthly$plot_id, input_monthly$year,
                                     input_monthly$month), ]
rownames(input_monthly) <- NULL

cat("Pre-filter dimensions:", nrow(input_monthly), "×", ncol(input_monthly), "\n")
cat("Pre-filter plots:     ", length(unique(input_monthly$plot_id)), "\n")


# ========================================================================================
# 2.1 Preliminary site_raw
# ========================================================================================

site_raw_prelim <- data.frame(
  plot_id = plot_data$plot_id,
  region  = plot_data$region
)


# ========================================================================================
# 3. Litter data quality checks (abbreviated)
# ========================================================================================

dir.create("./Data/data_quality_check", showWarnings = FALSE)

input_annual_diag <- aggregate(
  cbind(C_nwl_A, C_nwl_W, C_nwl_E, C_nwl_N,
        C_fwl_A, C_fwl_W, C_fwl_E, C_fwl_N,
        C_cwl_A, C_cwl_W, C_cwl_E, C_cwl_N) ~ plot_id + year,
  data = input_monthly, FUN = sum)

litter_diag_cols <- grep("^C_", names(input_annual_diag), value = TRUE)
input_annual_diag$total <- rowSums(input_annual_diag[, litter_diag_cols], na.rm = TRUE)

input_calib_diag <- input_annual_diag[input_annual_diag$plot_id %in% plot_data$plot_id, ]
input_calib_diag <- merge(input_calib_diag,
                          plot_data[, c("plot_id", "region", "species_name", "CODE")],
                          by = "plot_id")

mean_per_plot <- aggregate(total ~ plot_id, data = input_calib_diag, FUN = mean)

png("./Data/data_quality_check/01_hist_mean_total_litter.png", width = 800, height = 600)
hist(mean_per_plot$total, breaks = 40, col = "steelblue", border = "white",
     main = "Distribution of mean annual litter input across plots",
     xlab = "Mean annual total litter (Mg C ha⁻¹ yr⁻¹)", ylab = "Number of plots")
abline(v = median(mean_per_plot$total), col = "tomato", lwd = 2, lty = 2)
dev.off()


# ========================================================================================
# 4. Climate data extraction, PET, merge
# ========================================================================================

library(ncdf4)

nc              <- nc_open("./Data/Climate/nfi_plot_weather_data_1961_2025.nc")
nc_plot_ids_int <- as.integer(ncvar_get(nc, "plot_id"))
time_vals       <- ncvar_get(nc, "time")
time_dates      <- as.Date(time_vals, origin = "1961-01-01")

calib_ids  <- unique(site_raw_prelim$plot_id[site_raw_prelim$plot_id %in% nc_plot_ids_int])
nc_indices <- match(calib_ids, nc_plot_ids_int)

dup_in_nc <- sum(duplicated(nc_plot_ids_int))
if (dup_in_nc > 0) {
  message(sprintf("NOTE: %d duplicate plot_ids in netCDF — keeping first occurrence only",
                  dup_in_nc))
}

cat("Calibration plots matched in netCDF:", length(calib_ids), "\n")
cat("Unmatched calibration plots:        ",
    sum(!site_raw_prelim$plot_id %in% nc_plot_ids_int), "\n")

n_time <- length(time_dates); n_plots <- length(nc_indices)
Tavg  <- matrix(NA, nrow = n_time, ncol = n_plots)
Tmin  <- matrix(NA, nrow = n_time, ncol = n_plots)
Tmax  <- matrix(NA, nrow = n_time, ncol = n_plots)
Prec  <- matrix(NA, nrow = n_time, ncol = n_plots)
GlobR <- matrix(NA, nrow = n_time, ncol = n_plots)

cat("Extracting daily climate for", n_plots, "plots...\n")
for (i in seq_along(nc_indices)) {
  if (i %% 50 == 0) cat("  plot", i, "of", n_plots, "\n")
  idx <- nc_indices[i]
  Tavg[, i]  <- ncvar_get(nc, "Tavg",  start = c(1, idx), count = c(-1, 1))
  Tmin[, i]  <- ncvar_get(nc, "Tmin",  start = c(1, idx), count = c(-1, 1))
  Tmax[, i]  <- ncvar_get(nc, "Tmax",  start = c(1, idx), count = c(-1, 1))
  Prec[, i]  <- ncvar_get(nc, "Prec",  start = c(1, idx), count = c(-1, 1))
  GlobR[, i] <- ncvar_get(nc, "GlobR", start = c(1, idx), count = c(-1, 1))
}
nc_close(nc)
cat("Extraction complete.\n")

dT  <- Tmax - Tmin
dT[dT < 0] <- 0
PET <- 0.0023 * (GlobR / 1000) * sqrt(dT) * (Tavg + 17.8)
PET[PET < 0] <- 0

stopifnot(identical(dim(PET), dim(Tavg)))

time_df <- data.frame(date  = time_dates,
                      year  = as.integer(format(time_dates, "%Y")),
                      month = as.integer(format(time_dates, "%m")))

cat("Aggregating to monthly...\n")
climate_monthly <- do.call(rbind, lapply(seq_along(calib_ids), function(i) {
  pid       <- calib_ids[i]
  monthly_T <- aggregate(cbind(temp_air = Tavg[, i]),
                         by = list(year = time_df$year, month = time_df$month),
                         FUN = mean, na.rm = TRUE)
  monthly_T$precip  <- aggregate(Prec[, i], by = list(time_df$year, time_df$month),
                                 FUN = sum, na.rm = TRUE)$x
  monthly_T$evap    <- aggregate(PET[, i],  by = list(time_df$year, time_df$month),
                                 FUN = sum, na.rm = TRUE)$x
  monthly_T$plot_id <- pid
  monthly_T[, c("plot_id", "year", "month", "temp_air", "precip", "evap")]
}))

n_before <- nrow(climate_monthly)
climate_monthly <- climate_monthly[!duplicated(
  climate_monthly[, c("plot_id", "year", "month")]), ]
if (n_before != nrow(climate_monthly))
  message(sprintf("Removed %d duplicate rows from climate_monthly",
                  n_before - nrow(climate_monthly)))

cat("Monthly climate rows:", nrow(climate_monthly), "\n")
write.csv(climate_monthly, "./Data/model_inputs/climate_monthly.csv", row.names = FALSE)

input_monthly <- input_monthly[input_monthly$plot_id %in% calib_ids, ]
cat(sprintf("Filtered input_monthly to %d plots with climate (%d rows)\n",
            length(unique(input_monthly$plot_id)), nrow(input_monthly)))

input_monthly$temp_air <- NULL
input_monthly$precip   <- NULL
input_monthly$evap     <- NULL

input_raw_monthly <- merge(input_monthly,
                           climate_monthly[, c("plot_id", "year", "month",
                                               "temp_air", "precip", "evap")],
                           by = c("plot_id", "year", "month"), all.x = TRUE)

n_before <- nrow(input_raw_monthly)
input_raw_monthly <- input_raw_monthly[!duplicated(
  input_raw_monthly[, c("plot_id", "year", "month")]), ]
if (n_before != nrow(input_raw_monthly))
  message(sprintf("Removed %d duplicate rows from input_raw_monthly",
                  n_before - nrow(input_raw_monthly)))

input_raw_monthly <- input_raw_monthly[, c(
  "plot_id", "year", "month",
  "C_nwl_A", "C_nwl_W", "C_nwl_E", "C_nwl_N",
  "C_fwl_A", "C_fwl_W", "C_fwl_E", "C_fwl_N",
  "C_cwl_A", "C_cwl_W", "C_cwl_E", "C_cwl_N",
  "temp_air", "precip", "evap", "soc_obs_tCha", "soc_obs_tCha_sum",
  "soc_obs_year")]
input_raw_monthly <- input_raw_monthly[order(input_raw_monthly$plot_id,
                                             input_raw_monthly$year,
                                             input_raw_monthly$month), ]
rownames(input_raw_monthly) <- NULL

cat("\n--- Final monthly input table ---\n")
cat("Dimensions:    ", nrow(input_raw_monthly), "×", ncol(input_raw_monthly), "\n")
cat("Unique plots:  ", length(unique(input_raw_monthly$plot_id)), "\n")
cat("NAs temp_air:  ", sum(is.na(input_raw_monthly$temp_air)), "\n")
cat("SOC obs count (1m extrapolated): ", sum(!is.na(input_raw_monthly$soc_obs_tCha)), "\n")
cat("SOC obs count (raw sum):         ", sum(!is.na(input_raw_monthly$soc_obs_tCha_sum)), "\n")

write.csv(input_raw_monthly, "./Data/model_inputs/input_raw_monthly.csv", row.names = FALSE)


# ========================================================================================
# 4.1 Climate summary
# ========================================================================================

climate_annual <- aggregate(cbind(temp_air, precip) ~ plot_id + year,
                            data = climate_monthly, FUN = mean, na.rm = TRUE)
precip_annual <- aggregate(precip ~ plot_id + year, data = climate_monthly,
                           FUN = sum, na.rm = TRUE)
climate_annual$precip <- precip_annual$precip
temp_amp <- aggregate(temp_air ~ plot_id + year, data = climate_monthly,
                      FUN = function(x) (max(x) - min(x)) / 2)
climate_annual$temp_amplitude <- temp_amp$temp_air

cat("\n=== ANNUAL CLIMATE SUMMARY ===\n")
cat("Mean T (°C):\n"); print(summary(climate_annual$temp_air))
cat("Annual P (mm):\n"); print(summary(climate_annual$precip))
cat("T amplitude (°C):\n"); print(summary(climate_annual$temp_amplitude))


# ========================================================================================
# 4.2 Site-level derived covariates: climate metrics, Köppen class, litter quality
# ========================================================================================

cat("\n=== COMPUTING SITE-LEVEL COVARIATES ===\n")

# --- Monthly climate normals per plot (averaged across years) ---
T_norm <- aggregate(temp_air ~ plot_id + month, data = climate_monthly, FUN = mean)
P_norm <- aggregate(precip   ~ plot_id + month, data = climate_monthly, FUN = mean)
E_norm <- aggregate(evap     ~ plot_id + month, data = climate_monthly, FUN = mean)

T_wide <- reshape(T_norm, idvar = "plot_id", timevar = "month", direction = "wide")
P_wide <- reshape(P_norm, idvar = "plot_id", timevar = "month", direction = "wide")
E_wide <- reshape(E_norm, idvar = "plot_id", timevar = "month", direction = "wide")

T_wide <- T_wide[, c("plot_id", paste0("temp_air.", 1:12))]
P_wide <- P_wide[, c("plot_id", paste0("precip.",   1:12))]
E_wide <- E_wide[, c("plot_id", paste0("evap.",     1:12))]

# --- Climate variability metrics ---
days_per_month <- c(31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

climate_metrics_plot <- data.frame(
  plot_id            = T_wide$plot_id,
  mean_temp          = rowMeans(T_wide[, -1]),
  mean_precip_annual = rowSums(P_wide[, -1]),
  PET_annual         = rowSums(E_wide[, -1]),
  coldest_month_T    = apply(T_wide[, -1], 1, min),
  warmest_month_T    = apply(T_wide[, -1], 1, max),
  T_seasonality      = apply(T_wide[, -1], 1, sd),
  GDD5               = apply(T_wide[, -1], 1, function(t)
    sum(pmax(t - 5, 0) * days_per_month))
)
climate_metrics_plot$P_seasonality <- apply(P_wide[, -1], 1, function(p)
  sd(p) / mean(p))
climate_metrics_plot$aridity_index <- with(climate_metrics_plot,
                                           mean_precip_annual / PET_annual)

cat("\nClimate metrics summary:\n")
print(summary(climate_metrics_plot[, c("mean_temp", "GDD5", "coldest_month_T",
                                       "T_seasonality", "aridity_index")]))


# --- Köppen-Geiger classification (Finnish-relevant subset) ---
# Beck et al. (2018) simplified: ignores arid B and dry-season s/w subclasses.
# Expected codes for Finland: Dfb, Dfc, ET (rare Dfa/Cfb in southernmost plots).
classify_koppen <- function(monthly_T, monthly_P) {
  if (any(is.na(monthly_T)) || any(is.na(monthly_P))) return(NA_character_)
  if (length(monthly_T) != 12)                         return(NA_character_)
  T_cold <- min(monthly_T);  T_hot <- max(monthly_T);  n_warm <- sum(monthly_T >= 10)
  if (T_hot < 10)   return(if (T_hot < 0) "EF" else "ET")
  if (T_cold < -3)  {
    if (T_hot >= 22) return("Dfa")
    if (n_warm >= 4) return("Dfb")
    return("Dfc")
  }
  if (T_cold < 18) {
    if (T_hot >= 22) return("Cfa")
    if (n_warm >= 4) return("Cfb")
    return("Cfc")
  }
  NA_character_
}

climate_metrics_plot$koppen_class <- sapply(seq_len(nrow(T_wide)), function(i) {
  classify_koppen(as.numeric(T_wide[i, -1]), as.numeric(P_wide[i, -1]))
})

cat("\nKöppen-Geiger class distribution:\n")
print(table(climate_metrics_plot$koppen_class, useNA = "always"))


# --- Litter quality covariates per plot ---
litter_q_long <- aggregate(Cha ~ plot_id + size + compound,
                           data = inputs_matched[!is.na(inputs_matched$size), ],
                           FUN = sum)
litter_q_long$compound <- toupper(litter_q_long$compound)

plot_total <- aggregate(Cha ~ plot_id, data = litter_q_long, FUN = sum)
names(plot_total)[2] <- "litter_total"

awen_totals <- aggregate(Cha ~ plot_id + compound, data = litter_q_long, FUN = sum)
awen_wide   <- reshape(awen_totals, idvar = "plot_id",
                       timevar = "compound", direction = "wide")
names(awen_wide) <- sub("^Cha\\.", "litter_", names(awen_wide))

size_totals <- aggregate(Cha ~ plot_id + size, data = litter_q_long, FUN = sum)
size_wide   <- reshape(size_totals, idvar = "plot_id",
                       timevar = "size", direction = "wide")
names(size_wide) <- sub("^Cha\\.", "size_", names(size_wide))

# size_cwl (stump litter) is present directly in size_wide from the reshape.
# A previous placeholder (0.30 × fwl) has been removed.

litter_quality_plot <- merge(plot_total, awen_wide, by = "plot_id")
litter_quality_plot <- merge(litter_quality_plot, size_wide, by = "plot_id")

litter_quality_plot$litter_A_frac <- litter_quality_plot$litter_A / litter_quality_plot$litter_total
litter_quality_plot$litter_W_frac <- litter_quality_plot$litter_W / litter_quality_plot$litter_total
litter_quality_plot$litter_E_frac <- litter_quality_plot$litter_E / litter_quality_plot$litter_total
litter_quality_plot$litter_N_frac <- litter_quality_plot$litter_N / litter_quality_plot$litter_total

total_size <- with(litter_quality_plot, size_nwl + size_fwl + size_cwl)
litter_quality_plot$woody_share <- with(litter_quality_plot,
                                        (size_fwl + size_cwl) / total_size)

litter_quality_plot <- litter_quality_plot[, c(
  "plot_id", "litter_A_frac", "litter_W_frac", "litter_E_frac", "litter_N_frac",
  "woody_share")]

cat("\nLitter quality covariates (medians across plots):\n")
cat(sprintf("  AWEN fractions: A %.3f / W %.3f / E %.3f / N %.3f\n",
            median(litter_quality_plot$litter_A_frac, na.rm = TRUE),
            median(litter_quality_plot$litter_W_frac, na.rm = TRUE),
            median(litter_quality_plot$litter_E_frac, na.rm = TRUE),
            median(litter_quality_plot$litter_N_frac, na.rm = TRUE)))
cat(sprintf("  Woody share: %.3f [%.3f, %.3f]\n",
            median(litter_quality_plot$woody_share, na.rm = TRUE),
            min(litter_quality_plot$woody_share, na.rm = TRUE),
            max(litter_quality_plot$woody_share, na.rm = TRUE)))


# ========================================================================================
# 5. Build comprehensive site_raw
# ========================================================================================

site_raw <- data.frame(
  plot_id       = plot_data$plot_id,
  region        = plot_data$region,
  species_code  = plot_data$species,
  species_name  = plot_data$species_name,
  soil_code     = plot_data$CODE,
  soil_type     = plot_data$TEKSTI,
  shallow       = plot_data$shallow,
  n_soc_obs     = plot_data$n_soc_obs,
  mean_litter   = plot_data$Cha,
  mean_soc_Mgha = plot_data$C_Mgha,
  litter_cv     = plot_data$litter_cv
)

sp_cols  <- grep("^sp_frac_", names(plot_data), value = TRUE)
site_raw <- cbind(site_raw, plot_data[, sp_cols])

site_raw <- merge(site_raw,
                  coords_lookup[, c("plot_id", "x_ETRS", "y_ETRS",
                                    "lon_WGS84", "lat_WGS84")],
                  by = "plot_id", all.x = TRUE)

# Conifer share (Scots pine + Norway spruce)
sp_present <- intersect(c("sp_frac_1", "sp_frac_2"), names(site_raw))
site_raw$conifer_share <- if (length(sp_present) == 2) {
  site_raw$sp_frac_1 + site_raw$sp_frac_2
} else if (length(sp_present) == 1) {
  site_raw[[sp_present]]
} else {
  NA_real_
}

# Climate metrics + Köppen
site_raw <- merge(site_raw, climate_metrics_plot, by = "plot_id", all.x = TRUE)

# Temperature zone tertiles
temp_breaks      <- quantile(site_raw$mean_temp, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
site_raw$temp_zone <- cut(site_raw$mean_temp, breaks = temp_breaks,
                          labels = c("cold", "intermediate", "warm"),
                          include.lowest = TRUE)

# soil_depth from GTK proxy (coarse; profile_depth_cm from Biosoil is the measured value)
site_raw$soil_depth <- ifelse(site_raw$soil_code == "Tu", 100, 30)

# Merge all Biosoil-derived properties
site_raw <- merge(site_raw, biosoil_plot, by = "plot_id", all.x = TRUE)

# Clay: ClayContent only (clay_gtk fallback dropped — correlation 0.28 with Biosoil,
# median 10% vs 2.6%, fallback unreliable)
site_raw$clay <- site_raw$ClayContent

cat(sprintf("Clay: %d Biosoil measured, %d missing\n",
            sum(!is.na(site_raw$ClayContent)),
            sum(is.na(site_raw$clay))))

# --- Derived soil chemistry ---
site_raw$CN_ratio <- with(site_raw, ifelse(TotalNitrogen > 0,
                                           OrganicCarbon / TotalNitrogen,
                                           NA_real_))

exch_cols_for_cec <- c("ExchangeableCa", "ExchangeableMg", "ExchangeableK",
                       "ExchangeableNa", "ExchangeableAl", "FreeHAcidity")
site_raw$CEC <- rowSums(site_raw[, exch_cols_for_cec], na.rm = FALSE)
site_raw$base_saturation <- with(site_raw,
                                 (ExchangeableCa + ExchangeableMg + ExchangeableK + ExchangeableNa) / CEC)

cat(sprintf("Derived chemistry: CN_ratio median %.1f, CEC median %.2f, BS median %.2f\n",
            median(site_raw$CN_ratio,        na.rm = TRUE),
            median(site_raw$CEC,             na.rm = TRUE),
            median(site_raw$base_saturation, na.rm = TRUE)))

# --- Litter quality ---
site_raw <- merge(site_raw, litter_quality_plot, by = "plot_id", all.x = TRUE)

# --- NFI kuvio (stand-level) attributes ---
site_raw <- merge(site_raw, kuvio_summary, by = "plot_id", all.x = TRUE)

cat(sprintf("Stand age (1985) coverage in site_raw: %d non-NA of %d\n",
            sum(!is.na(site_raw$stand_age_85)), nrow(site_raw)))

# --- Peatland flag: KA 11-13 ---
site_raw$peatland <- !is.na(site_raw$KA) & site_raw$KA %in% c(11L, 12L, 13L)
cat(sprintf("Peatland plots (KA 11-13): %d\n", sum(site_raw$peatland, na.rm = TRUE)))

# --- Organic layer quality flags (derived in section 1.1) ---
# organic_missing: OFH never recorded in any year (structural gap)
# organic_zero:    OFH recorded but C_kgha == 0 in at least one year
# Both are mutually exclusive (verified); both excluded from calib_ready.
site_raw$organic_missing <- site_raw$plot_id %in% organic_missing_plots
site_raw$organic_zero    <- site_raw$plot_id %in% zero_organic_plots

cat(sprintf("Organic layer excluded — consistently absent: %d, zero C: %d\n",
            sum(site_raw$organic_missing, na.rm = TRUE),
            sum(site_raw$organic_zero,    na.rm = TRUE)))

# --- Implied MRT flag ---
site_raw$high_mrt <- site_raw$plot_id %in% high_mrt_plots
cat(sprintf("High implied MRT (>%dy or Inf) excluded: %d\n",
            MRT_EXCL_THRESHOLD, sum(site_raw$high_mrt, na.rm = TRUE)))

# --- calib_ready ---
plots_with_climate   <- unique(climate_monthly$plot_id[!is.na(climate_monthly$temp_air)])
plots_with_soc       <- unique(SOC_agg$plot_id)
site_raw$has_climate <- site_raw$plot_id %in% plots_with_climate
site_raw$has_soc     <- site_raw$plot_id %in% plots_with_soc

# --- SOC outlier flag (from the homogenized baseline; see section 1.3b) ---
# EXCLUSION RULE, applied and documented here rather than silently upstream:
#   a plot is excluded if its whole-profile stock exceeds 250 Mg/ha in ANY campaign year.
# Rationale: such stocks are not physically credible for a Finnish forest MINERAL soil
# (peatlands are already removed separately), and the affected plots read implausibly in
# every campaign or jump implausibly between campaigns -- i.e. they look like source-data
# errors, not real profiles. The flag itself is set in build_soc_homogenized.R; the
# DECISION to drop these plots from calibration is made here.
# Companion flag high_change (|consecutive rate| > 3 tC/ha/yr) is deliberately NOT an
# exclusion -- those plots are retained as genuine resampling noise.
site_raw$soc_outlier <- site_raw$plot_id %in% soc_outlier_plots
site_raw$high_change <- site_raw$plot_id %in% high_change_plots
cat(sprintf("SOC outlier (profile > 250 Mg/ha) excluded: %d  [plot_ids: %s]\n",
            sum(site_raw$soc_outlier, na.rm = TRUE),
            paste(sort(site_raw$plot_id[site_raw$soc_outlier]), collapse = ", ")))
cat(sprintf("High SOC change flagged but RETAINED: %d\n",
            sum(site_raw$high_change, na.rm = TRUE)))

# --- Constant-litter flag (see section 2) ---
# Litter identical to 4 dp for 1986-2024 => not a real series; see the block where
# const_litter_plots is derived. Excluded from calibration.
site_raw$const_litter <- site_raw$plot_id %in% const_litter_plots
cat(sprintf("Constant-litter plots excluded: %d\n", sum(site_raw$const_litter, na.rm = TRUE)))

site_raw$calib_ready <- site_raw$has_climate     &
  site_raw$has_soc          &
  !site_raw$peatland        &
  !site_raw$organic_missing &
  !site_raw$organic_zero    &
  !site_raw$high_mrt        &
  !site_raw$soc_outlier        &
  !site_raw$const_litter

cat("\n=== SITE_RAW SUMMARY ===\n")
cat("Total plots:               ", nrow(site_raw), "\n")
cat("With climate:              ", sum(site_raw$has_climate), "\n")
cat("With SOC:                  ", sum(site_raw$has_soc), "\n")
cat("Peatland excluded:         ", sum(site_raw$peatland,        na.rm = TRUE), "\n")
cat("Organic absent excluded:   ", sum(site_raw$organic_missing, na.rm = TRUE), "\n")
cat("Organic zero-C excluded:   ", sum(site_raw$organic_zero,    na.rm = TRUE), "\n")
cat("High implied MRT excluded: ", sum(site_raw$high_mrt,        na.rm = TRUE), "\n")
cat("SOC outlier excluded:      ", sum(site_raw$soc_outlier,     na.rm = TRUE), "\n")
cat("Calibration-ready:         ", sum(site_raw$calib_ready), "\n")
cat("Shallow:                   ", sum(site_raw$shallow), "\n")

cat("\nKA distribution (calib-ready plots):\n")
print(table(site_raw$KA[site_raw$calib_ready], useNA = "always"))
cat("kasvyo_syke distribution (calib-ready plots):\n")
print(table(site_raw$kasvyo_syke[site_raw$calib_ready], useNA = "always"))
cat("Köppen distribution (calib-ready plots):\n")
print(table(site_raw$koppen_class[site_raw$calib_ready], useNA = "always"))
cat("profile_depth_cm distribution (calib-ready plots):\n")
print(table(site_raw$profile_depth_cm[site_raw$calib_ready], useNA = "always"))


# =============================================================================
# HOLDOUT SPLIT (VALIDATION)
# =============================================================================
# Draws 20% of calib-ready plots once, at the data level.
# Stored as is_holdout column in site_raw.csv so all downstream scripts
# read the same fixed split with no seed dependency at runtime.

set.seed(42L)
calib_ids        <- as.character(site_raw$plot_id[site_raw$calib_ready])
n_holdout        <- floor(0.20 * length(calib_ids))
holdout_ids      <- sort(sample(calib_ids, n_holdout))
site_raw$is_holdout <- as.character(site_raw$plot_id) %in% holdout_ids &
  site_raw$calib_ready

cat(sprintf("\nHoldout split (20%% of calib-ready):\n"))
cat(sprintf("  Calibration plots: %d\n",
            sum(site_raw$calib_ready & !site_raw$is_holdout)))
cat(sprintf("  Holdout plots:     %d\n", sum(site_raw$is_holdout)))
cat(sprintf("  Total calib-ready: %d\n", sum(site_raw$calib_ready)))

write.csv(site_raw, "./Data/model_inputs/site_raw.csv", row.names = FALSE)

# ========================================================================================
# 6. Final sanity checks
# ========================================================================================

cat("\n=================================================================\n")
cat("  FINAL SANITY CHECKS\n")
cat("=================================================================\n")

chk_input   <- read.csv("./Data/model_inputs/input_raw_monthly.csv")
chk_climate <- read.csv("./Data/model_inputs/climate_monthly.csv")
chk_site    <- read.csv("./Data/model_inputs/site_raw.csv")

test_result <- function(label, pass, extra = "") {
  status <- if (pass) "PASS" else "FAIL"
  cat(sprintf("  [%s] %s%s\n", status, label,
              if (nchar(extra)) paste0(" — ", extra) else ""))
  pass
}

all_ok <- TRUE

# 1. No duplicate rows
dups   <- sum(duplicated(chk_input[, c("plot_id", "year", "month")]))
all_ok <- test_result("Input: no duplicate plot-year-month rows",
                      dups == 0, sprintf("%d duplicates", dups)) && all_ok

# 2. Complete climate
na_temp <- sum(is.na(chk_input$temp_air))
na_prec <- sum(is.na(chk_input$precip))
na_evap <- sum(is.na(chk_input$evap))
all_ok  <- test_result("Input: no NA in temp_air", na_temp == 0,
                       sprintf("%d NAs", na_temp)) && all_ok
all_ok  <- test_result("Input: no NA in precip",   na_prec == 0,
                       sprintf("%d NAs", na_prec)) && all_ok
all_ok  <- test_result("Input: no NA in evap",     na_evap == 0,
                       sprintf("%d NAs", na_evap)) && all_ok

# 3. Month counts per plot
rows_per_plot <- table(chk_input$plot_id)
all_ok <- test_result("Input: all plots have same month count",
                      length(unique(rows_per_plot)) == 1,
                      sprintf("range %d–%d",
                              min(rows_per_plot), max(rows_per_plot))) && all_ok

# 4. Non-negative litter
lit_cols <- grep("^C_", names(chk_input), value = TRUE)
neg_lit  <- sum(chk_input[, lit_cols] < 0, na.rm = TRUE)
all_ok   <- test_result("Input: no negative litter values",
                        neg_lit == 0, sprintf("%d negatives", neg_lit)) && all_ok

# 5. SOC obs only in June
soc_in_other_month <- sum(!is.na(chk_input$soc_obs_tCha) & chk_input$month != 6)
all_ok <- test_result("Input: SOC obs only in June",
                      soc_in_other_month == 0,
                      sprintf("%d obs outside June", soc_in_other_month)) && all_ok

# 6. SOC count matches expectation
soc_count_saved    <- sum(!is.na(chk_input$soc_obs_tCha))
soc_count_expected <- sum(SOC_agg$plot_id %in% unique(chk_input$plot_id))
all_ok <- test_result("Input: SOC obs count matches SOC_agg",
                      soc_count_saved == soc_count_expected,
                      sprintf("saved %d, expected %d",
                              soc_count_saved, soc_count_expected)) && all_ok

# 6b. soc_obs_tCha_sum present and matches soc_obs_tCha count
sum_count <- sum(!is.na(chk_input$soc_obs_tCha_sum))
all_ok <- test_result("Input: soc_obs_tCha_sum count matches soc_obs_tCha",
                      sum_count == soc_count_saved,
                      sprintf("sum count %d vs extrap count %d",
                              sum_count, soc_count_saved)) && all_ok

# 6c. Litter magnitude in expected boreal range (post-x10-removal)
ann_lit <- aggregate(rowSums(chk_input[, lit_cols]) ~ chk_input$plot_id + chk_input$year,
                     FUN = sum)
names(ann_lit)     <- c("plot_id", "year", "annual_litter")
plot_mean_lit      <- aggregate(annual_litter ~ plot_id, data = ann_lit, FUN = mean)
median_lit         <- median(plot_mean_lit$annual_litter)
all_ok <- test_result("Litter magnitude: median in boreal range 1.5-4.5 tC/ha/yr",
                      median_lit >= 1.0 && median_lit <= 5.5,
                      sprintf("median %.2f tC/ha/yr", median_lit)) && all_ok

# 7. Input plots all in site_raw
plots_input        <- unique(chk_input$plot_id)
plots_missing_site <- setdiff(plots_input, chk_site$plot_id)
all_ok <- test_result("Input plots all present in site_raw",
                      length(plots_missing_site) == 0,
                      sprintf("%d plots missing", length(plots_missing_site))) && all_ok

# 8. calib_ready plots in input_raw_monthly
calib_ready_ids     <- chk_site$plot_id[chk_site$calib_ready]
plots_missing_input <- setdiff(calib_ready_ids, plots_input)
all_ok <- test_result("All calib_ready plots present in input",
                      length(plots_missing_input) == 0,
                      sprintf("%d plots missing from input",
                              length(plots_missing_input))) && all_ok

# 9. Climate dedup
clim_dups <- sum(duplicated(chk_climate[, c("plot_id", "year", "month")]))
all_ok    <- test_result("Climate: no duplicate plot-year-month rows",
                         clim_dups == 0, sprintf("%d duplicates", clim_dups)) && all_ok

# 10. Climate plausibility
t_range <- range(chk_climate$temp_air, na.rm = TRUE)
p_range <- range(chk_climate$precip,   na.rm = TRUE)
all_ok  <- test_result("Climate: T within Finnish range (-40 to +30 °C)",
                       t_range[1] > -40 && t_range[2] < 30,
                       sprintf("actual: %.1f to %.1f",
                               t_range[1], t_range[2])) && all_ok
all_ok  <- test_result("Climate: P non-negative",
                       p_range[1] >= 0,
                       sprintf("min: %.1f", p_range[1])) && all_ok

# 11. site_raw completeness for calib_ready
site_ready <- chk_site[chk_site$calib_ready, ]
all_ok <- test_result("Site: region complete for calib_ready",
                      sum(is.na(site_ready$region)) == 0) && all_ok
all_ok <- test_result("Site: species complete for calib_ready",
                      sum(is.na(site_ready$species_code)) == 0) && all_ok
all_ok <- test_result("Site: soil_code complete for calib_ready",
                      sum(is.na(site_ready$soil_code)) == 0) && all_ok
all_ok <- test_result("Site: temp_zone complete for calib_ready",
                      sum(is.na(site_ready$temp_zone)) == 0) && all_ok

# 12. SOC plausibility (whole-profile target) — calib_ready plots only
# Bounds unchanged from the pre-baseline pipeline (5–400 tC/ha) and still comfortably met:
# observed calib-ready range is 11.4–244.8. Note the depth cap now lets genuinely thin
# soils (augering refusal at 10–20 cm) hold small honest stocks — the baseline minimum
# across ALL plots is 2.1 tC/ha — but no such plot is calibration-ready, so the 5 tC/ha
# floor stays meaningful rather than slack. The 400 ceiling is slack by design: the
# soc_outlier rule already removes anything above 250.
calib_ready_ids <- chk_site$plot_id[chk_site$calib_ready]
soc_obs_calib <- chk_input$soc_obs_tCha[
  !is.na(chk_input$soc_obs_tCha) &
    chk_input$plot_id %in% calib_ready_ids
]
soc_range <- range(soc_obs_calib)
all_ok <- test_result("SOC: calib-ready values within plausible range (5–400 tC/ha)",
                      soc_range[1] > 5 && soc_range[2] < 400,
                      sprintf("actual: %.1f to %.1f",
                              soc_range[1], soc_range[2])) && all_ok

# 12b. Whole-profile target >= measured 0-40 sum (direction check; deep layer only adds C)
soc_sum_vals <- chk_input$soc_obs_tCha_sum[!is.na(chk_input$soc_obs_tCha_sum)]
all_ok <- test_result("SOC: whole-profile median >= measured 0-40 sum median",
                      median(soc_obs_calib) >= median(soc_sum_vals),
                      sprintf("profile %.1f vs 0-40 sum %.1f tC/ha",
                              median(soc_obs_calib), median(soc_sum_vals))) && all_ok

# 13. No peatland plots in calib_ready
peat_in_calib <- sum(site_ready$peatland, na.rm = TRUE)
all_ok <- test_result("Site: no peatland plots (KA 11-13) in calib_ready",
                      peat_in_calib == 0,
                      sprintf("%d peatland plots still in calib_ready",
                              peat_in_calib)) && all_ok

# 13b. No organic-excluded plots in calib_ready
org_in_calib <- sum(site_ready$organic_missing | site_ready$organic_zero, na.rm = TRUE)
all_ok <- test_result("Site: no organic-excluded plots in calib_ready",
                      org_in_calib == 0,
                      sprintf("%d organic-excluded plots still in calib_ready",
                              org_in_calib)) && all_ok

# 13c. No high-MRT plots in calib_ready
high_mrt_in_calib <- sum(site_ready$high_mrt, na.rm = TRUE)
all_ok <- test_result("Site: no high-MRT plots in calib_ready",
                      high_mrt_in_calib == 0,
                      sprintf("%d high-MRT plots still in calib_ready",
                              high_mrt_in_calib)) && all_ok

# 13d. No SOC-outlier plots in calib_ready
soc_out_in_calib <- sum(site_ready$soc_outlier, na.rm = TRUE)
all_ok <- test_result("Site: no SOC-outlier plots (>250 Mg/ha) in calib_ready",
                      soc_out_in_calib == 0,
                      sprintf("%d SOC-outlier plots still in calib_ready",
                              soc_out_in_calib)) && all_ok

# 13e. Every calib_ready plot-year got a baseline profile value (no silent raw-sum fallback)
fallback_in_calib <- sum(!SOC_agg$soc_extrap_ok &
                           SOC_agg$plot_id %in% calib_ready_ids, na.rm = TRUE)
all_ok <- test_result("SOC: no calib-ready plot-year fell back to the raw 0-40 sum",
                      fallback_in_calib == 0,
                      sprintf("%d plot-years without a baseline profile value",
                              fallback_in_calib)) && all_ok

# 14. Existing covariates present in site_raw
for (v in c("KA", "kasvyo_syke", "profile_depth_cm", "CoarseFragments")) {
  all_ok <- test_result(sprintf("Site: %s present", v),
                        v %in% names(chk_site)) && all_ok
}

# 15. Derived chemistry covariates present and within plausible ranges
for (v in c("CN_ratio", "CEC", "base_saturation")) {
  all_ok <- test_result(sprintf("Site: %s present", v),
                        v %in% names(chk_site)) && all_ok
}
if (all(c("CN_ratio", "base_saturation") %in% names(chk_site))) {
  cn_med <- median(chk_site$CN_ratio[chk_site$calib_ready], na.rm = TRUE)
  all_ok <- test_result("Site: CN_ratio median in boreal range 15-30",
                        cn_med >= 15 && cn_med <= 30,
                        sprintf("median %.1f", cn_med)) && all_ok
  bs_med <- median(chk_site$base_saturation[chk_site$calib_ready], na.rm = TRUE)
  all_ok <- test_result("Site: base_saturation median in 0-1",
                        bs_med >= 0 && bs_med <= 1,
                        sprintf("median %.3f", bs_med)) && all_ok
}

# 16. Climate variability covariates present
for (v in c("GDD5", "coldest_month_T", "warmest_month_T",
            "T_seasonality", "P_seasonality", "aridity_index", "koppen_class")) {
  all_ok <- test_result(sprintf("Site: %s present", v),
                        v %in% names(chk_site)) && all_ok
}

# 17. Litter quality covariates present and AWEN fractions sum ~ 1
for (v in c("litter_A_frac", "litter_W_frac", "litter_E_frac", "litter_N_frac",
            "woody_share", "conifer_share")) {
  all_ok <- test_result(sprintf("Site: %s present", v),
                        v %in% names(chk_site)) && all_ok
}
if (all(c("litter_A_frac", "litter_W_frac", "litter_E_frac", "litter_N_frac")
        %in% names(chk_site))) {
  awen_sum <- with(chk_site[chk_site$calib_ready, ],
                   litter_A_frac + litter_W_frac + litter_E_frac + litter_N_frac)
  awen_ok  <- all(abs(awen_sum - 1) < 1e-6, na.rm = TRUE)
  all_ok   <- test_result("Site: AWEN fractions sum to 1 (per plot)",
                          awen_ok,
                          sprintf("range %.4f - %.4f",
                                  min(awen_sum, na.rm = TRUE),
                                  max(awen_sum, na.rm = TRUE))) && all_ok
}

# 18. clay_gtk should NO LONGER be present
all_ok <- test_result("Site: clay_gtk fallback removed",
                      !"clay_gtk" %in% names(chk_site)) && all_ok

# 19. NFI kuvio stand attributes present
kuvio_vars <- c("stand_age_85", "stand_age_2006_est", "basal_area_85",
                "mean_height_85_dm", "dev_class_85",
                "kasvup_tyyppi", "alaryhma", "ojitustilanne",
                "temp_sum_NFI", "elevation_m",
                "any_cut_85_95", "n_cuts_85_95", "any_trt_85_95",
                "soil_prep_pre85")
for (v in kuvio_vars) {
  all_ok <- test_result(sprintf("Site: %s present", v),
                        v %in% names(chk_site)) && all_ok
}

# 20. stand_age_85 plausibility
if ("stand_age_85" %in% names(chk_site)) {
  age_vals <- chk_site$stand_age_85[chk_site$calib_ready]
  age_vals <- age_vals[!is.na(age_vals)]
  if (length(age_vals) > 0) {
    age_range <- range(age_vals)
    all_ok <- test_result("stand_age_85: in plausible range (0-400 yr)",
                          age_range[1] >= 0 && age_range[2] <= 400,
                          sprintf("range %.0f-%.0f",
                                  age_range[1], age_range[2])) && all_ok
  }
}

# 21. n_cuts_85_95 in {0, 1, 2, 3}
if ("n_cuts_85_95" %in% names(chk_site)) {
  ncuts_vals <- chk_site$n_cuts_85_95[chk_site$calib_ready]
  all_ok <- test_result("n_cuts_85_95: values in 0:3",
                        all(ncuts_vals %in% 0:3, na.rm = TRUE)) && all_ok
}

# 22. Kuvio attribute coverage (informational)
if ("stand_age_85" %in% names(chk_site)) {
  age_na  <- sum(is.na(chk_site$stand_age_85[chk_site$calib_ready]))
  age_pct <- 100 * age_na / sum(chk_site$calib_ready)
  cat(sprintf("  [INFO] stand_age_85 NA in calib_ready: %d (%.1f%%)\n",
              age_na, age_pct))
}
if ("ojitustilanne" %in% names(chk_site)) {
  drainage_na  <- sum(is.na(chk_site$ojitustilanne[chk_site$calib_ready]) |
                        chk_site$ojitustilanne[chk_site$calib_ready] == "")
  drainage_pct <- 100 * drainage_na / sum(chk_site$calib_ready)
  cat(sprintf("  [INFO] ojitustilanne missing in calib_ready: %d (%.1f%%)\n",
              drainage_na, drainage_pct))
}

# 23. Bulk density coverage (informational)
if ("EstimatedBulkDensity" %in% names(chk_site)) {
  ebd_na  <- sum(is.na(chk_site$EstimatedBulkDensity[chk_site$calib_ready]))
  ebd_pct <- 100 * ebd_na / sum(chk_site$calib_ready)
  cat(sprintf("  [INFO] EstimatedBulkDensity NA in calib_ready: %d (%.1f%%)\n",
              ebd_na, ebd_pct))
}
if ("MeanBulkDensity" %in% names(chk_site)) {
  mbd_na  <- sum(is.na(chk_site$MeanBulkDensity[chk_site$calib_ready]))
  mbd_pct <- 100 * mbd_na / sum(chk_site$calib_ready)
  cat(sprintf(paste0("  [INFO] MeanBulkDensity NA in calib_ready: %d (%.1f%%)",
                     " -- prefer EstimatedBulkDensity for RF if coverage is better\n"),
              mbd_na, mbd_pct))
}

# --- Summary ---
cat(sprintf("\n  Plot counts:\n"))
cat(sprintf("    site_raw total:          %d\n", nrow(chk_site)))
cat(sprintf("    has_soc:                 %d\n", sum(chk_site$has_soc)))
cat(sprintf("    has_climate:             %d\n", sum(chk_site$has_climate)))
cat(sprintf("    peatland excluded:       %d\n", sum(chk_site$peatland,         na.rm = TRUE)))
cat(sprintf("    organic absent excluded: %d\n", sum(chk_site$organic_missing,  na.rm = TRUE)))
cat(sprintf("    organic zero-C excluded: %d\n", sum(chk_site$organic_zero,     na.rm = TRUE)))
cat(sprintf("    high MRT excluded:       %d\n", sum(chk_site$high_mrt,         na.rm = TRUE)))
cat(sprintf("    calib_ready:             %d\n", sum(chk_site$calib_ready)))
cat(sprintf("    in input_monthly:        %d\n", length(plots_input)))

cat(sprintf("\n  SOC observation years:\n"))
print(table(chk_input$year[!is.na(chk_input$soc_obs_tCha)]))

cat(sprintf("\n  SOC observations — 1m extrapolated (tC/ha):\n"))
cat(sprintf("    n:      %d\n",   length(soc_obs_calib)))
cat(sprintf("    min:    %.1f\n", min(soc_obs_calib)))
cat(sprintf("    median: %.1f\n", median(soc_obs_calib)))
cat(sprintf("    max:    %.1f\n", max(soc_obs_calib)))
cat(sprintf("    CV:     %.3f\n", sd(soc_obs_calib) / mean(soc_obs_calib)))

cat(sprintf("\n  SOC observations — raw sum 0-40cm + OFH (tC/ha):\n"))
cat(sprintf("    median: %.1f\n", median(soc_sum_vals)))
cat(sprintf("    max:    %.1f\n", max(soc_sum_vals)))

if (all_ok) {
  cat("\n  ✓ All checks passed. Data preparation complete.\n")
} else {
  cat("\n  ✗ One or more checks FAILED. Review output above before running calibration.\n")
}




### SOME PLOTS ARE BELONGING TO CLASSES (GTK) THAT SHOULD NOT BE THERE

# =============================================================================
# Investigation: Tu (peat) / Ve (water) / Ka (bedrock) plots in calib_ready
# =============================================================================

EXCLUDED_SOIL_CODES <- c("Tu", "Ve", "Ka")

# Subset: weird-class plots that are currently calib_ready
weird <- site_raw[site_raw$calib_ready &
                    site_raw$soil_code %in% EXCLUDED_SOIL_CODES, ]
cat(sprintf("Weird-class plots in calib_ready: %d (Tu=%d, Ve=%d, Ka=%d)\n",
            nrow(weird),
            sum(weird$soil_code == "Tu"),
            sum(weird$soil_code == "Ve"),
            sum(weird$soil_code == "Ka")))


# -----------------------------------------------------------------------------
# 1. Map: weird plots over Finland
# -----------------------------------------------------------------------------
# Finland boundary from rnaturalearth (already loaded in section 1.5).
# Background: all calib_ready mineral plots in light grey for context.
# Foreground: weird-class plots coloured by soil_code.

library(rnaturalearth)
finland <- ne_countries(country = "Finland", scale = "medium", returnclass = "sf")

# Plot coordinates already in site_raw (lon_WGS84, lat_WGS84)
mineral_calib <- site_raw[site_raw$calib_ready &
                            !site_raw$soil_code %in% EXCLUDED_SOIL_CODES, ]

class_cols <- c("Tu" = "#8c510a",   # peat — brown
                "Ve" = "#2166ac",   # water — blue
                "Ka" = "#999999")   # bedrock — grey

plot(st_geometry(finland), col = "white", border = "grey50",
     main = sprintf("Weird-class plots in calib_ready (n = %d)", nrow(weird)))
points(mineral_calib$lon_WGS84, mineral_calib$lat_WGS84,
       pch = 16, cex = 0.4, col = adjustcolor("grey70", 0.5))
points(weird$lon_WGS84, weird$lat_WGS84,
       pch = 16, cex = 1.0, col = class_cols[weird$soil_code])
legend("topright",
       legend = c(sprintf("Tu — peat (n=%d)",    sum(weird$soil_code == "Tu")),
                  sprintf("Ve — water (n=%d)",   sum(weird$soil_code == "Ve")),
                  sprintf("Ka — bedrock (n=%d)", sum(weird$soil_code == "Ka")),
                  "Other calib_ready (mineral)"),
       col = c(class_cols, adjustcolor("grey70", 0.5)),
       pch = 16, bty = "n", cex = 0.85)


# -----------------------------------------------------------------------------
# 2. Table of weird plots with SOC values
# -----------------------------------------------------------------------------
# One row per plot, with mean SOC across observation years.
# Columns chosen to give a quick read on plausibility: where they are,
# what GTK + Cajander say, how much SOC, OFH, and how many obs.

weird_table <- weird[, c("plot_id", "soil_code", "soil_type", "KA",
                         "kasvyo_syke", "lon_WGS84", "lat_WGS84",
                         "mean_soc_Mgha", "mean_litter",
                         "ofh_lower_cm", "ofh_weight_kgm2",
                         "n_soc_obs")]
# Sort by SOC descending — extremes first
weird_table <- weird_table[order(-weird_table$mean_soc_Mgha), ]
rownames(weird_table) <- NULL

cat("\n--- Weird-class plots, sorted by SOC ---\n")
print(weird_table, digits = 3)


# -----------------------------------------------------------------------------
# 3. Histogram of SOC in weird plots, vs mineral plots for context
# -----------------------------------------------------------------------------

soc_weird   <- weird$mean_soc_Mgha[is.finite(weird$mean_soc_Mgha)]
soc_mineral <- mineral_calib$mean_soc_Mgha[is.finite(mineral_calib$mean_soc_Mgha)]

# Common breaks so the two histograms are comparable
brks <- pretty(range(c(soc_weird, soc_mineral), na.rm = TRUE), n = 30)

par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))

hist(soc_mineral, breaks = brks, col = adjustcolor("grey60", 0.6), border = "white",
     main = sprintf("Calib_ready mineral plots (n = %d)", length(soc_mineral)),
     xlab = "Mean SOC (Mg C/ha)", freq = FALSE)
abline(v = median(soc_mineral), col = "tomato", lwd = 2, lty = 2)

hist(soc_weird, breaks = brks, col = adjustcolor("steelblue", 0.7), border = "white",
     main = sprintf("Weird-class plots Tu/Ve/Ka (n = %d)", length(soc_weird)),
     xlab = "Mean SOC (Mg C/ha)", freq = FALSE)
abline(v = median(soc_weird), col = "tomato", lwd = 2, lty = 2)
par(mfrow = c(1, 1))


# -----------------------------------------------------------------------------
# 4. EXTRA: SOC by soil_code, side-by-side boxplot
# -----------------------------------------------------------------------------
# Compares the weird-class SOC distributions against each other and against
# the bulk of mineral plots. Tells you whether Tu/Ve/Ka are systematically
# higher (consistent with peat/wetland C accumulation) or just noisier.

soc_by_class <- rbind(
  data.frame(soc = soc_mineral, class = "mineral (other)"),
  data.frame(soc = weird$mean_soc_Mgha[weird$soil_code == "Tu"], class = "Tu — peat"),
  data.frame(soc = weird$mean_soc_Mgha[weird$soil_code == "Ve"], class = "Ve — water"),
  data.frame(soc = weird$mean_soc_Mgha[weird$soil_code == "Ka"], class = "Ka — bedrock")
)
soc_by_class <- soc_by_class[is.finite(soc_by_class$soc), ]

boxplot(soc ~ class, data = soc_by_class,
        col = c(adjustcolor("grey70", 0.6), class_cols),
        main = "SOC distribution by soil class",
        xlab = "", ylab = "Mean SOC (Mg C/ha)",
        las = 1)
abline(h = median(soc_mineral), lty = 3, col = "grey40")


# -----------------------------------------------------------------------------
# 5. EXTRA: Cajander class (KA) cross-tab for the weird plots
# -----------------------------------------------------------------------------
# Catches the "what does the Cajander vegetation class say about these plots?"
# question. If Tu/Ve/Ka plots are in low-fertility KA classes (3-4 = VT, CT)
# they're probably real lean sites with thin/peaty soils. If they're in KA=1
# (OMT, rich fertility), the GTK and Cajander classifications disagree, which
# is itself diagnostic.

cat("\n--- Cajander class (KA) by weird soil_code ---\n")
print(table(weird$soil_code, weird$KA, useNA = "ifany"))

# And mean SOC by (soil_code × KA) to see if any combination is genuinely
# normal-looking and might be safe to keep
cat("\n--- Mean SOC (Mg C/ha) by (soil_code × KA), n = count ---\n")
agg <- aggregate(mean_soc_Mgha ~ soil_code + KA, data = weird,
                 FUN = function(x) c(n = length(x), median = median(x), max = max(x)))
print(agg)







# ========================================================================================
# METHODS & MATERIALS NOTES — Data Preparation
#
# This section assembles all decisions, sources, and filtering rationales applied
# in this script, organised for the M&M section of the manuscript. Numerical
# values (counts, ranges, medians) come from the most recent run of this script
# and should be re-extracted from the actual run output before submission.
# ========================================================================================
#
# DATA SOURCES
#
#   SOC observations
#     Finnish national soil C inventory, three waves: VMI8 (1985), Biosoil (2006)
#     and Komeetta (2024), all from one homogenized baseline.
#     File:   ./Data/SOC_homogeneized/soc_homogenized_layers.csv  (layer level)
#             ./Data/SOC_homogeneized/soc_homogenized_plot.csv    (plot level)
#     Built by ./Data/SOC_homogeneized/build_soc_homogenized.R from the LUKE
#     workbook Komeetta 150526hi--.xlsx (H. Ilvesniemi); stocks are bulk-density
#     and coarse-fragment corrected, and reproduce LUKE's official national
#     figures exactly (see "SOC CALIBRATION TARGET" below).
#     Layers measured: organic (OFH) + mineral 0-10 cm, 10-20 cm, 20-40 cm in the
#     2006/2024 protocol; 1985 uses 0-5 cm, 5-20 cm, 20-40 cm. Both are accepted.
#     Supersedes ./Data/SOC/soilC1985_2006.csv and the separate Komeetta extract
#     ./Data/Komeetta/Komeetta_mitatut hiilet.xlsx, neither of which is read now.
#
#   Litter inputs
#     Tree-derived litter from MUSTIKKA NFI permanent plots, allometric
#     equations applied to inventoried biomass.
#     File: ./Data/LitterData/tree_litter_per_site_year_awen_29.04.26.csv
#     Resolution: annual, per plot, per species, per AWEN compound, per size class
#     Coverage: 1985-2023 in the source; 2024 carried forward from 2023 (see below).
#     Note: a previous x10 multiplier (kg/m² → t/ha) was REMOVED from this script
#     after audit confirmed source values were already in t C/ha/yr (median 2.4
#     across calib_ready plots, within the boreal range 1.5-4.5 t C/ha/yr).
#
#     FIRST YEAR (1985) IS RECONSTRUCTED, NOT OBSERVED — see "LITTER INPUT SERIES:
#     THE RECONSTRUCTED FIRST YEAR" below. This is a reportable methods decision,
#     not an implementation detail.
#
#   Plot-key linkage
#     File: ./Data/soil_litter_site_key.csv
#     Maps three plot identifier conventions: koealatunnus_BIOSOIL (used
#     downstream), koealatunnus_MUSTIKKA (litter), koealatunnus_VANHA (NFI).
#
#   GTK soil classification
#     ArcGIS REST endpoint, Maapera 1:1M scale, layer 3
#     URL: https://gtkdata.gtk.fi/arcgis/rest/services/Hasu/maapera/MapServer/3
#     Used for descriptive soil_code field; NOT used as a calibration filter
#     (rationale below).
#
#   Biosoil mineral and organic layer properties
#     2006 European-wide soil monitoring campaign, Finnish subset
#     Files: ../../Datasets/BioSoil_maaperäaineisto_2006/{soil_data_main, plot_index}.csv
#     Mineral layers M01 (0-10 cm) + M12 (10-20 cm) averaged for topsoil
#     properties; OFH used for organic layer thickness and mass.
#     Provides Cajander fertility class (KA), SYKE biogeographic zones, soil
#     physical (texture, bulk density, coarse fragments) and chemical (pH,
#     exchangeable cations, total N, organic C) variables.
#
#   NFI permanent plot stand attributes
#     PysyvätKoealat 1985, 1990, 1995 measurement waves
#     Files: ./Data/PysyvätKoealat/{1985,1990,1995}/kuvio.csv
#     Stand age, basal area, dominant species, development class, drainage
#     status, treatment and cutting history, soil preparation history.
#     1985 wave provides full retrospective (cuts ≤10y, soil prep ≤30y);
#     1990 captures only new events since 1985; 1995 covers 1990-1995.
#     Combined coverage: ~1955-1995 for management indicators.
#
#   Climate
#     File: ./Data/Climate/nfi_plot_weather_data_1961_2025.nc
#     NetCDF with daily temperature (mean/min/max), precipitation, and global
#     radiation, 1961-2025, indexed by NFI plot_id.
#     PET computed via Hargreaves-Samani:
#       PET = 0.0023 × (GlobR/1000) × sqrt(Tmax-Tmin) × (Tavg+17.8)
#     Daily values aggregated to monthly: mean for temperature, sum for
#     precipitation and PET.
#
# ----------------------------------------------------------------------------------------
#
# LITTER INPUT SERIES: THE RECONSTRUCTED FIRST YEAR (1985)
#
#   (M&M-ready. Added 2026-08-04.)
#
#   What was found. The first year of the Tupek litter product is not usable. Source
#   median total litter is 0.060 t C/ha/yr in 1985 against 1.398 in 1986 — a factor of
#   23 — while the number of site-year records is identical in the two years (2805), so
#   this is not missing data: the values themselves collapse. The pattern is uniform
#   across months, plots and all three size classes.
#
#   Interpretation. This is the expected signature of litter estimated from between-
#   inventory biomass increments: the first year of such a series has no predecessor to
#   difference against and therefore carries no turnover signal. (To be confirmed with
#   the dataset author before this wording is used in the manuscript.)
#
#   Why it could not be left alone. 1985 is t0 — the year the transient pre-run ends and
#   the first SOC campaign is observed — so the artefact propagates into four places at
#   once. Three are litter anchors computed from the head of the series:
#     J_t0_mean    mean of the first 5 years   -> the 1917->1985 pre-run ENDPOINT
#     J_full_mean  mean of the whole series    -> the 1917 pre-run ANCHOR
#     J_total_mean mean of the first 20 years  -> feeds J_bar, the units bridge that
#                                                 converts the physical litter-flux
#                                                 window into the sigma_input multiplier
#   The fourth is the forward run's own first year, which coincides exactly with the
#   VMI8 observation. Uncorrected, J_t0_mean was 1.543 instead of 1.863 t C/ha/yr — the
#   spin-up terminated ~21% below the true 1985 litter level, biasing the initial carbon
#   state low in all six models and, through J_bar, shifting the physical flux bounds.
#
#   What was done. For each plot and each AWEN component, a linear trend is fitted over
#   1986-1990 and backcast one year to 1985 (negatives truncated at zero). Litter is
#   RISING through that window (population median 1.84 in 1986 to 2.17 in 1990), so a
#   plain 1986-1990 average would place 1985 ABOVE 1986 and contradict the growing-stock
#   history that the same pre-run uses for its input shape. The backcast instead puts
#   1985 just below 1986 (1.754 vs 1.838 t C/ha/yr), consistent with a still-accumulating
#   forest.
#
#   Robustness. The three defensible reconstructions agree closely — carry-back of 1986
#   gives J_t0_mean 1.890, the 1986-1990 mean gives 1.918, the backcast gives 1.863 —
#   a spread of ~3%, against the ~21% error being corrected. The result is therefore
#   insensitive to the choice among them; what matters is that the artefactual year is
#   not used. Dropping 1985 outright was considered and rejected: it would move t0 to
#   1986 and force all 441 VMI8 observations to be re-mapped, for the same numerical
#   answer (threshold statistic 0.826 reconstructed vs 0.857 dropped).
#
#   Reporting. 1985 litter is a RECONSTRUCTED value and should be described as such
#   wherever the litter history is presented.
#
# ----------------------------------------------------------------------------------------
#
# SOC CALIBRATION TARGET — HOMOGENIZED THREE-CAMPAIGN BASELINE
#
#   (M&M-ready. Revised 2026-08-04, replacing the in-script depth extrapolation
#   and the separate Komeetta ingest that preceded it.)
#
#   Source.
#     All three soil campaigns — VMI8 (1985), Biosoil (2006) and Komeetta (2024) —
#     are taken from a single LUKE workbook maintained by Hannu Ilvesniemi
#     (Komeetta 150526hi--.xlsx), rather than from campaign-specific files
#     processed along separate paths. The workbook carries LUKE's pre-computed
#     layer stocks, which apply both a bulk-density and a coarse-fragment
#     (stoniness) correction.
#
#   Why the source changed.
#     The previous SOC file (soilC1985_2006.csv, plus a separately-ingested
#     Komeetta 2024 extract) carried no coarse-fragment correction on the mineral
#     layers. Because stoniness is a volumetric property it inflates mineral
#     carbon multiplicatively while leaving the organic layer untouched — and that
#     is exactly the signature found: our stocks exceeded LUKE's official figures
#     by a near-constant factor of 1.60-1.66 across all three mineral layers in
#     both years, with an organic-layer ratio of 0.86. The median coarse-fragment
#     content of these plots is 43%, and 1/(1-0.43) = 1.76 accounts for the gap.
#     Two campaign-specific processing paths had also drifted apart: the Komeetta
#     ingest mislabelled its mineral layer codes as the 1985 protocol (201 = 0-5 cm)
#     when Komeetta follows the 2006 protocol (201 = 0-10 cm), corrupting the depth
#     model for that campaign alone.
#
#   Verification against the official inventory.
#     Applying LUKE's own cleaning rules to the workbook (one replicate per plot,
#     peat layers H01/H12 dropped, abandoned plots dropped, organic set to zero
#     where no organic-layer weight was recorded, litter folded into the organic
#     layer) reproduces the official national figures EXACTLY: organic + 0-40 cm,
#     region-weighted, 2006 = 59.1 and 2024 = 61.0 Mg/ha over 446 plots. The 1985
#     wave is held on the same corrected basis — the per-plot 2006/1985 stock ratio
#     is 1.17, whereas an uncorrected 1985 wave would give roughly 0.6.
#
#   Depth basis of the calibration target (soc_profile).
#     The models carry bulk SOC, so the target is a whole-profile stock:
#       organic layer (measured)
#         + mineral 0-40 cm (measured)
#         + mineral (deepest measured layer -> z_cap) (modelled)
#     The deep tail uses the same exponential depth model as before,
#     C(z) = C₀·exp(-λz), with one λ fitted per GTK soil class by pooling profiles
#     (C₀ concentrated out analytically, so the search is 1-D in λ) and C₀ then
#     re-solved per profile by no-intercept OLS against that profile's measured
#     layer stocks. Layer stocks are fitted in integral form,
#     Stock(z₁,z₂) = C₀/λ·(exp(-λz₁) - exp(-λz₂)), with no midpoint approximation.
#
#     The substantive change is the integration limit. The tail is integrated to
#     z_cap = min(100 cm, depth augering actually reached) rather than to a fixed
#     1 m. Depth reached is read from LAYER_LIMIT_INF: a value of 80 means the
#     target depth was reached and the profile is extended to 100 cm, whereas
#     10/20/40 records augering refusal and no deep extrapolation is applied.
#     Under the old fixed-1 m rule, 38 plots that hit refusal at 10-40 cm were each
#     assigned roughly 14 tC/ha of carbon below bedrock. The target is therefore a
#     per-plot VARIABLE depth that matches how each profile was actually measured,
#     rather than a fixed standardized depth.
#
#   Validation of the deep extrapolation.
#     Biosoil layer Krs 204 is a MEASURED 40-80 cm layer on 501 plots and is held
#     out of the fit, giving direct ground truth. Predicted 40-80 cm stocks track
#     it with median pred/obs = 0.96 and a bias of -2.5 Mg/ha: the extrapolation is
#     mildly conservative, and the concern that it might over-extrapolate rich
#     soils is not supported.
#
#   Exclusions applied to the SOC data.
#     Peat layers (CODE_LAYER H01/H12) are removed at source, independently of the
#     Cajander peatland filter applied later in this script. Plots whose 2024
#     profile is incomplete are dropped for 2024 rather than zero-filled (zero-
#     filling would deflate the 2024 mean); their 1985 and 2006 records are kept
#     and flagged. Plots with a whole-profile stock above 250 Mg/ha in any campaign
#     are excluded from calibration entirely (see calib_ready, criterion 8).
#
#   Region weighting.
#     Region is resolved to one value per plot through a documented provenance
#     chain (Biosoil design -> site_raw -> latitude), with the design value winning
#     the two conflicts found. Northern plots carry weight 3 and southern plots
#     weight 1, following the Biosoil sampling design; the same weights are assumed
#     for 1985.
#
#   Effect on the calibration target.
#     The target falls and, more importantly, becomes internally consistent across
#     campaigns. Median soc_obs_tCha by campaign is now 59.4 (1985), 66.7 (2006)
#     and 67.4 (2024) tC/ha — a gentle rise toward saturation. The previous
#     pipeline gave 63, 102 and 105: a 62% jump between 1985 and 2006 that was an
#     artefact of cross-campaign inconsistency, and a 2006 value contradicting the
#     official same-year figure by a factor of 1.7. The new-to-old ratio is ~1.0 in
#     1985 and ~0.69 in 2006/2024, i.e. the correction lands where the
#     inconsistency was.
#
#   Retained diagnostics.
#     soc_obs_tCha_sum retains the measured organic + 0-40 cm sum. It reproduces
#     the baseline's independently-computed 0-40 cm stock to 0.000 tC/ha, which is
#     the check that the layer table and the plot table describe the same data.
#     soc_extrap_ok records whether a baseline profile value was matched; no
#     calibration-ready plot-year falls back to the raw sum.
#
# ----------------------------------------------------------------------------------------
#
# CALIBRATION-READY PLOT FILTERING (calib_ready flag)
#
#   A plot is included in the calibration set if and only if all of the
#   following hold. The order below reflects the order of evaluation in the
#   script and the strength of the evidence supporting each criterion.
#
#   1. has_climate
#        Plot ID present in the climate netCDF with non-NA temperature
#        records. Excludes plots with no climate coverage.
#
#   2. has_soc
#        Plot ID present in SOC_agg with at least one non-zero observation.
#
#   3. NOT peatland (Cajander)
#        Cajander fertility class KA ∉ {11, 12, 13}. Cajander 11-13 are
#        peatland vegetation types (kasvyo); these plots are governed by
#        peatland C dynamics (waterlogging-driven accumulation) which Yasso
#        is not parameterised for.
#
#   4. NOT organic_missing
#        OFH layer never recorded in any year. Indicates a structural gap in
#        the organic layer record; a whole-profile stock cannot be assembled
#        without an OFH measurement.
#
#   5. NOT organic_zero
#        OFH recorded but C_kgha == 0 in at least one observation year.
#        Likely a measurement issue rather than absence; excluded out of
#        caution. Verified non-overlapping with organic_missing.
#
#   6. NOT zero-litter
#        Total litter input across all years and components > 0. Plots with
#        all-zero litter produce a steady-state SOC of zero (since
#        C* = -A⁻¹·b with b = 0), which propagates as -Inf in log-residuals
#        downstream. Caught after an earlier calibration run produced ~1% of
#        observations as +Inf residuals; added to the filter as a defensive
#        check.
#
#   7. NOT high_mrt  (implied MRT ≤ 100 years)
#        Implied MRT = mean(1m SOC) / mean(annual litter input). Plots with
#        MRT > 100y have litter inputs so low that the model cannot produce
#        positive SOC at any parameter value, causing -Inf likelihood at the
#        pre-MCMC sanity check. Threshold chosen from a clear natural gap in
#        the MRT distribution (see diagnostics_inputs.R plot 01_mrt_histogram).
#        Also catches four plots with MRT = Inf (zero litter after AWEN
#        mapping) that escape the raw zero-litter filter above, which operates
#        on inputs_matched before compound-fraction assignment.
#        Removed: 16 plots; several overlap with GTK Tu/Ve class, consistent
#        with peat-influenced C dynamics not supported by Yasso's mineral-soil
#        parameterisation.
#
#   8. NOT soc_outlier  (whole-profile stock ≤ 250 Mg/ha in every campaign)
#        A whole-profile stock above 250 Mg/ha is not credible for a Finnish
#        forest MINERAL soil, and peatlands are already removed by criterion 3.
#        The four affected plots read implausibly in every campaign or move
#        implausibly between campaigns — plot 33631 reads 290-417 Mg/ha in all
#        three waves, plot 31751 gains ~190 Mg/ha between 2006 and 2024 — which
#        is the signature of source-data error rather than of real profiles.
#        Removed: 4 plots (29232, 31751, 33631, 49571). The flag is computed in
#        build_soc_homogenized.R; the decision to exclude is taken in this script
#        so that it is visible alongside the other calibration filters.
#
#        Deliberately NOT an exclusion: the companion flag high_change
#        (|consecutive rate| > 3 tC/ha/yr, 23 plots) is carried into site_raw but
#        those plots are RETAINED. Large between-campaign changes at plot scale
#        are expected from resampling noise on a spatially heterogeneous soil and
#        are part of the signal the error model is meant to absorb; discarding
#        them would bias the calibration toward implausibly quiet plots.
#
##   9. NOT const_litter  (litter series is not a fixed repeated value)
#        Eight calibration-ready plots carry litter IDENTICAL to ~15 significant figures
#        for every year 1986-2024 (relative SD < 1e-6), constant in the source as well.
#        No stand produces identical litter for 39 consecutive years; this is the litter
#        model returning a fixed value where the underlying biomass series was static.
#        Same class as the artefactual 1985 first year -- plausible, positive, finite and
#        wrong -- and it survived every earlier check because nothing tested per-plot
#        constancy. Excluded rather than flagged: a constant series injects a spurious
#        "no temporal change" forcing into a study whose subject IS the temporal change.
#        Removed: 8 plots (64 exist dataset-wide; the rest already fail other filters).
#        Detected by rule, not hard-coded, so it survives a change in the litter product.
#
#        Borderline and deliberately NOT excluded: plots 39251 and 67631 take only three
#        distinct litter values across 39 years (relative SD ~4-7e-3). Suspiciously
#        quantised but not constant, so they fall outside the rule; noted here so the
#        decision is visible rather than implicit.
#
#   The full set of exclusions is reported in the section 6 sanity output
#   together with the final calib_ready count.
#
# ----------------------------------------------------------------------------------------
#
# GTK SOIL CLASS — REVIEWED, NOT USED AS A FILTER
#
#   GTK 1:1M soil classification was extracted via the ArcGIS REST API and
#   stored as soil_code. We considered using it to additionally exclude
#   plots classified as Tu (peat), Ve (water), or Ka (bedrock) — together
#   ~26% of calib_ready plots after the filters above. Inspection of these
#   "weird-class" plots showed:
#
#     - SOC distributions indistinguishable from mineral plots:
#       Tu median 59 Mg C/ha, Ve 59, Ka 70, mineral 62.
#       Both histograms span ~30-130 Mg C/ha with similar shape.
#
#     - Cajander vegetation surveys, conducted at plot scale, classified
#       these plots as mineral upland forest (KA = 1, OMT, rich mineral
#       upland) for 53/96 cases; none were classified as peatland (KA 11-13).
#
#     - The disagreement between the two classifications (GTK 1:1M soil map
#       vs Cajander field-surveyed vegetation) is consistent with the GTK
#       map's coarse scale: a 1:1,000,000 raster cell can easily misclassify
#       a 30 m forest plot as the surrounding water body or peatland that
#       dominates the kilometre-scale cell.
#
#   We therefore retained the GTK soil_code field for descriptive use and
#   for the per-class λ pooling in the deep-tail extrapolation (which now runs
#   upstream in build_soc_homogenized.R, where soil class is a useful covariate
#   for depth profile shape regardless of plot-level accuracy), but did NOT
#   exclude any GTK class from calib_ready.
#
#   Recommendation for the methods section: state that GTK class was reviewed
#   and not used as a calibration filter, with the reason — coarser than plot
#   resolution, and SOC distributions across GTK classes within calib_ready
#   were statistically indistinguishable.
#
# ----------------------------------------------------------------------------------------
#
# DERIVED COVARIATES (in site_raw, available for residual analysis)
#
#   Climate normals and variability
#     Computed per plot from monthly normals averaged across all available
#     years in the netCDF:
#       mean_temp           — annual mean air temperature
#       mean_precip_annual  — annual precipitation sum
#       PET_annual          — annual Hargreaves-Samani PET sum
#       coldest_month_T, warmest_month_T — extremes of monthly normals
#       T_seasonality       — SD of monthly normal temperatures
#       P_seasonality       — CV of monthly normal precipitation
#       GDD5                — growing degree days base 5 °C, weighted by
#                              calendar-month length (28.25 for Feb)
#       aridity_index       — mean_precip_annual / PET_annual
#       koppen_class        — Köppen-Geiger classification (Beck et al. 2018
#                              simplified, no arid B or s/w subclasses;
#                              expected codes for Finland: Dfb, Dfc, ET)
#       temp_zone           — tercile binning of mean_temp into cold /
#                              intermediate / warm
#
#   Soil chemistry (derived from Biosoil)
#     CN_ratio          — OrganicCarbon / TotalNitrogen (NA where N == 0)
#     CEC               — sum of ExchangeableCa, Mg, K, Na, Al + FreeHAcidity
#     base_saturation   — (Ca + Mg + K + Na) / CEC
#
#   Litter quality (derived from MUSTIKKA aggregates)
#     litter_{A,W,E,N}_frac — AWEN compound fractions of total plot litter
#     woody_share           — (fwl + cwl) / (nwl + fwl + cwl); cwl (stump
#                              litter) is present as a separate size class in
#                              the Tupek et al. dataset (stumps → cwl)
#     conifer_share         — sp_frac_1 (Scots pine) + sp_frac_2 (Norway
#                              spruce); fraction of litter from conifers
#
#   The CWL × 0.30 ratio is an assumption based on biomass partitioning in
#   Finnish boreal forests; document this explicitly in M&M as it directly
#   affects what enters the model as litter.
#
# ----------------------------------------------------------------------------------------
#
# COVARIATE QUALITY NOTES
#
#   For residual analysis covariate selection, prefer EstimatedBulkDensity
#   over MeanBulkDensity (the latter has 84% NA in calib_ready, audit data
#   not shown).
#
#   OrganicCarbon and OrganicMatter from Biosoil mineral layers are highly
#   collinear with the SOC observations themselves; should be excluded from
#   any residual model that aims to identify covariates of model error
#   (otherwise the analysis recovers the SOC signal trivially).
#
#   profile_depth_cm (measured Biosoil profile depth, max LowerDepthLimit
#   across all M-coded layers) is more informative than soil_depth (the GTK-
#   derived default of 30 cm mineral / 100 cm peat) for any analysis caring
#   about actual soil thickness.
#
# ----------------------------------------------------------------------------------------
#
# OUTPUTS
#
#   ./Data/model_inputs/input_raw_monthly.csv
#     Monthly time series, one row per (plot_id, year, month). Litter inputs
#     in t C/ha/month (annual values divided by 12), climate from netCDF,
#     SOC observations placed in June of each observation year. Filtered to
#     plots with full climate coverage; SOC observations are present for ALL
#     plots with SOC data, regardless of calib_ready status (calib_ready
#     gates filtering at calibration time, not in this table).
#
#   ./Data/model_inputs/climate_monthly.csv
#     Standalone monthly climate per plot, deduplicated.
#
#   ./Data/model_inputs/site_raw.csv
#     Comprehensive per-plot site attributes, stratification metadata, and
#     calib_ready flag. Used both as the calibration filter and as the
#     covariate table for residual analysis.
#
# ----------------------------------------------------------------------------------------
#
# SANITY CHECKS (section 6)
#
#   The script ends with ~25 automated checks covering: no duplicate rows,
#   no NA in climate fields, consistent month counts per plot, non-negative
#   litter, SOC observations placed only in June, SOC counts matching
#   between input_raw_monthly and SOC_agg, plausible SOC value ranges in
#   calib_ready plots (5-400 t C/ha), AWEN fractions summing to 1, climate
#   variables within Finnish ranges, presence of all derived covariates,
#   stand age plausibility, and exclusion of peatland and organic-flagged
#   plots from calib_ready.
#
#   The SOC plausibility check is restricted to calib_ready plots, whose
#   whole-profile stocks span 11.4-244.8 t C/ha. Plots outside the calibration
#   set range more widely (baseline minimum 2.1 t C/ha on a thin soil capped at
#   augering depth, maximum 415.5 t C/ha on plot 33631, which is excluded both
#   as soc_outlier and as organic_missing). Two checks were added with the
#   homogenized baseline: no soc_outlier plot survives into calib_ready, and no
#   calibration-ready plot-year silently falls back to the raw 0-40 cm sum
#   instead of receiving a baseline profile value.
#
# ========================================================================================