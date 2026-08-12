# =============================================================================
# build_soc_homogenized.R
# -----------------------------------------------------------------------------
# HIKET homogenized SOC baseline: ONE consistent soil-carbon source for all
# three Finnish campaigns (VMI8 1985, Biosoil 2006, Komeetta 2024), built from
# Hannu Ilvesniemi's LUKE workbook. This is intended to be THE baseline SOC
# for HIKET, replacing the earlier soilC1985_2006.csv + Komeetta_mitatut files
# (which over-counted mineral C ~1.6x by omitting the coarse-fragment/stoniness
# correction). The LUKE stocks here are already stoniness/bulk-density corrected
# and internally homogeneous across campaigns.
#
# SOURCE   Data/Komeetta/Hannu/Komeetta 150526hi--.xlsx
#   - sheet "BiSo"      : one row per plot-layer; 2006 & 2024 C stocks in kg/m2
#                          (cols FT/FU), litter in kg/ha (IM/IN). Cleaning recipe
#                          follows J. Heikkinen's BiosoilKomeetta15052026.R and
#                          reproduces his official means (2006=59.1, 2024=61.0).
#   - sheet "Data_1985" : one row per plot; 1985 organic + 0-40 cm stocks (kg/ha),
#                          1985 protocol layers (0-5 / 5-20 / 20-40 cm).
#
# DESIGN   (see README.md in this folder)
#   * Peat plots excluded.
#   * Depth: measured 0-40 cm AND an exponential extrapolation to 1 m are kept as
#     SEPARATE columns (soc_0_40, soc_1m). Extrapolation method ported verbatim
#     from Data_work.R (one decay rate lambda per GTK soil class; C0 per profile).
#   * Region: one consistent value per plot, resolved region.csv -> site_raw ->
#     latitude, with a region_source provenance column. weight = 3 (North) / 1 (South).
#   * All rows carry provenance flags; nothing is silently dropped.
#
# OUTPUTS (this folder)
#   soc_homogenized_layers.csv  long, original layer intervals (drop-in for soilC*.csv)
#   soc_homogenized_plot.csv    per plot x campaign: 0-40, 40-100, 1m, region, covariates
#   soc_homogenized.rds         bundle {layers, plot, meta} for the pipeline
#   plots/*.png                 diagnostic figures
#
# Run from repo root:  Rscript Data/SOC_homogeneized/build_soc_homogenized.R
# =============================================================================

suppressMessages({library(readxl); library(dplyr); library(tidyr)})
options(warn = 1)
set.seed(2025)

ROOT     <- getwd()
XLSX     <- file.path(ROOT, "Data/Komeetta/Hannu/Komeetta 150526hi--.xlsx")
REGION   <- file.path(ROOT, "Data/Komeetta/Juha/region.csv")
SITE_RAW <- file.path(ROOT, "Data/model_inputs/site_raw.csv")
SITE_KEY <- file.path(ROOT, "Data/soil_litter_site_key.csv")   # ';'-sep; carries northings for 1985-only plots
OUTDIR   <- file.path(ROOT, "Data/SOC_homogeneized")
PLOTDIR  <- file.path(OUTDIR, "plots")
dir.create(PLOTDIR, showWarnings = FALSE, recursive = TRUE)

NORTH_Y_THRESHOLD <- 7300000   # northing (m) boundary for the latitude fallback (ETRS/YKJ ~equal here)
SOC_OUTLIER_MAX   <- 250       # Mg/ha; soc_profile above this is physically implausible for boreal
                               #   mineral soil to 1 m -> soc_outlier = TRUE (excluded from calibration)
HIGH_CHANGE_RATE  <- 3         # tC/ha/yr; |consecutive-campaign change|/yr above this -> high_change
                               #   = TRUE (flagged for inspection but KEPT: likely resampling noise)
num <- function(x) suppressWarnings(as.numeric(x))            # quiet text->numeric coercion

# =============================================================================
# Helpers: exponential depth extrapolation (ported from Data_work.R, unchanged)
#   A single decay rate lambda is fitted per GTK soil class by concentrating out
#   each profile's surface concentration C0; the 40-100 cm stock is then the
#   integral of C0*exp(-lambda*z) over 40..100 cm.
# =============================================================================
layer_geom <- function(lyr) {
  b <- switch(lyr,
              "0-5cm" = c(0,5), "5-20cm" = c(5,20),
              "0-10cm" = c(0,10), "10-20cm" = c(10,20),
              "20-40cm" = c(20,40), NULL)
  if (is.null(b)) NULL else list(lower = b[1], upper = b[2])
}
layer_weights <- function(layer_names, lambda) {
  g  <- lapply(layer_names, layer_geom)
  z1 <- sapply(g, `[[`, "lower"); z2 <- sapply(g, `[[`, "upper")
  (1/lambda) * (exp(-lambda*z1) - exp(-lambda*z2))
}
fit_pooled_lambda <- function(class_dat) {
  g  <- lapply(class_dat$layer, layer_geom)
  z1 <- sapply(g, `[[`, "lower"); z2 <- sapply(g, `[[`, "upper")
  C  <- class_dat$C_kgha; py <- paste(class_dat$plot_id, class_dat$year)
  neg_rss <- function(lambda) {
    if (lambda <= 0) return(1e12)
    w <- (1/lambda) * (exp(-lambda*z1) - exp(-lambda*z2))
    rss <- 0
    for (p in unique(py)) {
      i <- py == p; C0 <- sum(C[i]*w[i]) / sum(w[i]^2)
      rss <- rss + sum((C[i] - C0*w[i])^2)
    }
    rss
  }
  optimise(neg_rss, c(1e-4, 2))$minimum
}
# Extrapolated stock from the deepest measured layer bottom (z_from) down to z_to, given
# lambda. Needs >=2 mineral layers (C>0). Returns 0 if z_to <= z_from (thin soil: nothing
# below the measured profile) so the extrapolation never extends past the plot's soil depth.
extrap_below <- function(min_dat, lambda, z_from, z_to) {
  if (nrow(min_dat) < 2) return(NA_real_)
  if (z_to <= z_from)   return(0)
  w  <- layer_weights(min_dat$layer, lambda)
  C0 <- sum(min_dat$C_kgha * w) / sum(w^2)
  C0 / lambda * (exp(-lambda*z_from) - exp(-lambda*z_to))
}

# =============================================================================
# 1. Read BiSo (2006 + 2024), by column index (row 3 = header, data from row 4)
# =============================================================================
# read the whole span as text (avoids per-column type-guess warnings), coerce explicitly
biso_span  <- suppressMessages(read_excel(
  XLSX, sheet = "BiSo", range = "A4:IN3031", col_names = FALSE, col_types = "text", .name_repair = "minimal"))
gc <- function(i) biso_span[[i]]   # grab column by spreadsheet index

biso <- tibble(
  VMI   = as.integer(num(gc(3))),   Krs   = as.integer(num(gc(4))),
  CODE_LAYER = as.character(gc(8)), REP = as.integer(num(gc(9))),
  LINF  = num(gc(11)),   # layer lower limit (cm); deepest mineral bottom = augering depth reached
  MAANAYTE = num(gc(140)),  # "maanaytteen ottovuosi": the FIRST campaign's soil sampling YEAR
  OLW   = num(gc(22)),  nhumus = num(gc(96)), notes = num(gc(100)),
  Ccont24 = num(gc(170)),
  C2006_kgm2 = num(gc(176)), C2024_kgm2 = num(gc(177)),
  lit2006 = num(gc(247)), lit2024 = num(gc(248)),
  # 2006 stand inventory (time-varying):
  basal_area = num(gc(92)), height = num(gc(93)), stand_age = num(gc(136)),
  # 2006-era management history:
  mgmt_op = num(gc(102)), mgmt_yrs_since = num(gc(103)), mgmt_residue = num(gc(104)),
  # static soil physics (Biosoil 2006, measured once):
  clay = num(gc(15)), silt = num(gc(16)), sand = num(gc(17)),
  texture_class = as.character(gc(18)), bd_est = num(gc(20)),
  coarse_frag = num(gc(21)), pH_CaCl2 = num(gc(23))
) |> filter(!is.na(VMI))

# --- provenance flags (Heikkinen cleaning) ---
peat_ids    <- biso |> filter(CODE_LAYER %in% c("H01","H12")) |> distinct(VMI) |> pull(VMI)
aband_ids   <- biso |> filter(notes == 0) |> distinct(VMI) |> pull(VMI)
unmeas_ids  <- biso |> group_by(VMI) |> summarise(u = all(is.na(Ccont24))) |> filter(u) |> pull(VMI)

# --- per-plot soil-depth proxy (augering depth reached) + measured 40-80 cm layer ---
# LAYER_LIMIT_INF of the deepest mineral layer: 80 = reached standard target (soil >= 80);
# 10/20/40 = augering refusal (bedrock/stones) => soil ~ that depth. Krs 204 = MEASURED 40-80 cm
# (2006 only; dropped from the homogeneous layer stack but kept here to validate/anchor the tail).
reached_depth <- biso |> filter(REP == 1, Krs %in% 201:204) |>
  group_by(plot_id = VMI) |> summarise(reached_cm = suppressWarnings(max(LINF, na.rm = TRUE)), .groups = "drop") |>
  mutate(reached_cm = ifelse(is.finite(reached_cm), reached_cm, NA_real_))
soc_40_80_meas <- biso |> filter(REP == 1, Krs == 204) |>
  transmute(plot_id = VMI, soc_40_80_meas = 1e4 * C2006_kgm2)

# Keep REP 1, drop 40-60 cm layer (204, not measured in 2024)
C <- biso |> filter(REP == 1, Krs != 204) |>
  mutate(
    # organic C set to 0 where there is no organic layer (no ORGANIC_LAYER_WEIGHT)
    C2006_kgm2 = ifelse(Krs == 101 & is.na(OLW), 0, C2006_kgm2),
    C2024_kgm2 = ifelse(Krs == 101 & is.na(OLW), 0, C2024_kgm2)
  )
# plots missing 2024 C in organic or topmost mineral layer (flag, do not hard-drop)
miss24_ids <- C |> filter(is.na(C2024_kgm2) & Krs %in% c(101,201)) |> distinct(VMI) |> pull(VMI)

C <- C |>
  replace_na(list(C2006_kgm2 = 0, C2024_kgm2 = 0, lit2006 = 0, lit2024 = 0)) |>
  mutate(
    layer  = recode(as.character(Krs), "101"="organic","201"="0-10cm","202"="10-20cm","203"="20-40cm"),
    # kg/ha: mineral = kg/m2 * 1e4; organic = kg/m2 * 1e4 + litter (kg/ha)  [Heikkinen]
    C2006  = ifelse(Krs == 101, 1e4*C2006_kgm2 + lit2006, 1e4*C2006_kgm2),
    C2024  = ifelse(Krs == 101, 1e4*C2024_kgm2 + lit2024, 1e4*C2024_kgm2)
  )

# long layer table for 2006 and 2024. Drop incomplete-2024 profiles (abandoned /
# unmeasured / missing organic|topmost) rather than zero-filling them, which would
# deflate the 2024 mean. Their valid 2006 rows are kept (flagged in plot_tab).
incomplete_2024 <- Reduce(union, list(aband_ids, unmeas_ids, miss24_ids))
layers_0624 <- bind_rows(
  C |> transmute(plot_id = VMI, year = 2006L, layer, C_kgha = C2006),
  C |> transmute(plot_id = VMI, year = 2024L, layer, C_kgha = C2024)
) |> filter(!is.na(layer), !(year == 2024L & plot_id %in% incomplete_2024))

# --- static soil covariates (plot-level; Biosoil-2006-measured, effectively constant) ---
soil_covars <- biso |> group_by(plot_id = VMI) |> summarise(
  clay = mean(clay, na.rm = TRUE), silt = mean(silt, na.rm = TRUE),
  sand = mean(sand, na.rm = TRUE), bd_est = mean(bd_est, na.rm = TRUE),
  coarse_frag = mean(coarse_frag, na.rm = TRUE), pH_CaCl2 = mean(pH_CaCl2, na.rm = TRUE),
  texture_class = { t <- na.omit(texture_class); if (length(t)) names(sort(table(t),decreasing=TRUE))[1] else NA_character_ },
  .groups = "drop") |> mutate(across(where(is.numeric), ~ ifelse(is.nan(.x), NA_real_, .x)))

# --- 2006 stand inventory + management (plot-level, from BiSo; year-tagged 2006) ---
stand_2006 <- biso |> group_by(plot_id = VMI) |> summarise(
  basal_area = first(na.omit(basal_area)), mean_height = first(na.omit(height)),
  stand_age = first(na.omit(stand_age)), mgmt_op = first(na.omit(mgmt_op)),
  mgmt_yrs_since = first(na.omit(mgmt_yrs_since)), mgmt_residue = first(na.omit(mgmt_residue)),
  .groups = "drop") |> mutate(across(where(is.numeric), ~ ifelse(is.nan(.x), NA_real_, .x)))

# =============================================================================
# 2. Read Data_1985 (organic + 0-5 / 5-20 / 20-40 cm), to long layer table
# =============================================================================
d85 <- suppressMessages(read_excel(XLSX, sheet = "Data_1985", .name_repair = "minimal")) |> as.data.frame()

# -----------------------------------------------------------------------------
# TREATMENT C (2026-08-12): the 1985 organic layer is OFH ONLY.
#
# LM (litter+moss) is a separately coded layer in the protocol -- see sheet
# Massat_2006, where Krs 100 = LM and Krs 101 = OFH. It was measured in 2006 and
# 2024 and folded into their organic layer (cols 247/248, Heikkinen's convention),
# and it is ABSENT from Data_1985. Evidence that the 1985 column is the OFH-only
# analogue: BiSo lays the three campaigns out in PARALLEL blocks with identical
# column names (1985 158-163, 2006 164-169, 2024 170-174); the 1985 block is
# byte-identical to this sheet (414/414 rows); and the 2006 C_kgha equals
# 1e4*C_kgm2 exactly, i.e. litter-free, since the litter is added on top of it.
# Corroborated independently by dry mass: 1985 43.7 vs 2006 OFH 50.0 vs OFH+LM 56.7
# (kg/ha x1e3) -- 1985 sits below the humus-only mass, nowhere near humus+litter.
#
# Leaving 1985 without a litter term is NOT the neutral option: it asserts the 1985
# litter layer weighed zero, which is certainly false. We add each plot's OWN 2006
# LM; where that is unavailable, the campaign median.
#
# NOT applied where the 1985 organic layer is zero or absent: those plots carry the
# organic_zero / organic_missing flags computed downstream in Data_work.R, and an
# imputed LM would silently make them calibration-ready.
# ⚠ REVERSE THIS if H. Ilvesniemi confirms that the 1985 value already includes LM.
# -----------------------------------------------------------------------------
ADD_1985_LM <- !identical(Sys.getenv("HIKET_ADD_1985_LM"), "0")   # env override, for A/B testing

lm_2006 <- biso |> filter(REP == 1, Krs == 101) |>
  transmute(plot_id = VMI, LM = lit2006) |>
  filter(is.finite(LM)) |> distinct(plot_id, .keep_all = TRUE)
LM_FALLBACK <- median(lm_2006$LM, na.rm = TRUE)
lm_lookup   <- setNames(lm_2006$LM, as.character(lm_2006$plot_id))

# the first campaign's per-plot sampling year (1986..1995), kept alongside `year`
samp_year_85 <- biso |> filter(REP == 1, MAANAYTE %in% 1980:2000) |>
  transmute(plot_id = VMI, samp_year = as.integer(MAANAYTE)) |>
  distinct(plot_id, .keep_all = TRUE)

.org85 <- d85$Corg_kgha
.lm85  <- unname(lm_lookup[as.character(d85$VMI)])
.lm85[!is.finite(.lm85)] <- LM_FALLBACK
.use   <- ADD_1985_LM & is.finite(.org85) & .org85 > 0        # never on zero/absent OFH
.org85_out <- ifelse(.use, .org85 + .lm85, .org85)
lm_added_1985 <- data.frame(plot_id = d85$VMI,
                            lm_added_1985 = ifelse(.use, .lm85, 0),
                            lm_imputed_1985 = .use & !(as.character(d85$VMI) %in% names(lm_lookup)))
message(sprintf("Treatment C: LM added to %d of %d 1985 organic rows (%d from the plot's own 2006 value, %d from the median %.3f Mg/ha); %d skipped (zero/absent OFH)",
                sum(.use), length(.org85),
                sum(.use & as.character(d85$VMI) %in% names(lm_lookup)),
                sum(.use & !(as.character(d85$VMI) %in% names(lm_lookup))),
                LM_FALLBACK/1000, sum(!.use)))

layers_1985 <- bind_rows(
  data.frame(plot_id = d85$VMI, year = 1985L, layer = "organic", C_kgha = .org85_out),
  data.frame(plot_id = d85$VMI, year = 1985L, layer = "0-5cm",   C_kgha = d85$C0_5_kgha),
  data.frame(plot_id = d85$VMI, year = 1985L, layer = "5-20cm",  C_kgha = d85$C5_20_kgha),
  data.frame(plot_id = d85$VMI, year = 1985L, layer = "20-40cm", C_kgha = d85$C20_40_kgha)
) |> filter(!is.na(plot_id))

# =============================================================================
# 3. Long table (all campaigns), original intervals + geometry, drop peat
# =============================================================================
site <- read.csv(SITE_RAW)                    # plot_id-keyed attributes (pipeline)
peat_all <- union(peat_ids, site$plot_id[site$peatland %in% TRUE])   # + 1985-only peat via site_raw

# site descriptors (per plot, Komeetta-era 2024): dev class (KEHLK), site fertility, dominant species
kohde <- suppressMessages(read_excel(XLSX, sheet = "Kohdekuvaukset", .name_repair = "minimal")) |>
  as.data.frame()
kohde <- kohde[, 1:5]; names(kohde) <- c("plot_id","site_main_class","site_fertility","dev_class_2024","dominant_species")
kohde[] <- lapply(kohde, function(x) suppressWarnings(as.numeric(x)))
kohde <- kohde[!is.na(kohde$plot_id), ]

# --- plot-level covariate table, EVERY column tagged with its measurement year/era ----------
#   Stand properties are time-varying, so each year gets its own column (wide) rather than one
#   ambiguous column. Height in metres (site_raw is dm -> /10; BiSo Keskipituus already m).
#     _1985 = NFI8 (site_raw)   _2006 = Biosoil (BiSo)   _2024 = Komeetta (Kohdekuvaukset)
#   Soil physics are Biosoil-2006 measurements (effectively static); site descriptors Komeetta-2024.
covars <- site |> transmute(plot_id,
    basal_area_1985 = basal_area_85, stand_age_1985 = stand_age_85,
    mean_height_1985_m = mean_height_85_dm / 10, dev_class_1985 = dev_class_85) |>
  full_join(stand_2006 |> transmute(plot_id,
    basal_area_2006 = basal_area, stand_age_2006 = stand_age, mean_height_2006_m = mean_height,
    mgmt_op_2006 = mgmt_op, mgmt_yrs_since_2006 = mgmt_yrs_since, mgmt_residue_2006 = mgmt_residue),
    by = "plot_id") |>
  full_join(kohde |> transmute(plot_id, dev_class_2024,
    site_main_class_2024 = site_main_class, site_fertility_2024 = site_fertility,
    dominant_species_2024 = dominant_species), by = "plot_id") |>
  full_join(soil_covars |> rename_with(~ paste0(.x, "_2006"), -plot_id), by = "plot_id")

layers <- bind_rows(layers_1985, layers_0624) |>
  filter(!plot_id %in% peat_all) |>
  mutate(
    campaign = recode(as.character(year), "1985"="VMI8","2006"="Biosoil","2024"="Komeetta"),
    depth_lower_cm = sapply(layer, function(l) { g<-layer_geom(l); if(is.null(g)) 0 else g$lower }),
    depth_upper_cm = sapply(layer, function(l) { g<-layer_geom(l); if(is.null(g)) NA else g$upper }),
    C_Mgha = C_kgha / 1000
  ) |>
  arrange(plot_id, year, depth_lower_cm)

# =============================================================================
# 4. Region resolution + design weights  (region.csv -> site_raw -> latitude)
# =============================================================================
region_csv <- read.csv(REGION)[, c("VMI","region")]  # 1 = South, 2 = North
site_reg   <- site[, c("plot_id","region","y_ETRS","x_ETRS","lon_WGS84","lat_WGS84","soil_code")]
key        <- read.csv(SITE_KEY, sep = ";")          # northing fallback for 1985-only plots
north_key  <- setNames(key$y, key$koealatunnus_BIOSOIL)   # YKJ northing keyed by VMI code

plots <- tibble(plot_id = sort(unique(layers$plot_id))) |>
  left_join(region_csv, by = c("plot_id" = "VMI")) |> rename(region_design = region) |>
  left_join(site_reg,  by = "plot_id") |> rename(region_site = region) |>
  mutate(
    northing_key = north_key[as.character(plot_id)],           # fallback northing (YKJ)
    region = case_when(
      !is.na(region_design)                          ~ region_design,
      !is.na(region_site)                            ~ region_site,
      !is.na(y_ETRS)       & y_ETRS       >= NORTH_Y_THRESHOLD ~ 2L,
      !is.na(y_ETRS)                                 ~ 1L,
      !is.na(northing_key) & northing_key >= NORTH_Y_THRESHOLD ~ 2L,
      !is.na(northing_key)                           ~ 1L,
      TRUE                                           ~ 1L        # last resort: South
    ),
    region_source = case_when(
      !is.na(region_design) ~ "biosoil_design",    # Juha's region.csv (design N/S)
      !is.na(region_site)   ~ "site_raw_region",    # our pipeline site_raw$region
      !is.na(y_ETRS)        ~ "latitude_ETRS",      # inferred from site_raw northing
      !is.na(northing_key)  ~ "latitude_sitekey",   # inferred from site-key northing
      TRUE                  ~ "default_South"        # no geodata anywhere
    ),
    region_conflict = !is.na(region_design) & !is.na(region_site) & region_design != region_site,
    region_label = ifelse(region == 2L, "North", "South"),
    weight = ifelse(region == 2L, 3L, 1L)
  )

# =============================================================================
# 5. Depth extrapolation to 1 m (fit lambda per GTK soil class, global fallback)
# =============================================================================
soil_lookup <- setNames(as.character(site$soil_code), site$plot_id)
minl <- layers |>
  filter(layer %in% c("0-5cm","5-20cm","10-20cm","0-10cm","20-40cm"), C_kgha > 0) |>
  mutate(soil_code = soil_lookup[as.character(plot_id)])

# per-class lambda from profiles with >=2 mineral layers
fit_dat <- minl |> group_by(plot_id, year) |> filter(n() >= 2) |> ungroup()
global_lambda <- fit_pooled_lambda(as.data.frame(fit_dat))
class_lambda  <- fit_dat |> filter(!is.na(soil_code)) |> group_by(soil_code) |>
  filter(length(unique(paste(plot_id,year))) >= 3) |>
  summarise(lambda = fit_pooled_lambda(pick(everything())), .groups = "drop")
lam_lookup <- setNames(class_lambda$lambda, class_lambda$soil_code)

# per profile: measured stock (organic + mineral to deepest measured layer, <=40 cm)
prof <- layers |> group_by(plot_id, year, campaign) |> summarise(
    organic         = sum(C_kgha[layer == "organic"], na.rm = TRUE),
    mineral_0_40    = sum(C_kgha[layer != "organic"], na.rm = TRUE),
    n_mineral       = sum(layer != "organic" & C_kgha > 0),
    .groups = "drop") |>
  mutate(soc_0_40 = organic + mineral_0_40)

# Per-plot integration cap at the soil depth: thin soils (augering refused at 10/20/40 cm)
# stop there; deep/unknown soils go to the full reference depth. The extrapolation runs from
# each profile's deepest measured layer down to z_cap, so it never extends past actual soil.
REF_DEPTH      <- 100
reached_lookup <- setNames(reached_depth$reached_cm, reached_depth$plot_id)
z_cap_of <- function(pid){ r <- reached_lookup[as.character(pid)]
  if (length(r) == 0 || is.na(r) || r >= 80) REF_DEPTH else r }

ext <- minl |> group_by(plot_id, year) |> group_modify(function(d, key){
    lam    <- lam_lookup[as.character(d$soil_code[1])]; if (is.na(lam)) lam <- global_lambda
    z_from <- max(sapply(d$layer, function(l) layer_geom(l)$upper))   # deepest measured bottom
    zc     <- z_cap_of(key$plot_id)
    tibble(lambda = lam, dm_bottom = z_from, z_cap = zc,
           soc_deep          = extrap_below(as.data.frame(d), lam, z_from, zc),
           soc_deep_uncapped = extrap_below(as.data.frame(d), lam, z_from, REF_DEPTH),
           soc_40_80_pred    = extrap_below(as.data.frame(d), lam, 40, 80))  # vs measured Krs 204
  }) |> ungroup()

prof <- prof |> left_join(ext, by = c("plot_id","year")) |>
  mutate(
    fit_ok            = !is.na(soc_deep),
    soc_deep          = ifelse(is.na(soc_deep), 0, soc_deep),
    soc_deep_uncapped = ifelse(is.na(soc_deep_uncapped), 0, soc_deep_uncapped),
    soc_profile       = soc_0_40 + soc_deep,           # depth-capped whole-profile (PRIMARY)
    soc_1m_uncapped   = soc_0_40 + soc_deep_uncapped   # old fixed-100 cm (for comparison)
  )

# =============================================================================
# 6. Per plot x campaign table: stocks (kg/ha AND Mg/ha) + region + covariates + flags
# =============================================================================
plot_tab <- prof |>
  left_join(plots |> select(plot_id, region, region_label, region_source, region_conflict,
                            weight, x_ETRS, y_ETRS, lon_WGS84, lat_WGS84, soil_code), by = "plot_id") |>
  left_join(covars, by = "plot_id") |>                          # year-tagged covariates (wide)
  left_join(reached_depth, by = "plot_id") |>                   # soil-depth proxy (augering reached)
  left_join(samp_year_85,  by = "plot_id") |>                   # true first-campaign year
  left_join(lm_added_1985, by = "plot_id") |>                   # treatment C provenance
  left_join(soc_40_80_meas, by = "plot_id") |>                  # measured 40-80 cm (2006, Krs 204)
  mutate(
    soil_depth_reached = reached_cm,
    # `year` stays the CAMPAIGN KEY (1985/2006/2024). `samp_year` is the true
    # sampling year: 1986-1995 for VMI8, equal to `year` for the later campaigns.
    # ⚠ Do not merge these two -- `year` keys sigma_infl and HIKET_DROP_CAMPAIGN.
    samp_year        = ifelse(year == 1985L, samp_year, year),
    lm_added_1985    = ifelse(year == 1985L, coalesce(lm_added_1985, 0), 0),
    lm_imputed_1985  = ifelse(year == 1985L, coalesce(lm_imputed_1985, FALSE), FALSE),
    soc_40_80_meas   = ifelse(year == 2006L, soc_40_80_meas, NA_real_),  # measured only in 2006
    is_peat          = FALSE,                          # peat already removed
    abandoned_2024   = plot_id %in% aband_ids,
    unmeasured_2024  = plot_id %in% unmeas_ids,
    missing_C2024    = plot_id %in% miss24_ids,
    # Heikkinen's paired 2006<->2024 subset (reproduces his official means):
    juha_subset      = !plot_id %in% union(aband_ids, union(unmeas_ids, miss24_ids)),
    across(c(organic, mineral_0_40, soc_0_40, soc_deep, soc_profile, soc_1m_uncapped,
             soc_40_80_meas, soc_40_80_pred),
           ~ .x/1000, .names = "{.col}_Mgha")
  ) |>
  arrange(plot_id, year)

# --- outlier flags (thresholds = SOC_OUTLIER_MAX / HIGH_CHANGE_RATE constants at top) ---
#   soc_outlier  (per plot-year): implausibly high stock -> EXCLUDE from calibration.
#   high_change  (per plot): a consecutive-campaign change exceeding the rate bound -> KEEP but flag
#                (likely resampling noise, since campaigns don't re-core the identical soil volume).
hc <- plot_tab |> filter(!is.na(soc_profile_Mgha)) |> arrange(plot_id, year) |>
  group_by(plot_id) |>
  summarise(high_change = length(year) >= 2 &&
              any(abs(diff(soc_profile_Mgha) / diff(year)) > HIGH_CHANGE_RATE), .groups = "drop")
plot_tab <- plot_tab |>
  left_join(hc, by = "plot_id") |>
  mutate(soc_outlier = soc_profile_Mgha > SOC_OUTLIER_MAX,
         high_change = ifelse(is.na(high_change), FALSE, high_change))

# =============================================================================
# 7. Validation: reproduce Heikkinen's official means EXACTLY (his exclusions:
#    peat via CODE_LAYER H01/H12 only, paired 2006<->2024 subset, org + 0-40).
#    (The baseline itself uses a stricter peat filter, so its means differ slightly.)
# =============================================================================
juha_drop <- Reduce(union, list(peat_ids, aband_ids, unmeas_ids, miss24_ids))
val <- layers_0624 |> filter(!plot_id %in% juha_drop) |>
  group_by(plot_id, year) |> summarise(soc_0_40 = sum(C_kgha), .groups = "drop") |>
  left_join(region_csv, by = c("plot_id" = "VMI")) |> mutate(w = ifelse(region == 2L, 3L, 1L)) |>
  group_by(year) |> summarise(wmean = weighted.mean(soc_0_40, w), n = n(), .groups = "drop")
tgt <- c(`2006` = 59055, `2024` = 61047)
cat("\n=== VALIDATION vs Heikkinen (target 2006 = 59055, 2024 = 61047 kg/ha) ===\n")
for (i in seq_len(nrow(val)))
  cat(sprintf("  %d: %.0f kg/ha (n=%d)  %s\n", val$year[i], val$wmean[i], val$n[i],
              ifelse(abs(val$wmean[i] - tgt[as.character(val$year[i])]) < 50, "MATCH", "MISMATCH")))

cat("\n=== region provenance ===\n"); print(table(plot_tab |> distinct(plot_id, region_source) |> pull(region_source)))
cat("region conflicts (design vs site_raw):", sum(plots$region_conflict), "\n")

cat(sprintf("\n=== outlier flags ===\n  soc_outlier (soc_profile > %d Mg/ha): %d plot-years, %d plots  [EXCLUDE from calibration]\n",
            SOC_OUTLIER_MAX, sum(plot_tab$soc_outlier, na.rm = TRUE),
            length(unique(plot_tab$plot_id[plot_tab$soc_outlier %in% TRUE]))))
cat(sprintf("  high_change (|rate| > %d tC/ha/yr): %d plots  [flagged, KEPT]\n",
            HIGH_CHANGE_RATE, length(unique(plot_tab$plot_id[plot_tab$high_change]))))

# =============================================================================
# 8. Write outputs
# =============================================================================
layers_out <- layers |> select(plot_id, year, campaign, layer, depth_lower_cm, depth_upper_cm, C_kgha, C_Mgha)
write.csv(layers_out, file.path(OUTDIR, "soc_homogenized_layers.csv"), row.names = FALSE)
write.csv(plot_tab,   file.path(OUTDIR, "soc_homogenized_plot.csv"),   row.names = FALSE)
saveRDS(list(layers = layers_out, plot = plot_tab,
             meta = list(source = basename(XLSX), built = Sys.time(),
                         global_lambda = global_lambda, class_lambda = class_lambda,
                         validation = val, north_y_threshold = NORTH_Y_THRESHOLD)),
        file.path(OUTDIR, "soc_homogenized.rds"))
cat(sprintf("\nWrote: %d layer rows, %d plot x campaign rows, %d unique plots\n",
            nrow(layers_out), nrow(plot_tab), length(unique(plot_tab$plot_id))))

# =============================================================================
# 9. Diagnostic plots
# =============================================================================
source(file.path(OUTDIR, "plots_soc_homogenized.R"), local = TRUE)
cat("\nDone. Outputs + plots in", OUTDIR, "\n")
