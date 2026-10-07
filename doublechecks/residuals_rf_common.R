# =============================================================================
# residuals_rf_common.R — shared by the residual-RF scripts in doublechecks/
#   prep()   NA handling and typing of the predictor table (see
#            residuals_rf_everything.R header for the rules)
#   LAB, lab(), DEV_VARS, DEV_LAB, PREFER   readable labels and the display
#            representative of near-identical predictor families
# =============================================================================

prep <- function(X) {
  X <- X[, colMeans(is.na(X) | X == "", na.rm = TRUE) <= 0.5, drop = FALSE]  # >50% missing out
  for (v in setdiff(names(X), "resid")) {
    x <- X[[v]]
    if (is.logical(x)) x <- as.character(x)
    if (is.numeric(x) && length(unique(na.omit(x))) <= 15) x <- as.character(x)
    if (is.character(x)) {
      x[is.na(x) | x == ""] <- "NA"; X[[v]] <- factor(x)
    } else {
      if (mean(is.na(x)) > 0.05) X[[paste0(v, "_missing")]] <- is.na(x)
      x[is.na(x)] <- median(x, na.rm = TRUE); X[[v]] <- x
    }
  }
  X <- X[, sapply(X, function(x) length(unique(x)) > 1)]          # constants out
  X
}

LAB <- c(
  "wb.bm85_stump" = "Stump biomass 1985", "wb.bm85_stem" = "Stem biomass 1985",
  "wb.bm85_roots" = "Root biomass 1985", "wb.bm85_branches" = "Branch biomass 1985",
  "wb.bm85_foliage" = "Foliage biomass 1985", "wb.ppa_kaikki_85" = "Basal area 1985 (m2/ha)",
  "wb.v_kaikki_85" = "Growing stock 1985 (m3/ha)", "cov.basal_area_85" = "Basal area 1985, stand record (m2/ha)",
  "k85.kuvion_pohja_1" = "Basal area 1985, stand record (m2/ha)",
  "cov.CoarseFragments" = "Coarse fragments (% vol.)", "cov.mean_litter" = "Mean litter input (tC/ha/yr)",
  "cov.dev_class_85" = "Development class 1985", "k85.keh_luokka" = "Development class 1985",
  "cov.woody_share" = "Woody share of litter", "cov.litter_N_frac" = "Litter N (non-soluble) fraction",
  "cov.litter_W_frac" = "Litter W (water-soluble) fraction", "cov.CEC" = "Cation exchange capacity",
  "cov.CN_ratio" = "Soil C:N", "cov.lat_WGS84" = "Latitude", "cov.mean_temp" = "Mean annual temperature",
  "wb.ppa_kaikki_mustikka" = "Basal area 2021-23 (m2/ha)", "m23.ppa_tot" = "Basal area 2021-23 (m2/ha)",
  "wb.ppa_ma_85" = "Scots pine basal area 1985 (m2/ha)", "wb.ppa_ku_85" = "Norway spruce basal area 1985 (m2/ha)",
  "wb.ppa_lehtip_85" = "Broadleaf basal area 1985 (m2/ha)",
  "m23.vall_kehityslk" = "Development class 2021-23", "k90.keh_luokka" = "Development class 1990",
  "k95.kehlk" = "Development class 1995",
  "k90.ehd_hakk" = "Proposed cutting 1990", "k85.ehd_hakk" = "Proposed cutting 1985", "k85.hakk_laatu" = "Cutting done before 1985 (type)",
  "m23.hakkuu1" = "Latest cutting by 2021-23 (type)", "wb.bm95_stump" = "Stump biomass 1995",
  "wb.bm23_stem" = "Stem biomass 2021-23", "wb.bm95_stem" = "Stem biomass 1995", "cov.sp_frac_1" = "Pine share of litter",
  "cov.sp_frac_2" = "Spruce share of litter", "cov.sp_frac_3" = "Birch share of litter",
  "tal.mort_share_85_95" = "Tree mortality 1985-95 (basal-area share)",
  "tal.mort_vol_share_85_95" = "Tree mortality 1985-95 (volume share)",
  "tal.cut_share_85_95" = "Trees cut 1985-95 (basal-area share)",
  "tal.deadwood_share_1985" = "Dead trees 1985 (share of tally)",
  "tal.deadwood_share_1990" = "Dead trees 1990 (share of tally)",
  "tal.deadwood_share_1995" = "Dead trees 1995 (share of tally)",
  "tal.site_index" = "Site index (height for age, log ratio)",
  "tal.stand_age_85_tally" = "Stand age 1985 (sample trees, yr)",
  "chem.bs06_ofh_CN" = "Organic-layer C:N 2006", "chem.mu23_CN" = "Organic-layer C:N 2021-23",
  "chem.ofh_CN_mean" = "Organic-layer C:N (mean 2006, 2021-23)",
  "chem.ofh_CN_logchange" = "Organic-layer C:N change 2006-2021/23 (log ratio)",
  "chem.ofh_N_logchange" = "Organic-layer N change 2006-2021/23 (log ratio)",
  "chem.bs06_ofh_TotalNitrogen" = "Organic-layer N 2006 (g/kg)", "chem.mu23_N" = "Organic-layer N 2021-23 (%)",
  "chem.bs06_ofh_pH.CaCl2." = "Organic-layer pH 2006", "chem.bs06_ofh_BS" = "Organic-layer base saturation 2006",
  "chem.bs06_ofh_ExchangeableCa" = "Organic-layer exch. Ca 2006",
  "chem.bs06_ofh_ExchangeableK" = "Organic-layer exch. K 2006",
  "chem.bs06_ofh_ExchangeableMg" = "Organic-layer exch. Mg 2006",
  "chem.bs06_ofh_ExchangeableAl" = "Organic-layer exch. Al 2006",
  "chem.bs06_m24_pH.CaCl2." = "Mineral 20-40 cm pH 2006", "chem.bs06_m24_BS" = "Mineral 20-40 cm base saturation 2006",
  "chem.mu23_fine_roots_per_soil" = "Fine roots per g organic soil 2021-23",
  "trend.ba_mean" = "Basal area, mean of four inventories (m2/ha)",
  "trend.ba_slope_per_decade" = "Basal area trend 1985-2022 (m2/ha per decade)",
  "trend.ba_change_85_95" = "Basal area change 1985-95 (m2/ha)",
  "trend.ba_change_95_22" = "Basal area change 1995-2022 (m2/ha)",
  "gis.twi16_point" = "Topographic wetness index (16 m)", "gis.twi16_mean3x3" = "Topographic wetness index (48 m mean)",
  "gis.dtw2_point" = "Depth to water, 2 ha (cm)", "gis.dtw2_mean10" = "Depth to water, 2 ha, 10 m mean (cm)",
  "gis.dtw2_wet10" = "Wet share within 10 m (DTW 2 ha < 1 m)", "gis.dtw050_point" = "Depth to water, 0.5 ha (cm)",
  "gis.dtw10_point" = "Depth to water, 10 ha (cm)",
  "gis.ndep_mean_2010_22" = "N deposition to forest 2010-22 (kg N/ha/yr)",
  "gis.ndep_slope_per_decade" = "N deposition trend 2010-22 (per decade)",
  "gis.ndep_oxn_mean" = "Oxidised N deposition (kg N/ha/yr)", "gis.ndep_rdn_mean" = "Reduced N deposition (kg N/ha/yr)",
  "gis.kaski_share" = "Slash-and-burn share of parish", "gis.kaski_1860" = "Slash-and-burn 1860 (class)",
  "gis.kaski_1913" = "Slash-and-burn 1913 (class)", "gis.pop1925_10km" = "Population within 10 km, 1925",
  "gis.pop1925_20km" = "Population within 20 km, 1925", "gis.pop1925_5km" = "Population within 5 km, 1925",
  "gis.dist_rail1925_km" = "Distance to 1925 railway (km)", "gis.state_forest1925" = "State forest in 1925",
  "gis.mki_n" = "Harvest declarations 1997-2026 (n)", "gis.mki_n_thin" = "Thinning declarations (n)",
  "gis.mki_n_regen" = "Regeneration-felling declarations (n)", "gis.mki_n_damage" = "Damage-driven cutting declarations (n)",
  "gis.kem_n" = "KEMERA completed works (n)", "gis.kem_n_tending" = "Young-stand tending, KEMERA (n)",
  "gis.kem_n_fertilise" = "Remedial fertilisation, KEMERA (n)")

# NFI development classes (same coding 1985-2023)
DEV_VARS <- c("cov.dev_class_85", "k85.keh_luokka", "k90.keh_luokka", "k95.kehlk", "m23.vall_kehityslk")
DEV_LAB  <- c(`0` = "Open", `1` = "Seed tree", `2` = "Small sapling", `3` = "Adv. sapling",
              `4` = "Young thinning", `5` = "Adv. thinning", `6` = "Mature",
              `7` = "Shelterwood", `9` = "Uneven-aged")

# Readable class labels for coded predictors drawn on a categorical axis
CUT85_LAB <- c(`0` = "None", `1` = "Young-stand tending", `2` = "Overstorey removal",
               `3` = "First thinning", `4` = "Other thinning", `5` = "Clearing",
               `6` = "Special cutting", `7` = "Regen. (planting)", `8` = "Regen. (natural)",
               `9` = "Clearing")
CODE_LABELS <- c(setNames(rep(list(DEV_LAB), length(DEV_VARS)), DEV_VARS),
                 list(k85.ehd_hakk = CUT85_LAB, k90.ehd_hakk = CUT85_LAB))

# Display representative for a family of near-identical predictors: 1985 basal
# area, volume and the five biomass components are one signal (rho > 0.9);
# basal area is the readable form, so it stands in for the family.
PREFER <- c(setNames(rep("wb.ppa_kaikki_85", 8),
                     c(paste0("wb.bm85_", c("stump", "stem", "roots", "branches", "foliage")),
                       "wb.v_kaikki_85", "k85.kuvion_pohja_1", "cov.basal_area_85")),
            setNames("k85.keh_luokka", "cov.dev_class_85"))
lab <- function(v) if (!is.na(LAB[v])) LAB[[v]] else
  sub("^k(85|90|95)\\.", "\\1 stand: ", sub("^m23\\.", "2021-23 stand: ",
      sub("^wb\\.", "", sub("^cov\\.", "", v))))
