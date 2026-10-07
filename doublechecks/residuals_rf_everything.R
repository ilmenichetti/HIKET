# =============================================================================
# residuals_rf_everything.R — every plausible predictor of the local residual
#
# WHY. The pipeline RF (run_residual_analysis.R) offers ~57 covariates and only
# the 1985 stand. The stand records hold much more (defoliation, damage, lichens,
# management done and proposed, structure, four inventories of basal area), and
# basal area's signal fades 1985 -> 2021-23. Offer everything that can remotely
# make sense and see what the data pick.
#
# WHAT. Target: mean log residual of the TARGET models, one value per plot
# (+ = model under-predicts). Default Yasso07 + Yasso15, the two models taken as
# the reference for the residual study; HIKET_RF_MODELS="SP1,TP2,..." overrides
# ("all" = the six). Predictors, by source:
#   cov    residuals_df plot covariates (climate, soil chemistry, litter, location)
#   k85/k90/k95  every stand-record column of the 1985/1990/1995 NFI (main stand)
#   m23    2021-23 stand record + dominant tree storey (MU21_23_*)
#   wb     tree-tally basal area, volume and biomass by year (kooste workbook)
#   tal    from the individual-tree tallies: mortality and cutting 1985-95, dead
#          wood 1985/90/95, site index, stand age (residuals_tally_predictors.R,
#          run it first)
#   chem   BioSoil 2006 organic-layer chemistry + MUSTIKKA 2021-23 C:N, N, fine
#          roots, levels and 2006->2021-23 trend (residuals_chemistry_predictors.R)
#   trend  basal area mean and trend over the four inventories
#   gis    GIS point extractions (GIS/README.md): wetness, N deposition, harvest
#          declarations, KEMERA works, 1925 history
#   und    understorey % cover 2023 by vegetation group (Data/Understorey)
# CIRCULAR predictors (measure soil carbon itself: organic-layer mass/thickness,
# soil organic C/N) are EXCLUDED from the main run; HIKET_RF_CIRC=1 adds a run
# with them.
# NA handling (ranger 0.15 has no na.learn): numeric -> median + missing flag,
# factor -> explicit "NA" level (in the NFI, "not recorded" often means "not
# applicable", e.g. open areas, so it is informative). Integer codes with
# <= 15 values are treated as factors.
# Importance, robust to collinearity: permutation (5 seeds), DROP-COLUMN for the
# top 25, and GROUP drop (refit without a whole theme). A pure-noise column
# gives the floor.
#
# Run from repo root:  Rscript doublechecks/residuals_rf_everything.R
# Writes doublechecks/figures/residuals_rf_everything_<target>_*.csv and
# residuals_rf_everything_<target>.rds (raw predictors + importance, read by
# residuals_rf_top_scatter.R).
# =============================================================================

suppressPackageStartupMessages({library(dplyr); library(readxl); library(ranger)})
source("manuscript/figures/run_ids.R")
source("doublechecks/residuals_rf_common.R")
cat("RUN_IDs:", paste(names(RID), RID, sep = "=", collapse = "  "), "\n\n")
OUT <- "doublechecks/figures"; dir.create(OUT, showWarnings = FALSE)
SEEDS <- 2025 + 0:4
TGT <- strsplit(Sys.getenv("HIKET_RF_MODELS", "Yasso07,Yasso15"), ",")[[1]]
if (identical(TGT, "all")) TGT <- FIG_MODELS
TAG <- if (setequal(TGT, FIG_MODELS)) "six" else paste(TGT, collapse = "_")
cat("Target = mean residual of:", paste(TGT, collapse = ", "), "\n")

# === 1. Target + residuals_df covariates =====================================
R <- bind_rows(lapply(TGT, function(m)
  as.data.frame(readRDS(sprintf(
    "Calibration_real_data_transient/runs/%s_posterior_predictive_%s.rds",
    m, RID[[m]]))$residuals_df) |> mutate(model = m)))
RM <- R |> group_by(plot_id, model) |>                 # per-model plot means
  summarise(r = mean(residual_log), .groups = "drop") |>
  tidyr::pivot_wider(names_from = model, values_from = r)
not_cov <- c("year", "soc_obs_tCha", "soc_mean", "soc_median", "soc_q025",
  "soc_q975", "soc_pp_q025", "soc_pp_q975", "log_obs", "log_hat_mean",
  "residual_log", "residual_abs", "is_first", "species_name", "soil_type",
  "peatland", "organic_missing", "organic_zero", "high_mrt", "has_climate",
  "has_soc", "soc_outlier", "high_change", "const_litter", "calib_ready",
  "is_holdout", "mean_soc_Mgha", "model")
E <- R |> group_by(plot_id) |>
  summarise(resid = mean(residual_log),
            across(-any_of(c(not_cov)), first), .groups = "drop")
names(E)[-(1:2)] <- paste0("cov.", names(E)[-(1:2)])

key <- read.csv("Data/soil_litter_site_key.csv", sep = ";")
key <- key[match(E$plot_id, key$koealatunnus_BIOSOIL), ]

# === 2. Stand records 1985/1990/1995 (main stand, vmikuvio 0) ================
admin <- c("koeala", "vmikuvio", "ed_kuvio", "kunta", "kunta_nro", "metsakeskus",
  "metsalautakunta", "mltk", "koko1", "koko3", "mitt_tapa", "mittaustapa",
  "ei_mita", "tiet_historia", "kuv_raj_et", "kuv_raj_su", "maisemaraj_suu",
  "maisemaraj_et", "reuna_1_et", "reuna_2_et")
kuv <- function(y, pre) {
  d <- read.csv(sprintf("Data/PysyvätKoealat/%d/kuvio.csv", y), stringsAsFactors = FALSE)
  d <- d[d$vmikuvio == 0, ]
  d <- d[match(key$koealatunnus_VANHA, d$koeala), setdiff(names(d), admin)]
  setNames(d, paste0(pre, ".", names(d)))
}
K <- cbind(kuv(1985, "k85"), kuv(1990, "k90"), kuv(1995, "k95"))

# === 3. 2021-23 stand record + dominant storey ===============================
m23 <- read.csv("Data/biomassatiedostot/23/MU21_23_metsikkokuvio.csv", sep = ";",
                stringsAsFactors = FALSE)
m23 <- m23[m23$vmikuvio == 0, ]
m23 <- m23[match(key$site, m23$koealatunnus),
           setdiff(names(m23), c("koealatunnus", "ryhma", "vmikuvio", "ed_kuvio",
                                 "koko1", "koko3", grep("^huom", names(m23), value = TRUE)))]
pj <- read.csv("Data/biomassatiedostot/23/MU21_23_puujakso.csv", sep = ";",
               stringsAsFactors = FALSE)
pj <- pj[pj$vmikuvio == 0, ]
for (v in c("ppa", "keskilpm", "keskipituus", "valtapituus", "jakson_ika",
            "paapl_osuus", "havu_lehti_osuus", "taimi_kok_lkm"))
  pj[[v]] <- suppressWarnings(as.numeric(pj[[v]]))
agg <- pj |> group_by(koealatunnus) |>
  summarise(n_storeys = n(), ppa_tot = sum(ppa, na.rm = TRUE), .groups = "drop")
dom <- pj[order(pj$koealatunnus, -replace(pj$ppa, is.na(pj$ppa), -1)), ]   # largest storey first
dom <- dom[!duplicated(dom$koealatunnus), c("koealatunnus", "vallitseva_pl",
  "paapl_osuus", "havu_lehti_osuus", "keskilpm", "keskipituus", "valtapituus",
  "jakson_ika", "taimi_kok_lkm", "tuhon_ilmiasu1", "tuhon_aiheuttaja1")]
pjs <- merge(agg, dom, by = "koealatunnus")
pjs <- pjs[match(key$site, pjs$koealatunnus), -1]
M <- setNames(cbind(m23, pjs), paste0("m23.", c(names(m23), names(pjs))))

# === 4. Workbook: basal area / volume by species & year; biomass by year =====
w <- suppressWarnings(read_excel("Data/biomassatiedostot/kasvillisuustk_kooste_13022026.xlsx"))
w <- w[match(E$plot_id, as.numeric(w$koealatunnus_BIOSOIL)), ]
num <- function(x) suppressWarnings(as.numeric(x))
W <- as.data.frame(lapply(w[, grep("^(ppa_|v_)", names(w))], num))
for (y in c("85", "90", "95", "23")) for (cmp in c("branches", "foliage", "roots", "stem", "stump")) {
  cols <- grep(sprintf("^%s_.*_%s$", y, cmp), names(w), value = TRUE)
  W[[sprintf("bm%s_%s", y, cmp)]] <- rowSums(sapply(w[cols], num), na.rm = TRUE)
}
W <- setNames(W, paste0("wb.", names(W)))

# === 4b. Tree-tally predictors =================================================
TL <- read.csv("doublechecks/figures/tally_predictors.csv")
TL <- TL[match(E$plot_id, TL$plot_id), setdiff(names(TL), c("plot_id", "n_si_trees", "n_live85"))]
TL <- setNames(TL, paste0("tal.", names(TL)))

# === 4c. Soil chemistry and fine roots (residuals_chemistry_predictors.R) =======
CH <- read.csv("doublechecks/figures/chemistry_predictors.csv")
CH <- setNames(CH[match(E$plot_id, CH$plot_id), -1], paste0("chem.", names(CH)[-1]))

# === 4c'. GIS point extractions (GIS/ scripts; Data/GIS_points) ===============
# Wetness (TWI 16 m, DTW 2 m), EMEP forest N deposition 2010-22 (mean + trend),
# Metsäkeskus harvest declarations 1997-2026, KEMERA completed works 2003-26,
# early-20th-century history (slash-and-burn, population, railways, state forests).
GS <- read.csv("Data/GIS_points/gis_point_predictors.csv")
GS <- setNames(GS[match(E$plot_id, GS$plot_id), -1], paste0("gis.", names(GS)[-1]))

# === 4c''. Understorey cover 2023 (Data/Understorey, LUKE) =====================
# % cover per vegetation group, mean over 1-4 quadrats. File codes plots by the
# old NFI code (koealatunnus_VANHA); mapped to our ids through the site key.
UF <- list.files("Data/Understorey", pattern = "\\.xlsx$", full.names = TRUE)[1]
UN <- as.data.frame(read_excel(UF, skip = 2, col_names = FALSE, .name_repair = "minimal"))
names(UN) <- c("vanha", "vmikuvio", "n_quadrats", "kangturv", "land_type", "region",
               "maaluokka", "herbs", "graminoids", "dwarf_shrubs", "dwarf_birch",
               "raspberry", "dwarf_shrubs_plus", "mosses", "lichens")[seq_len(ncol(UN))]
UN$plot_id <- key$koealatunnus_BIOSOIL[match(UN$vanha, key$koealatunnus_VANHA)]
UN$field_layer <- UN$herbs + UN$graminoids + UN$dwarf_shrubs_plus
UN <- UN[!is.na(UN$plot_id), c("plot_id", "n_quadrats", "land_type", "herbs",
         "graminoids", "dwarf_shrubs", "dwarf_birch", "raspberry", "dwarf_shrubs_plus",
         "field_layer", "mosses", "lichens")]
UN <- setNames(UN[match(E$plot_id, UN$plot_id), -1], paste0("und.", names(UN)[-1]))
cat("Understorey 2023 cover matched for", sum(!is.na(UN$und.herbs)), "of", nrow(E), "plots\n")

# === 4d. Basal-area LEVEL and TREND over the four inventories ===================
# Mean, and slope per decade from a per-plot regression on the inventory year.
ba_y  <- c(1985, 1990, 1995, 2022)
ba_m  <- as.matrix(W[, c("wb.ppa_kaikki_85", "wb.ppa_kaikki_90", "wb.ppa_kaikki_95",
                         "wb.ppa_kaikki_mustikka")])
BT <- data.frame(
  trend.ba_mean = rowMeans(ba_m, na.rm = TRUE),
  trend.ba_slope_per_decade = apply(ba_m, 1, function(b) {
    ok <- is.finite(b); if (sum(ok) < 3) return(NA)
    10 * coef(lm(b[ok] ~ ba_y[ok]))[2] }),
  trend.ba_change_85_95 = ba_m[, 3] - ba_m[, 1],
  trend.ba_change_95_22 = ba_m[, 4] - ba_m[, 3])
BT[] <- lapply(BT, function(v) replace(v, !is.finite(v), NA))

# === 5. Assemble, flag circular, impute =====================================
X <- cbind(E[, -1], K, M, W, TL, CH, BT, GS, UN)
X$noise.random <- { set.seed(1); rnorm(nrow(X)) }
circ_pat <- "^cov\\.(ofh_|OrganicCarbon|OrganicMatter|TotalNitrogen)|hum_paks|org_kerr?_paks"   # soil-carbon measures only
circular <- grep(circ_pat, names(X), value = TRUE)
cat("Circular (soil-carbon) predictors held out of the main run:\n ",
    paste(circular, collapse = ", "), "\n")

fit <- function(d, seed, imp = "none")
  ranger(resid ~ ., d, num.trees = 1000, importance = imp, seed = seed,
         respect.unordered.factors = "order")
oob <- function(d) mean(sapply(SEEDS, function(s) fit(d, s)$r.squared))

theme_of <- function(v) {
  s <- sub("_missing$", "", v)
  case_when(
    grepl("noise", s)                                           ~ "noise",
    grepl("^und\\.", s)                                        ~ "understorey cover 2023",
    grepl("mort_|deadwood_|lahopuu", s)                          ~ "mortality & dead wood",
    grepl("site_index", s)                                      ~ "site index",
    grepl("^chem\\.", s)                                         ~ "soil chemistry & roots",
    grepl("^gis\\.(twi|dtw)", s)                                 ~ "GIS: topographic wetness",
    grepl("^gis\\.ndep", s)                                      ~ "GIS: N deposition",
    grepl("^gis\\.mki", s)                                       ~ "GIS: harvest declarations",
    grepl("^gis\\.kem", s)                                       ~ "GIS: KEMERA works",
    grepl("^gis\\.(kaski|pop1925|dist_rail|state_forest)", s)    ~ "GIS: history 1925",
    grepl("^trend\\.ba", s)                                      ~ "basal area",
    grepl("kuvion_pohja|kuvion_ppa|ppa[0-9]|\\.ppa_|ppa_tot|basal_area|pohja_", s) ~ "basal area",
    grepl("^wb\\.v_|^wb\\.bm", s)                               ~ "volume & biomass",
    grepl("harsuunt", s)                                        ~ "defoliation",
    grepl("tuhon|tuho", s)                                      ~ "damage",
    grepl("jakal|_jak$", s)                                     ~ "epiphytic lichens",
    grepl("cut_share|hakk|mets_h|metsh|maanm|maanpar|aiemp|aiem_|ehd_|ojit|toimenpide|kuvion_hist|kuvion_historia|any_cut|n_cuts|any_trt|soil_prep|tehd|tehty|hoito|maanparannus|mp_aika", s) ~ "management",
    grepl("ika|age|keskipit|keskil|valtapit|runkol|kehlk|keh_luokka|dev_class|vall_|havu_le|sivu|paapl|mean_height|tukkik|mets_laatu|laat|laad|puujak|n_storeys|puusto|pensas|taimi|lahopuu|per_tapa|perust", s) ~ "stand structure & species",
    grepl("litter|woody|conifer|sp_frac|species", s)            ~ "litter & species mix",
    grepl("temp|precip|GDD|PET|month_T|seasonality|aridity|koppen|lamposumma|kasvyo|alavyo|suovyo|eliomk", s) ~ "climate & zone",
    grepl("lat_|lon_|x_ETRS|y_ETRS|region|korkeus|kork_merenp|elevation|topograf|kalt|maanp_muoto|luonnonolot", s) ~ "location & terrain",
    grepl(circ_pat, s)                                          ~ "CIRCULAR soil carbon",
    TRUE                                                        ~ "soil & site")
}

run <- function(X, label, n_drop = 25) {
  D <- prep(X); D <- D[!is.na(D$resid), ]
  cat(sprintf("\n========== %s: %d plots x %d predictors ==========\n",
              label, nrow(D), ncol(D) - 1))
  fits <- lapply(SEEDS, fit, d = D, imp = "permutation")
  r2   <- mean(sapply(fits, `[[`, "r.squared"))
  perm <- sort(rowMeans(sapply(fits, `[[`, "variable.importance")), decreasing = TRUE)
  cat(sprintf("OOB R² = %.3f   noise column permutation rank = %d of %d\n",
              r2, match("noise.random", names(perm)), length(perm)))
  top  <- names(perm)[1:n_drop]
  dcol <- sapply(top, function(v) r2 - oob(D[, names(D) != v]))
  vt <- data.frame(variable = top, theme = theme_of(top),
                   perm = round(perm[top], 4), drop_dR2 = round(dcol, 4))
  rownames(vt) <- NULL
  cat("\nTop", n_drop, "by permutation, with drop-column loss in OOB R²:\n"); print(vt)
  th  <- theme_of(setdiff(names(D), "resid"))
  grp <- sapply(setdiff(unique(th), "noise"), function(g) {
    keep <- c("resid", setdiff(names(D), "resid")[th != g]); r2 - oob(D[, keep])
  })
  gt <- data.frame(theme = names(grp), n_vars = as.integer(table(th)[names(grp)]),
                   drop_dR2 = round(grp, 4))
  gt <- gt[order(-gt$drop_dR2), ]; rownames(gt) <- NULL
  cat("\nGroup drop (refit without the whole theme):\n"); print(gt)
  tag <- paste0(TAG, "_", gsub("[^a-z]+", "_", tolower(label)))
  write.csv(vt, file.path(OUT, sprintf("residuals_rf_everything_%s_vars.csv", tag)), row.names = FALSE)
  write.csv(gt, file.path(OUT, sprintf("residuals_rf_everything_%s_groups.csv", tag)), row.names = FALSE)
  invisible(list(r2 = r2, perm = perm))
}

# === 6. Stand SIZE grouped by inventory year =================================
# Basal area, volume and biomass of one year are one signal in many forms; they
# substitute for each other, so only a joint drop measures it.
# HIKET_RF_SIZE_ONLY=1 runs this block alone (skips the two long runs above).
size_by_year <- function(D) {
  v  <- setdiff(names(D), "resid")
  sz <- function(y) v[grepl(sprintf("^wb\\.(ppa_|v_|bm)%s|^wb\\.(ppa|v)_.*_%s$", y, y), v) |
                      grepl(sprintf("^k%s\\.(kuvion_pohja|kuvion_ppa|pohja_|ppa[0-9])", y), v) |
                      (y == "85" & v %in% c("cov.basal_area_85", "cov.basal_area_85_missing"))]
  m23 <- v[grepl("^wb\\.(ppa_|v_).*mustikka$|^wb\\.(ppa_|v_).*_23$|^wb\\.bm23|^m23\\.ppa_tot", v)]
  G <- list(`1985` = sz("85"), `1990` = sz("90"), `1995` = sz("95"), `2021-23` = m23)
  G$`all years` <- unique(unlist(G))
  G$`all except 1985` <- setdiff(G$`all years`, G$`1985`)
  r2 <- oob(D)
  out <- data.frame(group = names(G), n_vars = lengths(G),
    drop_dR2 = round(sapply(G, function(g) r2 - oob(D[, setdiff(names(D), g)])), 4))
  rownames(out) <- NULL
  cat(sprintf("\nStand SIZE by inventory year (OOB R² full = %.3f):\n", r2)); print(out)
  write.csv(out, file.path(OUT, sprintf("residuals_rf_everything_%s_size_by_year.csv", TAG)), row.names = FALSE)
}

if (Sys.getenv("HIKET_RF_SIZE_ONLY") != "1") {
  main <- run(X[, setdiff(names(X), circular)], "main")
  if (Sys.getenv("HIKET_RF_CIRC") == "1") run(X, "with circular", n_drop = 10)
  saveRDS(list(target_models = TGT, run_ids = RID[TGT], plot_id = E$plot_id,
               X = X[, setdiff(names(X), circular)], per_model = RM[match(E$plot_id, RM$plot_id), ],
               r2 = main$r2, perm = main$perm),
          file.path(OUT, sprintf("residuals_rf_everything_%s.rds", TAG)))
}
D0 <- prep(X[, setdiff(names(X), circular)]); D0 <- D0[!is.na(D0$resid), ]
size_by_year(D0)
