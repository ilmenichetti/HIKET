# =============================================================================
# residuals_tally_predictors.R — plot predictors from the NFI tree tallies
#
# WHY. The residual follows the 1985 stand. Three candidate mechanisms need
# tree-level data the stand records do not carry:
#   MORTALITY  the Tupek litter product has no natural-mortality term, so dead
#              trees are a carbon input the models never see
#   DEAD WOOD  dead trees present on the plot at each survey
#   SITE INDEX productivity (height at a given age), to separate "productive
#              site" from "dense stand in 1985"
#   STAND AGE  real age from sample trees (the stand-record age is mostly NA)
#
# WHAT. Tallies (Data/PysyvätKoealat/<year>/lukupuu.csv), main stand
# (vmikuvio 0), trees followed by (koeala, puunro). On these permanent plots the
# tally trees are relascope-selected, so each tree stands for the same basal
# area: COUNT shares are BASAL-AREA shares.
# Dead classes: 1985 puuluokka 7/8 = natural removal (0-5 / >5 yr). 1990 data
# carry the SAME numeric coding (no letters, although the 1990 dictionary lists
# letters) — verified below by following 1985 trees. 1995: C-G = dead (standing
# or fallen, usable or not, incl. kelo). Cut: puulaji 0 in 1990/1995.
# Site index: national fit log(height) ~ species * log(total age) on dominant
# sample trees; plot index = mean residual of its dominant sample trees
# (> 0 = taller than expected for its age = more productive).
#
# Run from repo root:  Rscript doublechecks/residuals_tally_predictors.R
# Writes doublechecks/figures/tally_predictors.csv (one row per plot, keyed on
# koealatunnus_BIOSOIL = plot_id), read by residuals_rf_everything.R.
# =============================================================================

suppressPackageStartupMessages(library(dplyr))
key <- read.csv("Data/soil_litter_site_key.csv", sep = ";")

rd <- function(y) {
  t <- read.csv(sprintf("Data/PysyvätKoealat/%d/lukupuu.csv", y),
                colClasses = "character")
  t <- t[t$vmikuvio == "0", ]
  if (!"ptun" %in% names(t)) t$ptun <- t$koepuutunnus
  if (!"tilavuus" %in% names(t)) t$tilavuus <- if ("lask_tilavuus" %in% names(t)) t$lask_tilavuus else t$kuor_til
  t |> transmute(koeala, puunro, sp = puulaji, cls = puuluokka, latvus,
                 sample = ptun %in% c("2", "3", "4"),
                 h_dm = suppressWarnings(as.numeric(pituus)),
                 age = suppressWarnings(as.numeric(rinn_kork_ika)) +
                       suppressWarnings(as.numeric(ikalisays)),
                 vol = suppressWarnings(as.numeric(tilavuus)))
}
t85 <- rd(1985); t90 <- rd(1990); t95 <- rd(1995)

SPP   <- as.character(1:9)
dead85 <- function(cls) cls %in% c("7", "8")                 # 1985 & 1990 coding
dead95 <- function(cls) cls %in% c("C", "D", "E", "F", "G")
t85$status <- ifelse(dead85(t85$cls), "dead", ifelse(t85$sp %in% SPP, "live", "other"))
t90$status <- ifelse(t90$sp == "0", "cut", ifelse(dead85(t90$cls), "dead",
                     ifelse(t90$sp %in% SPP, "live", "other")))
t95$status <- ifelse(t95$sp == "0", "cut", ifelse(dead95(t95$cls), "dead",
                     ifelse(t95$sp %in% SPP, "live", "other")))

# --- verify the 1990 coding: 1990 "dead" trees should be 1985 live trees -----
v <- merge(t90[t90$status == "dead", c("koeala", "puunro")],
           t85[, c("koeala", "puunro", "status")], by = c("koeala", "puunro"))
cat("1990 trees coded dead (7/8): status of the same tree in 1985:\n"); print(table(v$status))

# --- follow 1985 live trees ---------------------------------------------------
fate <- t85[t85$status == "live", c("koeala", "puunro", "vol")] |>
  left_join(t90[, c("koeala", "puunro", "status")] |> rename(s90 = status), by = c("koeala", "puunro")) |>
  left_join(t95[, c("koeala", "puunro", "status")] |> rename(s95 = status), by = c("koeala", "puunro")) |>
  mutate(died = s90 %in% "dead" | s95 %in% "dead",
         cut  = !died & (s90 %in% "cut" | s95 %in% "cut"))
mort <- fate |> group_by(koeala) |>
  summarise(n_live85        = n(),
            mort_share_85_95 = mean(died),               # basal-area share that died
            mort_vol_share_85_95 = sum(vol[died], na.rm = TRUE) / sum(vol, na.rm = TRUE),
            cut_share_85_95  = mean(cut), .groups = "drop")

# --- dead wood present at each survey (share of tally trees that are dead) ----
dw <- function(t, y) t |> filter(status %in% c("live", "dead")) |> group_by(koeala) |>
  summarise(!!paste0("deadwood_share_", y) := mean(status == "dead"), .groups = "drop")

# --- stand age and site index from sample trees -------------------------------
st <- bind_rows(mutate(t85, year = 1985), mutate(t90, year = 1990), mutate(t95, year = 1995)) |>
  filter(status == "live", sample, latvus %in% c("B", "Y"), is.finite(h_dm), h_dm > 13,
         is.finite(age), age > 5, sp %in% SPP) |>
  mutate(spg = ifelse(sp %in% c("1", "2"), sp, "broadleaf"))
fit <- lm(log(h_dm) ~ spg * log(age), st)
cat(sprintf("\nSite-index fit: %d dominant sample trees, R² %.2f\n", nrow(st), summary(fit)$r.squared))
st$si <- resid(fit)
si  <- st |> group_by(koeala) |> summarise(site_index = mean(si), n_si_trees = n(), .groups = "drop")
age <- st |> group_by(koeala, year) |> summarise(a = median(age), .groups = "drop") |>
  mutate(a85 = a - (year - 1985)) |> group_by(koeala) |>
  summarise(stand_age_85_tally = a85[which.min(abs(year - 1985))], .groups = "drop")

# --- assemble ----------------------------------------------------------------
P <- data.frame(koeala = as.character(unique(c(t85$koeala, t90$koeala, t95$koeala)))) |>
  left_join(mort, by = "koeala") |> left_join(dw(t85, 1985), by = "koeala") |>
  left_join(dw(t90, 1990), by = "koeala") |> left_join(dw(t95, 1995), by = "koeala") |>
  left_join(si, by = "koeala") |> left_join(age, by = "koeala")
P$plot_id <- key$koealatunnus_BIOSOIL[match(as.integer(P$koeala), key$koealatunnus_VANHA)]
P <- P[!is.na(P$plot_id), ]; P$koeala <- NULL

source("manuscript/figures/run_ids.R")
ours <- unique(readRDS(sprintf("Calibration_real_data_transient/runs/Yasso15_posterior_predictive_%s.rds",
                                 RID[["Yasso15"]]))$residuals_df$plot_id)
cat("\nCoverage on our", length(ours), "plots:\n")
print(sapply(P[P$plot_id %in% ours, setdiff(names(P), "plot_id")], function(x) sum(is.finite(x))))
print(summary(P[P$plot_id %in% ours, c("mort_share_85_95", "cut_share_85_95", "deadwood_share_1985",
                                      "deadwood_share_1995", "site_index", "stand_age_85_tally")]))
write.csv(P, "doublechecks/figures/tally_predictors.csv", row.names = FALSE)
cat("\nWrote doublechecks/figures/tally_predictors.csv\n")
