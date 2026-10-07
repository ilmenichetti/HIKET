# =============================================================================
# residuals_chemistry_predictors.R — soil chemistry and fine roots per plot
#
# WHY. The residual RF already carries the BioSoil 2006 MINERAL 0-20 cm
# chemistry (cov.*). Two sources were never used:
#   BioSoil 2006 ORGANIC layer (OFH) and mineral 20-40 cm (M24) chemistry
#   MUSTIKKA 2021-23 organic-layer C and N and FINE ROOTS per g soil
#     (Data/MAAT_CN_mustikkaa_250205.xlsx, sheet dataForR)
#
# LEVEL + TREND. Where both dates exist (organic-layer C:N and N) the mean of the
# two and the log change 2006 -> 2021-23 are added.
#
# CIRCULARITY. Carbon CONCENTRATIONS are left out (they are half of the observed
# stock), and so is mineral-soil N concentration, which tracks mineral C almost
# one to one. Kept: C:N ratios, organic-layer N, pH, exchangeable cations, base
# saturation, fine roots.
#
# Run from repo root:  Rscript doublechecks/residuals_chemistry_predictors.R
# Writes doublechecks/figures/chemistry_predictors.csv (plot_id = BIOSOIL id),
# read by residuals_rf_everything.R.
# =============================================================================

suppressPackageStartupMessages({library(dplyr); library(readxl)})
key <- read.csv("Data/soil_litter_site_key.csv", sep = ";")
BS  <- "../../Datasets/BioSoil_maaperäaineisto_2006"

# --- BioSoil 2006: OFH and M24 ------------------------------------------------
idx <- read.csv(file.path(BS, "plot_index.csv"))
idx$plot_id <- key$koealatunnus_BIOSOIL[match(idx$KOEALA, key$koealatunnus_VANHA)]
d <- read.csv(file.path(BS, "soil_data_main.csv"))
d$plot_id <- idx$plot_id[match(d$PlotIndex, idx$PlotIndex)]
d <- d[!is.na(d$plot_id), ]
d$CN   <- with(d, ifelse(TotalNitrogen > 0, OrganicCarbon / TotalNitrogen, NA))
d$CEC  <- with(d, ExchangeableCa + ExchangeableMg + ExchangeableK + ExchangeableNa +
                  ExchangeableAl + FreeHAcidity)
d$BS   <- with(d, (ExchangeableCa + ExchangeableMg + ExchangeableK + ExchangeableNa) / CEC)
lay <- function(code, vars, pre) {
  x <- d[d$LayerCode == code, ] |> group_by(plot_id) |>
    summarise(across(all_of(vars), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
  x[vars] <- lapply(x[vars], function(v) replace(v, !is.finite(v), NA))
  setNames(x, c("plot_id", paste0(pre, vars)))
}
ofh <- lay("OFH", c("TotalNitrogen", "CN", "pH.CaCl2.", "ExchangeableCa", "ExchangeableMg",
                    "ExchangeableK", "ExchangeableAl", "BS"), "bs06_ofh_")
m24 <- lay("M24", c("CN", "pH.CaCl2.", "BS"), "bs06_m24_")

# --- MUSTIKKA 2021-23: organic-layer C:N, N, fine roots --------------------------
num <- function(x) suppressWarnings(as.numeric(x))
mu <- suppressWarnings(suppressMessages(read_excel("Data/MAAT_CN_mustikkaa_250205.xlsx", "dataForR")))
mu <- data.frame(site = num(mu$Mustikkakoodi), N = num(mu$Nka), C = num(mu$Cka),
                 roots = num(mu$JuuretPerMaa)) |>
  mutate(CN = ifelse(N > 0, C / N, NA), roots = ifelse(roots >= 0, roots, NA)) |>
  group_by(site) |>
  summarise(mu23_N = mean(N, na.rm = TRUE), mu23_CN = mean(CN, na.rm = TRUE),
            mu23_fine_roots_per_soil = mean(roots, na.rm = TRUE), .groups = "drop")
mu[] <- lapply(mu, function(v) replace(v, !is.finite(v), NA))
mu$plot_id <- key$koealatunnus_BIOSOIL[match(mu$site, key$koealatunnus_MUSTIKKA)]
mu <- mu[!is.na(mu$plot_id), setdiff(names(mu), "site")]

# --- assemble; LEVEL and TREND where two dates exist (organic layer 2006, 2021-23)
# Trend = log ratio 2021-23 / 2006, so a constant unit difference between the two
# laboratories (MUSTIKKA N is in %, BioSoil in g/kg) only shifts it.
P <- full_join(ofh, m24[, c("plot_id", "bs06_m24_pH.CaCl2.", "bs06_m24_BS")], by = "plot_id") |>
  full_join(mu, by = "plot_id") |>
  mutate(ofh_CN_mean      = rowMeans(cbind(bs06_ofh_CN, mu23_CN), na.rm = TRUE),
         ofh_CN_logchange = log(mu23_CN / bs06_ofh_CN),
         ofh_N_logchange  = log(mu23_N / bs06_ofh_TotalNitrogen))
P[] <- lapply(P, function(v) replace(v, !is.finite(v), NA))
source("manuscript/figures/run_ids.R")
ours <- unique(readRDS(sprintf("Calibration_real_data_transient/runs/Yasso15_posterior_predictive_%s.rds",
                               RID[["Yasso15"]]))$residuals_df$plot_id)
cat("Coverage on our", length(ours), "plots:\n")
print(sapply(P[P$plot_id %in% ours, -1], function(x) sum(is.finite(x))))
print(summary(P[P$plot_id %in% ours, c("bs06_ofh_CN", "bs06_ofh_TotalNitrogen", "mu23_CN",
                                      "mu23_N", "mu23_fine_roots_per_soil")]))
cat(sprintf("C:N organic layer 2006 vs 2021-23: Spearman %.2f\n",
            with(P[P$plot_id %in% ours, ], cor(bs06_ofh_CN, mu23_CN, use = "pair", method = "s"))))
write.csv(P, "doublechecks/figures/chemistry_predictors.csv", row.names = FALSE)
cat("Wrote doublechecks/figures/chemistry_predictors.csv\n")
