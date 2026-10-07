# =============================================================================
# understorey_litter.R — plot-level understorey litter from the 2023 NFI cover
#
# WHY. sigma_input exists to add the understorey litter our tree-only product
# (Tupek) leaves out. Its anchor rests on a NATIONAL figure: 50.6 / 66.6 gC/m2/yr
# (South / North; NID Table 6.4-3, L&H 2015 Table A5, after Muukkonen et al.
# 2006). With cover measured on our own plots we can compute the same quantity
# plot by plot, by the inventory's own method, and set it against J and the
# posterior sigma_input.
#
# METHOD (the NID's chain, NID 2024 §6.4.2):
#   cover -> above-ground biomass: Muukkonen, Mäkipää, Laiho, Minkkinen, Vasander
#     & Finér 2006, Silva Fennica 40(2):231-245, upland models, fixed effects:
#       bryophytes, lichens  y = x^2 / (b0 + b1 x)^2     (Table 3)
#       dwarf shrubs, herbs & grasses  y = b1 x          (Table 4)
#     y in g dry mass m-2, x = % cover. Pine and spruce forests separately;
#     broadleaved plots use the spruce models (as Liski et al. 2006 did).
#     Lichens have a pine model only; used for all.
#   biomass -> litter: Liski et al. 2006 Table I turnover rates (yr-1):
#     bryophytes 0.33, lichens 0.1, dwarf shrubs (above) 0.25,
#     herbs & grasses (above) 1.0. Carbon = 0.5 x dry mass.
#   ABOVE-GROUND ONLY: the cover models give no below-ground biomass.
#   Dwarf shrubs = varvut+ (incl. dwarf birch and raspberry).
#
# Data: Data/Understorey (2023 cover), codes = koealatunnus_VANHA, mapped to
# our plot_id via Data/soil_litter_site_key.csv.
# Run from repo root: Rscript doublechecks/understorey_litter.R
# Writes doublechecks/figures/understorey_litter.csv
# =============================================================================

suppressPackageStartupMessages({library(readxl)})
source("manuscript/figures/run_ids.R")

# === 1. Cover, mapped to our plot ids ========================================
UF <- list.files("Data/Understorey", pattern = "\\.xlsx$", full.names = TRUE)[1]
U  <- as.data.frame(read_excel(UF, skip = 2, col_names = FALSE, .name_repair = "minimal"))
names(U) <- c("vanha", "vmikuvio", "n_quadrats", "kangturv", "land_type", "NS",
              "maaluokka", "herbs", "graminoids", "dwarf_shrubs", "dwarf_birch",
              "raspberry", "dwarf_shrubs_plus", "mosses", "lichens")
key <- read.csv("Data/soil_litter_site_key.csv", sep = ";")
dup <- key$koealatunnus_BIOSOIL[duplicated(key$koealatunnus_BIOSOIL)]   # tracts 8347/8351/8355
U$plot_id <- key$koealatunnus_BIOSOIL[match(U$vanha, key$koealatunnus_VANHA)]
U <- U[!is.na(U$plot_id) & !U$plot_id %in% dup, ]

# === 2. Calibrated plots: species, region, tree litter =======================
lm_ <- readRDS(sprintf("Data/model_inputs/Yasso15_inputs_%s.rds", RID[["Yasso15"]]))
P <- do.call(rbind, lapply(lm_$plot_info, function(p) data.frame(
  plot_id = p$plot_id, region = p$region, species = p$species_code, J = p$mean_litter)))
D <- merge(P, U, by = "plot_id")
cat("Calibrated plots:", nrow(P), " with 2023 cover:", nrow(D), "\n")
cat("region (1=S?) vs file N/S:\n"); print(table(D$region, D$NS))

# === 3. Cover -> biomass (g DM m-2) -> litter (gC m-2 yr-1) ==================
nl <- function(x, b0, b1) x^2 / (b0 + b1 * x)^2
pine <- D$species == 1
D$B_bryo   <- ifelse(pine, nl(D$mosses, 4.3369, 0.0128), nl(D$mosses, 1.8304, 0.0482))
D$B_lichen <- nl(D$lichens, 1.1833, 0.0334)
D$B_shrub  <- D$dwarf_shrubs_plus * ifelse(pine, 2.1262, 1.3169)
D$B_herb   <- (D$herbs + D$graminoids) * ifelse(pine, 0.8416, 0.6552)
TURN <- c(bryo = 0.33, lichen = 0.1, shrub = 0.25, herb = 1.0)
for (g in names(TURN)) D[[paste0("L_", g)]] <- 0.5 * TURN[[g]] * D[[paste0("B_", g)]]
D$L_und <- rowSums(D[, paste0("L_", names(TURN))])           # gC m-2 yr-1
D$u <- D$L_und / 100                                         # tC ha-1 yr-1

# === 4. Against the national figure and against J ============================
cat("\nAbove-ground understorey litter, gC m-2 yr-1 (NID: South 50.6, North 66.6):\n")
print(round(do.call(rbind, tapply(D$L_und, D$NS, function(z)
  c(n = length(z), mean = mean(z), median = median(z), q10 = quantile(z, .1), q90 = quantile(z, .9)))), 1))
cat("\nShare by group (% of mean litter):\n")
print(round(100 * colMeans(D[, paste0("L_", names(TURN))]) / mean(D$L_und), 1))
cat("\nBiomass, g DM m-2 (mean): ", paste(names(TURN), round(colMeans(D[, paste0("B_", names(TURN))])), collapse = "  "), "\n")
cat(sprintf("\nTree litter J (mean over plots): %.3f   understorey u: %.3f tC/ha/yr\n", mean(D$J), mean(D$u)))
cat(sprintf("Implied multiplier (J+u)/J on the means: %.3f\n", (mean(D$J) + mean(D$u)) / mean(D$J)))
cat("Per-plot (J+u)/J:\n"); print(round(quantile((D$J + D$u) / D$J, c(.1, .25, .5, .75, .9)), 3))
cat("cor(u, J):", round(cor(D$u, D$J), 3), "\n")

# Posterior sigma_input of the current run, for comparison
for (m in FIG_MODELS) {
  f <- sprintf("Calibration_real_data_transient/runs/%s_posterior_%s.rds", m, RID[[m]])
  if (!file.exists(f)) next
  s <- readRDS(f)                                   # draws x params, PHYSICAL space
  if ("sigma_input" %in% colnames(s))
    cat(sprintf("  %-8s posterior sigma_input median %.2f\n", m, median(s[, "sigma_input"])))
}

# National biomass the NID itself uses (§ understorey biomass: field 782, bottom
# 1534 kg DM/ha), pushed through the same turnover: what the above-ground part
# of the 50.6 / 66.6 would be on the inventory's own numbers.
cat(sprintf("\nNID national biomass x Liski turnover (above-ground): %.1f gC m-2 yr-1\n",
            0.5 * (153.4 * 0.33 + 78.2 * mean(c(0.25, 1.0)))))

write.csv(D[, c("plot_id", "region", "NS", "species", "J", paste0("B_", names(TURN)),
                paste0("L_", names(TURN)), "L_und", "u")],
          "doublechecks/figures/understorey_litter.csv", row.names = FALSE)
