# =============================================================================
# subsoil_offset_mechanism.R — is the 1985->2006 mineral offset a CARBON signal
# or a multiplicative measurement/correction artefact (bulk density, stoniness,
# sampled volume)?
#
# The pattern to explain: over 1985->2006 the mineral 20-40 cm layer gains
# +20.5% while the 0-20 cm layer gains only +5.4%, and the humus layer is flat.
# Then over 2006->2024 the subsoil LOSES carbon. A real accumulation signal is
# driven by litter input and must be SURFACE-weighted; this one is
# DEPTH-weighted, which is backwards.
#
# Candidate mechanisms are all multiplicative on the stock:
#   stock = concentration x bulk density x thickness x (1 - coarse fraction)
# so an error in BD, in the stoniness correction, or in the sampled volume
# (Kramarenko 2012 sec. 4.4: mineral samples were mostly cut from the wall of a
# spade-dug pit, "not technically possible to represent the whole sampling
# depth evenly") rescales the layer rather than adding carbon to it.
#
# FOUR DISCRIMINATING TESTS
#   T1 multiplicative vs additive: is the offset proportional to the plot's own
#      stock (=> a rescaling) or a roughly constant carbon amount (=> a flux)?
#   T2 depth dependence: does the rescaling grow with depth, as stoniness and
#      spade-sampling difficulty do?
#   T3 stoniness: does the discrepancy scale with the plot's coarse fraction?
#      This is the specific fingerprint of a stoniness-correction difference.
#   T4 whole-profile coupling: if one plot-level factor (BD, volume) is wrong,
#      the 0-20 and 20-40 log-ratios should move TOGETHER within a plot.
#
# Run from repo root:  Rscript doublechecks/subsoil_offset_mechanism.R
# =============================================================================

suppressPackageStartupMessages({ library(dplyr); library(tidyr) })
options(width = 130)

ROOT   <- if (dir.exists("Data/SOC_homogeneized")) "." else ".."
SOCDIR <- file.path(ROOT, "Data/SOC_homogeneized")
layers <- read.csv(file.path(SOCDIR, "soc_homogenized_layers.csv"), stringsAsFactors = FALSE)
plotd  <- read.csv(file.path(SOCDIR, "soc_homogenized_plot.csv"),   stringsAsFactors = FALSE)

hdr <- function(x) cat("\n", strrep("=", 112), "\n", x, "\n", strrep("=", 112), "\n", sep = "")

g <- layers |>
  mutate(b = case_when(layer == "organic" ~ "organic",
                       layer %in% c("0-5cm","5-20cm","0-10cm","10-20cm") ~ "m0_20",
                       layer == "20-40cm" ~ "m20_40", TRUE ~ NA_character_)) |>
  filter(!is.na(b)) |>
  group_by(plot_id, campaign, b) |>
  summarise(C = sum(C_Mgha, na.rm = TRUE), n_src = n(), .groups = "drop") |>
  mutate(C = ifelse(b == "m0_20" & n_src < 2, NA_real_, C)) |> select(-n_src) |>
  pivot_wider(names_from = b, values_from = C) |>
  left_join(plotd |> select(plot_id, campaign, weight, coarse_frag_2006, bd_est_2006,
                            clay_2006, basal_area_1985, soc_outlier),
            by = c("plot_id", "campaign")) |>
  filter(!(soc_outlier %in% TRUE))

bal <- g |> filter(!is.na(organic), !is.na(m0_20), !is.na(m20_40)) |>
  count(plot_id) |> filter(n == 3) |> pull(plot_id)
W <- g |> filter(plot_id %in% bal) |>
  select(plot_id, campaign, m0_20, m20_40, coarse_frag_2006, bd_est_2006, clay_2006, basal_area_1985) |>
  pivot_wider(names_from = campaign, values_from = c(m0_20, m20_40)) |>
  distinct(plot_id, .keep_all = TRUE) |>
  filter(m20_40_VMI8 > 0, m20_40_Biosoil > 0, m0_20_VMI8 > 0, m0_20_Biosoil > 0) |>
  mutate(d_sub  = m20_40_Biosoil - m20_40_VMI8,
         lr_sub = log(m20_40_Biosoil / m20_40_VMI8),
         lr_top = log(m0_20_Biosoil  / m0_20_VMI8),
         mean_sub = (m20_40_Biosoil + m20_40_VMI8) / 2)
cat("plots:", nrow(W), "\n")

# --- T1 multiplicative vs additive -------------------------------------------
hdr("T1. IS THE OFFSET PROPORTIONAL TO THE STOCK (rescaling) OR CONSTANT (a flux)?")
f_mult <- lm(d_sub ~ mean_sub, data = W)     # slope>0 => bigger stocks shift more
f_add  <- lm(d_sub ~ 1, data = W)
cat(sprintf("  delta ~ level : slope = %+0.3f  [%+0.3f, %+0.3f]  p = %.2g\n",
            coef(f_mult)[2], confint(f_mult)[2,1], confint(f_mult)[2,2],
            summary(f_mult)$coefficients[2,4]))
cat(sprintf("  intercept     = %+0.3f Mg/ha   (a pure flux would put everything here)\n", coef(f_mult)[1]))
cat(sprintf("  R2 of the proportional term = %.3f ; AIC  level-model %.1f  vs  constant %.1f\n",
            summary(f_mult)$r.squared, AIC(f_mult), AIC(f_add)))
cat("  => a positive slope with a near-zero intercept is the signature of a RESCALING,\n")
cat("     not of carbon being added at a rate independent of how much is there.\n")

# --- T2 depth dependence ------------------------------------------------------
hdr("T2. DOES THE RESCALING GROW WITH DEPTH?")
wm <- function(x, w) sum(x * w, na.rm = TRUE) / sum(w[!is.na(x)], na.rm = TRUE)
cat(sprintf("  mineral  0-20 : median ratio 2006/1985 = %.3f\n", median(W$m0_20_Biosoil / W$m0_20_VMI8)))
cat(sprintf("  mineral 20-40 : median ratio 2006/1985 = %.3f\n", median(W$m20_40_Biosoil / W$m20_40_VMI8)))
cat(sprintf("  paired difference of log-ratios (sub - top) = %+0.3f  (wilcox p = %.2g)\n",
            median(W$lr_sub - W$lr_top), wilcox.test(W$lr_sub, W$lr_top, paired = TRUE)$p.value))
cat("  Stoniness and spade-sampling error both INCREASE with depth; litter input does not.\n")

# --- T3 stoniness fingerprint -------------------------------------------------
hdr("T3. DOES THE DISCREPANCY SCALE WITH STONINESS?  (the specific fingerprint)")
for (v in c("coarse_frag_2006", "bd_est_2006", "clay_2006", "basal_area_1985")) {
  ok <- is.finite(W[[v]]) & is.finite(W$lr_sub)
  if (sum(ok) < 30) { cat(sprintf("  %-18s too few plots\n", v)); next }
  ct <- cor.test(W[[v]][ok], W$lr_sub[ok], method = "spearman", exact = FALSE)
  cat(sprintf("  log-ratio(20-40) vs %-18s  rho = %+0.3f  p = %.3g   (n=%d)\n",
              v, ct$estimate, ct$p.value, sum(ok)))
}
cat("\n  coarse_frag / bd are the MEASUREMENT-correction terms; clay and basal area are\n")
cat("  the biological controls. A stoniness artefact loads on the first pair only.\n")

# --- T4 whole-profile coupling ------------------------------------------------
hdr("T4. DO THE TWO MINERAL LAYERS MOVE TOGETHER WITHIN A PLOT?")
ct <- cor.test(W$lr_top, W$lr_sub)
cat(sprintf("  corr( log-ratio 0-20 , log-ratio 20-40 ) = %+0.3f  [%+0.3f, %+0.3f]  p = %.2g\n",
            ct$estimate, ct$conf.int[1], ct$conf.int[2], ct$p.value))
cat("  A single plot-level factor (bulk density, sampled volume, a mis-set stoniness)\n")
cat("  would rescale BOTH mineral layers together => strong positive correlation.\n")
cat("  Layer-specific causes (spade tilt within one pit, a per-layer C% equation)\n")
cat("  would not.\n")

hdr("READ-OUT")
cat("The four tests together say whether to treat the mineral offset as a plot-level\n")
cat("rescaling (correctable in principle, if the correction terms can be recovered)\n")
cat("or as an irreducible campaign-level uncertainty to be carried in the error model.\n")
