# =============================================================================
# organic_mineral_boundary.R — where does the apparent SOC change actually sit,
# once the litter definition is taken out, and is it a boundary artefact?
#
# Motivation. Removing the litter term leaves the humus layer with almost no
# 1985->2006 trend (+3% C). So does the accumulation disappear, or does it move
# somewhere else? And is any of it just the ORGANIC/MINERAL boundary shifting?
#
# The boundary between the organic layer and mineral soil is set BY EYE in the
# field (Kramarenko 2012 sec. 4.4; Tamminen 1999 found a quarter of NFI crews'
# organic thicknesses differed >2 cm from the soil-research crews', against a
# mean thickness of 4.3 cm). If the boundary is placed deeper in one campaign,
# carbon is reassigned from "mineral 0-20" to "organic" with no carbon moving
# at all. Kramarenko saw a single plot at organic -381 / mineral 0-20 +283
# g m-2 yr-1 and attributed it to exactly this.
#
# THE TEST. Per plot, correlate the organic-layer change against the mineral
# 0-20 cm change over the same interval.
#   process       -> both bands respond to the same drivers => corr >= 0
#   boundary shift-> one gains what the other loses         => corr << 0
# Also check the DILUTION signature: a deeper boundary adds mass at lower C%.
#
# Run from repo root:  Rscript doublechecks/organic_mineral_boundary.R
# =============================================================================

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readxl) })
options(width = 130)

ROOT   <- if (dir.exists("Data/SOC_homogeneized")) "." else ".."
SOCDIR <- file.path(ROOT, "Data/SOC_homogeneized")
OUTDIR <- file.path(ROOT, "doublechecks/soc_depth"); dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

layers <- read.csv(file.path(SOCDIR, "soc_homogenized_layers.csv"), stringsAsFactors = FALSE)
plotd  <- read.csv(file.path(SOCDIR, "soc_homogenized_plot.csv"),   stringsAsFactors = FALSE)
CAMPS  <- c("VMI8", "Biosoil", "Komeetta")

grid <- layers |>
  mutate(band = case_when(layer == "organic" ~ "organic",
                          layer %in% c("0-5cm","5-20cm","0-10cm","10-20cm") ~ "m0_20",
                          layer == "20-40cm" ~ "m20_40", TRUE ~ NA_character_)) |>
  filter(!is.na(band)) |>
  group_by(plot_id, campaign, band) |>
  summarise(C = sum(C_Mgha, na.rm = TRUE), n_src = n(), .groups = "drop") |>
  mutate(C = ifelse(band == "m0_20" & n_src < 2, NA_real_, C)) |> select(-n_src) |>
  pivot_wider(names_from = band, values_from = C) |>
  left_join(plotd |> select(plot_id, campaign, weight, soc_deep_Mgha, soc_profile_Mgha, soc_outlier),
            by = c("plot_id", "campaign")) |>
  filter(!(soc_outlier %in% TRUE))

bal_ids <- grid |> filter(!is.na(organic), !is.na(m0_20), !is.na(m20_40)) |>
  count(plot_id) |> filter(n == 3) |> pull(plot_id)
w <- grid |> filter(plot_id %in% bal_ids) |>
  select(plot_id, campaign, weight, organic, m0_20, m20_40, soc_deep_Mgha, soc_profile_Mgha) |>
  pivot_wider(names_from = campaign, values_from = c(organic, m0_20, m20_40, soc_deep_Mgha, soc_profile_Mgha))
cat("balanced panel:", length(bal_ids), "plots\n")
wm <- function(x, wt) sum(x * wt, na.rm = TRUE) / sum(wt[!is.na(x)], na.rm = TRUE)

# --- litter term, to build the like-for-like organic layer -------------------
XLSX <- file.path(ROOT, "Data/Komeetta/Hannu/Komeetta 150526hi--.xlsx")
sp <- suppressMessages(read_excel(XLSX, sheet = "BiSo", range = "A4:IN3031", col_names = FALSE,
                                  col_types = "text", .name_repair = "minimal"))
nn <- function(i) suppressWarnings(as.numeric(sp[[i]]))
lit <- tibble(plot_id = as.integer(nn(3)), Krs = as.integer(nn(4)), REP = as.integer(nn(9)),
              lit06 = nn(247) / 1000, lit24 = nn(248) / 1000) |>
  filter(!is.na(plot_id), REP == 1, Krs == 101) |> select(plot_id, lit06, lit24)
w <- w |> left_join(lit, by = "plot_id") |>
  mutate(lit06 = coalesce(lit06, 0), lit24 = coalesce(lit24, 0),
         humus_VMI8 = organic_VMI8, humus_Biosoil = organic_Biosoil - lit06,
         humus_Komeetta = organic_Komeetta - lit24)

# =============================================================================
# 1. WHERE THE CHANGE SITS, like-for-like (humus, not organic+litter)
# =============================================================================
cat("\n=== 1. LIKE-FOR-LIKE CHANGE BY BAND (weighted, Mg C/ha) ===\n")
bands <- c(humus = "humus", `mineral 0-20` = "m0_20", `mineral 20-40` = "m20_40",
           `deep tail` = "soc_deep_Mgha")
tab <- lapply(names(bands), function(nm) {
  v <- bands[[nm]]
  d1 <- w[[paste0(v, "_Biosoil")]]  - w[[paste0(v, "_VMI8")]]
  d2 <- w[[paste0(v, "_Komeetta")]] - w[[paste0(v, "_Biosoil")]]
  data.frame(band = nm, d_85_06 = wm(d1, w$weight), d_06_24 = wm(d2, w$weight))
}) |> bind_rows()
tab <- rbind(tab, data.frame(band = "TOTAL (like-for-like)",
                             d_85_06 = sum(tab$d_85_06), d_06_24 = sum(tab$d_06_24)))
tab$rate_85_06 <- tab$d_85_06 / 21; tab$rate_06_24 <- tab$d_06_24 / 18
print(as.data.frame(tab |> mutate(across(where(is.numeric), \(x) round(x, 3)))), row.names = FALSE)
cat("\n  => the accumulation does NOT vanish with the litter term removed. It RELOCATES:\n")
cat("     the humus layer goes flat, and essentially all of the 1985->2006 sink ends up\n")
cat("     in the mineral soil — the compartment with the documented method problem.\n")

# =============================================================================
# 2. IS IT THE ORGANIC/MINERAL BOUNDARY MOVING?
# =============================================================================
cat("\n=== 2. BOUNDARY TEST: corr(delta humus, delta mineral 0-20) per plot ===\n")
cat("    process >= 0   |   boundary reallocation << 0\n\n")
bt <- function(a, b, lab) {
  dh <- w[[paste0("humus_", b)]] - w[[paste0("humus_", a)]]
  dm <- w[[paste0("m0_20_", b)]] - w[[paste0("m0_20_", a)]]
  ok <- is.finite(dh) & is.finite(dm)
  ct <- cor.test(dh[ok], dm[ok])
  cat(sprintf("  %-18s n=%3d  r = %+0.3f  [%+0.3f,%+0.3f]  p = %.2g   spearman %+0.3f\n",
              lab, sum(ok), ct$estimate, ct$conf.int[1], ct$conf.int[2], ct$p.value,
              cor(dh[ok], dm[ok], method = "spearman")))
  # how much of each band's variance is the pure trade-off?
  cat(sprintf("  %-18s   sd(humus)=%.2f sd(min0_20)=%.2f  sd(SUM)=%.2f  <- if it were a pure\n",
              "", sd(dh[ok]), sd(dm[ok]), sd(dh[ok] + dm[ok])))
  cat(sprintf("  %-18s   trade-off, sd(SUM) would collapse below both (quadrature = %.2f)\n",
              "", sqrt(sd(dh[ok])^2 + sd(dm[ok])^2)))
  invisible(NULL)
}
bt("VMI8", "Biosoil", "1985 -> 2006"); bt("Biosoil", "Komeetta", "2006 -> 2024")

# control: a band pair with NO shared boundary should not show the effect
dh <- w$humus_Komeetta - w$humus_Biosoil; dd <- w$m20_40_Komeetta - w$m20_40_Biosoil
ok <- is.finite(dh) & is.finite(dd)
cat(sprintf("\n  CONTROL humus vs mineral 20-40 (no shared boundary), 2006->2024: r = %+0.3f\n",
            cor(dh[ok], dd[ok])))

# =============================================================================
# 3. DILUTION SIGNATURE — a deeper boundary adds mass at lower C%
# =============================================================================
cat("\n=== 3. DILUTION: does the organic layer gain mass while losing C%? ===\n")
b2 <- tibble(plot_id = as.integer(nn(3)), Krs = as.integer(nn(4)), REP = as.integer(nn(9)),
             OLW = nn(22), Cpct85 = nn(159), Cpct06 = nn(165)) |>
  filter(!is.na(plot_id), REP == 1, Krs == 101)
d85 <- suppressMessages(read_excel(XLSX, sheet = "Data_1985", .name_repair = "minimal")) |> as.data.frame()
a85 <- tibble(plot_id = as.integer(d85$VMI), C85 = as.numeric(d85$Corg_kgha),
              Cp85 = as.numeric(d85$`Corg%`)) |> filter(!is.na(plot_id))
mm <- inner_join(a85, b2, by = "plot_id") |>
  mutate(mass85 = C85 / (Cp85 / 100), mass06 = OLW * 1e4) |>
  filter(is.finite(mass85), is.finite(mass06), mass85 > 0, mass06 > 0, is.finite(Cpct06))
cat(sprintf("  1985 -> 2006 humus:  mass %+5.1f %%   C%% %+5.2f pp   carbon %+5.1f %%   (n=%d)\n",
            100 * (median(mm$mass06) / median(mm$mass85) - 1),
            median(mm$Cpct06) - median(mm$Cp85),
            100 * (median(mm$mass06 * mm$Cpct06 / 100) / median(mm$C85) - 1), nrow(mm)))
dm <- mm$mass06 - mm$mass85; dc <- mm$Cpct06 - mm$Cp85
cat(sprintf("  per plot, corr(delta mass, delta C%%) = %+0.3f  (p = %.2g)\n",
            cor(dm, dc), cor.test(dm, dc)$p.value))
cat("\n  !! DO NOT READ THIS AS DILUTION. The test is confounded BY CONSTRUCTION: the 1985\n")
cat("  mass is not measured, it is derived as C85/(C%85). An over-estimated 1985 C% lowers\n")
cat("  mass85 and raises delta-mass while lowering delta-C%, i.e. it manufactures exactly\n")
cat("  this negative correlation. With only a stock and a concentration for 1985, mass and\n")
cat("  concentration cannot be separated, so signature 3 is INCONCLUSIVE.\n")
cat("\n  The same confound attacks the +14.5% mass gain itself: Kramarenko sec. 4.4 records\n")
cat("  that campaign 1's C% was largely LOI-regression PREDICTED. A C% biased high by the\n")
cat("  observed ~2 pp (~5% relative) biases mass85 low by ~5% and the 1985 organic STOCK\n")
cat("  high by ~5% (~0.9 Mg/ha) if the mass was the measured quantity.\n")
cat("  Robust statement that survives either way: the humus layer's 1985->2006 change is\n")
cat("  within about +/-1 Mg/ha over 21 yr, i.e. |rate| <= 0.05 Mg/ha/yr against a\n")
cat("  like-for-like total of +0.209. It is small however the C% question resolves.\n")

png(file.path(OUTDIR, "03_boundary_tradeoff.png"), width = 1400, height = 700, res = 130)
par(mfrow = c(1, 2), mar = c(4.3, 4.3, 3, 1))
for (p in list(c("VMI8","Biosoil"), c("Biosoil","Komeetta"))) {
  dh <- w[[paste0("humus_", p[2])]] - w[[paste0("humus_", p[1])]]
  dmn <- w[[paste0("m0_20_", p[2])]] - w[[paste0("m0_20_", p[1])]]
  ok <- is.finite(dh) & is.finite(dmn)
  plot(dh[ok], dmn[ok], pch = 16, col = "#2E6F9E88", cex = .6,
       xlab = expression(Delta * " humus (Mg C/ha)"), ylab = expression(Delta * " mineral 0-20 (Mg C/ha)"),
       main = sprintf("%s -> %s   r = %+0.2f", p[1], p[2], cor(dh[ok], dmn[ok])))
  abline(h = 0, v = 0, col = "grey70"); abline(lm(dmn[ok] ~ dh[ok]), col = "#C1553B", lwd = 2)
  abline(0, -1, lty = 2)
  legend("topright", c("fit", "pure trade-off (slope -1)"), col = c("#C1553B", "black"),
         lty = c(1, 2), lwd = c(2, 1), bty = "n", cex = .75)
}
dev.off()
cat("\nfigure:", file.path(OUTDIR, "03_boundary_tradeoff.png"), "\n")
