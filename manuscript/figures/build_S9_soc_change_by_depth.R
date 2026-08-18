setwd("/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling")
# Appendix: WHERE the apparent SOC change sits by depth, and how much the
# unresolved litter definition moves the observed trend.
#
# ONE claim per panel.
#   (a),(b) The change by layer over each inter-campaign interval. The humus
#           layer (OFH) is flat over 1985->2006 and DECLINES over 2006->2024,
#           while the mineral soil carries the entire apparent sink and the
#           subsoil changes sign between the two intervals. This holds whichever
#           litter treatment is chosen, and it is what motivates treating the
#           mineral bands as the weak link (Kramarenko 2012 sec. 4.4: mineral
#           samples are spade-cut, not volumetric).
#   (c)     What the litter treatment does to the trend, as a bracket.
#
# LAYERS ARE KEPT SEPARATE rather than lumped into one "organic" band, so the
# litter contribution is visible instead of being folded into the humus bar:
#   OFH = humus (Krs 101)     LM = litter+moss (Krs 100, a separately coded layer)
#
# TREATMENTS. The 1985 sheet carries no LM; 2006 and 2024 do.
#   A  as built        1985 without LM, 2006/24 with  -- internally INCONSISTENT
#   B  LM removed      from all three campaigns
#   C  LM added to 1985 (that plot's own 2006 value)  -- the working default:
#      it repairs only the campaign that is broken and keeps the official LUKE
#      stock basis for 2006/24. Panels (a),(b) are on basis C.
# Under C the 1985->2006 LM change is zero BY CONSTRUCTION, and is drawn hollow
# to say so.
#
# Balanced panel (measured in all three campaigns), region-weighted, documented
# outliers dropped.

source("manuscript/figures/model_palette.R")
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readxl) })

SOCDIR <- "Data/SOC_homogeneized"
XLSX   <- "Data/Komeetta/Hannu/Komeetta 150526hi--.xlsx"
layers <- read.csv(file.path(SOCDIR, "soc_homogenized_layers.csv"), stringsAsFactors = FALSE)
plotd  <- read.csv(file.path(SOCDIR, "soc_homogenized_plot.csv"),   stringsAsFactors = FALSE)
CA <- CAMPAIGN_ORDER

meta <- plotd |> select(plot_id, campaign, weight, soc_deep_Mgha, soc_outlier, lm_added_1985)
band <- layers |>
  mutate(b = case_when(layer == "organic" ~ "organic",
                       layer %in% c("0-5cm","5-20cm","0-10cm","10-20cm") ~ "m0_20",
                       layer == "20-40cm" ~ "m20_40", TRUE ~ NA_character_)) |>
  filter(!is.na(b)) |>
  group_by(plot_id, campaign, b) |>
  summarise(C = sum(C_Mgha, na.rm = TRUE), n_src = n(), .groups = "drop") |>
  mutate(C = ifelse(b == "m0_20" & n_src < 2, NA_real_, C)) |> select(-n_src) |>
  pivot_wider(names_from = b, values_from = C) |>
  left_join(meta, by = c("plot_id", "campaign")) |> filter(!(soc_outlier %in% TRUE))
bal_ids <- band |> filter(!is.na(organic), !is.na(m0_20), !is.na(m20_40)) |>
  count(plot_id) |> filter(n == 3) |> pull(plot_id)
band <- band |> filter(plot_id %in% bal_ids)

sp <- suppressMessages(read_excel(XLSX, sheet = "BiSo", range = "A4:IN3031", col_names = FALSE,
                                  col_types = "text", .name_repair = "minimal"))
nn  <- function(i) suppressWarnings(as.numeric(sp[[i]]))
lit <- tibble(plot_id = as.integer(nn(3)), Krs = as.integer(nn(4)), REP = as.integer(nn(9)),
              Biosoil = nn(247) / 1000, Komeetta = nn(248) / 1000) |>
  filter(!is.na(plot_id), REP == 1, Krs == 101) |> select(plot_id, Biosoil, Komeetta) |>
  pivot_longer(-plot_id, names_to = "campaign", values_to = "LM")
band <- band |> left_join(lit, by = c("plot_id", "campaign")) |>
  # ⚠ 2026-08-12: the baseline now applies TREATMENT C — the 1985 organic layer
# already contains the imputed LM. `lm_added_1985` (kg/ha, 0 for 2006/2024) records
# how much, so OFH is recovered by subtracting it. Without this the 1985 "humus"
# would silently include the litter and the whole comparison would invert.
  mutate(LM = coalesce(LM, coalesce(lm_added_1985, 0) / 1000),
         OFH = organic - LM)

W <- band |> select(plot_id, campaign, weight, OFH, LM, m0_20, m20_40, soc_deep_Mgha) |>
  pivot_wider(names_from = campaign, values_from = c(OFH, LM, m0_20, m20_40, soc_deep_Mgha))
W$LM_VMI8_C <- W$LM_Biosoil        # treatment C: carry the plot's own 2006 LM back to 1985
wm <- function(x, w) sum(x * w, na.rm = TRUE) / sum(w[!is.na(x)], na.rm = TRUE)
cat(sprintf("balanced panel: %d plots\n", nrow(W)))

# True mean interval, not the nominal 39 yr: the first campaign was sampled
# 1986-1995 (mean ~1989), so 39 overstates it by ~11%. samp_year is in the baseline.
.sy <- plotd$samp_year[plotd$campaign == "VMI8"]; .sy <- .sy[is.finite(.sy)]
YRS_85_24 <- 2024 - mean(.sy)
YRS_85_06 <- 2006 - mean(.sy)   # 17.1 yr, NOT the nominal 21
YRS_06_24 <- 18
VMI8_YR   <- round(mean(.sy))   # ~1989; F4 labels the campaign the same way

BANDS <- c("OFH", "LM", "m0_20", "m20_40", "soc_deep_Mgha")
BLAB  <- c("humus\nOFH", "litter\nLM", "mineral\n0-20 cm", "mineral\n20-40 cm", "modelled\ntail")
BCOL  <- c("#6E8B74", "#C9B458", "#4B7FA8", "#2E4E63", "#9AA0A6")
d1 <- sapply(BANDS, function(v) {
  a <- if (v == "LM") "LM_VMI8_C" else paste0(v, "_VMI8")
  wm(W[[paste0(v, "_Biosoil")]] - W[[a]], W$weight) })
d2 <- sapply(BANDS, function(v) wm(W[[paste0(v, "_Komeetta")]] - W[[paste0(v, "_Biosoil")]], W$weight))

# treatment levels for panel (c)
lv <- function(orgfun) sapply(CA, function(cc)
  wm(orgfun(cc) + W[[paste0("m0_20_", cc)]] + W[[paste0("m20_40_", cc)]] +
     W[[paste0("soc_deep_Mgha_", cc)]], W$weight))
levA <- lv(function(cc) W[[paste0("OFH_", cc)]] + W[[paste0("LM_", cc)]])
levB <- lv(function(cc) W[[paste0("OFH_", cc)]])
levC <- lv(function(cc) W[[paste0("OFH_", cc)]] +
             (if (cc == "VMI8") W$LM_VMI8_C else W[[paste0("LM_", cc)]]))
YRS  <- c(VMI8_YR, 2006, 2024)   # first campaign at its MEAN sampling year (~1989), not 1985
TR   <- list(A = levA, B = levB, C = levC)
TLAB <- c(A = "A  as built (inconsistent)", B = "B  litter removed from all",
          C = "C  litter added to 1985  [default]")
TCOL <- c(A = "#B0B0B0", B = "#7FA8C4", C = "#1F4E79")
TLTY <- c(A = 2, B = 1, C = 1)

# =============================================================================
png("manuscript/figures/S9_soc_change_by_depth.png",
    width = 12, height = 4.8, units = "in", res = 200)
par(mfrow = c(1, 3), mar = c(4.6, 4.5, 3.9, 1.0), mgp = c(2.5, 0.7, 0), las = 1)
sub <- function(txt) mtext(txt, side = 3, line = 0.35, cex = 0.62, col = "grey35")

bandpanel <- function(d, yrs, ttl, note, hollow = NULL) {
  yl <- range(0, d) + c(-0.45, 0.65)
  fill <- BCOL; if (!is.null(hollow)) fill[hollow] <- "white"
  bp <- barplot(d, col = fill, border = BCOL, names.arg = BLAB, ylim = yl,
                ylab = sprintf("Change over %.0f yr (Mg C/ha)", yrs), main = ttl, cex.names = 0.8)
  abline(h = 0, col = "grey30")
  text(bp, d, sprintf("%+.2f", d), pos = ifelse(d >= 0, 3, 1), cex = 0.72, offset = 0.28, font = 2)
  text(mean(bp[3:5]), yl[2] * 0.88, sprintf("total %+.2f", sum(d)), cex = 0.78, font = 3, col = "grey25")
  sub(note); invisible(bp)
}
bp1 <- bandpanel(d1, YRS_85_06, sprintf("(a)  VMI8 (mean %d) -> 2006, by layer", VMI8_YR),
                 "humus flat; the mineral soil carries the whole sink", hollow = 2)
text(bp1[2], 0, "zero by\nconstruction", pos = 3, cex = 0.58, col = "grey45", font = 3, offset = 1.9)
bandpanel(d2, YRS_06_24, "(b)  2006 -> 2024, by layer",
          "humus now DECLINES; the subsoil reverses sign")

## (c) the litter treatment as a bracket ---------------------------------------
plot(NA, xlim = c(1981, 2029), ylim = range(unlist(TR)) + c(-2.2, 3.4),
     xlab = "Year", ylab = "Mean whole-profile SOC (Mg C/ha)",
     main = "(c)  Litter treatment effect")
for (k in names(TR)) {
  y <- TR[[k]]
  lines(YRS, y, col = TCOL[k], lwd = 2.6, lty = TLTY[k])
  points(YRS, y, pch = 21, bg = "white", col = TCOL[k], cex = 1.1, lwd = 2)
  text(2024.8, y[3], sprintf(" %.1f", y[3]), adj = 0, cex = 0.72, col = TCOL[k], font = 2, xpd = NA)
  text(1984.2, y[1], sprintf("%.1f ", y[1]), adj = 1, cex = 0.72, col = TCOL[k], font = 2, xpd = NA)
}
legend("topleft", lwd = 2.6, lty = TLTY[names(TLAB)], col = TCOL[names(TLAB)],
       legend = sprintf("%s   %+.3f Mg/ha/yr", TLAB, sapply(TR, function(y) (y[3] - y[1]) / YRS_85_24)),
       bty = "n", cex = 0.72)
sub("A and C coincide after the first campaign; A and B share it -- each pair differs in one campaign only")

dev.off()
cat("Wrote manuscript/figures/S9_soc_change_by_depth.png\n")
for (k in names(TR)) cat(sprintf("  %s: %s | 85-24 %+.3f | 06-24 %+.3f\n", k,
  paste(sprintf("%.2f", TR[[k]]), collapse = " "), (TR[[k]][3]-TR[[k]][1])/YRS_85_24, (TR[[k]][3]-TR[[k]][2])/18))
