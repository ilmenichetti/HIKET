setwd("/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling")
# Appendix: how repeatable is each depth band between campaigns, and is any of
# the apparent change just the organic/mineral boundary moving?
#
# Rationale. Subsoil carbon is a near-static plot property on a 20-40 yr horizon,
# so the between-campaign correlation of a band is an upper bound on how well
# that band is measured. 2006->2024 is the reference pair: a band whose 1985
# pairing is much weaker than its 2006->2024 pairing is failing at the 1985 end.
# Panel (c) tests the other failure mode -- the organic/mineral boundary is set
# by eye in the field, so carbon can be reassigned between the organic layer and
# mineral 0-20 with no carbon moving at all (Kramarenko 2012 sec. 4.4).
#
# Bands are like-for-like (litter removed from 2006/2024, which the 1985 sheet
# lacks). Balanced panel, documented outliers dropped.

source("manuscript/figures/model_palette.R")
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readxl) })

SOCDIR <- "Data/SOC_homogeneized"
XLSX   <- "Data/Komeetta/Hannu/Komeetta 150526hi--.xlsx"
layers <- read.csv(file.path(SOCDIR, "soc_homogenized_layers.csv"), stringsAsFactors = FALSE)
plotd  <- read.csv(file.path(SOCDIR, "soc_homogenized_plot.csv"),   stringsAsFactors = FALSE)

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
  pivot_longer(-plot_id, names_to = "campaign", values_to = "litter")
band <- band |> left_join(lit, by = c("plot_id", "campaign")) |>
  # ⚠ 2026-08-12: the baseline applies TREATMENT C, so the 1985 organic layer already
  # contains the imputed LM; lm_added_1985 (kg/ha, 0 elsewhere) recovers OFH.
  mutate(litter = coalesce(litter, coalesce(lm_added_1985, 0) / 1000),
         humus = organic - litter,
         profile = organic + m0_20 + m20_40 + soc_deep_Mgha)

W <- band |> select(plot_id, campaign, humus, m0_20, m20_40, profile) |>
  pivot_wider(names_from = campaign, values_from = c(humus, m0_20, m20_40, profile))
cat(sprintf("balanced panel: %d plots\n", nrow(W)))

BANDS <- c("humus", "m0_20", "m20_40", "profile")
BLAB  <- c("humus\nlayer", "mineral\n0-20 cm", "mineral\n20-40 cm", "whole\nprofile")
PAIRS <- list(c("VMI8", "Biosoil"), c("Biosoil", "Komeetta"))
PLAB  <- c("1985 vs 2006", "2006 vs 2024")
PCOL  <- c("1985 vs 2006" = unname(CAMPAIGN_COL["VMI8"]),
           "2006 vs 2024" = unname(CAMPAIGN_COL["Komeetta"]))

rr <- sapply(PAIRS, function(p) sapply(BANDS, function(v) {
  x <- W[[paste0(v, "_", p[1])]]; y <- W[[paste0(v, "_", p[2])]]
  ok <- is.finite(x) & is.finite(y); cor(x[ok], y[ok])
}))
colnames(rr) <- PLAB

# =============================================================================
# --- Bland-Altman: difference vs average, per campaign pair -------------------
# Neither campaign is a gold standard, so plot the DIFFERENCE against the
# AVERAGE of the two (not against one of them -- the difference contains -x, so
# regressing on x alone manufactures a negative slope out of measurement error).
#   intercept != 0  => constant bias
#   slope     != 0  => PROPORTIONAL bias, the fingerprint of a scale/calibration
#                      problem (density, volume) rather than an added amount
ba <- function(x, y) {
  ok <- is.finite(x) & is.finite(y) & x > 0 & y > 0
  d <- y[ok] - x[ok]; m <- (x[ok] + y[ok]) / 2; f <- lm(d ~ m)
  list(d = d, m = m, fit = f, bias = mean(d), sd = sd(d),
       slope = coef(f)[2], p = summary(f)$coefficients[2, 4],
       loa = mean(d) + c(-1.96, 1.96) * sd(d))
}
BA <- list(); for (i in seq_along(PAIRS)) {
  pr <- PAIRS[[i]]
  BA[[PLAB[i]]] <- ba(W[[paste0("m20_40_", pr[1])]], W[[paste0("m20_40_", pr[2])]])
}
cat("\nBland-Altman, mineral 20-40 cm:\n")
for (k in names(BA)) cat(sprintf("  %-14s bias %+5.2f  slope %+0.3f (p=%.2g)  LoA %+.1f to %+.1f Mg/ha\n",
  k, BA[[k]]$bias, BA[[k]]$slope, BA[[k]]$p, BA[[k]]$loa[1], BA[[k]]$loa[2]))
cat("  (mineral 0-20 for contrast: 1985->2006 slope +0.165, 2006->2024 slope +0.339 --\n")
cat("   proportional bias is present in EVERY pair, so it is a property of the mineral\n")
cat("   measurement in general, not of 1985 alone.)\n")

png("manuscript/figures/S10_soc_campaign_reliability.png",
    width = 10.5, height = 8.8, units = "in", res = 200)
par(mfrow = c(2, 2), mar = c(4.6, 4.5, 3.9, 1.0), mgp = c(2.5, 0.7, 0), las = 1)
sub <- function(txt) mtext(txt, side = 3, line = 0.35, cex = 0.62, col = "grey35")

## (a) repeatability by band ---------------------------------------------------
bp <- barplot(t(rr), beside = TRUE, col = PCOL[PLAB], border = "grey30",
              names.arg = BLAB, ylim = c(0, 1), cex.names = 0.82,
              ylab = "Between-campaign correlation (Pearson r)",
              main = "(a)  Repeatability by depth band")
text(as.vector(bp), as.vector(t(rr)), sprintf("%.2f", as.vector(t(rr))),
     pos = 3, cex = 0.68, offset = 0.25)
legend("topright", fill = PCOL[PLAB], legend = PLAB, bty = "n", cex = 0.76, border = "grey30")
sub("2006 vs 2024 is the reference: both ends measured with one protocol")

## (b) subsoil repeatability scatter -------------------------------------------
lim <- c(0, 55)
plot(NA, xlim = lim, ylim = lim, xlab = "Mineral 20-40 cm, earlier campaign (Mg C/ha)",
     ylab = "Mineral 20-40 cm, later campaign (Mg C/ha)",
     main = "(b)  Subsoil, plot by plot")
abline(0, 1, lty = 2, col = "grey40")
for (i in seq_along(PAIRS)) {
  p <- PAIRS[[i]]
  points(W[[paste0("m20_40_", p[1])]], W[[paste0("m20_40_", p[2])]],
         pch = 16, col = adjustcolor(PCOL[PLAB[i]], 0.55), cex = 0.62)
}
legend("topleft", pch = 16, col = PCOL[PLAB], bty = "n", cex = 0.76,
       legend = sprintf("%s  (r = %.2f)", PLAB, rr["m20_40", ]))
sub("1985 is no noisier here than 2006 -- the subsoil problem is a level, not scatter")

## (c) Bland-Altman of the subsoil ---------------------------------------------
xr <- range(unlist(lapply(BA, `[[`, "m"))); yr <- c(-26, 26)
plot(NA, xlim = xr, ylim = yr,
     xlab = "Mean of the two campaigns, 20-40 cm (Mg C/ha)",
     ylab = "Difference, later - earlier (Mg C/ha)",
     main = "(c)  Do the campaigns agree? (Bland-Altman)")
abline(h = 0, col = "grey40", lwd = 1.4)
for (i in seq_along(BA)) {
  k <- names(BA)[i]; b <- BA[[k]]; cl <- PCOL[k]
  abline(h = b$loa, col = adjustcolor(cl, 0.45), lty = 3)
  points(b$m, b$d, pch = 16, col = adjustcolor(cl, 0.45), cex = 0.6)
  xx <- seq(xr[1], xr[2], length.out = 50)
  lines(xx, predict(b$fit, data.frame(m = xx)), col = cl, lwd = 2.6)
}
legend("topleft", bty = "n", cex = 0.72, lwd = 2.6, col = PCOL[PLAB],
       legend = sprintf("%s   slope %+.2f%s", PLAB,
                        sapply(BA, `[[`, "slope"),
                        ifelse(sapply(BA, `[[`, "p") < 0.05, "*", "")))
sub("dotted = 95% limits of agreement; a tilted line = PROPORTIONAL bias, i.e. a scale error")

## (d) organic/mineral boundary trade-off --------------------------------------
dh <- lapply(PAIRS, function(p) W[[paste0("humus_", p[2])]] - W[[paste0("humus_", p[1])]])
dm <- lapply(PAIRS, function(p) W[[paste0("m0_20_", p[2])]] - W[[paste0("m0_20_", p[1])]])
rb <- sapply(seq_along(PAIRS), function(i) {
  ok <- is.finite(dh[[i]]) & is.finite(dm[[i]]); cor(dh[[i]][ok], dm[[i]][ok]) })
rg <- c(-32, 32)
plot(NA, xlim = rg, ylim = rg, xlab = expression(Delta * " humus layer (Mg C/ha)"),
     ylab = expression(Delta * " mineral 0-20 cm (Mg C/ha)"),
     main = "(d)  Is the boundary moving?")
abline(h = 0, v = 0, col = "grey80"); abline(0, -1, lty = 2, col = "grey40")
for (i in seq_along(PAIRS)) {
  ok <- is.finite(dh[[i]]) & is.finite(dm[[i]])
  points(dh[[i]][ok], dm[[i]][ok], pch = 16, col = adjustcolor(PCOL[PLAB[i]], 0.5), cex = 0.62)
  abline(lm(dm[[i]][ok] ~ dh[[i]][ok]), col = PCOL[PLAB[i]], lwd = 2.4)
}
legend("topright", bty = "n", cex = 0.72, lwd = c(2.4, 2.4, 1), lty = c(1, 1, 2),
       col = c(PCOL[PLAB], "grey40"),
       legend = c(sprintf("%s  (r = %+.2f)", PLAB, rb), "pure trade-off (slope -1)"))
sub("a strong negative slope would mean carbon was reassigned, not gained")

dev.off()
cat("Wrote manuscript/figures/S10_soc_campaign_reliability.png\n")
print(round(rr, 3)); cat(sprintf("boundary r: %s\n", paste(sprintf("%s %+.3f", PLAB, rb), collapse = "  ")))
