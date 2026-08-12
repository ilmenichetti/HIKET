setwd("/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling")
# =============================================================================
# Bland-Altman agreement between the three SOC campaigns, layer by layer.
# Internal report, August 2026.
#
# WHAT THIS ASKS. Do two campaigns, measuring the same plot, agree well enough
# for their difference to be read as carbon change? Neither campaign is a gold
# standard, so each panel plots
#       DIFFERENCE (later - earlier)   against   AVERAGE of the two
# rather than one against the other. Regressing the difference on the EARLIER
# value alone would manufacture a negative slope, because the difference
# contains -earlier and measurement error then drives both axes in opposite
# directions; the average splits that error between the axes.
#
# HOW TO READ A PANEL
#   solid horizontal line, offset from 0 -> CONSTANT bias: one campaign reads
#       systematically higher by a fixed amount.
#   dashed band (mean +/- 1.96 SD)       -> LIMITS OF AGREEMENT: the range in
#       which 95% of plot-level differences fall. Compare it with the layer's
#       own stock (printed): if the band is as wide as the stock, the pair
#       cannot resolve change at the plot level at all.
#   TILTED regression line               -> PROPORTIONAL bias: the disagreement
#       grows with the amount present. That is the fingerprint of a SCALE error
#       (bulk density, sampled volume, a calibration factor), because those
#       multiply the value, whereas carbon arriving from litter above would add
#       a roughly constant amount regardless of what is already there.
#
# LAYERS. Litter (LM) is excluded from the humus row: it exists for 2006/2024
# but not for 1985, so including it would compare different quantities. The
# whole-profile row is the actual calibration target (soc_profile).
#
# Run from repo root:
#   Rscript Reporting/August2026_internal_report/build_bland_altman_campaigns.R
# =============================================================================

source("manuscript/figures/model_palette.R")     # CAMPAIGN_COL, for pair colouring
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readxl) })

SOCDIR <- "Data/SOC_homogeneized"
XLSX   <- "Data/Komeetta/Hannu/Komeetta 150526hi--.xlsx"
OUT    <- "Reporting/August2026_internal_report"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

layers <- read.csv(file.path(SOCDIR, "soc_homogenized_layers.csv"), stringsAsFactors = FALSE)
plotd  <- read.csv(file.path(SOCDIR, "soc_homogenized_plot.csv"),   stringsAsFactors = FALSE)

band <- layers |>
  mutate(b = case_when(layer == "organic" ~ "organic",
                       layer %in% c("0-5cm","5-20cm","0-10cm","10-20cm") ~ "m0_20",
                       layer == "20-40cm" ~ "m20_40", TRUE ~ NA_character_)) |>
  filter(!is.na(b)) |>
  group_by(plot_id, campaign, b) |>
  summarise(C = sum(C_Mgha, na.rm = TRUE), n_src = n(), .groups = "drop") |>
  mutate(C = ifelse(b == "m0_20" & n_src < 2, NA_real_, C)) |> select(-n_src) |>
  pivot_wider(names_from = b, values_from = C) |>
  left_join(plotd |> select(plot_id, campaign, soc_deep_Mgha, soc_outlier),
            by = c("plot_id", "campaign")) |>
  filter(!(soc_outlier %in% TRUE))
bal <- band |> filter(!is.na(organic), !is.na(m0_20), !is.na(m20_40)) |>
  count(plot_id) |> filter(n == 3) |> pull(plot_id)
band <- band |> filter(plot_id %in% bal)

# strip the litter term so the humus row compares like with like
sp <- suppressMessages(read_excel(XLSX, sheet = "BiSo", range = "A4:IN3031", col_names = FALSE,
                                  col_types = "text", .name_repair = "minimal"))
nn  <- function(i) suppressWarnings(as.numeric(sp[[i]]))
lit <- tibble(plot_id = as.integer(nn(3)), Krs = as.integer(nn(4)), REP = as.integer(nn(9)),
              Biosoil = nn(247) / 1000, Komeetta = nn(248) / 1000) |>
  filter(!is.na(plot_id), REP == 1, Krs == 101) |> select(plot_id, Biosoil, Komeetta) |>
  pivot_longer(-plot_id, names_to = "campaign", values_to = "LM")
band <- band |> left_join(lit, by = c("plot_id", "campaign")) |>
  mutate(LM = coalesce(LM, 0), humus = organic - LM,
         profile = humus + m0_20 + m20_40 + soc_deep_Mgha)

W <- band |> select(plot_id, campaign, humus, m0_20, m20_40, profile) |>
  pivot_wider(names_from = campaign, values_from = c(humus, m0_20, m20_40, profile))
cat(sprintf("balanced panel: %d plots\n\n", nrow(W)))

ROWS <- list(c("humus",   "Humus layer (OFH)"),
             c("m0_20",   "Mineral 0-20 cm"),
             c("m20_40",  "Mineral 20-40 cm"),
             c("profile", "Whole profile  [the calibration target]"))
PAIRS <- list(c("VMI8", "Biosoil", "1985 → 2006", 21),
              c("Biosoil", "Komeetta", "2006 → 2024", 18))
PCOL <- c(unname(CAMPAIGN_COL["VMI8"]), unname(CAMPAIGN_COL["Komeetta"]))

ba <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  d <- y[ok] - x[ok]; m <- (x[ok] + y[ok]) / 2; f <- lm(d ~ m)
  list(d = d, m = m, fit = f, n = sum(ok), bias = mean(d), sdd = sd(d),
       level = mean(m), slope = coef(f)[2], p = summary(f)$coefficients[2, 4],
       loa = mean(d) + c(-1.96, 1.96) * sd(d))
}

# =============================================================================
png(file.path(OUT, "bland_altman_campaigns.png"),
    width = 10, height = 13.4, units = "in", res = 200)
par(mfcol = c(4, 2), mar = c(4.2, 4.6, 3.4, 1.2), mgp = c(2.6, 0.7, 0), las = 1,
    oma = c(0, 0, 3.0, 0))

res <- list()
for (j in seq_along(PAIRS)) {
  pr <- PAIRS[[j]]; col <- PCOL[j]
  for (i in seq_along(ROWS)) {
    v <- ROWS[[i]][1]; lab <- ROWS[[i]][2]
    b <- ba(W[[paste0(v, "_", pr[1])]], W[[paste0(v, "_", pr[2])]])
    res[[paste(v, pr[3])]] <- b
    yl <- max(abs(c(b$loa, b$d))) * 1.06
    plot(b$m, b$d, pch = 16, col = adjustcolor(col, 0.42), cex = 0.62,
         xlab = "Mean of the two campaigns (Mg C/ha)",
         ylab = "Difference, later − earlier (Mg C/ha)",
         main = sprintf("%s\n%s", lab, pr[3]), cex.main = 1.0,
         ylim = c(-yl, yl))
    abline(h = 0, col = "grey55", lwd = 1.3)
    abline(h = b$bias, col = col, lwd = 2, lty = 1)
    abline(h = b$loa, col = col, lwd = 1.2, lty = 2)
    xx <- seq(min(b$m), max(b$m), length.out = 60)
    lines(xx, predict(b$fit, data.frame(m = xx)), col = "grey15", lwd = 2.4)
    # furniture labels
    usr <- par("usr")
    text(usr[2], b$bias, sprintf(" bias %+.2f", b$bias), adj = c(1, -0.45),
         cex = 0.66, col = col, font = 2)
    text(usr[2], b$loa[2], sprintf(" +1.96 SD  %+.1f", b$loa[2]), adj = c(1, -0.35),
         cex = 0.6, col = col)
    text(usr[2], b$loa[1], sprintf(" −1.96 SD  %+.1f", b$loa[1]), adj = c(1, 1.25),
         cex = 0.6, col = col)
    legend("topleft", bty = "n", cex = 0.68,
           legend = c(sprintf("slope %+.3f%s", b$slope, ifelse(b$p < 0.05, "  (p<0.05)", " (n.s.)")),
                      sprintf("agreement band %.0f%% of the mean stock",
                              100 * diff(b$loa) / b$level),
                      sprintf("n = %d", b$n)))
  }
}
mtext("Bland–Altman agreement between SOC campaigns, by layer   |   tilted grey line = proportional bias (a scale error);   dashed = 95% limits of agreement",
      outer = TRUE, side = 3, line = 0.6, cex = 0.78, font = 2)
dev.off()
cat("Wrote", file.path(OUT, "bland_altman_campaigns.png"), "\n\n")

# --- console summary ---------------------------------------------------------
cat(sprintf("%-38s %6s %8s %8s %9s %22s\n",
            "layer / pair", "n", "bias", "slope", "p(slope)", "limits of agreement"))
for (k in names(res)) { b <- res[[k]]
  cat(sprintf("%-38s %6d %+8.2f %+8.3f %9.2g   %+7.1f to %+7.1f  (%.0f%% of stock)\n",
              k, b$n, b$bias, b$slope, b$p, b$loa[1], b$loa[2],
              100 * diff(b$loa) / b$level)) }
cat("\nRead-out: proportional bias (non-zero slope) appears in every mineral pair, so it is\n")
cat("a property of the mineral measurement rather than of one campaign. The limits of\n")
cat("agreement are the operative number: where the band is comparable to the stock itself,\n")
cat("that layer cannot resolve plot-level change in that interval.\n")

# =============================================================================
# Single-panel version — the subsoil, both campaign pairs overlaid.
# Same construction as the panels above; this is the S10(c) design standalone,
# because overlaying the two pairs is what makes the contrast legible: the
# 1985-2006 line is steeply tilted, the 2006-2024 line nearly flat, and both
# agreement bands are enormous relative to the stock.
# =============================================================================
png(file.path(OUT, "bland_altman_subsoil.png"),
    width = 7.6, height = 5.8, units = "in", res = 200)
par(mar = c(4.4, 4.6, 4.0, 1.2), mgp = c(2.6, 0.7, 0), las = 1)

B <- list(`1985 vs 2006`  = res[["m20_40 1985 → 2006"]],
          `2006 vs 2024`  = res[["m20_40 2006 → 2024"]])
xr <- range(unlist(lapply(B, `[[`, "m"))); yr <- c(-26, 26)
plot(NA, xlim = xr, ylim = yr,
     xlab = "Mean of the two campaigns, 20-40 cm (Mg C/ha)",
     ylab = "Difference, later − earlier (Mg C/ha)",
     main = "Do the campaigns agree on the subsoil?")
abline(h = 0, col = "grey40", lwd = 1.4)
for (i in seq_along(B)) {
  b <- B[[i]]; cl <- PCOL[i]
  abline(h = b$loa, col = adjustcolor(cl, 0.45), lty = 3)
  points(b$m, b$d, pch = 16, col = adjustcolor(cl, 0.45), cex = 0.6)
  xx <- seq(xr[1], xr[2], length.out = 60)
  lines(xx, predict(b$fit, data.frame(m = xx)), col = cl, lwd = 2.6)
}
legend("topleft", bty = "n", cex = 0.78, lwd = 2.6, col = PCOL,
       legend = sprintf("%s   slope %+.2f%s", names(B),
                        sapply(B, `[[`, "slope"),
                        ifelse(sapply(B, `[[`, "p") < 0.05, "*", "")))
mtext("dotted = 95% limits of agreement;  a tilted line = PROPORTIONAL bias, i.e. a scale error",
      side = 3, line = 0.35, cex = 0.68, col = "grey35")
dev.off()
cat("Wrote", file.path(OUT, "bland_altman_subsoil.png"), "\n")
