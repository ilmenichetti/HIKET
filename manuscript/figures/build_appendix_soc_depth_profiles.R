setwd("/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling")
# Appendix: mean SOC DEPTH PROFILE of the three campaigns.
# The three campaigns do NOT share layer intervals (VMI8: organic/0-5/5-20/20-40;
# Biosoil & Komeetta: organic/0-10/10-20/20-40), so the panels are built two ways:
#   (a) each campaign at its OWN native resolution, as a concentration-equivalent
#       (stock / layer thickness) step profile -- the honest depth picture;
#   (b) cumulative stock with depth, including the modelled sub-40 cm tail;
#   (c) the per-plot decay rate implied by the two common mineral windows,
#       lambda = -log(C(20-40)/C(0-20))/20, which is interval-free and so is the
#       one shape statistic that compares the campaigns without harmonisation.
# ONE CONSISTENT BASIS THROUGHOUT: humus (OFH) + mineral. The LITTER layer (LM,
# a separately coded layer in the protocol -- see Massat_2006, Krs 100 = LM vs
# Krs 101 = OFH) is EXCLUDED, because it was measured in 2006/2024 but is absent
# from the 1985 sheet. Including it would put a definitional step between 1985
# and the later campaigns inside a figure about depth distribution. Its size is
# printed in panel (a) so nothing is hidden; what it does to the trend is the
# subject of appendix_soc_change_by_depth.
# Balanced panel (plots measured in all three campaigns), region-weighted.

source("manuscript/figures/model_palette.R")   # shared palettes (CAMPAIGN_COL)
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readxl) })

SOCDIR <- "Data/SOC_homogeneized"
layers <- read.csv(file.path(SOCDIR, "soc_homogenized_layers.csv"), stringsAsFactors = FALSE)
plotd  <- read.csv(file.path(SOCDIR, "soc_homogenized_plot.csv"),   stringsAsFactors = FALSE)
CA <- CAMPAIGN_ORDER

# --- balanced panel, documented outliers dropped ------------------------------
meta <- plotd |> select(plot_id, campaign, weight, soc_deep_Mgha, soc_profile_Mgha,
                        z_cap, soc_outlier)
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
cat(sprintf("balanced panel: %d plots\n", length(bal_ids)))

# --- split the organic layer into OFH (humus) and LM (litter) -----------------
XLSX <- "Data/Komeetta/Hannu/Komeetta 150526hi--.xlsx"
sp <- suppressMessages(read_excel(XLSX, sheet = "BiSo", range = "A4:IN3031", col_names = FALSE,
                                  col_types = "text", .name_repair = "minimal"))
nn  <- function(i) suppressWarnings(as.numeric(sp[[i]]))
lit <- tibble(plot_id = as.integer(nn(3)), Krs = as.integer(nn(4)), REP = as.integer(nn(9)),
              Biosoil = nn(247) / 1000, Komeetta = nn(248) / 1000) |>
  filter(!is.na(plot_id), REP == 1, Krs == 101) |> select(plot_id, Biosoil, Komeetta) |>
  pivot_longer(-plot_id, names_to = "campaign", values_to = "LM")
band <- band |> left_join(lit, by = c("plot_id", "campaign")) |>
  mutate(LM = coalesce(LM, 0), OFH = organic - LM)   # 1985: LM absent => OFH = organic

wm <- function(x, w) sum(x * w, na.rm = TRUE) / sum(w[!is.na(x)], na.rm = TRUE)

# --- native-resolution layer means (panel a) ---------------------------------
GEOM <- list("0-5cm" = c(0,5), "5-20cm" = c(5,20), "0-10cm" = c(0,10),
             "10-20cm" = c(10,20), "20-40cm" = c(20,40))
nat <- layers |> filter(plot_id %in% bal_ids, layer %in% names(GEOM)) |>
  left_join(meta |> select(plot_id, campaign, weight), by = c("plot_id","campaign")) |>
  group_by(campaign, layer) |>
  summarise(C = wm(C_Mgha, weight), n = n(), .groups = "drop") |>
  mutate(top = sapply(layer, function(l) GEOM[[l]][1]),
         bot = sapply(layer, function(l) GEOM[[l]][2]),
         conc = C / (bot - top))                      # Mg C /ha /cm

org <- band |> group_by(campaign) |>
  summarise(OFH = wm(OFH, weight), LM = wm(LM, weight), .groups = "drop")

# --- cumulative stock with depth (panel b) -----------------------------------
cum <- band |> group_by(campaign) |>
  summarise(organic = wm(OFH, weight), m0_20 = wm(m0_20, weight),
            m20_40 = wm(m20_40, weight), deep = wm(soc_deep_Mgha, weight),
            zcap = wm(pmin(100, z_cap), weight), .groups = "drop")

# --- interval-free shape statistic (panel c) ---------------------------------
shape <- band |> filter(m0_20 > 0, m20_40 > 0) |>
  mutate(lambda = -log(m20_40 / m0_20) / 20)
lam_ref <- median(plotd$lambda, na.rm = TRUE)   # pooled GTK-class lambda for the deep tail

# =============================================================================
png("manuscript/figures/appendix_soc_depth_profiles.png",
    width = 12, height = 4.8, units = "in", res = 200)
par(mfrow = c(1, 3), mar = c(4.0, 4.5, 3.9, 1.0), mgp = c(2.5, 0.7, 0), las = 1)
sub <- function(txt) mtext(txt, side = 3, line = 0.35, cex = 0.62, col = "grey35")

## (a) concentration-equivalent depth profile, native layers -------------------
# The organic layer has no measured thickness here, so it is NOT drawn on the
# concentration axis (that would invite reading a bar length as a concentration).
# Its stock is printed instead, in the shaded band above the mineral surface.
ORG_H <- 11
xmax  <- max(nat$conc) * 1.18
plot(NA, xlim = c(0, xmax), ylim = c(40, -ORG_H), xaxs = "i",
     xlab = "SOC concentration (Mg C/ha per cm)", ylab = "Depth (cm)",
     main = "(a)  Mean depth profile, native layers", axes = FALSE)
axis(1); axis(2, at = seq(0, 40, 10)); box()
rect(0, -ORG_H, xmax, 0, col = "grey95", border = NA)
abline(h = 0, col = "grey40", lwd = 1.2)
text(xmax * 0.03, -ORG_H * 0.84, "humus layer OFH, stock (Mg C/ha):",
     adj = 0, cex = 0.68, col = "grey35", font = 3)
for (i in seq_along(CA)) {
  cc <- CA[i]
  o  <- org$OFH[org$campaign == cc]; l <- org$LM[org$campaign == cc]
  text(xmax * (0.06 + 0.30 * (i - 1)), -ORG_H * 0.46, sprintf("%.1f", o),
       adj = 0, cex = 0.86, col = CAMPAIGN_COL[cc], font = 2)
  text(xmax * (0.06 + 0.30 * (i - 1)), -ORG_H * 0.14,
       if (l > 0) sprintf("+%.1f LM", l) else "LM absent",
       adj = 0, cex = 0.62, col = CAMPAIGN_COL[cc], font = 3)
  # mineral layers: step profile at native resolution
  d  <- nat |> filter(campaign == cc) |> arrange(top)
  lines(as.vector(rbind(d$conc, d$conc)), as.vector(rbind(d$top, d$bot)),
        col = CAMPAIGN_COL[cc], lwd = 2.6)
  segments(d$conc[-nrow(d)], d$bot[-nrow(d)], d$conc[-1], d$bot[-nrow(d)],
           col = CAMPAIGN_COL[cc], lwd = 2.6)
}
legend("bottomright", bty = "n", cex = 0.78, lwd = 2.6, col = CAMPAIGN_COL[CA],
       legend = CAMPAIGN_LAB[CA], inset = c(0, 0.02))
sub("mineral curves at native layers; LM (litter) excluded everywhere -- see text")

## (b) cumulative stock with depth --------------------------------------------
ymax <- max(cum$organic + cum$m0_20 + cum$m20_40 + cum$deep) * 1.08
plot(NA, xlim = c(0, 108), ylim = c(0, ymax), xlab = "Depth (cm)",
     ylab = "Cumulative SOC, humus + mineral (Mg C/ha)",
     main = "(b)  Cumulative stock with depth")
rect(40, 0, 108, ymax, col = adjustcolor("grey70", 0.16), border = NA)
text(74, ymax * 0.10, "modelled tail\n(not measured)", cex = 0.68, col = "grey35", font = 3)
for (cc in CA) {
  r  <- cum[cum$campaign == cc, ]
  zz <- c(0, 0, 20, 40, r$zcap)
  cs <- c(0, r$organic, r$organic + r$m0_20, r$organic + r$m0_20 + r$m20_40,
          r$organic + r$m0_20 + r$m20_40 + r$deep)
  lines(zz, cs, col = CAMPAIGN_COL[cc], lwd = 2.6)
  points(0, cs[2], pch = 16, col = CAMPAIGN_COL[cc], cex = 0.9)
  text(r$zcap, tail(cs, 1), sprintf(" %.1f", tail(cs, 1)), adj = 0,
       cex = 0.72, col = CAMPAIGN_COL[cc], font = 2, xpd = NA)
}
legend("topleft", bty = "n", cex = 0.78, lwd = 2.6, col = CAMPAIGN_COL[CA],
       legend = CAMPAIGN_LAB[CA], inset = c(0.02, 0.02))
sub("dots at depth 0 = the humus layer, which sits above the mineral surface")

## (c) interval-free shape statistic ------------------------------------------
bx  <- split(shape$lambda, factor(shape$campaign, levels = CA))
med <- sapply(bx, median)
ylim <- range(unlist(lapply(bx, quantile, c(0.02, 0.98))))
ylim[2] <- ylim[2] + diff(ylim) * 0.10
boxplot(bx, col = adjustcolor(CAMPAIGN_COL[CA], 0.75), border = "grey25",
        outline = FALSE, names = c("1985", "2006", "2024"), ylim = ylim,
        ylab = expression("Depth decay rate " * lambda * "  (cm"^-1 * ")"),
        main = "(c)  Mineral profile shape, per plot")
abline(h = lam_ref, lty = 2, col = "grey30")
legend("topleft", lty = 2, col = "grey30", bty = "n", cex = 0.68, seg.len = 1.6,
       legend = expression("pooled class " * lambda * ", used for the deep tail"))
ups <- sapply(bx, quantile, 0.75)
text(seq_along(med), ups, sprintf("%.3f", med), pos = 3, cex = 0.72, offset = 0.30, font = 2)
sub(expression(lambda == -log(C["20-40"] / C["0-20"]) / 20 * ";  larger = steeper decline"))

dev.off()
cat("Wrote manuscript/figures/appendix_soc_depth_profiles.png\n")
cat(sprintf("  OFH: %s\n", paste(sprintf("%s %.2f", CA, org$OFH[match(CA, org$campaign)]), collapse = "  ")))
cat(sprintf("  lambda medians: %s\n", paste(sprintf("%s %.4f", CA, med), collapse = "  ")))
