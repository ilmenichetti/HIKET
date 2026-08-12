# =============================================================================
# litter_layer_vs_input.R — is the measured LITTER LAYER (LM) growth real?
#
# THE DISCREPANCY. Between 2006 and 2024 the measured litter layer grew +27%
# (3.34 -> 4.25 Mg C/ha) while the Tupek litter INPUT series fell ~13% over the
# same window. A litter layer is a fast pool fed directly by input; it should
# track input, not oppose it.
#
# TWO READINGS, and they make opposite per-plot predictions.
#   H1 REAL          the layer genuinely grew   => plots whose INPUT rose should
#                    be the plots whose LM rose:  corr(dLM, dJ) > 0
#   H2 RECLASSIFIED  the LM/OFH split was drawn deeper in 2024, moving material
#                    out of humus into litter   => corr(dLM, dOFH) << 0, and no
#                    relation to input
#
# A cross-sectional check backs it up: if LM is a real litter pool in quasi-
# steady state, LM ~ J * tau, so LM should correlate with J ACROSS plots in each
# campaign separately. A reclassification artefact need not.
#
# Run from repo root:  Rscript doublechecks/litter_layer_vs_input.R
# =============================================================================

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readxl) })
options(width = 125)

ROOT <- if (dir.exists("Data/SOC_homogeneized")) "." else ".."
XLSX <- file.path(ROOT, "Data/Komeetta/Hannu/Komeetta 150526hi--.xlsx")
hdr  <- function(x) cat("\n", strrep("=", 108), "\n", x, "\n", strrep("=", 108), "\n", sep = "")

# --- measured organic sub-layers, 2006 and 2024 -------------------------------
L <- read.csv(file.path(ROOT, "Data/SOC_homogeneized/soc_homogenized_layers.csv"))
P <- read.csv(file.path(ROOT, "Data/SOC_homogeneized/soc_homogenized_plot.csv"))
sp <- suppressMessages(read_excel(XLSX, sheet = "BiSo", range = "A4:IN3031", col_names = FALSE,
                                  col_types = "text", .name_repair = "minimal"))
nn  <- function(i) suppressWarnings(as.numeric(sp[[i]]))
lit <- tibble(plot_id = as.integer(nn(3)), Krs = as.integer(nn(4)), REP = as.integer(nn(9)),
              Biosoil = nn(247) / 1000, Komeetta = nn(248) / 1000) |>
  filter(!is.na(plot_id), REP == 1, Krs == 101) |> select(plot_id, Biosoil, Komeetta) |>
  pivot_longer(-plot_id, names_to = "campaign", values_to = "LM")

org <- L |> filter(layer == "organic") |> select(plot_id, campaign, organic = C_Mgha) |>
  left_join(P |> select(plot_id, campaign, weight, soc_outlier), by = c("plot_id", "campaign")) |>
  filter(!(soc_outlier %in% TRUE)) |>
  left_join(lit, by = c("plot_id", "campaign")) |>
  filter(campaign %in% c("Biosoil", "Komeetta"), is.finite(LM)) |>
  mutate(OFH = organic - LM) |>
  select(plot_id, campaign, LM, OFH, weight) |>
  pivot_wider(names_from = campaign, values_from = c(LM, OFH, weight)) |>
  filter(is.finite(LM_Biosoil), is.finite(LM_Komeetta), is.finite(OFH_Biosoil), is.finite(OFH_Komeetta)) |>
  mutate(dLM = LM_Komeetta - LM_Biosoil, dOFH = OFH_Komeetta - OFH_Biosoil)

# --- litter input series, per plot -------------------------------------------
f <- sort(list.files(file.path(ROOT, "Data/model_inputs"), pattern = "^SP1_inputs_.*\\.rds$",
                     full.names = TRUE), decreasing = TRUE)[1]
cat("input bundle:", basename(f), "\n")
inp <- bind_rows(readRDS(f)$inputs_by_plot)
cat("input years:", paste(range(inp$year), collapse = "-"), " plots:", n_distinct(inp$plot_id), "\n")

# 5-year means centred on each campaign, to damp interannual noise
W06 <- 2004:2008; W24 <- intersect(2020:2024, unique(inp$year))
J <- inp |> filter(year %in% c(W06, W24)) |>
  mutate(win = ifelse(year %in% W06, "J06", "J24")) |>
  group_by(plot_id, win) |> summarise(J = mean(J_total, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(names_from = win, values_from = J) |>
  mutate(dJ = J24 - J06, rJ = J24 / J06)

d <- inner_join(org, J, by = "plot_id") |> filter(is.finite(dJ), is.finite(dLM))
cat("plots with both the layer pair and the input series:", nrow(d), "\n")

wm <- function(x, w) sum(x * w, na.rm = TRUE) / sum(w[!is.na(x)], na.rm = TRUE)
hdr("AGGREGATE: do the layer and its input move the same way?")
cat(sprintf("  litter INPUT  %s -> %s : %.3f -> %.3f tC/ha/yr  (%+.1f%%)\n",
            paste(range(W06), collapse = "-"), paste(range(W24), collapse = "-"),
            wm(d$J06, d$weight_Biosoil), wm(d$J24, d$weight_Biosoil),
            100 * (wm(d$J24, d$weight_Biosoil) / wm(d$J06, d$weight_Biosoil) - 1)))
cat(sprintf("  litter LAYER  2006 -> 2024      : %.3f -> %.3f tC/ha      (%+.1f%%)\n",
            wm(d$LM_Biosoil, d$weight_Biosoil), wm(d$LM_Komeetta, d$weight_Biosoil),
            100 * (wm(d$LM_Komeetta, d$weight_Biosoil) / wm(d$LM_Biosoil, d$weight_Biosoil) - 1)))
cat(sprintf("  humus  OFH    2006 -> 2024      : %.2f -> %.2f tC/ha      (%+.1f%%)\n",
            wm(d$OFH_Biosoil, d$weight_Biosoil), wm(d$OFH_Komeetta, d$weight_Biosoil),
            100 * (wm(d$OFH_Komeetta, d$weight_Biosoil) / wm(d$OFH_Biosoil, d$weight_Biosoil) - 1)))

hdr("H1 (real growth):  does dLM track dJ across plots?")
ct1 <- cor.test(d$dLM, d$dJ)
cat(sprintf("  corr(dLM, dJ)   = %+0.3f  [%+0.3f, %+0.3f]  p = %.3g   spearman %+0.3f\n",
            ct1$estimate, ct1$conf.int[1], ct1$conf.int[2], ct1$p.value,
            cor(d$dLM, d$dJ, method = "spearman")))

hdr("H2 (reclassification):  does dLM come out of dOFH?")
ct2 <- cor.test(d$dLM, d$dOFH)
f2  <- lm(dLM ~ dOFH, data = d)
cat(sprintf("  corr(dLM, dOFH) = %+0.3f  [%+0.3f, %+0.3f]  p = %.3g   slope %+0.3f\n",
            ct2$estimate, ct2$conf.int[1], ct2$conf.int[2], ct2$p.value, coef(f2)[2]))
cat("  (a pure LM<->OFH transfer would give correlation -1 and slope -1)\n")

hdr("CROSS-SECTIONAL: is LM a real litter pool at all?  (LM ~ J * tau)")
for (cc in c("Biosoil", "Komeetta")) {
  lmv <- d[[paste0("LM_", cc)]]; jv <- if (cc == "Biosoil") d$J06 else d$J24
  cat(sprintf("  %-9s corr(LM, J) = %+0.3f  (p = %.2g)\n", cc,
              cor(lmv, jv), cor.test(lmv, jv)$p.value))
}
cat("  A litter layer in quasi-steady state with its input should show a clear positive\n")
cat("  relation here. If it does not, LM is not behaving like a litter pool.\n")

hdr("VERDICT  (result recorded 2026-08-12)")
cat("BOTH hypotheses fail at the plot level:\n")
cat("  H1 real growth      corr(dLM, dJ)   = +0.09 (p=0.07)  -- marginal, effectively null\n")
cat("  H2 reclassification corr(dLM, dOFH) = -0.04 (p=0.45), slope -0.01 -- null\n")
cat("  and cross-sectionally corr(LM, J) = +0.03 / +0.01 in the two campaigns: the measured\n")
cat("  litter layer has NO relation to modelled litter input in either year.\n\n")
cat("WHAT THAT LEAVES. The aggregate contradiction is real (input -4.9% on this window,\n")
cat("layer +27.7%), but no PLOT-VARYING mechanism explains it. A uniform campaign-level\n")
cat("shift in what was called 'litter' in 2024 would produce exactly this signature --\n")
cat("a large aggregate change with no plot-level correlate -- because plot-level\n")
cat("correlations are blind to an offset applied equally to every plot. That is Hannu's\n")
cat("protocol concern, and it is the surviving explanation.\n\n")
cat("CAVEAT, and it is a real one. The null corr(LM, J) cannot distinguish 'LM is not\n")
cat("behaving like a litter pool' from 'the Tupek series carries little plot-level\n")
cat("signal'. The latter is independently known (litter does not order SOC across plots;\n")
cat("model R^2 ~ 0). So H1 is not so much rejected as untestable with this predictor.\n")
cat("Do not report this as evidence that the litter layer is fake.\n")
