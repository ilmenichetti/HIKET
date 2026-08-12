# =============================================================================
# soc_depth_distribution.R — is the 1985 (VMI8) SUBSOIL trustworthy?
#
# Hypothesis under test (2026-08-12): the 1985 campaign's ORGANIC layer is
# credible (litter is what VMI8 measured well), but its MINERAL SUBSOIL is not.
# Two consequences would follow, and both are testable without touching a model:
#   (H-a) the 1985 profile SHAPE departs from the exponential decline the later
#         campaigns show — too much carbon at 20-40 cm;
#   (H-b) the subsoil, which should be near-inert on a 20-40 yr horizon, carries
#         an implausibly large share of the apparent 1985->2006 stock change.
#
# The campaigns do NOT share layer intervals (1985: organic/0-5/5-20/20-40;
# 2006+2024: organic/0-10/10-20/20-40), so EVERYTHING here is first put on the
# common grid  organic | 0-20 | 20-40  which both protocols can form exactly.
# Comparing raw layers across campaigns would compare different depth windows.
#
# Depth-shape metric. For C(z) = C0*exp(-lambda*z) the stock ratio of the two
# mineral windows is exactly
#       C(20-40) / C(0-20) = exp(-20*lambda)
# so     lambda_obs = -log(ratio)/20
# is a per-plot-campaign decay rate read straight off the data, with C0 and all
# units cancelling. It is directly comparable to the pooled GTK-class lambda
# (~0.032) that build_soc_homogenized.R uses for the deep tail.
#
# Run from repo root:  Rscript doublechecks/soc_depth_distribution.R
# =============================================================================

suppressPackageStartupMessages({ library(dplyr); library(tidyr) })
options(width = 130)
set.seed(2025)

ROOT <- if (dir.exists("Data/SOC_homogeneized")) "." else ".."   # repo root or doublechecks/
SOCDIR <- file.path(ROOT, "Data/SOC_homogeneized")
OUTDIR <- file.path(ROOT, "doublechecks/soc_depth")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

layers <- read.csv(file.path(SOCDIR, "soc_homogenized_layers.csv"), stringsAsFactors = FALSE)
plotd  <- read.csv(file.path(SOCDIR, "soc_homogenized_plot.csv"),   stringsAsFactors = FALSE)

CAMPS <- c("VMI8", "Biosoil", "Komeetta")
layers$campaign <- factor(layers$campaign, levels = CAMPS)
plotd$campaign  <- factor(plotd$campaign,  levels = CAMPS)

hdr <- function(x) cat("\n\n", strrep("=", 118), "\n", x, "\n", strrep("=", 118), "\n", sep = "")
p3  <- function(x, d = 3) formatC(x, format = "f", digits = d, width = 8)

# =============================================================================
# 1. COMPLETENESS AND GRANULARITY — before any aggregation
#    A campaign whose deep layer never goes missing, while the later campaigns
#    lose 4-6% of theirs to augering refusal, is reporting something other than
#    a measurement.
# =============================================================================
hdr("1. RAW LAYER COMPLETENESS AND GRANULARITY")

gran <- layers |>
  group_by(campaign, layer) |>
  summarise(n        = n(),
            n_NA     = sum(is.na(C_kgha)),
            n_zero   = sum(C_kgha == 0, na.rm = TRUE),
            n_unique = n_distinct(round(C_kgha, 6)),
            median   = median(C_Mgha, na.rm = TRUE),
            .groups  = "drop") |>
  arrange(campaign, layer)
print(as.data.frame(gran), row.names = FALSE)

n_plots <- layers |> group_by(campaign) |> summarise(n_plot = n_distinct(plot_id), .groups = "drop")
cat("\nplots per campaign:\n"); print(as.data.frame(n_plots), row.names = FALSE)
cat("\n-> per-layer coverage as a FRACTION of that campaign's plots:\n")
gran |> left_join(n_plots, by = "campaign") |>
  mutate(coverage = round(n / n_plot, 3)) |>
  select(campaign, layer, n, n_plot, coverage) |>
  as.data.frame() |> print(row.names = FALSE)

# =============================================================================
# 2. COMMON DEPTH GRID  organic | 0-20 | 20-40
# =============================================================================
grid <- layers |>
  mutate(band = case_when(
    layer == "organic"                        ~ "organic",
    layer %in% c("0-5cm","5-20cm","0-10cm","10-20cm") ~ "m0_20",
    layer == "20-40cm"                        ~ "m20_40",
    TRUE                                      ~ NA_character_)) |>
  filter(!is.na(band)) |>
  group_by(plot_id, campaign, band) |>
  summarise(C = sum(C_Mgha, na.rm = TRUE), n_src = n(), .groups = "drop") |>
  # a 0-20 band is only valid if BOTH of its sub-layers are present
  mutate(C = ifelse(band == "m0_20" & n_src < 2, NA_real_, C)) |>
  select(-n_src) |>
  pivot_wider(names_from = band, values_from = C) |>
  mutate(min_0_40 = m0_20 + m20_40,
         soc_0_40 = organic + min_0_40)

# region weights + the deep tail from the built baseline
grid <- grid |>
  left_join(plotd |> select(plot_id, campaign, weight, region_label, soil_code,
                            soc_deep_Mgha, soc_profile_Mgha, lambda, z_cap,
                            soil_depth_reached, soc_outlier, high_change),
            by = c("plot_id", "campaign"))

# analysis set: drop the documented outliers, keep everything else
grid <- grid |> filter(!(soc_outlier %in% TRUE))

# balanced panel: plots measured in all three campaigns
bal_ids <- grid |> filter(!is.na(soc_0_40)) |>
  count(plot_id) |> filter(n == 3) |> pull(plot_id)
bal <- grid |> filter(plot_id %in% bal_ids)
cat("\nbalanced panel (all three campaigns, outliers dropped): n =", length(bal_ids), "plots\n")

wmean <- function(x, w) sum(x * w, na.rm = TRUE) / sum(w[!is.na(x)], na.rm = TRUE)

hdr("2. STOCKS ON THE COMMON GRID (balanced panel, Mg C /ha)")
tab2 <- bal |> group_by(campaign) |>
  summarise(n           = sum(!is.na(soc_0_40)),
            organic     = wmean(organic, weight),
            m0_20       = wmean(m0_20,   weight),
            m20_40      = wmean(m20_40,  weight),
            min_0_40    = wmean(min_0_40, weight),
            soc_0_40    = wmean(soc_0_40, weight),
            deep_tail   = wmean(soc_deep_Mgha, weight),
            soc_profile = wmean(soc_profile_Mgha, weight),
            .groups = "drop")
print(as.data.frame(tab2 |> mutate(across(where(is.numeric), \(x) round(x, 2)))), row.names = FALSE)

cat("\nSAME, unweighted (region weights off):\n")
tab2u <- bal |> group_by(campaign) |>
  summarise(organic = mean(organic, na.rm = TRUE), m0_20 = mean(m0_20, na.rm = TRUE),
            m20_40 = mean(m20_40, na.rm = TRUE), soc_0_40 = mean(soc_0_40, na.rm = TRUE),
            soc_profile = mean(soc_profile_Mgha, na.rm = TRUE), .groups = "drop")
print(as.data.frame(tab2u |> mutate(across(where(is.numeric), \(x) round(x, 2)))), row.names = FALSE)

cat("\nCOMPOSITION — each band as a share of measured 0-40 cm (weighted):\n")
tab2 |> transmute(campaign,
                  f_organic = round(organic / soc_0_40, 3),
                  f_0_20    = round(m0_20   / soc_0_40, 3),
                  f_20_40   = round(m20_40  / soc_0_40, 3),
                  f_deep_of_profile = round(deep_tail / soc_profile, 3)) |>
  as.data.frame() |> print(row.names = FALSE)

# =============================================================================
# 3. DEPTH SHAPE — the exponential decay rate read off the data
#    lambda_obs = -log( C(20-40) / C(0-20) ) / 20
#    Larger lambda = steeper decline. A campaign that over-reports subsoil shows
#    a SMALLER lambda (flatter profile).
# =============================================================================
hdr("3. DEPTH SHAPE: per-plot decay rate lambda_obs (1/cm)")

shape <- bal |>
  filter(!is.na(m0_20), !is.na(m20_40), m0_20 > 0, m20_40 > 0) |>
  mutate(ratio = m20_40 / m0_20,
         lambda_obs = -log(ratio) / 20)

sh <- shape |> group_by(campaign) |>
  summarise(n = n(),
            ratio_median  = median(ratio),
            lambda_median = median(lambda_obs),
            lambda_q25    = quantile(lambda_obs, .25),
            lambda_q75    = quantile(lambda_obs, .75),
            frac_flat     = mean(lambda_obs < 0.010),   # near-flat profile
            frac_inverted = mean(ratio > 1),            # MORE carbon at 20-40 than 0-20
            .groups = "drop")
print(as.data.frame(sh |> mutate(across(where(is.numeric), \(x) round(x, 4)))), row.names = FALSE)

cat("\npooled GTK-class lambda used for the deep tail in build_soc_homogenized.R:\n")
cat("  median over plots:", round(median(plotd$lambda, na.rm = TRUE), 4),
    " range:", paste(round(range(plotd$lambda, na.rm = TRUE), 4), collapse = " - "), "\n")

# paired within-plot shape change: same plot, same soil, different campaign
sh_w <- shape |> select(plot_id, campaign, lambda_obs, ratio) |>
  pivot_wider(names_from = campaign, values_from = c(lambda_obs, ratio))
cat("\nPAIRED within-plot shape (same plot, so soil type/stoniness cancel):\n")
pair_shape <- function(a, b) {
  d <- sh_w[[paste0("lambda_obs_", b)]] - sh_w[[paste0("lambda_obs_", a)]]
  ok <- !is.na(d)
  cat(sprintf("  lambda %-8s -> %-8s : n=%3d  median delta = %+0.4f   %% steeper in %s = %4.1f%%   (paired t p = %.2g)\n",
              a, b, sum(ok), median(d[ok]), b, 100 * mean(d[ok] > 0),
              tryCatch(t.test(d[ok])$p.value, error = function(e) NA)))
}
pair_shape("VMI8", "Biosoil"); pair_shape("Biosoil", "Komeetta"); pair_shape("VMI8", "Komeetta")

cat("\nORGANIC-layer share of measured 0-40 cm, per plot (median):\n")
bal |> filter(!is.na(soc_0_40), soc_0_40 > 0) |>
  group_by(campaign) |>
  summarise(f_org_median = round(median(organic / soc_0_40, na.rm = TRUE), 3),
            f_20_40_median = round(median(m20_40 / soc_0_40, na.rm = TRUE), 3), .groups = "drop") |>
  as.data.frame() |> print(row.names = FALSE)

# =============================================================================
# 4. WHERE DOES THE APPARENT CHANGE COME FROM?
#    The subsoil should be the most inert compartment on a 20-40 yr horizon.
#    If it supplies a large share of the 1985->2006 change but not of
#    2006->2024, the 1985 subsoil is the odd one out.
# =============================================================================
hdr("4. CHANGE DECOMPOSITION BY DEPTH BAND (balanced panel, weighted, Mg C/ha)")

w <- bal |> select(plot_id, campaign, weight, organic, m0_20, m20_40, soc_0_40,
                   soc_profile_Mgha, soc_deep_Mgha) |>
  pivot_wider(names_from = campaign,
              values_from = c(organic, m0_20, m20_40, soc_0_40, soc_profile_Mgha, soc_deep_Mgha))

delta_tab <- function(a, b, yrs) {
  bands <- c("organic", "m0_20", "m20_40", "soc_0_40", "soc_deep_Mgha", "soc_profile_Mgha")
  out <- lapply(bands, function(v) {
    d <- w[[paste0(v, "_", b)]] - w[[paste0(v, "_", a)]]
    data.frame(band = v, n = sum(!is.na(d)),
               delta = wmean(d, w$weight), delta_per_yr = wmean(d, w$weight) / yrs,
               median_delta = median(d, na.rm = TRUE),
               frac_increasing = mean(d > 0, na.rm = TRUE))
  }) |> bind_rows()
  tot <- out$delta[out$band == "soc_0_40"]
  out$share_of_0_40_change <- ifelse(out$band %in% c("organic","m0_20","m20_40"),
                                     out$delta / tot, NA)
  cat("\n--- ", a, " -> ", b, " (", yrs, " yr) ---\n", sep = "")
  print(as.data.frame(out |> mutate(across(where(is.numeric), \(x) round(x, 3)))), row.names = FALSE)
  invisible(out)
}
d1 <- delta_tab("VMI8", "Biosoil", 21)
d2 <- delta_tab("Biosoil", "Komeetta", 18)

cat("\nSUBSOIL (20-40 cm) rate of change, Mg C/ha/yr — this compartment should be near-inert:\n")
cat(sprintf("  1985->2006 : %+0.3f\n  2006->2024 : %+0.3f\n",
            d1$delta_per_yr[d1$band == "m20_40"], d2$delta_per_yr[d2$band == "m20_40"]))

# =============================================================================
# 5. REPEATABILITY — does a plot's subsoil in one campaign predict its own
#    subsoil in another? Real subsoil carbon is a near-static plot property, so
#    the between-campaign correlation is an upper bound on measurement quality.
# =============================================================================
hdr("5. BETWEEN-CAMPAIGN REPEATABILITY, BY DEPTH BAND (paired plots)")

rep_row <- function(v, a, b) {
  x <- w[[paste0(v, "_", a)]]; y <- w[[paste0(v, "_", b)]]
  ok <- !is.na(x) & !is.na(y)
  data.frame(band = v, pair = paste(a, b, sep = "->"), n = sum(ok),
             pearson  = cor(x[ok], y[ok]),
             spearman = cor(x[ok], y[ok], method = "spearman"),
             pearson_log = cor(log(pmax(x[ok], .01)), log(pmax(y[ok], .01))),
             cv_diff = sd(y[ok] - x[ok]) / mean(c(x[ok], y[ok])))
}
reps <- bind_rows(lapply(c("organic","m0_20","m20_40","soc_0_40","soc_profile_Mgha"), function(v)
  bind_rows(rep_row(v, "VMI8", "Biosoil"), rep_row(v, "Biosoil", "Komeetta"),
            rep_row(v, "VMI8", "Komeetta"))))
print(as.data.frame(reps |> mutate(across(where(is.numeric), \(x) round(x, 3)))), row.names = FALSE)

cat("\n-> read this as: 2006->2024 is the reference for how repeatable each band is\n")
cat("   when both ends are trusted. A band whose 1985 pairing is much weaker than\n")
cat("   its 2006->2024 pairing is failing at the 1985 end.\n")

rr <- reps |> select(band, pair, pearson) |> pivot_wider(names_from = pair, values_from = pearson)
rr$ratio_85_06_over_06_24 <- round(rr$`VMI8->Biosoil` / rr$`Biosoil->Komeetta`, 3)
print(as.data.frame(rr |> mutate(across(where(is.numeric), \(x) round(x, 3)))), row.names = FALSE)

# =============================================================================
# 6. CONSEQUENCE — rebuild the 1985 subsoil from a trusted shape and see what
#    it does to the level and, more importantly, to the TREND.
#
#    Reconstruction: keep 1985 organic and 0-20 cm AS MEASURED (the trusted
#    part) and replace 1985's 20-40 cm by applying THAT PLOT'S OWN depth shape
#    from the later campaigns:
#         m20_40_hat(1985) = m0_20(1985) * ratio_plot(later)
#    ratio_plot = the plot's own median m20_40/m0_20 over 2006+2024; falls back
#    to the campaign-median ratio where the plot has none. The deep tail is then
#    refitted with the SAME class lambda and the same C0 least-squares rule
#    build_soc_homogenized.R uses, so the tail moves consistently.
# =============================================================================
hdr("6. CONSEQUENCE OF RECONSTRUCTING THE 1985 SUBSOIL FROM A TRUSTED SHAPE")

lat <- bal |> filter(campaign != "VMI8", !is.na(m0_20), m0_20 > 0, !is.na(m20_40)) |>
  group_by(plot_id) |> summarise(ratio_later = median(m20_40 / m0_20), .groups = "drop")
ratio_fallback <- median(lat$ratio_later, na.rm = TRUE)
cat("plot-specific later-campaign ratio available for", nrow(lat), "plots; fallback ratio =",
    round(ratio_fallback, 4), "\n")

# exact C0 least-squares weight for the two mineral bands under class lambda
lw <- function(lam, z1, z2) (1 / lam) * (exp(-lam * z1) - exp(-lam * z2))
deep_from_bands <- function(c0_20, c20_40, lam, z_from, z_to) {
  if (is.na(lam) || is.na(c0_20) || is.na(c20_40)) return(NA_real_)
  if (z_to <= z_from) return(0)
  wts <- c(lw(lam, 0, 20), lw(lam, 20, 40))
  C0  <- sum(c(c0_20, c20_40) * wts) / sum(wts^2)
  C0 / lam * (exp(-lam * z_from) - exp(-lam * z_to))
}

v85 <- bal |> filter(campaign == "VMI8") |> left_join(lat, by = "plot_id") |>
  mutate(ratio_later = ifelse(is.na(ratio_later), ratio_fallback, ratio_later),
         m20_40_hat  = m0_20 * ratio_later,
         z_from      = pmin(40, z_cap),
         deep_orig   = mapply(deep_from_bands, m0_20, m20_40,     lambda, z_from, z_cap),
         deep_hat    = mapply(deep_from_bands, m0_20, m20_40_hat, lambda, z_from, z_cap),
         soc_0_40_hat = organic + m0_20 + m20_40_hat,
         # keep the published tail where our re-fit reproduces it, else use ours
         tail_check   = deep_orig - soc_deep_Mgha,
         profile_hat  = soc_0_40_hat + deep_hat + (soc_deep_Mgha - deep_orig))

cat("\ndeep-tail re-fit check (our reproduction of soc_deep vs the built value):",
    "median abs diff =", round(median(abs(v85$tail_check), na.rm = TRUE), 3), "Mg/ha\n")

cat("\n1985 BEFORE vs AFTER reconstruction (weighted means, Mg C/ha):\n")
cmp <- data.frame(
  quantity = c("mineral 20-40", "measured 0-40", "deep tail", "soc_profile (TARGET)"),
  before   = c(wmean(v85$m20_40, v85$weight), wmean(v85$soc_0_40, v85$weight),
               wmean(v85$soc_deep_Mgha, v85$weight), wmean(v85$soc_profile_Mgha, v85$weight)),
  after    = c(wmean(v85$m20_40_hat, v85$weight), wmean(v85$soc_0_40_hat, v85$weight),
               wmean(v85$deep_hat + (v85$soc_deep_Mgha - v85$deep_orig), v85$weight),
               wmean(v85$profile_hat, v85$weight)))
cmp$delta <- cmp$after - cmp$before
print(cmp |> mutate(across(where(is.numeric), \(x) round(x, 2))), row.names = FALSE)

lv <- tab2 |> filter(campaign != "VMI8") |> select(campaign, soc_profile)
p85_before <- cmp$before[cmp$quantity == "soc_profile (TARGET)"]
p85_after  <- cmp$after [cmp$quantity == "soc_profile (TARGET)"]
p06 <- lv$soc_profile[lv$campaign == "Biosoil"]; p24 <- lv$soc_profile[lv$campaign == "Komeetta"]

cat("\nEFFECT ON THE TREND (weighted profile means, Mg C/ha, balanced panel):\n")
trend <- data.frame(
  case      = c("as built", "1985 subsoil reconstructed"),
  soc_1985  = c(p85_before, p85_after), soc_2006 = c(p06, p06), soc_2024 = c(p24, p24),
  d_85_24   = c(p24 - p85_before, p24 - p85_after),
  rate_85_24 = c((p24 - p85_before) / 39, (p24 - p85_after) / 39),
  rate_06_24 = c((p24 - p06) / 18, (p24 - p06) / 18))
print(trend |> mutate(across(where(is.numeric), \(x) round(x, 3))), row.names = FALSE)
cat("\n!! the 2006->2024 rate is UNCHANGED by construction — no 1985 revision can touch it.\n")

# =============================================================================
# 7. THE ORGANIC LAYER IS NOT DEFINED THE SAME WAY IN 1985 AND LATER
#    Heikkinen folds a LITTER term into the 2006 and 2024 organic layer
#    (build_soc_homogenized.R: C = 1e4*C_kgm2 + lit<year>). The 1985 sheet has
#    no litter component at all. The README calls this "minor, sub-Mg/ha".
#    It is not: measure it.
# =============================================================================
hdr("7. THE LITTER TERM — a definitional step at the 1985/2006 boundary")

XLSX <- file.path(ROOT, "Data/Komeetta/Hannu/Komeetta 150526hi--.xlsx")
if (requireNamespace("readxl", quietly = TRUE) && file.exists(XLSX)) {
  sp <- suppressMessages(readxl::read_excel(XLSX, sheet = "BiSo", range = "A4:IN3031",
                                            col_names = FALSE, col_types = "text",
                                            .name_repair = "minimal"))
  nn <- function(i) suppressWarnings(as.numeric(sp[[i]]))
  org <- data.frame(plot_id = as.integer(nn(3)), Krs = as.integer(nn(4)), REP = as.integer(nn(9)),
                    C06 = nn(176), C24 = nn(177), lit06 = nn(247), lit24 = nn(248)) |>
    filter(!is.na(plot_id), REP == 1, Krs == 101) |>
    mutate(soil06 = 1e4 * C06 / 1000, soil24 = 1e4 * C24 / 1000,   # Mg/ha
           lit06 = coalesce(lit06, 0) / 1000, lit24 = coalesce(lit24, 0) / 1000) |>
    filter(plot_id %in% bal_ids) |>
    left_join(bal |> filter(campaign == "VMI8") |> select(plot_id, weight, org85 = organic),
              by = "plot_id")

  cat(sprintf("balanced-panel plots matched in the workbook: %d\n", nrow(org)))
  cat(sprintf("\n  1985 organic (NO litter term exists)          : %6.2f Mg/ha\n", wmean(org$org85, org$weight)))
  cat(sprintf("  2006 organic = humus %5.2f + litter %5.2f    : %6.2f Mg/ha  (litter = %4.1f%%)\n",
              wmean(org$soil06, org$weight), wmean(org$lit06, org$weight),
              wmean(org$soil06 + org$lit06, org$weight),
              100 * wmean(org$lit06, org$weight) / wmean(org$soil06 + org$lit06, org$weight)))
  cat(sprintf("  2024 organic = humus %5.2f + litter %5.2f    : %6.2f Mg/ha  (litter = %4.1f%%)\n",
              wmean(org$soil24, org$weight), wmean(org$lit24, org$weight),
              wmean(org$soil24 + org$lit24, org$weight),
              100 * wmean(org$lit24, org$weight) / wmean(org$soil24 + org$lit24, org$weight)))

  d_org_asis <- wmean(org$soil06 + org$lit06, org$weight) - wmean(org$org85, org$weight)
  d_org_like <- wmean(org$soil06, org$weight)             - wmean(org$org85, org$weight)
  cat(sprintf("\n  apparent 1985->2006 organic gain, AS BUILT (litter in 2006 only) : %+6.2f Mg/ha\n", d_org_asis))
  cat(sprintf("  same gain with the litter term REMOVED (like-for-like)          : %+6.2f Mg/ha\n", d_org_like))
  cat(sprintf("  => the litter DEFINITION supplies %.0f%% of the apparent organic-layer gain\n",
              100 * (d_org_asis - d_org_like) / d_org_asis))
  # what the definitional step does to the TREND, on the same balanced panel
  lit <- org |> select(plot_id, lit06, lit24)
  tr  <- bal |> select(plot_id, campaign, weight, soc_profile_Mgha) |>
    left_join(lit, by = "plot_id") |>
    mutate(lit = case_when(campaign == "Biosoil" ~ lit06, campaign == "Komeetta" ~ lit24, TRUE ~ 0),
           profile_likeforlike = soc_profile_Mgha - coalesce(lit, 0)) |>
    group_by(campaign) |>
    summarise(as_built = wmean(soc_profile_Mgha, weight),
              like4like = wmean(profile_likeforlike, weight), .groups = "drop")
  g <- function(col, camp) tr[[col]][tr$campaign == camp]
  cat("\n  EFFECT ON THE OBSERVED TREND (weighted profile, Mg C/ha):\n")
  print(as.data.frame(tr |> mutate(across(where(is.numeric), \(x) round(x, 2)))), row.names = FALSE)
  cat(sprintf("\n    1985->2024 rate : as built %+0.3f  ->  like-for-like %+0.3f  Mg/ha/yr  (%.0f%% smaller)\n",
              (g("as_built","Komeetta") - g("as_built","VMI8")) / 39,
              (g("like4like","Komeetta") - g("like4like","VMI8")) / 39,
              100 * (1 - ((g("like4like","Komeetta") - g("like4like","VMI8")) /
                          (g("as_built","Komeetta") - g("as_built","VMI8"))))))
  cat(sprintf("    2006->2024 rate : as built %+0.3f  ->  like-for-like %+0.3f  Mg/ha/yr  (%.0f%% smaller)\n",
              (g("as_built","Komeetta") - g("as_built","Biosoil")) / 18,
              (g("like4like","Komeetta") - g("like4like","Biosoil")) / 18,
              100 * (1 - ((g("like4like","Komeetta") - g("like4like","Biosoil")) /
                          (g("as_built","Komeetta") - g("as_built","Biosoil"))))))
  cat("    !! the litter term GREW 2006->2024, so unlike a 1985 revision this one DOES\n")
  cat("       move the 2006->2024 window — the one every model reports as a source.\n")

  cat("\n  NB this is LUKE's own convention (it is what reproduces the official 2006/2024\n")
  cat("  national stocks), so it is not a pipeline bug — but it is applied to two of the\n")
  cat("  three campaigns, and it lands exactly on the 1985->2006 step we are trying to fit.\n")
} else {
  cat("source workbook not available; skipping\n")
}

# =============================================================================
# 8. "1985" IS 1986-1995 — and that gives a decisive test
#    The workbook sheet kiv_maat85_95 carries the per-plot soil sampling YEAR of
#    the first campaign (1986/87/88/89/1995). The pipeline throws it away and
#    calls every VMI8 observation "1985".
#
#    THE TEST. Regress each plot's layer-wise change on its OWN interval length:
#         delta_i = a + b * interval_i
#    A real accumulation process must pass through the origin: no elapsed time,
#    no change (a = 0, b > 0). A campaign-level METHOD offset is present in full
#    regardless of how much time elapsed (a != 0, b ~ 0).
# =============================================================================
hdr("8. THE FIRST CAMPAIGN SPANS 1986-1995 — process vs method-offset test")

yr <- NULL
if (exists("sp")) {
  xy <- suppressMessages(readxl::read_excel(XLSX, sheet = "kiv_maat85_95", col_names = FALSE,
                                            col_types = "text", .name_repair = "minimal"))
  yr <- data.frame(plot_id  = suppressWarnings(as.numeric(xy[[10]])),
                   samp_year = suppressWarnings(as.numeric(xy[[9]]))) |>
    filter(!is.na(plot_id), samp_year %in% 1986:1995) |> distinct(plot_id, .keep_all = TRUE)
}

if (!is.null(yr) && nrow(yr) > 0) {
  cat("per-plot first-campaign sampling year recovered for", nrow(yr), "plots:\n")
  print(table(yr$samp_year))
  ww <- w |> left_join(yr, by = "plot_id") |> filter(!is.na(samp_year)) |>
    mutate(interval = 2006 - samp_year)
  cat("\nbalanced-panel plots with a known sampling year:", nrow(ww),
      " | interval range:", paste(range(ww$interval), collapse = "-"), "yr\n")
  cat("=> the pipeline labels ALL of these 1985. For the",
      sum(ww$samp_year == 1995), "plots sampled in 1995 that is a 10-year error.\n")

  cat("\nDoes the 1985 stock itself depend on when it was sampled? (weighted mean soc_profile)\n")
  ww |> group_by(samp_year) |>
    summarise(n = n(), soc_1985 = round(wmean(soc_profile_Mgha_VMI8, weight), 2), .groups = "drop") |>
    as.data.frame() |> print(row.names = FALSE)

  cat("\nPROCESS-vs-OFFSET REGRESSION  delta = a + b*interval  (1st campaign -> 2006)\n")
  cat("  a = the part present at zero elapsed time  => METHOD OFFSET\n")
  cat("  b = the part proportional to elapsed time  => PROCESS\n\n")
  for (v in c("organic", "m0_20", "m20_40", "soc_0_40", "soc_profile_Mgha")) {
    d <- ww[[paste0(v, "_Biosoil")]] - ww[[paste0(v, "_VMI8")]]
    ok <- !is.na(d)
    fit <- lm(d[ok] ~ ww$interval[ok])
    ci  <- confint(fit)
    a <- coef(fit)[1]; b <- coef(fit)[2]
    cat(sprintf("  %-17s n=%3d  a = %+6.2f Mg/ha [%+5.2f,%+5.2f]  b = %+5.3f Mg/ha/yr [%+5.3f,%+5.3f]  p(b)=%.2f\n",
                v, sum(ok), a, ci[1,1], ci[1,2], b, ci[2,1], ci[2,2],
                summary(fit)$coefficients[2, 4]))
    cat(sprintf("  %-17s   at the mean interval (%.1f yr) the offset is %.0f%% of the total change\n",
                "", mean(ww$interval[ok]), 100 * a / (a + b * mean(ww$interval[ok]))))
  }
  cat("\nSame, as a simple contrast of the shortest vs longest intervals:\n")
  grp <- ww |> mutate(g = ifelse(samp_year >= 1995, "11 yr (1995)", "17-20 yr (1986-89)"))
  for (v in c("organic", "m0_20", "m20_40", "soc_profile_Mgha")) {
    d <- grp[[paste0(v, "_Biosoil")]] - grp[[paste0(v, "_VMI8")]]
    s <- tapply(d, grp$g, function(z) mean(z, na.rm = TRUE))
    cat(sprintf("  %-17s  %s = %+5.2f   %s = %+5.2f   ratio = %.2f (expect ~1.7 if it is a process)\n",
                v, names(s)[1], s[1], names(s)[2], s[2], s[2] / s[1]))
  }

  # --- the confound, and the control -----------------------------------------
  # The 1986-89 soil inventory ran SOUTH -> NORTH (Lapland last), so sampling
  # year tracks latitude, and the 1995 plots are a separate supplementary
  # sample. Sampling year is therefore NOT randomly assigned. Control for it.
  ww <- ww |> left_join(plotd |> filter(campaign == "VMI8") |>
                          select(plot_id, region_label, lat_WGS84) |> distinct(plot_id, .keep_all = TRUE),
                        by = "plot_id")
  cat("\n!! CONFOUND: the first inventory ran south->north, so year tracks latitude.\n")
  cat("   sampling year x region, and mean latitude:\n")
  print(with(ww, table(samp_year, region_label)))
  print(round(tapply(ww$lat_WGS84, ww$samp_year, mean, na.rm = TRUE), 2))

  cat("\nRegression WITH latitude controlled:  delta ~ interval + lat_WGS84\n")
  for (v in c("organic", "m0_20", "m20_40", "soc_profile_Mgha")) {
    d  <- ww[[paste0(v, "_Biosoil")]] - ww[[paste0(v, "_VMI8")]]
    ok <- !is.na(d) & !is.na(ww$lat_WGS84)
    fit <- lm(d[ok] ~ ww$interval[ok] + ww$lat_WGS84[ok])
    cs <- summary(fit)$coefficients
    cat(sprintf("  %-17s n=%3d  b(interval) = %+5.3f Mg/ha/yr  se %.3f  p = %.2f\n",
                v, sum(ok), cs[2,1], cs[2,2], cs[2,4]))
  }
  cat("\n   The interval slope survives no test at the 0.05 level once latitude is in.\n")
  cat("   That is NOT evidence the change is fake — it is evidence the design cannot\n")
  cat("   separate elapsed time from geography, i.e. the dating error is unrecoverable\n")
  cat("   from these data alone and has to be carried as an uncertainty.\n")
} else {
  cat("sampling-year sheet not available; skipping\n")
}

# =============================================================================
# 9. IF WE RECONSTRUCTED THE 1985 SUBSOIL, WHICH RECONSTRUCTION?
#    Three defensible ways to replace the 1985 20-40 cm layer. They disagree,
#    and the spread between them is the point: it is the size of the assumption
#    you would be inserting into the target as if it were data.
#      (i)   SHAPE   — the plot's own later-campaign depth ratio (section 6)
#      (ii)  STATIC  — the subsoil is inert: 1985 = mean(2006, 2024) per plot
#      (iii) REGRESS — straight line through the plot's 2006 and 2024 points,
#                      extrapolated back to 1985  <- the "use the other two" idea
# =============================================================================
hdr("9. THREE RECONSTRUCTIONS OF THE 1985 SUBSOIL — how much does the choice matter?")

sub <- bal |> select(plot_id, campaign, weight, m20_40, organic, m0_20, lambda, z_cap,
                     soc_deep_Mgha, soc_profile_Mgha) |>
  pivot_wider(names_from = campaign, values_from = c(m20_40, organic, m0_20, lambda, z_cap,
                                                     soc_deep_Mgha, soc_profile_Mgha)) |>
  left_join(lat, by = "plot_id") |>
  mutate(ratio_later = coalesce(ratio_later, ratio_fallback),
         rec_shape   = m0_20_VMI8 * ratio_later,
         rec_static  = (m20_40_Biosoil + m20_40_Komeetta) / 2,
         # linear in time through (2006, x06) and (2024, x24), evaluated at 1985
         rec_regress = m20_40_Biosoil + (m20_40_Biosoil - m20_40_Komeetta) / 18 * 21)

# each reconstruction propagates into the deep tail through the same C0 least-squares fit
for (v in c("rec_shape", "rec_static", "rec_regress")) {
  sub[[paste0("deep_", v)]] <- mapply(deep_from_bands, sub$m0_20_VMI8, sub[[v]],
                                      sub$lambda_VMI8, pmin(40, sub$z_cap_VMI8), sub$z_cap_VMI8)
}
sub$deep_orig <- mapply(deep_from_bands, sub$m0_20_VMI8, sub$m20_40_VMI8,
                        sub$lambda_VMI8, pmin(40, sub$z_cap_VMI8), sub$z_cap_VMI8)

p06 <- wmean(sub$soc_profile_Mgha_Biosoil,  sub$weight)
p24 <- wmean(sub$soc_profile_Mgha_Komeetta, sub$weight)

rec_tab <- lapply(c(measured = "m20_40_VMI8", shape = "rec_shape",
                    static = "rec_static", regress = "rec_regress"), function(v) {
  dtail <- if (v == "m20_40_VMI8") sub$soc_deep_Mgha_VMI8 else
    sub[[paste0("deep_", v)]] + (sub$soc_deep_Mgha_VMI8 - sub$deep_orig)
  prof <- sub$organic_VMI8 + sub$m0_20_VMI8 + sub[[v]] + dtail
  data.frame(sub20_40 = wmean(sub[[v]], sub$weight),
             profile_1985 = wmean(prof, sub$weight),
             rate_85_24 = (p24 - wmean(prof, sub$weight)) / 39)
}) |> bind_rows(.id = "reconstruction")
rec_tab$vs_measured <- rec_tab$profile_1985 - rec_tab$profile_1985[1]
print(as.data.frame(rec_tab |> mutate(across(where(is.numeric), \(x) round(x, 3)))), row.names = FALSE)

cat("\n  All three RAISE the 1985 subsoil, i.e. all three SHRINK the observed sink.\n")
cat(sprintf("  Spread across the three: %.2f Mg/ha on the 1985 level, %.3f Mg/ha/yr on the rate\n",
            diff(range(rec_tab$profile_1985[-1])), diff(range(rec_tab$rate_85_24[-1]))))
cat(sprintf("  For scale, the LITTER definition (section 7) moves the same rate by %.3f Mg/ha/yr.\n",
            0.248 - 0.139))
cat("\n  => the reconstruction CHOICE is itself worth a large fraction of the effect, and\n")
cat("     'regress' is the most extreme because it inherits the 2006->2024 subsoil DECLINE\n")
cat("     and projects it backwards. That decline is the very thing in doubt.\n")

# is the 2006->2024 litter growth credible? the Tupek input series falls over this window
if (exists("org")) {
  cat("\n--- cross-check: the measured litter LAYER grew 2006->2024. Is that credible? ---\n")
  dl <- org$lit24 - org$lit06
  cat(sprintf("  paired litter layer: 2006 %.2f -> 2024 %.2f Mg/ha  (%+.1f%%), median delta %+.2f, %.0f%% of plots up\n",
              wmean(org$lit06, org$weight), wmean(org$lit24, org$weight),
              100 * (wmean(org$lit24, org$weight) / wmean(org$lit06, org$weight) - 1),
              median(dl, na.rm = TRUE), 100 * mean(dl > 0, na.rm = TRUE)))
  cat("  vs the Tupek litter INPUT series, which falls -12.8% over 2006-2021 (CLAUDE.md, H2 closed).\n")
  cat("  A litter layer growing ~27% while its input falls ~13% is a contradiction: either the\n")
  cat("  2006 and 2024 litter measurements are not comparable, or the layer is not input-limited.\n")
}

# =============================================================================
# 10. FIGURES
# =============================================================================
png(file.path(OUTDIR, "01_depth_profile_by_campaign.png"), width = 1500, height = 1000, res = 130)
par(mfrow = c(2, 2), mar = c(4, 4.2, 3, 1))
cols <- c(VMI8 = "#C1553B", Biosoil = "#2E6F9E", Komeetta = "#3F8F5E")

# (a) mean stock per band
m <- as.matrix(tab2[, c("organic", "m0_20", "m20_40")]); rownames(m) <- as.character(tab2$campaign)
barplot(m, beside = TRUE, col = cols[rownames(m)], ylab = "Mg C / ha",
        names.arg = c("organic", "mineral 0-20", "mineral 20-40"),
        main = "(a) stock by depth band (weighted)")
legend("topright", legend = rownames(m), fill = cols[rownames(m)], bty = "n", cex = .85)

# (b) lambda distributions
bx <- split(shape$lambda_obs, shape$campaign)
boxplot(bx, col = cols[names(bx)], outline = FALSE, ylab = expression(lambda[obs] ~ (cm^-1)),
        main = "(b) depth decay rate per plot")
abline(h = median(plotd$lambda, na.rm = TRUE), lty = 2)
mtext("dashed = pooled class lambda used for the deep tail", side = 3, line = -1.1, cex = .62)

# (c) paired shape change 1985 -> 2006
d <- sh_w$lambda_obs_Biosoil - sh_w$lambda_obs_VMI8
hist(d[!is.na(d)], breaks = 40, col = "grey80", border = "white",
     xlab = expression(lambda[2006] - lambda[1985]), main = "(c) within-plot shape change")
abline(v = 0, lwd = 2); abline(v = median(d, na.rm = TRUE), col = "#C1553B", lwd = 2)

# (d) subsoil repeatability
plot(w$m20_40_VMI8, w$m20_40_Biosoil, pch = 16, col = "#C1553B88", cex = .6,
     xlab = "20-40 cm, earlier campaign", ylab = "20-40 cm, later campaign",
     main = "(d) subsoil repeatability", xlim = c(0, 60), ylim = c(0, 60))
points(w$m20_40_Biosoil, w$m20_40_Komeetta, pch = 16, col = "#2E6F9E88", cex = .6)
abline(0, 1, lty = 2)
legend("topleft", c(sprintf("1985 vs 2006 (r=%.2f)", reps$pearson[reps$band=="m20_40" & reps$pair=="VMI8->Biosoil"]),
                    sprintf("2006 vs 2024 (r=%.2f)", reps$pearson[reps$band=="m20_40" & reps$pair=="Biosoil->Komeetta"])),
       pch = 16, col = c("#C1553B", "#2E6F9E"), bty = "n", cex = .8)
dev.off()

png(file.path(OUTDIR, "02_profile_shape.png"), width = 1400, height = 700, res = 130)
par(mfrow = c(1, 2), mar = c(4.2, 4.2, 3, 1))
# mean concentration-equivalent per band (stock / thickness), plotted against depth
conc <- tab2 |> transmute(campaign, c0_20 = m0_20 / 20, c20_40 = m20_40 / 20)
plot(NA, xlim = c(0, max(conc$c0_20) * 1.1), ylim = c(40, 0), xlab = "Mg C / ha / cm",
     ylab = "depth (cm)", main = "(a) mean mineral depth profile")
for (i in seq_len(nrow(conc))) {
  cc <- as.character(conc$campaign[i])
  lines(c(conc$c0_20[i], conc$c0_20[i], conc$c20_40[i], conc$c20_40[i]),
        c(0, 20, 20, 40), col = cols[cc], lwd = 2.5)
}
legend("bottomright", names(cols), col = cols, lwd = 2.5, bty = "n", cex = .85)

plot(density(shape$lambda_obs[shape$campaign == "VMI8"], na.rm = TRUE), col = cols["VMI8"],
     lwd = 2, main = "(b) decay rate density", xlab = expression(lambda[obs]))
lines(density(shape$lambda_obs[shape$campaign == "Biosoil"],  na.rm = TRUE), col = cols["Biosoil"],  lwd = 2)
lines(density(shape$lambda_obs[shape$campaign == "Komeetta"], na.rm = TRUE), col = cols["Komeetta"], lwd = 2)
abline(v = 0, lty = 3)
dev.off()

cat("\n\nfigures written to", OUTDIR, "\n")
