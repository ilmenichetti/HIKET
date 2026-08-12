# =============================================================================
# observed_soc_basis.R — ONE canonical table of the observed SOC levels and
# trends, on every basis in circulation, so that no number is quoted without
# knowing which one it is.
#
# WHY. "The observed 2006-2024 trend" currently appears as +0.209, +0.139 and
# +0.111 in different documents. All three are correct; they differ in how plots
# are averaged and which plots are included. Mixing them inside one comparison
# is the error, not any single value.
#
# THE AXES
#   WEIGHTING  unweighted  = the average PLOT in our sample
#              weighted    = the average HECTARE of Finnish forest. The network
#                            samples the North at 1/3 the southern density, so
#                            the design weights are 3 (North) and 1 (South).
#                            LUKE's official national stocks are weighted.
#   PLOT SET   paired      = plots measured at BOTH ends of the interval. The only
#                            valid basis for a CHANGE: differencing means over
#                            different plot subsets is not an estimate of change.
#              balanced    = plots measured in all THREE campaigns (one common set
#                            for every interval, so intervals are comparable).
#              all         = every plot with a value in that campaign. Fine for a
#                            LEVEL, invalid for a change.
#   DEPTH      soc_0_40    = organic + MEASURED mineral 0-40 cm. What LUKE's official
#                            stocks and the Heikkinen/Ilvesniemi paper report.
#              soc_profile = soc_0_40 + the MODELLED tail from the deepest measured
#                            layer to z_cap. This is HIKET's calibration target.
#   OUTLIERS   soc_outlier plot-years are excluded throughout, as the pipeline does.
#
# WHICH TO USE
#   * comparing against MODEL trajectories -> UNWEIGHTED (the models are summarised
#     as mean(soc_mean) across plots, unweighted), PAIRED.
#   * any statement about FINLAND -> WEIGHTED.
#   * comparing intervals with each other -> BALANCED, so the plot set is constant.
#
# NB the denominators here are NOMINAL (39 / 21 / 18 yr). The first campaign is
# really 1986-1995, mean ~1989, so the true 1985->2024 interval is ~34.7 yr; that
# correction is separate and not applied here (see NEXT_SESSION.md sec. 0c).
#
# Run from repo root:  Rscript doublechecks/observed_soc_basis.R
# =============================================================================

suppressPackageStartupMessages({ library(dplyr); library(tidyr) })
options(width = 132)

ROOT <- if (dir.exists("Data/SOC_homogeneized")) "." else ".."
P <- read.csv(file.path(ROOT, "Data/SOC_homogeneized/soc_homogenized_plot.csv"))
DEPTH <- Sys.getenv("BASIS_DEPTH", "soc_profile_Mgha")   # or soc_0_40_Mgha
P <- P |> filter(!(soc_outlier %in% TRUE), is.finite(.data[[DEPTH]])) |>
  select(plot_id, campaign, weight, soc = all_of(DEPTH))
cat("DEPTH BASIS:", DEPTH, "\n")

W <- P |> pivot_wider(names_from = campaign, values_from = c(soc, weight))
wt <- coalesce(W$weight_VMI8, W$weight_Biosoil, W$weight_Komeetta)
set_paired  <- function(a, b) is.finite(W[[paste0("soc_", a)]]) & is.finite(W[[paste0("soc_", b)]])
set_balanced <- is.finite(W$soc_VMI8) & is.finite(W$soc_Biosoil) & is.finite(W$soc_Komeetta)

mn <- function(x, keep, weighted) {
  x <- x[keep]; w <- wt[keep]; ok <- is.finite(x)
  if (weighted) sum(x[ok] * w[ok]) / sum(w[ok]) else mean(x[ok])
}

hdr <- function(s) cat("\n", strrep("=", 118), "\n", s, "\n", strrep("=", 118), "\n", sep = "")

# =============================================================================
hdr("A. CAMPAIGN LEVELS  (Mg C/ha, whole profile = soc_profile_Mgha)")
lev <- lapply(c(FALSE, TRUE), function(wgt) {
  data.frame(weighting = ifelse(wgt, "weighted", "unweighted"),
             set = c("all available", "balanced (3 campaigns)"),
             n_1985 = c(sum(is.finite(W$soc_VMI8)), sum(set_balanced)),
             y1985 = c(mn(W$soc_VMI8, rep(TRUE, nrow(W)), wgt), mn(W$soc_VMI8, set_balanced, wgt)),
             y2006 = c(mn(W$soc_Biosoil, rep(TRUE, nrow(W)), wgt), mn(W$soc_Biosoil, set_balanced, wgt)),
             y2024 = c(mn(W$soc_Komeetta, rep(TRUE, nrow(W)), wgt), mn(W$soc_Komeetta, set_balanced, wgt)))
}) |> bind_rows()
print(lev |> mutate(across(where(is.numeric), \(x) round(x, 2))), row.names = FALSE)

# =============================================================================
hdr("B. TRENDS  (Mg C/ha/yr, nominal denominators)")
rows <- list(c("VMI8","Biosoil",21), c("Biosoil","Komeetta",18), c("VMI8","Komeetta",39))
out <- list()
for (r in rows) {
  a <- r[1]; b <- r[2]; yrs <- as.numeric(r[3])
  lab <- sprintf("%s->%s", c(VMI8="1985",Biosoil="2006",Komeetta="2024")[a],
                           c(VMI8="1985",Biosoil="2006",Komeetta="2024")[b])
  for (wgt in c(FALSE, TRUE)) for (st in c("paired", "balanced")) {
    keep <- if (st == "paired") set_paired(a, b) else set_balanced
    d <- W[[paste0("soc_", b)]] - W[[paste0("soc_", a)]]
    out[[length(out) + 1]] <- data.frame(
      interval = lab, set = st, weighting = ifelse(wgt, "weighted", "unweighted"),
      n = sum(keep & is.finite(d)), rate = mn(d, keep, wgt) / yrs)
  }
  # the invalid one, shown so it is recognisable when met in the wild
  for (wgt in c(FALSE, TRUE)) {
    la <- mn(W[[paste0("soc_", a)]], rep(TRUE, nrow(W)), wgt)
    lb <- mn(W[[paste0("soc_", b)]], rep(TRUE, nrow(W)), wgt)
    out[[length(out) + 1]] <- data.frame(
      interval = lab, set = "UNPAIRED (invalid)", weighting = ifelse(wgt, "weighted", "unweighted"),
      n = NA_integer_, rate = (lb - la) / yrs)
  }
}
tb <- bind_rows(out) |> mutate(rate = round(rate, 4))
print(tb |> pivot_wider(names_from = weighting, values_from = c(n, rate)) |>
        select(interval, set, n = n_unweighted, unweighted = rate_unweighted, weighted = rate_weighted),
      row.names = FALSE)

# =============================================================================
hdr("C. THE OTHER DEPTH BASIS  (re-run with BASIS_DEPTH=soc_0_40_Mgha for the full table)")
cat("  soc_0_40    = organic + measured mineral 0-40 cm  <- LUKE official / the Hannu-Juha paper\n")
cat("  soc_profile = the above + the modelled deep tail   <- HIKET calibration target\n")
cat("  The tail is ~16% of the profile, so the two bases give different LEVELS and,\n")
cat("  because the tail is refitted per campaign, different TRENDS.\n")

hdr("D. WHERE THE NUMBERS IN THE DOCUMENTS COME FROM")
cat("  +0.209  = 2006->2024, paired, UNWEIGHTED   (HIKET_next_session, M&M working doc)\n")
cat("  +0.312  = 1985->2024, paired, UNWEIGHTED   (same documents; consistent with the above)\n")
cat("  +0.111  = 2006->2024, LUKE official weighted set (59.06 -> 61.05 Mg/ha, n=446)\n")
cat("\n  => the manuscript pair (+0.209 / +0.312) is INTERNALLY CONSISTENT and is the right\n")
cat("     basis for model comparison, because the model trajectories are unweighted plot\n")
cat("     means. The clash is only with LUKE's weighted national figures, which must not\n")
cat("     appear in the same sentence as a model-observation comparison.\n")

hdr("E. RULE")
cat("  model comparison  -> paired, UNWEIGHTED     national statement -> WEIGHTED\n")
cat("  comparing intervals with each other -> BALANCED (constant plot set)\n")
cat("  never difference means over different plot subsets (the 'UNPAIRED' rows above)
  state the DEPTH basis whenever a level is quoted: 0-40 measured, or whole profile\n")
