# === Shared OBSERVED-SOC basis (source this wherever observations meet models) ===
# ONE place, so that every model-vs-observation figure compares the same population.
#
# WHY THIS EXISTS. The observed SOC level and trend depend on three choices that
# together span a 4x range for the same data (see doublechecks/observed_soc_basis.R):
#   weighting  unweighted = the average PLOT in our sample
#              weighted   = the average HECTARE of Finland (North sampled at 1/3
#                           density => design weight 3). LUKE's official stocks.
#   plot set   paired/balanced vs "every plot with a value that year"
#   depth      measured 0-40 cm vs whole profile incl. the modelled deep tail
#
# THE RULE APPLIED HERE. Model trajectories are UNWEIGHTED cross-plot means, so the
# observations they are drawn against must be unweighted too -- and on the SAME
# plots. We use the BALANCED set (plots observed in all three campaigns), because a
# model curve is a single line and cannot change population from year to year: with
# the balanced set, every point on the curve and every observed marker describe one
# fixed group of plots.
#
# WHAT IT FIXED (2026-08-12). F2/F3/F4 drew a model curve averaged over ALL plots
# against observed means taken over per-campaign subsets (404/456/409). On the
# companion appendix figure the same mismatch moved Yasso15's 1985 level 63.3 -> 65.1
# and its 1985-2024 rate +0.422 -> +0.379, inflating the apparent model-observation
# gap by ~40%. Levels here move less (~2-3%) but in the same direction.
#
# ⚠ NOT for national statements. Anything describing Finland must be region-WEIGHTED
# and on the measured 0-40 cm basis; that is a different number and does not belong
# on a model comparison plot.

# Plots observed in all three campaigns, from an input bundle's obs_meta.
balanced_plots <- function(om, n_campaigns = 3L) {
  keep <- vapply(om, function(z) length(z$soc_obs) >= n_campaigns, logical(1))
  as.integer(names(om)[keep])
}

# Observed campaign means (+/- 95% CI of the mean), restricted to `plots`.
# Returns data.frame(year, m, lo, hi, n).
obs_campaign_means <- function(om, plots = NULL) {
  nm <- names(om)
  if (!is.null(plots)) nm <- nm[as.integer(nm) %in% plots]
  obs <- do.call(rbind, lapply(nm, function(p) {
    z <- om[[p]]; if (!length(z$soc_obs)) return(NULL)
    data.frame(year = 1984L + z$idx, soc = z$soc_obs)
  }))
  a <- aggregate(soc ~ year, obs, function(x)
    c(m = mean(x), lo = mean(x) - 1.96*sd(x)/sqrt(length(x)),
      hi = mean(x) + 1.96*sd(x)/sqrt(length(x)), n = length(x)))
  data.frame(year = a$year, m = a$soc[, "m"], lo = a$soc[, "lo"],
             hi = a$soc[, "hi"], n = a$soc[, "n"])
}

# One-line provenance string for figure captions / logs.
basis_note <- function(plots)
  sprintf("basis: balanced plot set (n = %d), unweighted, whole profile", length(plots))
