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

# --- CAMPAIGN collapse (added 2026-08-17) ------------------------------------
# obs_campaign_means() aggregates by TRUE observation year, and since the 1985
# dating correction the first campaign spans 1986/87/88/89/95 -- so it returns
# SEVEN rows, not three. Every figure that draws "three campaign markers" and
# indexed them as year == 1985 broke on that (F4 died with a zero-length
# subscript; the label vector c("VMI8","Biosoil","Komeetta") silently recycled).
#
# THE RULE. A campaign is ONE marker, placed at the MEAN TRUE YEAR of its own
# observations -- ~1989 for VMI8, not 1985. That is the honest x-position: it is
# when the campaign actually happened, and it makes every rate denominator right
# without further correction (1985->2024 becomes 35.1 yr, not 39).
# Use campaign INDEX (1/2/3) to address them, never a literal year.
campaign_of <- function(y) ifelse(y <= 2000L, 1L, ifelse(y <= 2015L, 2L, 3L))

obs_campaigns <- function(om, plots = NULL) {
  nm <- names(om)
  if (!is.null(plots)) nm <- nm[as.integer(nm) %in% plots]
  obs <- do.call(rbind, lapply(nm, function(p) {
    z <- om[[p]]; if (!length(z$soc_obs)) return(NULL)
    data.frame(year = 1984L + z$idx, soc = z$soc_obs)
  }))
  obs$camp <- campaign_of(obs$year)
  do.call(rbind, lapply(sort(unique(obs$camp)), function(k) {
    x <- obs[obs$camp == k, ]
    se <- sd(x$soc) / sqrt(nrow(x))
    data.frame(camp = k, year = mean(x$year), m = mean(x$soc),
               lo = mean(x$soc) - 1.96 * se, hi = mean(x$soc) + 1.96 * se,
               n = nrow(x))
  }))
}

CAMPAIGN_LABELS <- c("VMI8", "Biosoil", "Komeetta")

# --- WHY COLLAPSING IS SAFE, AND WHEN IT WOULD STOP BEING (checked 2026-08-17) ---
# The CALIBRATION indexes every observation on its own obs_year, so the likelihood
# compares each plot to the model in the year that plot was actually sampled. A figure
# that collapses VMI8 to one marker is only consistent with that if the two agree.
# MEASURED, balanced set: mean-over-plots of the model AT EACH PLOT'S OWN YEAR vs the
# model curve at the campaign mean year differ by <= 0.6 tC/ha in all six models
# (SP1 0.6, TP2 0.1, TP3 0.1, Y07 0.2, Y15 0.0, Y20 0.4), against a model-observation
# gap of -9 to -12. So the marker is not distorting the comparison.
# ⚠ CONDITIONAL, not general: it holds because the plot-count-weighted mean of the VMI8
# sampling years is ~1989 -- exactly where the marker sits -- and the trajectory is
# near-linear across 1986-1995, making the discrepancy second-order in the curvature.
# Re-check if the marker is ever placed anywhere but the weighted mean year, or if a
# model's curve becomes strongly bent over that decade.
#
# ⚠ DO NOT plot the VMI8 sampling-year groups as separate markers. They differ by 18.3
# tC/ha at first visit (p<0.001) -- but by 37.2 in 2006 and 34.2 in 2024, years when ALL
# of those plots were measured simultaneously. The between-group differences are a PLOT
# SET (geographic) effect, not a temporal one: VMI8 was worked through region by region.
# Drawn as a time series they read as a 16 tC/ha rise-and-fall inside the campaign that
# does not exist. If the sub-year structure must be shown, plot each group as an ANOMALY
# FROM ITS OWN 2006 VALUE, which differences the plot-set effect out.

# Per-sampling-year means WITHIN one campaign (VMI8 is the only one with sub-structure).
# Drawn in a distinct lighter BLUE, as UNCONNECTED points, precisely because each is a
# different SUBSET of plots -- see the warning above. Never join them with a line: the
# apparent trend between them is geographic, not temporal.
SUBYEAR_COL <- "#5B9BD5"

obs_subyears <- function(om, plots = NULL, camp = 1L) {
  nm <- names(om)
  if (!is.null(plots)) nm <- nm[as.integer(nm) %in% plots]
  obs <- do.call(rbind, lapply(nm, function(p) {
    z <- om[[p]]; if (!length(z$soc_obs)) return(NULL)
    data.frame(year = 1984L + z$idx, soc = z$soc_obs)
  }))
  obs <- obs[campaign_of(obs$year) == camp, ]
  do.call(rbind, lapply(sort(unique(obs$year)), function(y) {
    x <- obs$soc[obs$year == y]; se <- sd(x) / sqrt(length(x))
    data.frame(year = y, m = mean(x), lo = mean(x) - 1.96 * se,
               hi = mean(x) + 1.96 * se, n = length(x))
  }))
}

# One-line provenance string for figure captions / logs.
basis_note <- function(plots)
  sprintf("basis: balanced plot set (n = %d), unweighted, whole profile", length(plots))
