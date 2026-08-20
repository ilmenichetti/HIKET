# =============================================================================
# preinit_input_shape.R  (C3, 2026-07-16)
#
# Growing-stock-derived shape for the transient pre-run litter input.
#
# The 1917->1985 pre-run previously ramped litter LINEARLY between the calibrated
# endpoints J_1917 and J_1985. But litter scales with standing biomass, and the
# Finnish growing stock was NOT linear over that period: a depleted, near-stationary
# base (~1400-1500 M m3) into the ~1970s, then a sustained rise (Korhonen et al.
# 2024, Silva Fennica 58(5) art. 24045, Fig 10a). This replaces the linear ramp
# with a shape that follows that trajectory.
#
# Returns a length-n_pre vector in [0,1]: 0 at preinit_year, 1 at t0_year, tracking
# the NFI total growing stock in between (constant extrapolation before NFI1). The
# calibrated endpoints J_1917/J_1985 are UNCHANGED; only the shape between moves.
#
# Source of truth: Data/forest_history/nfi_growing_stock.csv (endpoints exact from
# the paper; NFI2-NFI10 digitized from Fig 10a). Same file feeds figure F10b.
# =============================================================================

growing_stock_preinit_shape <- function(preinit_year = 1917L, t0_year = 1985L,
                                        n_pre = 68L,
                                        csv = "Data/forest_history/nfi_growing_stock.csv") {
  # Ablation switch (doublechecks/run_ablation.R): HIKET_PREINIT_LINEAR=1 returns
  # the OLD linear ramp, so C3 can be turned off for all six models from one place
  # without touching the run scripts. Unset => growing-stock shape, as in production.
  if (identical(Sys.getenv("HIKET_PREINIT_LINEAR"), "1")) {
    message("[C3 ABLATION] preinit shape forced LINEAR (HIKET_PREINIT_LINEAR=1)")
    return((seq_len(n_pre) - 1L) / (n_pre - 1L))
  }
  gs    <- read.csv(csv)
  years <- seq(preinit_year, t0_year, length.out = n_pre)
  # linear interpolation of the NFI series; rule=2 holds the NFI1 level flat back
  # to 1917 (pre-inventory, near-stationary) and would hold the last NFI forward.
  gs_at <- approx(gs$midpoint_year, gs$total_stock_Mm3, xout = years, rule = 2)$y
  shape <- (gs_at - gs_at[1L]) / (gs_at[n_pre] - gs_at[1L])
  # Clamp to [0,1]: growing stock dips ~2% below the 1917 level in the 1930s (well
  # within the +-50 M m3 digitization error), which a raw shape would turn into
  # litter BELOW the calibrated J_1917 -> near-zero/negative J for low sigma_init.
  # Clamping keeps J in [J_1917, J_1985] > 0 while preserving the near-stationary
  # base and the delayed post-1970 rise (the load-bearing feature).
  shape <- pmin(pmax(shape, 0), 1)
  shape[1L]     <- 0                     # pin endpoints exactly
  shape[n_pre]  <- 1
  shape
}


# =============================================================================
# liski_preinit_shape()  (2026-08-20)
#
# Pre-run shape from the RECONSTRUCTED INPUT TO SOIL of Liski et al. (2006),
# Ann. For. Sci. 63:687-697 -- Finland, same period, same model family, same NFI
# source; A. Lehtonen and M. Peltoniemi are coauthors.
#
# WHY THIS REPLACES THE GROWING-STOCK SHAPE. Standing volume is the wrong driver
# twice over: it is a STOCK, not an input, and it is a NATIONAL TOTAL, not per
# hectare. Liski et al. reconstructed the input directly, INCLUDING HARVEST
# RESIDUES -- the term that buffers the trajectory. Their own S4.4: large harvests
# decrease tree carbon but temporarily INCREASE litter and soil carbon, "as a
# result of these contrasting effects, the compounded carbon balance ... was less
# variable than that of any of the components alone". Our growing-stock shape has
# the sign backwards through the heavy-cutting decades: it says input sat at its
# 1917 minimum through 1922-1970, while the reconstruction peaks around 1960.
#
# COMPOSITION (decision 2026-08-20). Uses total_input_tree_basis = tree litter +
# harvest residues + natural mortality, matching the post-1985 Tupek product
# exactly (understorey EXCLUDED). Ground vegetation is deliberately left out: it
# is absorbed by sigma_input, as it is after 1985. Including it would make
# sigma_input mean one thing before the join and another after -- a discontinuity
# in a fitted parameter that no diagnostic here would reveal.
#
# AREA (decision 2026-08-20). Divided to a per-hectare basis using the upland soil
# area implied by Liski's own numbers (stock / density: 13.90 Mha 1922, 15.22 Mha
# 2004), with the growth placed in 1965-1980 following Korhonen et al. 2024
# ("since the mid 1960s ... most of this change was before the 1980s").
# ⚠ Constant area is NOT the neutral alternative -- their numbers falsify it.
# ⚠ Only the TIMING matters here: the shape is normalised, so sigma_init absorbs
#   the amplitude. Path is sensitive to the timing (mean 0.330-0.487); amplitude
#   is not (+14.6..+24.5%).
#
# ⚠ Residual bias, stated not fixed: the tree terms cover all forest land while we
#   divide by the UPLAND area, so the implied rise is an UPPER bound on the true
#   upland per-hectare rise. It runs against our own argument, so it is the
#   conservative direction.
#
# Data + provenance + validation: Data/liski2006_litter_input/README.md
# Contract is unchanged: length-68 vector on [0,1], 0 at 1917, 1 at 1985.
# =============================================================================
liski_preinit_shape <- function(preinit_year = 1917L, t0_year = 1985L, n_pre = 68L,
                                csv  = "Data/liski2006_litter_input/liski2006_fig5_input_to_soil.csv",
                                acsv = "Data/liski2006_litter_input/liski2006_upland_soil_area.csv") {
  d  <- read.csv(csv)
  ar <- read.csv(acsv)
  years <- seq(preinit_year, t0_year, length.out = n_pre)
  # input to soil, tree basis; rule=2 holds the 1922 value back to 1917 (pre-series)
  inp <- approx(d$year, d$total_input_tree_basis_TgC_yr, xout = years, rule = 2)$y
  # upland soil area: flat, then all growth 1965-1980, then flat (Korhonen dating)
  a0 <- ar$area_Mha[which.min(ar$year)]; a1 <- ar$area_Mha[which.max(ar$year)]
  a  <- rep(a0, n_pre)
  g  <- years >= 1965 & years <= 1980
  if (any(g)) a[g] <- seq(a0, a1, length.out = sum(g))
  a[years > 1980] <- a1
  per_ha <- inp / a
  shape  <- (per_ha - per_ha[1L]) / (per_ha[n_pre] - per_ha[1L])
  # Clamp as for the growing-stock shape: the per-hectare series is NOT monotone
  # (it peaks ~1960 on the residue pulse), so raw values can exceed [0,1].
  shape <- pmin(pmax(shape, 0), 1)
  shape[1L] <- 0; shape[n_pre] <- 1
  shape
}

# Dispatcher used by the six run_*_transient_calibration.R scripts.
#   HIKET_PREINIT_SHAPE = liski (default) | growing_stock | linear
# The older HIKET_PREINIT_LINEAR=1 still forces the linear ramp (C3 ablation).
preinit_shape <- function(...) {
  mode <- Sys.getenv("HIKET_PREINIT_SHAPE", "liski")
  if (identical(Sys.getenv("HIKET_PREINIT_LINEAR"), "1")) mode <- "linear"
  out <- switch(mode,
    liski         = liski_preinit_shape(...),
    growing_stock = growing_stock_preinit_shape(...),
    linear        = { n <- 68L; (seq_len(n) - 1L) / (n - 1L) },
    stop("preinit_shape: unknown HIKET_PREINIT_SHAPE '", mode, "'", call. = FALSE))
  message(sprintf("[PREINIT SHAPE] %s | mean %.3f | 1950 %.3f | 1970 %.3f",
                  mode, mean(out), out[34], out[54]))
  out
}
