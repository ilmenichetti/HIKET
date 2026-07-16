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
