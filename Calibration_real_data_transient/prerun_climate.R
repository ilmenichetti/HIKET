# =============================================================================
# prerun_climate.R -- the climate of the 1917-1984 pre-run (2026-10-12)
# -----------------------------------------------------------------------------
# HIKET_PRERUN_CLIM_YEARS=<y0>-<y1>, e.g. 1961-1984, sets the window whose MEAN
# climate drives the pre-run (on Roihu: SINGULARITYENV_HIKET_PRERUN_CLIM_YEARS).
# Unset => the first STEADY_STATE_YEARS of the model climate (1985-2004), i.e.
# production unchanged.
#
# WHY. The pre-run applies ONE constant climate modifier (xi_for_ss) to all of
# 1917-1984. Taken from 1985-2004 it carries warming the pre-run predates (plot
# means: 1961-84 2.13 C, 1985-2004 2.67 C), and the attribution's climate part
# is then measured against years the pre-run never saw (it sums to ~0 until
# ~2004). A pre-1985 window gives ONE long-term mean, used for the pre-run AND
# as the attribution's reference climate (decided 2026-10-09).
#
# SOURCE. Data/model_inputs/climate_monthly.csv: the same FMI 10 km extraction
# the model climate comes from (identical values from 1985 on), 1961-2025. Its
# pre-1985 years were simply dropped by the merge onto the litter rows.
#
# HOW. Every call site that took clim[seq_len(n_ss), ] calls
# hiket_ss_clim(pid, clim, n_ss). The frames travel in the input bundle
# (climate_ss_by_plot), so the predictive stage and the attribution replay a
# posterior with ITS OWN pre-run climate whatever the shell says; old bundles
# lack the field => old behaviour. litter_means$precip_mean (Yasso15/20 pre-run
# leaching term, inactive at leac = 0) is replaced with the same window.
#
# NOT CHANGED: the litter side of the pre-run (1985-89 endpoint, Liski shape),
# the forward climate, and the doublechecks/ scripts that slice clim themselves.
# =============================================================================

.hiket_prerun_years <- local({
  s <- Sys.getenv("HIKET_PRERUN_CLIM_YEARS", "")
  if (!nzchar(s)) NULL else {
    y <- as.integer(strsplit(s, "-", fixed = TRUE)[[1]])
    if (length(y) != 2L || anyNA(y) || y[1] > y[2] || y[1] < 1961L)
      stop("HIKET_PRERUN_CLIM_YEARS must be '<y0>-<y1>' within 1961-2025, got '", s, "'")
    # 1985-2004 is allowed ON PURPOSE: it must reproduce production exactly (the
    # identity test in doublechecks/prerun_climate_test.R). Any post-1984 year is
    # otherwise climate the pre-run predates.
    if (y[2] >= 1985L) message("[PRERUN CLIMATE] NOTE: window reaches past 1984 (test use only)")
    y
  }
})
if (!is.null(.hiket_prerun_years))
  message(sprintf("[PRERUN CLIMATE] ON -- pre-run climate = %d-%d mean (not the first steady-state years)",
                  .hiket_prerun_years[1], .hiket_prerun_years[2]))

# Registry read by hiket_ss_clim(); set once per script, inherited by mclapply forks.
.hiket_prerun_env <- new.env()
.hiket_prerun_env$clim <- NULL

# Per-plot pre-run climate frames in the model's own format, or NULL when off.
#   mapper: map_climate_yasso07 (annual: SP1/TP2/TP3/Yasso07/Yasso15) or
#           map_climate_monthly (monthly: Yasso20).
hiket_prerun_climate <- function(plots, mapper,
                                 file = "./Data/model_inputs/climate_monthly.csv") {
  if (is.null(.hiket_prerun_years)) return(NULL)
  y <- .hiket_prerun_years
  cm <- read.csv(file)
  cm <- cm[cm$year >= y[1] & cm$year <= y[2] &
           as.character(cm$plot_id) %in% plots, ]
  cm <- cm[order(cm$plot_id, cm$year, cm$month), ]
  # Complete record or stop: a gap would silently shift the mean.
  n_need <- 12L * (y[2] - y[1] + 1L)
  n_have <- table(factor(as.character(cm$plot_id), levels = plots))
  bad <- names(n_have)[n_have != n_need]
  if (length(bad))
    stop(sprintf("[PRERUN CLIMATE] %d plot(s) lack a complete %d-%d record, e.g. %s",
                 length(bad), y[1], y[2], paste(head(bad, 5), collapse = ", ")))
  if (anyNA(cm[, c("temp_air", "precip")]))
    stop("[PRERUN CLIMATE] NA in temp_air/precip over the pre-run window")
  mapped <- mapper(cm)
  out <- split(mapped, as.character(mapped$plot_id))[plots]
  out <- lapply(out, function(d) { rownames(d) <- NULL; d })
  message(sprintf("[PRERUN CLIMATE] %d plots x %d yr | mean T %.2f C",
                  length(out), y[2] - y[1] + 1L,
                  mean(vapply(out, function(d) mean(d[[intersect(c("temp_mean", "temp_air"), names(d))[1]]]),
                              numeric(1)))))
  out
}

hiket_set_prerun_clim <- function(x, plots = NULL) {
  if (!is.null(x) && !is.null(plots)) {
    miss <- setdiff(plots, names(x))
    if (length(miss))
      stop(sprintf("[PRERUN CLIMATE] bundle lacks %d plot(s), e.g. %s",
                   length(miss), paste(head(miss, 5), collapse = ", ")))
  }
  .hiket_prerun_env$clim <- x
  invisible(x)
}

# The climate frame the pre-run mean is computed from.
hiket_ss_clim <- function(pid, clim, n_ss) {
  pc <- .hiket_prerun_env$clim
  if (is.null(pc)) return(clim[seq_len(n_ss), , drop = FALSE])
  d <- pc[[as.character(pid)]]
  if (is.null(d)) stop("[PRERUN CLIMATE] no pre-run climate for plot ", pid)
  d
}

# Replace litter_means$precip_mean (Yasso15/20) with the pre-run window's mean
# annual total. No-op when off, or for models whose litter_means carry none.
hiket_prerun_precip <- function(litter_means) {
  pc <- .hiket_prerun_env$clim
  if (is.null(pc)) return(litter_means)
  for (pid in names(litter_means)) {
    if (is.null(litter_means[[pid]]$precip_mean)) next
    d <- pc[[pid]]
    litter_means[[pid]]$precip_mean <- if ("month" %in% names(d))
      sum(tapply(d$precip, d$month, mean, na.rm = TRUE))   # Yasso20 convention
    else mean(d$precip)                                    # Yasso15 convention
  }
  litter_means
}

# Label for metadata / logs.
hiket_prerun_label <- function(x = .hiket_prerun_env$clim) {
  if (is.null(x)) "first_steady_state_years" else {
    yr <- range(x[[1]]$year); sprintf("%d-%d", yr[1], yr[2])
  }
}
