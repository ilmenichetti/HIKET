# =============================================================================
# INPUT MAPPING FUNCTIONS (3-class AWEN format)
# =============================================================================
#
# Unified input format uses three litter classes with AWEN fractions:
#   C_nwl_A/W/E/N  -- non-woody litter (foliage + fine roots + understorey)
#   C_fwl_A/W/E/N  -- fine woody litter (branches)
#   C_cwl_A/W/E/N  -- coarse woody litter (stems + coarse roots)
#
# Climate and input mapping are fully separated. Each model has:
#   map_climate_<model>()  -- climate forcing at model-native resolution
#   map_inputs_<model>()   -- litter carbon inputs at model-native resolution
#
# Yasso20 and RothC share map_climate_monthly() for climate.
# RothC additionally needs clay, depth and soil_cover in the climate data frame;
# these come from site_df and are ignored by Yasso20.
#
# All output data frames carry plot_id as the first column.
#
# Model-specific resolutions:
#   Yasso07/15:  annual climate (mean T, amplitude, precip)
#                annual litter (nwl/fwl/cwl AWEN)
#   Yasso20:     monthly climate via map_climate_monthly()
#                annual litter (nwl/fwl/cwl AWEN)
#   RothC:       monthly climate via map_climate_monthly(site_df=)
#                monthly litter (C_DPM, C_RPM only)
#
# =============================================================================


# =============================================================================
# HELPERS
# =============================================================================

get_id_cols <- function(df) {
  intersect(names(df), c("plot_id", "plot_ID", "site", "site_id"))
}

#' Sum all four AWEN fractions for a litter class.
#' class should be one of "nwl", "fwl", "cwl".
sum_awen_class <- function(df, class) {
  cols <- paste0("C_", class, "_", c("A", "W", "E", "N"))
  cols <- cols[cols %in% names(df)]
  if (length(cols) == 0) return(rep(0, nrow(df)))
  rowSums(df[, cols, drop = FALSE], na.rm = TRUE)
}

.rn <- function(df) { row.names(df) <- NULL; df }

.by_plot <- function(input_df, id_cols, FUN, ...) {
  plots  <- split(input_df, input_df[[id_cols[1]]], drop = TRUE)
  result <- do.call(rbind, lapply(plots, FUN, id_cols = id_cols, ...))
  .rn(result)
}

# Look up a per-plot scalar from site_df; returns default if not found
.site_val <- function(pid, site_df, col, default, warn_msg = NULL) {
  if (!is.null(site_df) && col %in% names(site_df)) {
    id_site <- get_id_cols(site_df)
    if (length(id_site) > 0) {
      row <- site_df[site_df[[id_site[1]]] == pid, ]
      if (nrow(row) > 0) return(as.numeric(row[[col]][1]))
    }
  }
  if (!is.null(warn_msg)) warning(warn_msg)
  default
}

# Check that all required litter columns exist in a data frame
.check_litter_cols <- function(df, context = "") {
  required <- paste0("C_", rep(c("nwl","fwl","cwl"), each=4), "_",
                     rep(c("A","W","E","N"), times=3))
  missing <- setdiff(required, names(df))
  if (length(missing) > 0)
    stop(context, "Missing litter columns: ", paste(missing, collapse=", "),
         call. = FALSE)
  invisible(TRUE)
}


# =============================================================================
# MONTHLY CLIMATE  (shared by Yasso20 and RothC)
# =============================================================================

#' Monthly climate data frame.
#'
#' Columns always present: plot_id, year, month, temp_air, precip, evap
#' Columns added when site_df is supplied (used by RothC, ignored by Yasso20):
#'   clay       -- clay content [%]
#'   depth      -- soil depth [cm]  (from site_df$soil_depth)
#'   soil_cover -- plant cover (1 = vegetated, 0 = bare); always 1 for forest
#'
#' @param input_df  Monthly data frame with temp_air, precip, evap columns.
#' @param site_df   Optional site data frame with plot_id, clay, soil_depth.
map_climate_monthly <- function(input_df, site_df = NULL) {
  if (!"month" %in% names(input_df))
    stop("map_climate_monthly requires monthly input data", call. = FALSE)

  id_cols <- get_id_cols(input_df)

  base <- .rn(cbind(
    input_df[, id_cols, drop = FALSE],
    data.frame(
      year     = input_df$year,
      month    = input_df$month,
      temp_air = input_df$temp_air,
      precip   = input_df$precip,
      evap     = if ("evap" %in% names(input_df)) input_df$evap
                 else input_df$precip * 0.5
    )
  ))

  # Add site properties for RothC if site_df supplied
  if (!is.null(site_df)) {
    pid_col <- id_cols[1]
    base$clay  <- vapply(base[[pid_col]], function(pid)
      .site_val(pid, site_df, "clay",       default = 20,
                warn_msg = paste("No clay for plot", pid, "-- using 20%")),
      numeric(1))
    base$depth <- vapply(base[[pid_col]], function(pid)
      .site_val(pid, site_df, "soil_depth", default = 23,
                warn_msg = paste("No soil_depth for plot", pid, "-- using 23 cm")),
      numeric(1))
    base$soil_cover <- if ("soil_cover" %in% names(input_df))
      as.integer(input_df$soil_cover)
    else
      1L
  }

  base
}


# =============================================================================
# YASSO07 / YASSO15
# =============================================================================

.climate_yasso07_single <- function(plot_df, id_cols) {
  meta <- plot_df[1, id_cols, drop = FALSE]
  if ("month" %in% names(plot_df)) {
    years <- sort(unique(plot_df$year))
    do.call(rbind, lapply(years, function(y) {
      d <- plot_df[plot_df$year == y, ]
      .rn(cbind(meta, data.frame(
        year           = y,
        temp_mean      = mean(d$temp_air, na.rm = TRUE),
        temp_amplitude = (max(d$temp_air, na.rm = TRUE) -
                          min(d$temp_air, na.rm = TRUE)) / 2,
        precip         = sum(d$precip, na.rm = TRUE)
      )))
    }))
  } else {
    .rn(cbind(
      plot_df[, id_cols, drop = FALSE],
      data.frame(
        year           = plot_df$year,
        temp_mean      = plot_df$temp_air,
        temp_amplitude = if ("temp_amplitude" %in% names(plot_df))
                           plot_df$temp_amplitude else 12,
        precip         = plot_df$precip
      )
    ))
  }
}

.inputs_yasso07_single <- function(plot_df, id_cols) {
  meta <- plot_df[1, id_cols, drop = FALSE]
  if ("month" %in% names(plot_df)) {
    years <- sort(unique(plot_df$year))
    do.call(rbind, lapply(years, function(y) {
      d <- plot_df[plot_df$year == y, ]
      .rn(cbind(meta, data.frame(
        year  = y,
        nwl_A = sum(d$C_nwl_A, na.rm = TRUE),
        nwl_W = sum(d$C_nwl_W, na.rm = TRUE),
        nwl_E = sum(d$C_nwl_E, na.rm = TRUE),
        nwl_N = sum(d$C_nwl_N, na.rm = TRUE),
        fwl_A = sum(d$C_fwl_A, na.rm = TRUE),
        fwl_W = sum(d$C_fwl_W, na.rm = TRUE),
        fwl_E = sum(d$C_fwl_E, na.rm = TRUE),
        fwl_N = sum(d$C_fwl_N, na.rm = TRUE),
        cwl_A = sum(d$C_cwl_A, na.rm = TRUE),
        cwl_W = sum(d$C_cwl_W, na.rm = TRUE),
        cwl_E = sum(d$C_cwl_E, na.rm = TRUE),
        cwl_N = sum(d$C_cwl_N, na.rm = TRUE)
      )))
    }))
  } else {
    .rn(cbind(
      plot_df[, id_cols, drop = FALSE],
      data.frame(
        year  = plot_df$year,
        nwl_A = plot_df$C_nwl_A,
        nwl_W = plot_df$C_nwl_W,
        nwl_E = plot_df$C_nwl_E,
        nwl_N = plot_df$C_nwl_N,
        fwl_A = plot_df$C_fwl_A,
        fwl_W = plot_df$C_fwl_W,
        fwl_E = plot_df$C_fwl_E,
        fwl_N = plot_df$C_fwl_N,
        cwl_A = plot_df$C_cwl_A,
        cwl_W = plot_df$C_cwl_W,
        cwl_E = plot_df$C_cwl_E,
        cwl_N = plot_df$C_cwl_N
      )
    ))
  }
}

#' Annual climate for Yasso07/15: plot_id, year, temp_mean, temp_amplitude, precip
map_climate_yasso07 <- function(input_df) {
  id_cols <- get_id_cols(input_df)
  if (length(id_cols) == 0) return(.climate_yasso07_single(input_df, id_cols))
  .by_plot(input_df, id_cols, .climate_yasso07_single)
}

#' Annual litter inputs for Yasso07/15: plot_id, year, nwl/fwl/cwl AWEN
map_inputs_yasso07 <- function(input_df) {
  .check_litter_cols(input_df, "map_inputs_yasso07: ")
  id_cols <- get_id_cols(input_df)
  if (length(id_cols) == 0) return(.inputs_yasso07_single(input_df, id_cols))
  .by_plot(input_df, id_cols, .inputs_yasso07_single)
}

map_climate_yasso15 <- map_climate_yasso07
map_inputs_yasso15  <- map_inputs_yasso07


# =============================================================================
# YASSO20  (climate: map_climate_monthly)
# =============================================================================

.inputs_yasso20_single <- function(plot_df, id_cols) {
  meta  <- plot_df[1, id_cols, drop = FALSE]
  years <- sort(unique(plot_df$year))
  sum_ann <- function(col) {
    as.numeric(tapply(plot_df[[col]], plot_df$year, sum, na.rm = TRUE)[as.character(years)])
  }
  .rn(cbind(
    meta[rep(1, length(years)), , drop = FALSE],
    data.frame(
      year  = years,
      nwl_A = sum_ann("C_nwl_A"),
      nwl_W = sum_ann("C_nwl_W"),
      nwl_E = sum_ann("C_nwl_E"),
      nwl_N = sum_ann("C_nwl_N"),
      fwl_A = sum_ann("C_fwl_A"),
      fwl_W = sum_ann("C_fwl_W"),
      fwl_E = sum_ann("C_fwl_E"),
      fwl_N = sum_ann("C_fwl_N"),
      cwl_A = sum_ann("C_cwl_A"),
      cwl_W = sum_ann("C_cwl_W"),
      cwl_E = sum_ann("C_cwl_E"),
      cwl_N = sum_ann("C_cwl_N")
    )
  ))
}

#' Annual litter inputs for Yasso20: plot_id, year, nwl/fwl/cwl AWEN
#' Climate: use map_climate_monthly()
map_inputs_yasso20 <- function(input_df) {
  .check_litter_cols(input_df, "map_inputs_yasso20: ")
  id_cols <- get_id_cols(input_df)
  if (length(id_cols) == 0) return(.inputs_yasso20_single(input_df, id_cols))
  .by_plot(input_df, id_cols, .inputs_yasso20_single)
}


# =============================================================================
# ROTHC  (climate: map_climate_monthly(site_df=))
# =============================================================================
#
# DPM/RPM split is based on litter class chemistry:
#   W + E fractions are labile  -> DPM
#   A + N fractions are recalcitrant -> RPM
#
# A small correction (+0.1) shifts the raw labile fraction upward to account
# for microbial processing, then the result is clamped to [0.15, 0.50].
# This is applied uniformly across all three litter classes.
# =============================================================================

.inputs_rothc_single <- function(plot_df, id_cols) {
  # Labile (W+E) and recalcitrant (A+N) summed across all three classes
  C_labile <- rowSums(plot_df[, c("C_nwl_W","C_nwl_E",
                                   "C_fwl_W","C_fwl_E",
                                   "C_cwl_W","C_cwl_E"), drop = FALSE],
                      na.rm = TRUE)
  C_recalc <- rowSums(plot_df[, c("C_nwl_A","C_nwl_N",
                                   "C_fwl_A","C_fwl_N",
                                   "C_cwl_A","C_cwl_N"), drop = FALSE],
                      na.rm = TRUE)
  C_total  <- C_labile + C_recalc
  dpm_frac <- C_labile / pmax(C_total, 1e-6)
  dpm_frac <- pmax(0.15, pmin(0.50, dpm_frac + 0.10))
  .rn(cbind(
    plot_df[, id_cols, drop = FALSE],
    data.frame(
      year  = plot_df$year,
      month = plot_df$month,
      C_DPM = C_total * dpm_frac,
      C_RPM = C_total * (1 - dpm_frac)
    )
  ))
}

#' Monthly litter inputs for RothC: plot_id, year, month, C_DPM, C_RPM
#'
#' Clay, depth and soil_cover belong in the climate data frame produced by
#' map_climate_monthly(site_df=), not here.
#'
#' @param input_df  Monthly data frame with C_nwl/fwl/cwl_A/W/E/N columns.
map_inputs_rothc <- function(input_df) {
  if (!"month" %in% names(input_df))
    stop("RothC requires monthly input data", call. = FALSE)
  .check_litter_cols(input_df, "map_inputs_rothc: ")
  id_cols <- get_id_cols(input_df)
  if (length(id_cols) == 0) return(.inputs_rothc_single(input_df, id_cols))
  .by_plot(input_df, id_cols, .inputs_rothc_single)
}


# =============================================================================
# VALIDATION HELPER
# =============================================================================

#' Validate 3-class AWEN input data
#'
#' Checks that all required litter and climate columns are present.
#' @param input_df   Input data frame.
#' @param verbose    If TRUE, prints found columns.
validate_input <- function(input_df, verbose = FALSE) {
  required_climate <- c("year", "temp_air", "precip")
  missing_climate  <- setdiff(required_climate, names(input_df))
  if (length(missing_climate) > 0)
    stop("Missing required climate columns: ", paste(missing_climate, collapse = ", "))

  .check_litter_cols(input_df)

  if (verbose) {
    litter_cols <- grep("^C_(nwl|fwl|cwl)_[AWEN]$", names(input_df), value = TRUE)
    cat("Litter columns found:", paste(litter_cols, collapse = ", "), "\n")
    cat("Total litter C:",
        sum(rowSums(input_df[, litter_cols], na.rm = TRUE)), "t C ha-1\n")
  }
  invisible(TRUE)
}


# =============================================================================
# SYNTHETIC DATA GENERATOR  (for testing)
# =============================================================================

#' Generate a minimal synthetic input data frame for testing.
#'
#' @param n_plots   Number of plots.
#' @param n_years   Number of years per plot.
#' @param monthly   If TRUE, generate monthly rows (12 per year).
#' @param seed      Random seed.
make_synthetic_input <- function(n_plots = 2, n_years = 5,
                                 monthly = TRUE, seed = 42) {
  set.seed(seed)
  rows <- list()
  for (p in seq_len(n_plots)) {
    pid <- paste0("plot_", p)
    for (y in 2000 + seq_len(n_years)) {
      months <- if (monthly) 1:12 else 1L
      for (m in months) {
        # Simple seasonal temperature curve
        t_air <- 5 + 10 * sin((m - 3) * pi / 6) + rnorm(1, 0, 0.5)
        rows[[length(rows) + 1]] <- data.frame(
          plot_id  = pid,
          year     = y,
          month    = m,
          temp_air = t_air,
          precip   = runif(1, 30, 80),
          evap     = runif(1, 10, 40),
          # NWL: most of the litter
          C_nwl_A  = runif(1, 0.005, 0.020),
          C_nwl_W  = runif(1, 0.002, 0.008),
          C_nwl_E  = runif(1, 0.001, 0.004),
          C_nwl_N  = runif(1, 0.003, 0.012),
          # FWL: smaller amounts
          C_fwl_A  = runif(1, 0.001, 0.005),
          C_fwl_W  = runif(1, 0.001, 0.003),
          C_fwl_E  = runif(1, 0.000, 0.002),
          C_fwl_N  = runif(1, 0.001, 0.004),
          # CWL: smallest
          C_cwl_A  = runif(1, 0.001, 0.004),
          C_cwl_W  = runif(1, 0.000, 0.002),
          C_cwl_E  = runif(1, 0.000, 0.001),
          C_cwl_N  = runif(1, 0.001, 0.003),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  do.call(rbind, rows)
}
