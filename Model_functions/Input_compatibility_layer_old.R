# =============================================================================
# INPUT MAPPING FUNCTIONS (Cohort-based format)
# =============================================================================
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
#   Q-model:     annual climate (mean T)
#                annual litter (5 cohorts, total C)
#
# =============================================================================


# =============================================================================
# HELPERS
# =============================================================================

get_id_cols <- function(df) {
  intersect(names(df), c("plot_id", "plot_ID", "site", "site_id"))
}

sum_awen <- function(df, cohort) {
  cols <- paste0("C_", cohort, "_", c("A", "W", "E", "N"))
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
    # soil_cover: always 1 for forest (no fallow); override per-row if column exists
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
        nwl_A = sum(d$C_foliage_A, na.rm=TRUE) + sum(d$C_understorey_A, na.rm=TRUE) + sum(d$C_fine_roots_A, na.rm=TRUE),
        nwl_W = sum(d$C_foliage_W, na.rm=TRUE) + sum(d$C_understorey_W, na.rm=TRUE) + sum(d$C_fine_roots_W, na.rm=TRUE),
        nwl_E = sum(d$C_foliage_E, na.rm=TRUE) + sum(d$C_understorey_E, na.rm=TRUE) + sum(d$C_fine_roots_E, na.rm=TRUE),
        nwl_N = sum(d$C_foliage_N, na.rm=TRUE) + sum(d$C_understorey_N, na.rm=TRUE) + sum(d$C_fine_roots_N, na.rm=TRUE),
        fwl_A = sum(d$C_branches_A, na.rm=TRUE),
        fwl_W = sum(d$C_branches_W, na.rm=TRUE),
        fwl_E = sum(d$C_branches_E, na.rm=TRUE),
        fwl_N = sum(d$C_branches_N, na.rm=TRUE),
        cwl_A = sum(d$C_stems_A, na.rm=TRUE) + sum(d$C_coarse_roots_A, na.rm=TRUE),
        cwl_W = sum(d$C_stems_W, na.rm=TRUE) + sum(d$C_coarse_roots_W, na.rm=TRUE),
        cwl_E = sum(d$C_stems_E, na.rm=TRUE) + sum(d$C_coarse_roots_E, na.rm=TRUE),
        cwl_N = sum(d$C_stems_N, na.rm=TRUE) + sum(d$C_coarse_roots_N, na.rm=TRUE)
      )))
    }))
  } else {
    .rn(cbind(
      plot_df[, id_cols, drop = FALSE],
      data.frame(
        year  = plot_df$year,
        nwl_A = rowSums(plot_df[, c("C_foliage_A","C_understorey_A","C_fine_roots_A")],  na.rm=TRUE),
        nwl_W = rowSums(plot_df[, c("C_foliage_W","C_understorey_W","C_fine_roots_W")],  na.rm=TRUE),
        nwl_E = rowSums(plot_df[, c("C_foliage_E","C_understorey_E","C_fine_roots_E")],  na.rm=TRUE),
        nwl_N = rowSums(plot_df[, c("C_foliage_N","C_understorey_N","C_fine_roots_N")],  na.rm=TRUE),
        fwl_A = plot_df$C_branches_A,
        fwl_W = plot_df$C_branches_W,
        fwl_E = plot_df$C_branches_E,
        fwl_N = plot_df$C_branches_N,
        cwl_A = rowSums(plot_df[, c("C_stems_A","C_coarse_roots_A")], na.rm=TRUE),
        cwl_W = rowSums(plot_df[, c("C_stems_W","C_coarse_roots_W")], na.rm=TRUE),
        cwl_E = rowSums(plot_df[, c("C_stems_E","C_coarse_roots_E")], na.rm=TRUE),
        cwl_N = rowSums(plot_df[, c("C_stems_N","C_coarse_roots_N")], na.rm=TRUE)
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
  sum_ann <- function(cols) {
    x <- rowSums(plot_df[, cols, drop = FALSE], na.rm = TRUE)
    as.numeric(tapply(x, plot_df$year, sum, na.rm = TRUE)[as.character(years)])
  }
  .rn(cbind(
    meta[rep(1, length(years)), , drop = FALSE],
    data.frame(
      year  = years,
      nwl_A = sum_ann(c("C_foliage_A","C_understorey_A","C_fine_roots_A")),
      nwl_W = sum_ann(c("C_foliage_W","C_understorey_W","C_fine_roots_W")),
      nwl_E = sum_ann(c("C_foliage_E","C_understorey_E","C_fine_roots_E")),
      nwl_N = sum_ann(c("C_foliage_N","C_understorey_N","C_fine_roots_N")),
      fwl_A = sum_ann("C_branches_A"),
      fwl_W = sum_ann("C_branches_W"),
      fwl_E = sum_ann("C_branches_E"),
      fwl_N = sum_ann("C_branches_N"),
      cwl_A = sum_ann(c("C_stems_A","C_coarse_roots_A")),
      cwl_W = sum_ann(c("C_stems_W","C_coarse_roots_W")),
      cwl_E = sum_ann(c("C_stems_E","C_coarse_roots_E")),
      cwl_N = sum_ann(c("C_stems_N","C_coarse_roots_N"))
    )
  ))
}

#' Annual litter inputs for Yasso20: plot_id, year, nwl/fwl/cwl AWEN
#' Climate: use map_climate_monthly()
map_inputs_yasso20 <- function(input_df) {
  id_cols <- get_id_cols(input_df)
  if (length(id_cols) == 0) return(.inputs_yasso20_single(input_df, id_cols))
  .by_plot(input_df, id_cols, .inputs_yasso20_single)
}


# =============================================================================
# ROTHC  (climate: map_climate_monthly(site_df=))
# =============================================================================

.inputs_rothc_single <- function(plot_df, id_cols) {
  cohorts  <- c("foliage","branches","stems","fine_roots","coarse_roots","understorey")
  C_labile <- rep(0, nrow(plot_df))
  C_recalc <- rep(0, nrow(plot_df))
  for (coh in cohorts) {
    for (fr in c("W","E")) {
      col <- paste0("C_", coh, "_", fr)
      if (col %in% names(plot_df)) C_labile <- C_labile + plot_df[[col]]
    }
    for (fr in c("A","N")) {
      col <- paste0("C_", coh, "_", fr)
      if (col %in% names(plot_df)) C_recalc <- C_recalc + plot_df[[col]]
    }
  }
  C_total  <- C_labile + C_recalc
  dpm_frac <- C_labile / pmax(C_total, 1e-6)
  dpm_frac <- pmax(0.15, pmin(0.5, dpm_frac + 0.1))
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
#' Clay, depth and soil_cover are NOT included here -- they belong in the
#' climate data frame produced by map_climate_monthly(site_df=).
#'
#' @param input_df  Monthly data frame with AWEN cohort columns.
map_inputs_rothc <- function(input_df) {
  if (!"month" %in% names(input_df))
    stop("RothC requires monthly input data", call. = FALSE)
  id_cols <- get_id_cols(input_df)
  if (length(id_cols) == 0) return(.inputs_rothc_single(input_df, id_cols))
  .by_plot(input_df, id_cols, .inputs_rothc_single)
}


# =============================================================================
# Q-MODEL
# =============================================================================

.climate_q_model_single <- function(plot_df, id_cols) {
  meta <- plot_df[1, id_cols, drop = FALSE]
  if ("month" %in% names(plot_df)) {
    years <- sort(unique(plot_df$year))
    do.call(rbind, lapply(years, function(y) {
      d <- plot_df[plot_df$year == y, ]
      .rn(cbind(meta, data.frame(year = y, temp_mean = mean(d$temp_air, na.rm=TRUE))))
    }))
  } else {
    .rn(cbind(
      plot_df[, id_cols, drop = FALSE],
      data.frame(year = plot_df$year, temp_mean = plot_df$temp_air)
    ))
  }
}

.inputs_q_model_single <- function(plot_df, id_cols) {
  meta <- plot_df[1, id_cols, drop = FALSE]
  if ("month" %in% names(plot_df)) {
    years <- sort(unique(plot_df$year))
    do.call(rbind, lapply(years, function(y) {
      d <- plot_df[plot_df$year == y, ]
      .rn(cbind(meta, data.frame(
        year          = y,
        C_needles     = sum(sum_awen(d, "foliage")),
        C_branches    = sum(sum_awen(d, "branches")),
        C_stems       = sum(sum_awen(d, "stems")),
        C_fine_roots  = sum(sum_awen(d, "fine_roots")) + sum(sum_awen(d, "coarse_roots")),
        C_understorey = sum(sum_awen(d, "understorey"))
      )))
    }))
  } else {
    .rn(cbind(
      plot_df[, id_cols, drop = FALSE],
      data.frame(
        year          = plot_df$year,
        C_needles     = sum_awen(plot_df, "foliage"),
        C_branches    = sum_awen(plot_df, "branches"),
        C_stems       = sum_awen(plot_df, "stems"),
        C_fine_roots  = sum_awen(plot_df, "fine_roots") + sum_awen(plot_df, "coarse_roots"),
        C_understorey = sum_awen(plot_df, "understorey")
      )
    ))
  }
}

#' Annual climate for Q-model: plot_id, year, temp_mean
map_climate_q_model <- function(input_df) {
  id_cols <- get_id_cols(input_df)
  if (length(id_cols) == 0) return(.climate_q_model_single(input_df, id_cols))
  .by_plot(input_df, id_cols, .climate_q_model_single)
}

#' Annual litter inputs for Q-model: plot_id, year, C_needles/branches/stems/fine_roots/understorey
map_inputs_q_model <- function(input_df) {
  id_cols <- get_id_cols(input_df)
  if (length(id_cols) == 0) return(.inputs_q_model_single(input_df, id_cols))
  .by_plot(input_df, id_cols, .inputs_q_model_single)
}


# =============================================================================
# VALIDATION HELPER
# =============================================================================

#' Validate cohort-based input data
validate_input <- function(input_df, verbose = FALSE) {
  required_climate <- c("year", "temp_air", "precip")
  missing_climate  <- setdiff(required_climate, names(input_df))
  if (length(missing_climate) > 0)
    stop("Missing required climate columns: ", paste(missing_climate, collapse=", "))
  c_cols <- grep("^C_.*_[AWEN]$", names(input_df), value=TRUE)
  if (length(c_cols) == 0)
    stop("No litter input columns found (expected C_cohort_A/W/E/N format)")
  cohorts_found <- unique(gsub("^C_(.*)_[AWEN]$", "\\1", c_cols))
  if (verbose) cat("Found cohorts:", paste(cohorts_found, collapse=", "), "\n")
  invisible(TRUE)
}