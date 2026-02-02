# =============================================================================
# INPUT MAPPING FUNCTIONS (Cohort-based format)
# =============================================================================
#
# Maps from unified input format (AWEN by cohort + covariates) to 
# model-specific formats.
#
# Input data structure:
#   - C_foliage_A, C_foliage_W, C_foliage_E, C_foliage_N
#   - C_branches_A, C_branches_W, C_branches_E, C_branches_N
#   - C_stems_A, C_stems_W, C_stems_E, C_stems_N
#   - C_fine_roots_A, C_fine_roots_W, C_fine_roots_E, C_fine_roots_N
#   - C_coarse_roots_A, C_coarse_roots_W, C_coarse_roots_E, C_coarse_roots_N
#   - C_understorey_A, C_understorey_W, C_understorey_E, C_understorey_N
#   - temp_air, precip, evap, year, month
#
# Timestep handling:
#   - Yasso07:  Annual (aggregates monthly → annual), 3 litter types
#   - RothC:    Monthly (uses directly)
#   - Q-model:  Annual (aggregates monthly → annual), 5 litter types
#
# =============================================================================

library(dplyr)
library(tidyr)

# -----------------------------------------------------------------------------
# HELPER FUNCTIONS
# -----------------------------------------------------------------------------

#' Sum AWEN fractions for a cohort
#' @param df Data frame with C_cohort_A/W/E/N columns
#' @param cohort Name of cohort (foliage, branches, stems, fine_roots, coarse_roots, understorey)
sum_awen <- function(df, cohort) {
  cols <- paste0("C_", cohort, "_", c("A", "W", "E", "N"))
  cols <- cols[cols %in% names(df)]
  if (length(cols) == 0) return(rep(0, nrow(df)))
  rowSums(df[, cols, drop = FALSE], na.rm = TRUE)
}

#' Get AWEN vector for a cohort (for a single row or aggregated)
#' @param df Data frame
#' @param cohort Cohort name
get_awen <- function(df, cohort) {
  cols <- paste0("C_", cohort, "_", c("A", "W", "E", "N"))
  if (all(cols %in% names(df))) {
    c(sum(df[[cols[1]]], na.rm = TRUE),
      sum(df[[cols[2]]], na.rm = TRUE),
      sum(df[[cols[3]]], na.rm = TRUE),
      sum(df[[cols[4]]], na.rm = TRUE))
  } else {
    c(0, 0, 0, 0)
  }
}


# -----------------------------------------------------------------------------
# YASSO07 MAPPING (expects annual output with 3 litter types)
# 
# USED ALSO FOR YASSO15!
#
# -----------------------------------------------------------------------------

#' Map unified input to Yasso07 format
#'
#' Yasso07 needs 3 litter types based on wood size:
#'   - nwl: non-woody (foliage + understorey + fine_roots)
#'   - fwl: fine woody (branches)
#'   - cwl: coarse woody (stems + coarse_roots)
#'
#' Aggregates monthly to annual if needed.
#' 
#' @param input_df Monthly or annual input data with C_cohort_AWEN columns
#' @return Data frame with annual data: year, temp_mean, temp_amplitude, precip, nwl_A/W/E/N, fwl_A/W/E/N, cwl_A/W/E/N
map_input_to_yasso07 <- function(input_df) {
  
  is_monthly <- "month" %in% names(input_df)
  
  if (is_monthly) {
    years <- sort(unique(input_df$year))
    
    out_list <- lapply(years, function(y) {
      year_data <- input_df[input_df$year == y, ]
      
      # Climate: mean temp, amplitude, total precip
      temps <- year_data$temp_air
      temp_mean <- mean(temps, na.rm = TRUE)
      temp_amplitude <- (max(temps, na.rm = TRUE) - min(temps, na.rm = TRUE)) / 2
      precip <- sum(year_data$precip, na.rm = TRUE)
      
      # NWL: foliage + understorey + fine_roots (sum across months)
      nwl_A <- sum(year_data$C_foliage_A, na.rm = TRUE) + 
        sum(year_data$C_understorey_A, na.rm = TRUE) + 
        sum(year_data$C_fine_roots_A, na.rm = TRUE)
      nwl_W <- sum(year_data$C_foliage_W, na.rm = TRUE) + 
        sum(year_data$C_understorey_W, na.rm = TRUE) + 
        sum(year_data$C_fine_roots_W, na.rm = TRUE)
      nwl_E <- sum(year_data$C_foliage_E, na.rm = TRUE) + 
        sum(year_data$C_understorey_E, na.rm = TRUE) + 
        sum(year_data$C_fine_roots_E, na.rm = TRUE)
      nwl_N <- sum(year_data$C_foliage_N, na.rm = TRUE) + 
        sum(year_data$C_understorey_N, na.rm = TRUE) + 
        sum(year_data$C_fine_roots_N, na.rm = TRUE)
      
      # FWL: branches only
      fwl_A <- sum(year_data$C_branches_A, na.rm = TRUE)
      fwl_W <- sum(year_data$C_branches_W, na.rm = TRUE)
      fwl_E <- sum(year_data$C_branches_E, na.rm = TRUE)
      fwl_N <- sum(year_data$C_branches_N, na.rm = TRUE)
      
      # CWL: stems + coarse_roots
      cwl_A <- sum(year_data$C_stems_A, na.rm = TRUE) + 
        sum(year_data$C_coarse_roots_A, na.rm = TRUE)
      cwl_W <- sum(year_data$C_stems_W, na.rm = TRUE) + 
        sum(year_data$C_coarse_roots_W, na.rm = TRUE)
      cwl_E <- sum(year_data$C_stems_E, na.rm = TRUE) + 
        sum(year_data$C_coarse_roots_E, na.rm = TRUE)
      cwl_N <- sum(year_data$C_stems_N, na.rm = TRUE) + 
        sum(year_data$C_coarse_roots_N, na.rm = TRUE)
      
      data.frame(
        year = y,
        temp_mean = temp_mean,
        temp_amplitude = temp_amplitude,
        precip = precip,
        nwl_A = nwl_A, nwl_W = nwl_W, nwl_E = nwl_E, nwl_N = nwl_N,
        fwl_A = fwl_A, fwl_W = fwl_W, fwl_E = fwl_E, fwl_N = fwl_N,
        cwl_A = cwl_A, cwl_W = cwl_W, cwl_E = cwl_E, cwl_N = cwl_N
      )
    })
    
    out_df <- do.call(rbind, out_list)
    
  } else {
    # Annual input - just reorganize
    out_df <- data.frame(
      year = input_df$year,
      temp_mean = input_df$temp_air,
      temp_amplitude = if ("temp_amplitude" %in% names(input_df)) input_df$temp_amplitude else 12,
      precip = input_df$precip
    )
    
    # NWL
    out_df$nwl_A <- rowSums(input_df[, c("C_foliage_A", "C_understorey_A", "C_fine_roots_A"), drop = FALSE], na.rm = TRUE)
    out_df$nwl_W <- rowSums(input_df[, c("C_foliage_W", "C_understorey_W", "C_fine_roots_W"), drop = FALSE], na.rm = TRUE)
    out_df$nwl_E <- rowSums(input_df[, c("C_foliage_E", "C_understorey_E", "C_fine_roots_E"), drop = FALSE], na.rm = TRUE)
    out_df$nwl_N <- rowSums(input_df[, c("C_foliage_N", "C_understorey_N", "C_fine_roots_N"), drop = FALSE], na.rm = TRUE)
    
    # FWL
    out_df$fwl_A <- input_df$C_branches_A
    out_df$fwl_W <- input_df$C_branches_W
    out_df$fwl_E <- input_df$C_branches_E
    out_df$fwl_N <- input_df$C_branches_N
    
    # CWL
    out_df$cwl_A <- rowSums(input_df[, c("C_stems_A", "C_coarse_roots_A"), drop = FALSE], na.rm = TRUE)
    out_df$cwl_W <- rowSums(input_df[, c("C_stems_W", "C_coarse_roots_W"), drop = FALSE], na.rm = TRUE)
    out_df$cwl_E <- rowSums(input_df[, c("C_stems_E", "C_coarse_roots_E"), drop = FALSE], na.rm = TRUE)
    out_df$cwl_N <- rowSums(input_df[, c("C_stems_N", "C_coarse_roots_N"), drop = FALSE], na.rm = TRUE)
  }
  
  # Add totals for reference
  out_df$C_nwl <- out_df$nwl_A + out_df$nwl_W + out_df$nwl_E + out_df$nwl_N
  out_df$C_fwl <- out_df$fwl_A + out_df$fwl_W + out_df$fwl_E + out_df$fwl_N
  out_df$C_cwl <- out_df$cwl_A + out_df$cwl_W + out_df$cwl_E + out_df$cwl_N
  out_df$C_total <- out_df$C_nwl + out_df$C_fwl + out_df$C_cwl
  
  out_df
}


# -----------------------------------------------------------------------------
# ROTHC MAPPING (requires monthly input)
# -----------------------------------------------------------------------------

#' Map unified input to RothC format
#'
#' RothC uses monthly timestep natively.
#' Splits litter into DPM (decomposable) and RPM (resistant) based on AWEN:
#'   - W + E (water/ethanol soluble) → more labile → higher DPM fraction
#'   - A + N (acid soluble + non-soluble) → more recalcitrant → higher RPM fraction
#'
#' @param input_df Monthly time series with C_cohort_AWEN columns
#' @param site_row Optional: single row with clay content
#' @return Data frame with: year, month, C_DPM, C_RPM, temp_air, precip, evap, clay, soil_cover
map_input_to_rothc <- function(input_df, site_row = NULL) {
  
  if (!"month" %in% names(input_df)) {
    stop("RothC requires monthly input data")
  }
  
  # Get clay from site_row or default
  clay <- if (!is.null(site_row) && "clay" %in% names(site_row)) {
    as.numeric(site_row$clay)
  } else {
    warning("No clay value provided, using default 20%")
    20
  }
  
  # Calculate total C from all cohorts
  cohorts <- c("foliage", "branches", "stems", "fine_roots", "coarse_roots", "understorey")
  
  C_total <- rep(0, nrow(input_df))
  C_labile <- rep(0, nrow(input_df))  # W + E fractions
  C_recalc <- rep(0, nrow(input_df))  # A + N fractions
  
  for (coh in cohorts) {
    A_col <- paste0("C_", coh, "_A")
    W_col <- paste0("C_", coh, "_W")
    E_col <- paste0("C_", coh, "_E")
    N_col <- paste0("C_", coh, "_N")
    
    if (all(c(A_col, W_col, E_col, N_col) %in% names(input_df))) {
      C_total <- C_total + input_df[[A_col]] + input_df[[W_col]] + 
        input_df[[E_col]] + input_df[[N_col]]
      C_labile <- C_labile + input_df[[W_col]] + input_df[[E_col]]
      C_recalc <- C_recalc + input_df[[A_col]] + input_df[[N_col]]
    }
  }
  
  # Calculate DPM fraction based on labile content
  # DPM/RPM ratio typically 1.44 for grassland, 0.25 for forest
  # We use AWEN to estimate: more W+E = more DPM
  dpm_frac <- C_labile / pmax(C_total, 1e-6)
  # Constrain to reasonable range and add offset for forest
  dpm_frac <- pmax(0.15, pmin(0.5, dpm_frac + 0.1))
  
  data.frame(
    year = input_df$year,
    month = input_df$month,
    C_DPM = C_total * dpm_frac,
    C_RPM = C_total * (1 - dpm_frac),
    temp_air = input_df$temp_air,
    precip = input_df$precip,
    evap = if ("evap" %in% names(input_df)) input_df$evap else input_df$precip * 0.5,
    clay = clay,
    soil_cover = if ("soil_cover" %in% names(input_df)) input_df$soil_cover else 1
  )
}


# -----------------------------------------------------------------------------
# Q-MODEL MAPPING (expects annual output with 5 litter types)
# -----------------------------------------------------------------------------

#' Map unified input to Q-model format
#'
#' Q-model uses 5 litter types:
#'   - C_needles: from foliage
#'   - C_branches: from branches
#'   - C_stems: from stems
#'   - C_fine_roots: from fine_roots + coarse_roots
#'   - C_understorey: from understorey
#'
#' Aggregates monthly to annual if needed.
#'
#' @param input_df Monthly or annual input data
#' @return Data frame with annual: year, temp_mean, C_needles, C_branches, C_stems, C_fine_roots, C_understorey
map_input_to_q_model <- function(input_df) {
  
  is_monthly <- "month" %in% names(input_df)
  
  if (is_monthly) {
    years <- sort(unique(input_df$year))
    
    out_list <- lapply(years, function(y) {
      year_data <- input_df[input_df$year == y, ]
      
      # Climate
      temp_mean <- mean(year_data$temp_air, na.rm = TRUE)
      
      # Sum cohorts across months
      C_needles <- sum(sum_awen(year_data, "foliage"))
      C_branches <- sum(sum_awen(year_data, "branches"))
      C_stems <- sum(sum_awen(year_data, "stems"))
      C_fine_roots <- sum(sum_awen(year_data, "fine_roots")) + 
        sum(sum_awen(year_data, "coarse_roots"))
      C_understorey <- sum(sum_awen(year_data, "understorey"))
      
      data.frame(
        year = y,
        temp_mean = temp_mean,
        C_needles = C_needles,
        C_branches = C_branches,
        C_stems = C_stems,
        C_fine_roots = C_fine_roots,
        C_understorey = C_understorey
      )
    })
    
    out_df <- do.call(rbind, out_list)
    
  } else {
    # Annual input
    out_df <- data.frame(
      year = input_df$year,
      temp_mean = input_df$temp_air,
      C_needles = sum_awen(input_df, "foliage"),
      C_branches = sum_awen(input_df, "branches"),
      C_stems = sum_awen(input_df, "stems"),
      C_fine_roots = sum_awen(input_df, "fine_roots") + sum_awen(input_df, "coarse_roots"),
      C_understorey = sum_awen(input_df, "understorey")
    )
  }
  
  # Add total for reference
  out_df$C_total <- out_df$C_needles + out_df$C_branches + out_df$C_stems + 
    out_df$C_fine_roots + out_df$C_understorey
  
  out_df
}


# -----------------------------------------------------------------------------
# VALIDATION HELPER
# -----------------------------------------------------------------------------

#' Validate cohort-based input data
validate_input <- function(input_df, verbose = FALSE) {
  
  # Required climate columns
  required_climate <- c("year", "temp_air", "precip")
  missing_climate <- setdiff(required_climate, names(input_df))
  if (length(missing_climate) > 0) {
    stop(paste("Missing required climate columns:", paste(missing_climate, collapse = ", ")))
  }
  
  # Check for at least some litter columns
  c_cols <- grep("^C_.*_[AWEN]$", names(input_df), value = TRUE)
  if (length(c_cols) == 0) {
    stop("No litter input columns found (expected C_cohort_A/W/E/N format)")
  }
  
  # Check cohorts present
  cohorts_found <- unique(gsub("^C_(.*)_[AWEN]$", "\\1", c_cols))
  if (verbose){cat("Found cohorts:", paste(cohorts_found, collapse = ", "), "\n")}
  
  invisible(TRUE)
}



# # -----------------------------------------------------------------------------
# # YASSO20 MAPPING (STRICT MONTHLY INPUT REQUIRED)
# # -----------------------------------------------------------------------------
# 
# map_input_to_yasso20 <- function(input_df) {
#   
#   # ---- Required monthly forcing columns ----
#   req_monthly <- c("year", "month", "temp_air", "precip")
#   
#   # ---- Required litter columns (monthly inputs; will be summed to annual) ----
#   req_litter <- c(
#     "C_foliage_A","C_foliage_W","C_foliage_E","C_foliage_N",
#     "C_branches_A","C_branches_W","C_branches_E","C_branches_N",
#     "C_stems_A","C_stems_W","C_stems_E","C_stems_N",
#     "C_fine_roots_A","C_fine_roots_W","C_fine_roots_E","C_fine_roots_N",
#     "C_coarse_roots_A","C_coarse_roots_W","C_coarse_roots_E","C_coarse_roots_N",
#     "C_understorey_A","C_understorey_W","C_understorey_E","C_understorey_N"
#   )
#   
#   # ---- Identify plot/site id column if present (optional but recommended) ----
#   id_cols <- intersect(names(input_df), c("plot_id", "plot_ID", "site", "site_id"))
#   
#   # ---- Validate structure: must be monthly ----
#   if (!("month" %in% names(input_df))) {
#     stop("Yasso20 mapping requires MONTHLY input: missing column `month`.", call. = FALSE)
#   }
#   
#   missing_cols <- setdiff(c(req_monthly, req_litter), names(input_df))
#   if (length(missing_cols) > 0) {
#     stop(
#       paste0(
#         "Yasso20 mapping requires monthly forcing and litter columns.\nMissing columns: ",
#         paste(missing_cols, collapse = ", ")
#       ),
#       call. = FALSE
#     )
#   }
#   
#   # ---- Validate month range ----
#   bad_month <- input_df$month[!(input_df$month %in% 1:12)]
#   if (length(bad_month) > 0) {
#     stop("`month` must be integers in 1..12 for Yasso20 mapping.", call. = FALSE)
#   }
#   
#   # ---- Validate complete months per plot-year (no approximation allowed) ----
#   # We require exactly 12 months per (id, year).
#   key_vars <- c(id_cols, "year")
#   df_keyed <- input_df
#   
#   # helper to count unique months per group
#   month_counts <- df_keyed |>
#     dplyr::group_by(dplyr::across(dplyr::all_of(key_vars))) |>
#     dplyr::summarise(n_months = dplyr::n_distinct(.data$month), .groups = "drop")
#   
#   incomplete <- month_counts |> dplyr::filter(.data$n_months != 12)
#   if (nrow(incomplete) > 0) {
#     # show up to a few offending groups
#     preview <- utils::capture.output(utils::head(incomplete, 10))
#     stop(
#       paste0(
#         "Yasso20 mapping requires complete monthly forcing (12 months) for each plot-year.\n",
#         "Found plot-years with != 12 months. First rows:\n",
#         paste(preview, collapse = "\n")
#       ),
#       call. = FALSE
#     )
#   }
#   
#   # ---- Annual litter aggregation (sum of monthly inputs) ----
#   annual_litter <- input_df |>
#     dplyr::group_by(dplyr::across(dplyr::all_of(key_vars))) |>
#     dplyr::summarise(
#       # annual precip total (mm/yr)
#       precip = sum(.data$precip, na.rm = TRUE),
#       
#       # nwl: foliage + understorey + fine_roots
#       nwl_A = sum(.data$C_foliage_A + .data$C_understorey_A + .data$C_fine_roots_A, na.rm = TRUE),
#       nwl_W = sum(.data$C_foliage_W + .data$C_understorey_W + .data$C_fine_roots_W, na.rm = TRUE),
#       nwl_E = sum(.data$C_foliage_E + .data$C_understorey_E + .data$C_fine_roots_E, na.rm = TRUE),
#       nwl_N = sum(.data$C_foliage_N + .data$C_understorey_N + .data$C_fine_roots_N, na.rm = TRUE),
#       
#       # fwl: branches
#       fwl_A = sum(.data$C_branches_A, na.rm = TRUE),
#       fwl_W = sum(.data$C_branches_W, na.rm = TRUE),
#       fwl_E = sum(.data$C_branches_E, na.rm = TRUE),
#       fwl_N = sum(.data$C_branches_N, na.rm = TRUE),
#       
#       # cwl: stems + coarse_roots
#       cwl_A = sum(.data$C_stems_A + .data$C_coarse_roots_A, na.rm = TRUE),
#       cwl_W = sum(.data$C_stems_W + .data$C_coarse_roots_W, na.rm = TRUE),
#       cwl_E = sum(.data$C_stems_E + .data$C_coarse_roots_E, na.rm = TRUE),
#       cwl_N = sum(.data$C_stems_N + .data$C_coarse_roots_N, na.rm = TRUE),
#       
#       .groups = "drop"
#     )
#   
#   # ---- Monthly temperature to wide format T_1..T_12 (must exist) ----
#   annual_T <- input_df |>
#     dplyr::select(dplyr::all_of(c(key_vars, "month", "temp_air"))) |>
#     dplyr::distinct() |>
#     dplyr::mutate(month = as.integer(.data$month)) |>
#     tidyr::pivot_wider(names_from = month, values_from = temp_air, names_prefix = "T_")
#   
#   # Enforce presence of all months after pivot (should hold due to checks above)
#   missing_T <- setdiff(paste0("T_", 1:12), names(annual_T))
#   if (length(missing_T) > 0) {
#     stop(
#       paste0("Internal error: missing temperature columns after reshaping: ",
#              paste(missing_T, collapse = ", ")),
#       call. = FALSE
#     )
#   }
#   
#   annual_T <- annual_T |>
#     dplyr::select(dplyr::all_of(c(key_vars, paste0("T_", 1:12))))
#   
#   # ---- Combine ----
#   out_df <- annual_litter |>
#     dplyr::left_join(annual_T, by = key_vars) |>
#     dplyr::arrange(dplyr::across(dplyr::all_of(key_vars)))
#   
#   out_df
# }

# -----------------------------------------------------------------------------
# YASSO20 MAPPING (BASE R; strict monthly input required)
# -----------------------------------------------------------------------------
map_input_to_yasso20 <- function(input_df) {
  
  req_monthly <- c("year", "month", "temp_air", "precip")
  req_litter  <- c(
    "C_foliage_A","C_foliage_W","C_foliage_E","C_foliage_N",
    "C_branches_A","C_branches_W","C_branches_E","C_branches_N",
    "C_stems_A","C_stems_W","C_stems_E","C_stems_N",
    "C_fine_roots_A","C_fine_roots_W","C_fine_roots_E","C_fine_roots_N",
    "C_coarse_roots_A","C_coarse_roots_W","C_coarse_roots_E","C_coarse_roots_N",
    "C_understorey_A","C_understorey_W","C_understorey_E","C_understorey_N"
  )
  
  id_cols <- intersect(names(input_df), c("plot_id", "plot_ID", "site", "site_id"))
  key_vars <- c(id_cols, "year")
  
  missing_cols <- setdiff(c(req_monthly, req_litter), names(input_df))
  if (length(missing_cols) > 0) {
    stop("Missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  
  input_df$month <- as.integer(input_df$month)
  if (any(!input_df$month %in% 1:12)) stop("month must be 1..12", call. = FALSE)
  if (any(!is.finite(input_df$temp_air))) stop("temp_air has NA/Inf", call. = FALSE)
  
  # group factor (canonical order via levels)
  grp <- interaction(input_df[, key_vars, drop = FALSE], drop = TRUE, lex.order = TRUE)
  grp <- factor(grp)
  grp_lvls <- levels(grp)
  
  # require 12 months per group and no duplicates
  n_months <- tapply(input_df$month, grp, function(m) length(unique(m)))
  if (any(n_months != 12)) stop("Some plot-years missing months", call. = FALSE)
  if (any(duplicated(paste(grp, input_df$month, sep = "|")))) {
    stop("Duplicate rows for same plot-year-month", call. = FALSE)
  }
  
  # helper: annual sum of rowSums(cols) by group
  sum_cols_by_grp <- function(cols) {
    x <- rowSums(input_df[, cols, drop = FALSE], na.rm = TRUE)
    tapply(x, grp, sum, na.rm = TRUE)
  }
  
  precip_ann <- tapply(input_df$precip, grp, sum, na.rm = TRUE)
  
  nwl_A <- sum_cols_by_grp(c("C_foliage_A","C_understorey_A","C_fine_roots_A"))
  nwl_W <- sum_cols_by_grp(c("C_foliage_W","C_understorey_W","C_fine_roots_W"))
  nwl_E <- sum_cols_by_grp(c("C_foliage_E","C_understorey_E","C_fine_roots_E"))
  nwl_N <- sum_cols_by_grp(c("C_foliage_N","C_understorey_N","C_fine_roots_N"))
  
  fwl_A <- sum_cols_by_grp("C_branches_A")
  fwl_W <- sum_cols_by_grp("C_branches_W")
  fwl_E <- sum_cols_by_grp("C_branches_E")
  fwl_N <- sum_cols_by_grp("C_branches_N")
  
  cwl_A <- sum_cols_by_grp(c("C_stems_A","C_coarse_roots_A"))
  cwl_W <- sum_cols_by_grp(c("C_stems_W","C_coarse_roots_W"))
  cwl_E <- sum_cols_by_grp(c("C_stems_E","C_coarse_roots_E"))
  cwl_N <- sum_cols_by_grp(c("C_stems_N","C_coarse_roots_N"))
  
  # Monthly temperatures wide matrix [group x 12]
  Tmat <- matrix(NA_real_, nrow = length(grp_lvls), ncol = 12,
                 dimnames = list(grp_lvls, paste0("T_", 1:12)))
  Tmat[cbind(as.integer(grp), input_df$month)] <- input_df$temp_air
  if (anyNA(Tmat)) stop("Missing T_1..T_12 after reshaping", call. = FALSE)
  
  # Keys (first row per group), then optional stable sort by keys
  first_idx <- tapply(seq_len(nrow(input_df)), grp, function(ii) ii[1])
  keys_unique <- input_df[as.integer(first_idx[grp_lvls]), key_vars, drop = FALSE]
  ord <- do.call(order, keys_unique)
  keys_unique <- keys_unique[ord, , drop = FALSE]
  grp_lvls <- grp_lvls[ord]
  
  # Assemble output
  out_df <- data.frame(
    keys_unique,
    precip = as.numeric(precip_ann[grp_lvls]),
    nwl_A = as.numeric(nwl_A[grp_lvls]), nwl_W = as.numeric(nwl_W[grp_lvls]),
    nwl_E = as.numeric(nwl_E[grp_lvls]), nwl_N = as.numeric(nwl_N[grp_lvls]),
    fwl_A = as.numeric(fwl_A[grp_lvls]), fwl_W = as.numeric(fwl_W[grp_lvls]),
    fwl_E = as.numeric(fwl_E[grp_lvls]), fwl_N = as.numeric(fwl_N[grp_lvls]),
    cwl_A = as.numeric(cwl_A[grp_lvls]), cwl_W = as.numeric(cwl_W[grp_lvls]),
    cwl_E = as.numeric(cwl_E[grp_lvls]), cwl_N = as.numeric(cwl_N[grp_lvls]),
    as.data.frame(Tmat[grp_lvls, , drop = FALSE], check.names = FALSE)
  )
  rownames(out_df) <- NULL
  out_df
}
