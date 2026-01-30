# =============================================================================
# TESTING
# =============================================================================
#
# This script is for testing the whole SOC model comparison framework.
# It generates synthetic data, runs the models and compares results.
#
# =============================================================================

# Load required libraries
library(dplyr)
library(tidyr)

#load multi-model wrapper functions
source("./Models_functions/multi_model.R")


# =============================================================================
# STEP 1: Generate Synthetic Input Data
# =============================================================================

# Generate monthly input data for 2 plots, 2000-2020
generate_input_template <- function(plot_ids = c("PLOT_001", "PLOT_002"),
                                    years = 2000:2020) {
  
  grid <- expand.grid(
    plot_id = plot_ids,
    year = years,
    month = 1:12,
    stringsAsFactors = FALSE
  )
  grid <- grid[order(grid$plot_id, grid$year, grid$month), ]
  n <- nrow(grid)
  
  # AWEN proportions by cohort
  awen_props <- list(
    foliage      = c(A = 0.52, W = 0.10, E = 0.02, N = 0.36),
    branches     = c(A = 0.50, W = 0.05, E = 0.03, N = 0.42),
    stems        = c(A = 0.45, W = 0.02, E = 0.01, N = 0.52),
    fine_roots   = c(A = 0.44, W = 0.05, E = 0.01, N = 0.50),
    coarse_roots = c(A = 0.45, W = 0.03, E = 0.01, N = 0.51),
    understorey  = c(A = 0.50, W = 0.15, E = 0.03, N = 0.32)
  )
  
  # Annual litter totals (Mg C/ha/yr)
  litter_totals <- list(
    foliage = 1.2, branches = 0.4, stems = 0.1,
    fine_roots = 0.8, coarse_roots = 0.3, understorey = 0.3
  )
  
  cohorts <- names(litter_totals)
  
  # Generate litter columns (monthly = annual/12)
  for (coh in cohorts) {
    total <- litter_totals[[coh]] / 12
    props <- awen_props[[coh]]
    
    for (frac in c("A", "W", "E", "N")) {
      col_name <- paste0("C_", coh, "_", frac)
      base_value <- total * props[frac]
      
      # Variation by plot and year
      plot_effect <- ifelse(grid$plot_id == plot_ids[1], 1.0, 1.1)
      year_effect <- 1 + 0.1 * sin(2 * pi * (grid$year - min(years)) / 10)
      
      grid[[col_name]] <- round(base_value * plot_effect * year_effect, 6)
    }
  }
  
  # Climate
  base_temp <- ifelse(grid$plot_id == plot_ids[1], 3.5, 2.0)
  grid$temp_air <- round(base_temp + 12 * sin(2 * pi * (grid$month - 4) / 12) + rnorm(n, 0, 1), 1)
  
  monthly_precip <- c(40, 35, 35, 40, 45, 55, 70, 75, 60, 55, 50, 45)
  grid$precip <- round(monthly_precip[grid$month] * 
                         ifelse(grid$plot_id == plot_ids[1], 1.0, 1.1) * 
                         (1 + rnorm(n, 0, 0.2)), 1)
  
  grid$evap <- round(pmax(5, 20 + 40 * sin(2 * pi * (grid$month - 4) / 12) + rnorm(n, 0, 5)), 1)
  
  grid
}

# Generate site data
generate_site_template <- function(plot_ids = c("PLOT_001")) {  # Only one plot ID
  data.frame(
    plot_id = plot_ids,
    pine_prop = 0.3,      # Single value
    spruce_prop = 0.6,
    birch_prop = 0.1,
    clay = 18,
    soil_depth = 23,
    stringsAsFactors = FALSE
  )
}


# Generate data
input_data <- generate_input_template(plot_ids = "PLOT_001", years = 2000:2020)
site_data <- generate_site_template(plot_ids = "PLOT_001")

cat(sprintf("Generated input data: %d rows (months)\n", nrow(input_data)))
cat(sprintf("Time span: %d-%d\n", min(input_data$year), max(input_data$year)))


# -----------------------------------------------------
# Climate of the scenario
# -----------------------------------------------------

boxplot(input_data$temp_air ~ input_data$month,
        main = "Synthetic monthly temperature (°C)",
        xlab = "Month", ylab = "Temp (°C)")

boxplot(input_data$precip ~ input_data$month,
        main = "Synthetic precipitation (mm)",
        xlab = "Month", ylab = "Precipitation (mm)")

boxplot(input_data$evap ~ input_data$month,
        main = "Synthetic evaporation (mm)",
        xlab = "Month", ylab = "ET (mm)")


# -----------------------------------------------------
# Monthly litter totals by source (A + W + E + N) of the scenario
# -----------------------------------------------------

# Helper to compute total litter per cohort
total_litter <- function(df, cohort) {
  df[[paste0("C_", cohort, "_A")]] +
    df[[paste0("C_", cohort, "_W")]] +
    df[[paste0("C_", cohort, "_E")]] +
    df[[paste0("C_", cohort, "_N")]]
}

cohorts <- c("foliage", "branches", "stems",
             "fine_roots", "coarse_roots", "understorey")

par(mfrow = c(3, 2), mar = c(4, 4, 3, 1))

for (coh in cohorts) {
  boxplot(
    total_litter(input_data, coh) ~ input_data$month,
    main = paste("Synthetic monthly litter –", coh),
    xlab = "Month",
    ylab = "Litter input (Mg C/ha/month)"
  )
}

par(mfrow = c(1, 1))


#all cohorts together
cols <- c("forestgreen", "sienna", "grey40",
          "orange", "brown", "darkolivegreen")

boxplot(
  lapply(cohorts, function(coh) total_litter(input_data, coh)),
  names = cohorts,
  las = 2,
  col = cols,
  ylab = "Mean monthly litter (Mg C/ha/month)",
  main = "Relative magnitude of litter inputs by cohort"
)


# =============================================================================
# STEP 2: Run Models
# =============================================================================

# Run all three models
results <- run_soc_models(
  input_df = input_data,
  site_row = site_data,
  models = c("yasso07", "yasso15", "yasso20", "rothc", "q_model"),
  spinup_years = 200  # Shorter spinup for demonstration
)

cat("\nModels completed successfully!\n")



# =============================================================================
# STEP 3: Examine Results
# =============================================================================

print(head(results$summary, 10))
cat("\n")
print(tail(results$summary, 10))



# =============================================================================
# ADDITIONAL ANALYSES (ALL MODELS + 2 FIGURES)
# =============================================================================

# ---- Helper: annual input (Mg C/ha/yr) for each model ----
get_annual_input <- function(model_name, input_data, site_data) {
  
  if (model_name %in% c("yasso07", "yasso15")) {
    mapped <- map_input_to_yasso07(input_data)
    return(mean(mapped$C_nwl + mapped$C_fwl + mapped$C_cwl, na.rm = TRUE))
  }
  
  if (model_name == "yasso20") {
    mapped <- map_input_to_yasso20(input_data)  # strict monthly forcing required
    nwl <- mapped$nwl_A + mapped$nwl_W + mapped$nwl_E + mapped$nwl_N
    fwl <- mapped$fwl_A + mapped$fwl_W + mapped$fwl_E + mapped$fwl_N
    cwl <- mapped$cwl_A + mapped$cwl_W + mapped$cwl_E + mapped$cwl_N
    return(mean(nwl + fwl + cwl, na.rm = TRUE))
  }
  
  if (model_name == "q_model") {
    mapped <- map_input_to_q_model(input_data)
    return(mean(mapped$C_needles + mapped$C_branches + mapped$C_stems +
                  mapped$C_fine_roots + mapped$C_understorey, na.rm = TRUE))
  }
  
  if (model_name == "rothc") {
    mapped <- map_input_to_rothc(input_data, site_data)
    annual <- aggregate(mapped$C_DPM + mapped$C_RPM, by = list(year = mapped$year), FUN = sum)
    return(mean(annual$x, na.rm = TRUE))
  }
  
  stop(sprintf("Unknown model for annual input: %s", model_name))
}

# ---- 1) Turnover times for ALL models in the summary ----
model_cols <- setdiff(names(results$summary), "year")

turnover_df <- data.frame(
  model = model_cols,
  mean_soc = NA_real_,
  annual_input = NA_real_,
  turnover_years = NA_real_
)

for (m in model_cols) {
  avg_soc <- mean(results$summary[[m]], na.rm = TRUE)
  ann_in <- get_annual_input(m, input_data, site_data)
  
  turnover_df[turnover_df$model == m, "mean_soc"] <- avg_soc
  turnover_df[turnover_df$model == m, "annual_input"] <- ann_in
  turnover_df[turnover_df$model == m, "turnover_years"] <- avg_soc / ann_in
}

cat("\nTurnover times (mean SOC / mean annual input):\n")
print(turnover_df)

# ---- 2) Divergence across models per year ----
soc_divergence <- apply(results$summary[, model_cols, drop = FALSE], 1, function(x) {
  (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)) / mean(x, na.rm = TRUE) * 100
})
cat(sprintf("\nModel divergence across years (range/mean):\n  Mean: %.1f%%\n  Range: %.1f%% - %.1f%%\n",
            mean(soc_divergence, na.rm = TRUE),
            min(soc_divergence, na.rm = TRUE),
            max(soc_divergence, na.rm = TRUE)))

# =============================================================================
# FIGURE 1: SOC time series comparison (all models)
# =============================================================================

png("./Figures/syntetic_data_soc_model_comparison_C_timeseries.png", width = 800, height = 600)

years <- results$summary$year

ylim_all <- range(results$summary[, model_cols, drop = FALSE], na.rm = TRUE)

plot(years, results$summary[[model_cols[1]]],
     type = "l", lwd = 2,
     xlab = "Year", ylab = "SOC (Mg C/ha)",
     ylim = ylim_all,
     main = "SOC model comparison (annual)")

if (length(model_cols) > 1) {
  for (i in 2:length(model_cols)) {
    lines(years, results$summary[[model_cols[i]]], lwd = 2, lty = i)
  }
}

legend("topleft",
       legend = model_cols,
       lty = seq_along(model_cols),
       lwd = 2,
       bty = "n")

dev.off()


# =============================================================================
# SET OF FIGUREs 2: Pairwise scatter vs a reference model (all models vs ref)
# =============================================================================

model_cols <- setdiff(names(results$summary), "year")

# Choose which models to cycle as references
ref_models <- c("yasso07", "yasso15", "yasso20", "rothc")
ref_models <- ref_models[ref_models %in% model_cols]

for (ref_model in ref_models) {
  
  # One PNG per reference model
  fn <- sprintf("./Figures/synthetic_data_soc_scatter_ref_%s.png", ref_model)
  png(fn, width = 1200, height = 800)
  
  others <- setdiff(model_cols, ref_model)
  x <- results$summary[[ref_model]]
  
  op <- par(mfrow = c(2, 2), mar = c(4, 4, 2.5, 1))
  
  for (m in head(others, 4)) {
    
    y <- results$summary[[m]]
    ok <- is.finite(x) & is.finite(y)
    
    plot(x[ok], y[ok],
         xlab = paste(ref_model, "SOC"),
         ylab = paste(m, "SOC"),
         main = paste(m, "vs", ref_model))
    
    abline(0, 1, lty = 2)
    
    rmse <- sqrt(mean((y[ok] - x[ok])^2))
    bias <- mean(y[ok] - x[ok])
    corv <- suppressWarnings(cor(x[ok], y[ok]))
    
    mtext(sprintf("r = %.2f   RMSE = %.2f   bias = %.2f",
                  corv, rmse, bias),
          side = 3, line = -1.2, cex = 0.8)
  }
  
  mtext(paste("SOC comparison – reference:", ref_model),
        outer = TRUE, line = -1, cex = 1.2)
  
  par(op)
  dev.off()
  
  cat("Saved:", fn, "\n")
}

# =============================================================================
# WRITE THE TEMPLATE WITH THE INPUT DATA
# =============================================================================
write.csv(input_data, "synthetic_input_data_template.csv", row.names = FALSE)
write.csv(site_data, "synthetic_site_data_template.csv", row.names = FALSE)





# =============================================================================
# TESTING SUITE (testthat) — HEAVILY COMMENTED
# =============================================================================
#
# Goal: turn your SOC comparison framework into something you can trust.
# These tests are designed to catch the most common sources of
# "models disagree a lot":
#
# 1) Input mapping errors (wrong fractions / missing cohorts / wrong units)
# 2) Non-finite or negative values introduced by mapping or model step
# 3) Spinup not actually reaching equilibrium (so "turnover" is meaningless)
# 4) Climate aggregation bugs (monthly vs annual, nonlinear modifiers)
# 5) Sanity physics: more input => more SOC; warmer => less SOC (when inputs fixed)
# 6) Fraction sensitivity: labile-only should store less than recalcitrant-only
#
# Note: tests here are intentionally "framework-level".
# They do NOT assume that different models should match exactly.
# They assume:
#   - bookkeeping is consistent (mass, sign, finiteness)
#   - qualitative behavior is consistent under controlled scenarios
#   - spinup is long enough to compare equilibrium-like states
#
# =============================================================================

if (!requireNamespace("testthat", quietly = TRUE)) {
  install.packages("testthat")
}
library(testthat)

# -----------------------------------------------------------------------------
# HELPER 1: "all finite and non-negative"
# -----------------------------------------------------------------------------
# Why:
#   In SOC models, negative pools or NaNs almost always indicate a bug:
#     - bad mapping (subtraction, NA propagation)
#     - numerical instability
#     - unit mismatch pushing rates insane
#     - climate modifier producing negative rates
#
# When to relax:
#   Some intermediate fluxes may be negative depending on conventions,
#   but stocks/pools and inputs should never be negative in this framework.
expect_all_finite_nonneg <- function(x, name = deparse(substitute(x))) {
  expect_true(all(is.finite(x)), info = paste("Non-finite values in", name))
  expect_true(all(x >= 0), info = paste("Negative values in", name))
}

# -----------------------------------------------------------------------------
# HELPER 2: Annual total litter input from the synthetic template
# -----------------------------------------------------------------------------
# What it means:
#   Your template has many columns like C_foliage_A, C_foliage_W, ...
#   which are monthly litter inputs. Summing over all these columns gives
#   total monthly litter input (Mg C/ha/month in your synthetic generation).
#
# Then:
#   - We sum per year (12 months) => annual litter (Mg C/ha/yr)
#   - Then average across years => mean annual input.
#
# Why:
#   This is the "ground truth" litter input from the template that mapping
#   functions should preserve (unless mapping intentionally discards/creates C).
annual_total_input_from_template <- function(input_df) {
  litter_cols <- grep("^C_", names(input_df), value = TRUE)
  annual <- input_df |>
    dplyr::group_by(year) |>
    dplyr::summarise(
      total = sum(dplyr::across(dplyr::all_of(litter_cols)), na.rm = TRUE),
      .groups = "drop"
    )
  mean(annual$total, na.rm = TRUE)
}

# -----------------------------------------------------------------------------
# HELPER 3: Make climate constant
# -----------------------------------------------------------------------------
# Why:
#   Many SOC models apply nonlinear climate modifiers (exp, Q10, moisture curves).
#   With seasonal cycles + noise, you can get surprises from aggregation.
#
# A constant-climate scenario is a "unit test world":
#   - removes Jensen’s inequality issues (mean f(T) != f(mean T))
#   - makes spinup behave more predictably
#   - makes cross-model comparisons more interpretable
make_constant_climate <- function(input_df,
                                  temp = 5,
                                  precip = 50,
                                  evap = 20) {
  input_df$temp_air <- temp
  input_df$precip   <- precip
  input_df$evap     <- evap
  input_df
}

# -----------------------------------------------------------------------------
# HELPER 4: Scale ALL litter inputs
# -----------------------------------------------------------------------------
# Why:
#   Monotonicity test: if you double the litter input, equilibrium SOC
#   should increase (all else equal). If not, something is wrong:
#     - mapping drops or clips litter
#     - model treats litter as something else
#     - spinup too short and you're seeing transient weirdness
scale_litter <- function(input_df, factor = 1) {
  litter_cols <- grep("^C_", names(input_df), value = TRUE)
  input_df[litter_cols] <- lapply(input_df[litter_cols], function(z) z * factor)
  input_df
}

# -----------------------------------------------------------------------------
# HELPER 5: Labile-only scenario
# -----------------------------------------------------------------------------
# What:
#   For each cohort, move all AWEN litter into A (acid hydrolysable),
#   set W/E/N to zero.
#
# Why:
#   A fraction is typically fastest decomposing => should reduce SOC storage.
#   This test checks that:
#     - your AWEN columns are actually being used in mapping
#     - fraction allocation affects results in expected direction
#
# Caveat:
#   Some models (e.g. RothC) don't use AWEN directly; their mapping may
#   collapse fractions. So this is a "soft" expectation.
make_labile_only <- function(input_df) {
  cohorts <- unique(gsub("^C_(.*)_[AWEN]$", "\\1",
                         grep("^C_.*_[AWEN]$", names(input_df), value = TRUE)))
  for (coh in cohorts) {
    A <- paste0("C_", coh, "_A")
    W <- paste0("C_", coh, "_W")
    E <- paste0("C_", coh, "_E")
    N <- paste0("C_", coh, "_N")
    if (all(c(A,W,E,N) %in% names(input_df))) {
      tot <- input_df[[A]] + input_df[[W]] + input_df[[E]] + input_df[[N]]
      input_df[[A]] <- tot
      input_df[[W]] <- 0
      input_df[[E]] <- 0
      input_df[[N]] <- 0
    }
  }
  input_df
}

# -----------------------------------------------------------------------------
# HELPER 6: Recalcitrant-only scenario
# -----------------------------------------------------------------------------
# Same idea as labile-only, but push everything into N (non-soluble),
# typically slower decomposing => should increase SOC storage.
make_recalcitrant_only <- function(input_df) {
  cohorts <- unique(gsub("^C_(.*)_[AWEN]$", "\\1",
                         grep("^C_.*_[AWEN]$", names(input_df), value = TRUE)))
  for (coh in cohorts) {
    A <- paste0("C_", coh, "_A")
    W <- paste0("C_", coh, "_W")
    E <- paste0("C_", coh, "_E")
    N <- paste0("C_", coh, "_N")
    if (all(c(A,W,E,N) %in% names(input_df))) {
      tot <- input_df[[A]] + input_df[[W]] + input_df[[E]] + input_df[[N]]
      input_df[[N]] <- tot
      input_df[[A]] <- 0
      input_df[[W]] <- 0
      input_df[[E]] <- 0
    }
  }
  input_df
}

# -----------------------------------------------------------------------------
# HELPER 7: Tail-mean SOC
# -----------------------------------------------------------------------------
# Why:
#   "Mean SOC over entire run" can be dominated by transient approach to equilibrium.
#   For equilibrium comparisons, use the tail period (e.g. last 10 years).
tail_mean_soc <- function(results, model, n_years = 10) {
  df <- results$summary
  yrs <- sort(unique(df$year))
  last_years <- tail(yrs, n_years)
  mean(df[df$year %in% last_years, model], na.rm = TRUE)
}

# -----------------------------------------------------------------------------
# HELPER 8: Tail slope of SOC (quasi steady state check)
# -----------------------------------------------------------------------------
# Why:
#   If SOC is still trending strongly at the end, you are NOT at equilibrium.
#   Then turnover time = SOC/input is not interpretable and models won't
#   be comparable.
#
# We fit SOC ~ year on the last n_years and extract the slope.
# Expected: near 0 (tolerance depends on model + units).
tail_soc_slope <- function(results, model, n_years = 20) {
  df <- results$summary
  yrs <- sort(unique(df$year))
  last_years <- tail(yrs, n_years)
  sub <- df[df$year %in% last_years, c("year", model)]
  sub <- sub[is.finite(sub[[model]]), ]
  if (nrow(sub) < 5) return(NA_real_)
  coef(lm(sub[[model]] ~ sub$year))[2]
}

# -----------------------------------------------------------------------------
# HELPER 9: Safe model runner wrapper
# -----------------------------------------------------------------------------
# Why:
#   Keeps tests short, and if run_soc_models throws, testthat will show context.
run_models_safe <- function(input_df, site_row,
                            models = c("yasso07","yasso15","yasso20","rothc","q_model"),
                            spinup_years = 200) {
  run_soc_models(input_df = input_df,
                 site_row  = site_row,
                 models    = models,
                 spinup_years = spinup_years)
}

# =============================================================================
# TEST 1: synthetic template structure + basic climate sanity
# =============================================================================
test_that("Synthetic input template has required columns and sane climate", {
  
  # Structural: the minimum set required by your wrappers/mappers
  expect_true(all(c("plot_id","year","month","temp_air","precip","evap") %in% names(input_data)))
  
  # Time grid: should be years * 12 for one plot (your generation is single plot here)
  expect_equal(nrow(input_data), length(unique(input_data$year)) * 12)
  
  # Precip should be non-negative and finite
  expect_all_finite_nonneg(input_data$precip, "precip")
  
  # Evap should be finite; you clip to >=5 so it should be positive
  expect_true(all(is.finite(input_data$evap)))
  expect_true(all(input_data$evap > 0))
})

# =============================================================================
# TEST 2: litter columns exist, are finite, and produce a positive annual input
# =============================================================================
test_that("Raw litter totals are stable, finite and non-negative", {
  
  litter_cols <- grep("^C_", names(input_data), value = TRUE)
  expect_true(length(litter_cols) > 0)
  
  # Each litter column should be a non-negative monthly input
  for (cc in litter_cols) {
    expect_all_finite_nonneg(input_data[[cc]], cc)
  }
  
  # Annual total litter input should be positive
  ann_in <- annual_total_input_from_template(input_data)
  expect_true(is.finite(ann_in))
  expect_true(ann_in > 0)
})

# =============================================================================
# TEST 3: mapping functions preserve total litter mass (within tolerance)
# =============================================================================
test_that("Mapping functions preserve total litter input (within tolerance)", {
  
  # Ground truth from the template
  ann_template <- annual_total_input_from_template(input_data)
  
  # ---- Yasso07 / Yasso15 mapping ----
  # Expectation: mapping to nwl/fwl/cwl partitions the litter,
  # but should conserve total C (aside from numerical rounding).
  if (exists("map_input_to_yasso07")) {
    mapped <- map_input_to_yasso07(input_data)
    
    # Must have these columns
    expect_true(all(c("C_nwl","C_fwl","C_cwl") %in% names(mapped)))
    
    # No negative or non-finite mapped litter
    expect_all_finite_nonneg(mapped$C_nwl, "mapped$C_nwl")
    expect_all_finite_nonneg(mapped$C_fwl, "mapped$C_fwl")
    expect_all_finite_nonneg(mapped$C_cwl, "mapped$C_cwl")
    
    # Compare mean annual input
    ann_mapped <- mean(mapped$C_nwl + mapped$C_fwl + mapped$C_cwl, na.rm = TRUE)
    
    # 2% tolerance: should be tight because this is simple partitioning.
    expect_true(abs(ann_mapped - ann_template) / ann_template < 0.02,
                info = sprintf(
                  "Yasso07 mapping total differs: mapped=%.4f template=%.4f",
                  ann_mapped, ann_template
                ))
  }
  
  # ---- Yasso20 mapping ----
  # Expectation: mapping returns monthly AWEN for nwl/fwl/cwl, still mass-conserving.
  if (exists("map_input_to_yasso20")) {
    mapped <- map_input_to_yasso20(input_data)
    
    needed <- c("nwl_A","nwl_W","nwl_E","nwl_N",
                "fwl_A","fwl_W","fwl_E","fwl_N",
                "cwl_A","cwl_W","cwl_E","cwl_N")
    expect_true(all(needed %in% names(mapped)))
    
    total <- mapped$nwl_A + mapped$nwl_W + mapped$nwl_E + mapped$nwl_N +
      mapped$fwl_A + mapped$fwl_W + mapped$fwl_E + mapped$fwl_N +
      mapped$cwl_A + mapped$cwl_W + mapped$cwl_E + mapped$cwl_N
    
    expect_all_finite_nonneg(total, "yasso20 mapped total litter")
    
    ann_mapped <- mean(total, na.rm = TRUE)
    
    # 2% tolerance: should be mostly partitioning/relabeling + rounding.
    expect_true(abs(ann_mapped - ann_template) / ann_template < 0.02,
                info = sprintf(
                  "Yasso20 mapping total differs: mapped=%.4f template=%.4f",
                  ann_mapped, ann_template
                ))
  }
  
  # ---- RothC mapping ----
  # Expectation: RothC input mapping splits litter into DPM and RPM (and maybe
  # routes some to HUM depending on design). In your get_annual_input() you
  # use DPM+RPM, so here we check that DPM+RPM captures ~all litter.
  if (exists("map_input_to_rothc")) {
    mapped <- map_input_to_rothc(input_data, site_data)
    expect_true(all(c("year","C_DPM","C_RPM") %in% names(mapped)))
    
    annual <- aggregate(mapped$C_DPM + mapped$C_RPM, by = list(year = mapped$year), FUN = sum)
    ann_mapped <- mean(annual$x, na.rm = TRUE)
    
    # 5% tolerance: allow a bit more leeway because some RothC mappings
    # may allocate parts differently or do minor transformations.
    expect_true(abs(ann_mapped - ann_template) / ann_template < 0.05,
                info = sprintf(
                  "RothC mapping total differs: mapped=%.4f template=%.4f",
                  ann_mapped, ann_template
                ))
  }
})

# =============================================================================
# TEST 4: model outputs are finite and non-negative
# =============================================================================
test_that("Model runs produce finite, non-negative annual SOC series", {
  
  res <- run_models_safe(input_data, site_data, spinup_years = 200)
  
  expect_true(is.list(res))
  expect_true("summary" %in% names(res))
  
  model_cols <- setdiff(names(res$summary), "year")
  expect_true(length(model_cols) >= 2)
  
  # SOC stock should never be negative or non-finite
  for (m in model_cols) {
    expect_true(all(is.finite(res$summary[[m]])), info = paste("Non-finite SOC in", m))
    expect_true(all(res$summary[[m]] >= 0), info = paste("Negative SOC in", m))
  }
})

# =============================================================================
# TEST 5: quasi-steady state under constant climate (spinup check)
# =============================================================================
test_that("Spinup achieves quasi-steady state (small slope) in constant climate", {
  
  # Constant climate reduces seasonal nonlinearities and should allow
  # clean equilibration.
  input_const <- make_constant_climate(input_data, temp = 5, precip = 50, evap = 20)
  
  # Use longer spinup here; this test is about equilibration itself.
  res <- run_models_safe(input_const, site_data, spinup_years = 500)
  
  model_cols <- setdiff(names(res$summary), "year")
  slopes <- sapply(model_cols, function(m) tail_soc_slope(res, m, n_years = 20))
  
  # Threshold:
  #   0.5 Mg C/ha/year is a loose "still drifting" cutoff.
  # Tune after seeing typical slopes in your system.
  expect_true(all(is.na(slopes) | abs(slopes) < 0.5),
              info = paste("Large SOC slopes:", paste(names(slopes), round(slopes, 3), collapse = "; ")))
})

# =============================================================================
# TEST 6: monotonicity — doubling litter increases equilibrium SOC
# =============================================================================
test_that("Monotonicity: doubling litter increases equilibrium SOC (constant climate)", {
  
  input_const <- make_constant_climate(input_data, temp = 5, precip = 50, evap = 20)
  
  res1 <- run_models_safe(input_const, site_data, spinup_years = 500)
  res2 <- run_models_safe(scale_litter(input_const, 2), site_data, spinup_years = 500)
  
  model_cols <- intersect(setdiff(names(res1$summary), "year"),
                          setdiff(names(res2$summary), "year"))
  
  for (m in model_cols) {
    soc1 <- tail_mean_soc(res1, m, n_years = 10)
    soc2 <- tail_mean_soc(res2, m, n_years = 10)
    
    # If this fails:
    #   - mapping is broken (inputs not actually doubled)
    #   - model has internal feedbacks altering input (unlikely here)
    #   - spinup not long enough / slope not near zero
    expect_true(soc2 > soc1,
                info = sprintf("%s: doubled input did not increase SOC (%.2f -> %.2f)", m, soc1, soc2))
  }
})

# =============================================================================
# TEST 7: fraction sensitivity — labile-only stores less than recalcitrant-only
# =============================================================================
test_that("Labile-only yields lower SOC than recalcitrant-only (constant climate)", {
  
  input_const <- make_constant_climate(input_data, temp = 5, precip = 50, evap = 20)
  
  res_labile <- run_models_safe(make_labile_only(input_const), site_data, spinup_years = 500)
  res_recalc <- run_models_safe(make_recalcitrant_only(input_const), site_data, spinup_years = 500)
  
  model_cols <- intersect(setdiff(names(res_labile$summary), "year"),
                          setdiff(names(res_recalc$summary), "year"))
  
  for (m in model_cols) {
    soc_l <- tail_mean_soc(res_labile, m, n_years = 10)
    soc_r <- tail_mean_soc(res_recalc, m, n_years = 10)
    
    # Soft expectation:
    #   Under a simple decomposition logic, recalcitrant inputs increase storage.
    # If this fails:
    #   - your mapping collapses AWEN (so fractions don't matter)
    #   - the model has different conceptual pools
    #   - or there is a sign/unit bug in fraction routing
    expect_true(soc_l <= soc_r,
                info = sprintf("%s: labile-only SOC not <= recalcitrant-only SOC (%.2f vs %.2f)", m, soc_l, soc_r))
  }
})


# =============================================================================
# TEST 7 (REVISED): AWEN manipulation has a noticeable effect (fractions used)
# =============================================================================
test_that("AWEN manipulation changes equilibrium SOC (fractions are used)", {
  
  input_const <- make_constant_climate(input_data, temp = 5, precip = 50, evap = 20)
  
  res_labile <- run_models_safe(make_labile_only(input_const), site_data,
                                models = c("yasso07","yasso15","yasso20"), spinup_years = 500)
  res_recalc <- run_models_safe(make_recalcitrant_only(input_const), site_data,
                                models = c("yasso07","yasso15","yasso20"), spinup_years = 500)
  
  model_cols <- intersect(setdiff(names(res_labile$summary), "year"),
                          setdiff(names(res_recalc$summary), "year"))
  
  for (m in model_cols) {
    soc_l <- tail_mean_soc(res_labile, m, n_years = 10)
    soc_r <- tail_mean_soc(res_recalc, m, n_years = 10)
    
    # This test checks that changing fractions actually changes outcomes.
    # Direction depends on parameterization, so we only require "non-trivial change".
    rel_change <- abs(soc_l - soc_r) / mean(c(soc_l, soc_r))
    expect_true(rel_change > 0.05,
                info = sprintf("%s: SOC barely changed when AWEN was altered (labile=%.2f recalc=%.2f)",
                               m, soc_l, soc_r))
  }
})



# =============================================================================
# TEST 8: temperature sensitivity sanity — warming lowers SOC (inputs fixed)
# =============================================================================
test_that("Temperature sanity: warmer should not increase SOC when litter is fixed", {
  
  # We isolate temperature only.
  # Precip/evap kept constant to avoid moisture function coupling effects.
  input_cold <- make_constant_climate(input_data, temp = 0, precip = 50, evap = 20)
  input_warm <- make_constant_climate(input_data, temp = 8, precip = 50, evap = 20)
  
  res_cold <- run_models_safe(input_cold, site_data, spinup_years = 500)
  res_warm <- run_models_safe(input_warm, site_data, spinup_years = 500)
  
  model_cols <- intersect(setdiff(names(res_cold$summary), "year"),
                          setdiff(names(res_warm$summary), "year"))
  
  for (m in model_cols) {
    soc_c <- tail_mean_soc(res_cold, m, n_years = 10)
    soc_w <- tail_mean_soc(res_warm, m, n_years = 10)
    
    # Expectation:
    #   With fixed litter inputs, warmer => faster decomposition => lower SOC.
    #
    # If this fails:
    #   - temperature units bug (°C vs K)
    #   - modifier inverted (increasing with lower temp)
    #   - moisture modifier tied to evap in a way you didn't isolate
    expect_true(soc_w <= soc_c,
                info = sprintf("%s: warm SOC not <= cold SOC (cold=%.2f warm=%.2f)", m, soc_c, soc_w))
  }
})

# =============================================================================
# TEST 9: annual input computation returns positive finite values
# =============================================================================
test_that("Annual-input calculator returns positive finite values per model", {
  
  model_cols <- c("yasso07","yasso15","yasso20","rothc","q_model")
  
  for (m in model_cols) {
    ann_in <- get_annual_input(m, input_data, site_data)
    
    # If this fails:
    #   - mapping returns NA/NaN due to missing columns
    #   - get_annual_input uses monthly mean when it should annual sum (or vice versa)
    expect_true(is.finite(ann_in), info = paste("Annual input not finite for", m))
    expect_true(ann_in > 0, info = paste("Annual input not positive for", m))
  }
})

cat("\n✅ Finished running testthat suite.\n")




# =============================================================================
# FOCUS ON YASSO20 DISCREPANCIES
# =============================================================================

# Runs:
#   S1) Original forcing (your current synthetic dataset)
#   S2) No temperature seasonality (monthly temps flattened to annual mean per year)
#   S3) Fully constant climate (temp, precip, evap fixed constants)
#
# For each scenario:
#   - run yasso07, yasso15, yasso20
#   - compute tail-mean SOC over last N years
#   - compute ratios y20/y07 and y15/y07
#
# Paste this AFTER you have input_data, site_data, and multi_model.R sourced.

library(dplyr)

# -----------------------------
# Helpers
# -----------------------------

# Flatten temp seasonality but preserve annual mean per year
# (monthly temp_air replaced by the within-year mean)
remove_temp_seasonality <- function(input_df) {
  input_df %>%
    group_by(plot_id, year) %>%
    mutate(temp_air = mean(temp_air, na.rm = TRUE)) %>%
    ungroup()
}

# Fully constant forcing (maximally comparable)
make_fully_constant_forcing <- function(input_df, temp = 5, precip = 50, evap = 20) {
  input_df$temp_air <- temp
  input_df$precip   <- precip
  input_df$evap     <- evap
  input_df
}

# Tail-mean SOC from your run_soc_models output
tail_mean_soc <- function(results, model, n_years = 10) {
  df <- results$summary
  yrs <- sort(unique(df$year))
  last_years <- tail(yrs, n_years)
  mean(df[df$year %in% last_years, model], na.rm = TRUE)
}

# Convenience: run only Yasso models
run_yasso <- function(input_df, site_data, spinup_years = 1000) {
  run_soc_models(
    input_df = input_df,
    site_row = site_data,
    models = c("yasso07", "yasso15", "yasso20"),
    spinup_years = spinup_years
  )
}

# Summarise a run into a single row
summarise_yasso_run <- function(res, scenario, n_tail_years = 10) {
  y07 <- tail_mean_soc(res, "yasso07", n_years = n_tail_years)
  y15 <- tail_mean_soc(res, "yasso15", n_years = n_tail_years)
  y20 <- tail_mean_soc(res, "yasso20", n_years = n_tail_years)
  
  data.frame(
    scenario = scenario,
    yasso07 = y07,
    yasso15 = y15,
    yasso20 = y20,
    ratio_y15_y07 = y15 / y07,
    ratio_y20_y07 = y20 / y07,
    stringsAsFactors = FALSE
  )
}

# -----------------------------
# Build scenarios
# -----------------------------

inputs <- list(
  S1_original = input_data,
  S2_no_temp_seasonality = remove_temp_seasonality(input_data),
  S3_fully_constant = make_fully_constant_forcing(input_data, temp = 5, precip = 50, evap = 20)
)

# -----------------------------
# Run scenarios + collect results
# -----------------------------
# Notes:
# - spinup_years is set to 1000 to reduce "still drifting" effects.
#   If this is slow, drop to 500, but then check tail slopes as you did earlier.

out <- lapply(names(inputs), function(scn) {
  cat("\n--- Running", scn, "---\n")
  res <- run_yasso(inputs[[scn]], site_data, spinup_years = 1000)
  summarise_yasso_run(res, scenario = scn, n_tail_years = 10)
})

yasso_scenario_summary <- bind_rows(out)

print(yasso_scenario_summary)

# -----------------------------
# Quick visualization (base R)
# -----------------------------
#png("./Figures/yasso_scenarios_tailmean_soc.png", width = 900, height = 600)
mat <- as.matrix(yasso_scenario_summary[, c("yasso07","yasso15","yasso20")])
rownames(mat) <- yasso_scenario_summary$scenario
barplot(t(mat),
        beside = TRUE,
        legend.text = colnames(mat),
        args.legend = list(bty = "n"),
        ylab = "Tail-mean SOC (Mg C/ha)",
        main = "Yasso tail-mean SOC across forcing scenarios")
#dev.off()

cat("\nSaved figure: ./Figures/yasso_scenarios_tailmean_soc.png\n")
