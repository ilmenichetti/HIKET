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



# =============================================================================
# UPDATED INPUT TEMPLATE (4 plots, 1990-2020, SOC observations)
# =============================================================================

generate_input_template <- function(plot_ids = c("PLOT_001", "PLOT_002", "PLOT_003", "PLOT_004"),
                                    years = 1990:2020,
                                    soc_meas_years = c(1997, 2006),
                                    soc_meas_month = 7) {
  
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
  
  # Plot-to-plot multiplier (deterministic, but slightly different per plot)
  plot_mult <- setNames(seq(0.95, 1.15, length.out = length(plot_ids)), plot_ids)
  
  # Generate litter columns (monthly = annual/12)
  for (coh in cohorts) {
    total <- litter_totals[[coh]] / 12
    props <- awen_props[[coh]]
    
    for (frac in c("A", "W", "E", "N")) {
      col_name <- paste0("C_", coh, "_", frac)
      base_value <- total * props[frac]
      
      # Variation by plot and year
      plot_effect <- plot_mult[grid$plot_id]
      year_effect <- 1 + 0.1 * sin(2 * pi * (grid$year - min(years)) / 10)
      
      grid[[col_name]] <- round(base_value * plot_effect * year_effect, 6)
    }
  }
  
  # Climate (slightly different baseline per plot)
  base_temp_by_plot <- setNames(c(3.5, 2.5, 2.0, 1.5)[seq_along(plot_ids)], plot_ids)
  grid$temp_air <- round(
    base_temp_by_plot[grid$plot_id] +
      12 * sin(2 * pi * (grid$month - 4) / 12) +
      rnorm(n, 0, 1),
    1
  )
  
  monthly_precip <- c(40, 35, 35, 40, 45, 55, 70, 75, 60, 55, 50, 45)
  precip_mult <- setNames(c(1.00, 1.05, 1.10, 0.95)[seq_along(plot_ids)], plot_ids)
  
  grid$precip <- round(
    monthly_precip[grid$month] *
      precip_mult[grid$plot_id] *
      (1 + rnorm(n, 0, 0.2)),
    1
  )
  
  grid$evap <- round(
    pmax(5, 20 + 40 * sin(2 * pi * (grid$month - 4) / 12) + rnorm(n, 0, 5)),
    1
  )
  
  # ---------------------------------------------------------------------------
  # Measured total SOC stock (t C / ha): only one point in 1997 and one in 2006
  # ---------------------------------------------------------------------------
  grid$soc_obs_tCha <- NA_real_
  
  # Reasonable Finnish forest mineral soil SOC stock ballpark (0-30 cm-ish):
  # ~60–140 t C/ha (wide, but plausible). We'll generate per-plot baselines
  # and add small measurement noise + a modest change between 1997 and 2006.
  set.seed(42)  # remove or change if you want different synthetic SOC each run
  
  soc_base <- setNames(runif(length(plot_ids), min = 70, max = 130), plot_ids)
  
  for (pid in plot_ids) {
    for (yy in soc_meas_years) {
      idx <- which(grid$plot_id == pid & grid$year == yy & grid$month == soc_meas_month)
      if (length(idx) == 1) {
        # modest trend over time + measurement noise
        # (could be slightly increasing or decreasing)
        trend_per_year <- rnorm(1, mean = 0.15, sd = 0.25)  # t C/ha/year (small)
        expected <- soc_base[pid] + trend_per_year * (yy - min(soc_meas_years))
        grid$soc_obs_tCha[idx] <- round(expected + rnorm(1, 0, 3), 1)  # ±3 t/ha noise
      }
    }
  }
  
  grid
}

# =============================================================================
# UPDATED SITE TEMPLATE (4 plots)
# =============================================================================
generate_site_template <- function(plot_ids = c("PLOT_001", "PLOT_002", "PLOT_003", "PLOT_004")) {
  
  # simple, plausible variation across plots
  n <- length(plot_ids)
  df <- data.frame(
    plot_id = plot_ids,
    pine_prop   = c(0.30, 0.20, 0.40, 0.25)[seq_len(n)],
    spruce_prop = c(0.60, 0.70, 0.50, 0.65)[seq_len(n)],
    birch_prop  = c(0.10, 0.10, 0.10, 0.10)[seq_len(n)],
    clay        = c(18, 12, 25, 15)[seq_len(n)],
    soil_depth  = c(23, 18, 30, 20)[seq_len(n)],
    stringsAsFactors = FALSE
  )
  
  # ensure proportions sum to 1 (in case you tweak above later)
  s <- df$pine_prop + df$spruce_prop + df$birch_prop
  df$pine_prop <- df$pine_prop / s
  df$spruce_prop <- df$spruce_prop / s
  df$birch_prop <- df$birch_prop / s
  
  df
}

# =============================================================================
# USE IT
# =============================================================================
input_data <- generate_input_template(
  plot_ids = c("PLOT_001", "PLOT_002", "PLOT_003", "PLOT_004"),
  years = 1990:2020,
  soc_meas_years = c(1997, 2006),
  soc_meas_month = 7
)

site_data <- generate_site_template(
  plot_ids = c("PLOT_001", "PLOT_002", "PLOT_003", "PLOT_004")
)

cat(sprintf("Generated input data: %d rows (months)\n", nrow(input_data)))
cat(sprintf("Time span: %d-%d\n", min(input_data$year), max(input_data$year)))
cat("SOC observations (non-NA):\n")
print(input_data %>%
        dplyr::filter(!is.na(soc_obs_tCha)) %>%
        dplyr::select(plot_id, year, month, soc_obs_tCha) %>%
        dplyr::arrange(plot_id, year))


# =============================================================================
# TEST THE WRAPPER FOR CALIBRATION
# =============================================================================

# First, source the new wrapper
source("./Models_functions/multi_model_calibration_ready.R")

# Test 1: Run with DEFAULT parameters (should match your original fast wrapper)
timed_default <- system.time(
  invisible(capture.output(
    results_default <- run_soc_models_multipplot_calibration(
      input_df = input_data,
      site_df = site_data,
      models = c("yasso07", "yasso15", "yasso20", "rothc", "q_model"),
      params_list = NULL,  # NULL means use defaults
      spinup_years = 100
    )
  ))
)

cat("Default parameters timing:\n")
print(timed_default)
cat("\nFirst few rows of summary:\n")
print(head(results_default$summary))


# Test 2: Run with CUSTOM Q-model parameters (to verify parameter passing works)
custom_q_params <- Q_MODEL_DEFAULT_PARAMS
custom_q_params$u01 <- 0.020  # Increase temperature sensitivity (default is 0.0157)

timed_custom <- system.time(
  invisible(capture.output(
    results_custom <- run_soc_models_multipplot_calibration(
      input_df = input_data,
      site_df = site_data,
      models = "q_model",  # Only run Q-model for speed
      params_list = list(q_model = custom_q_params),  # Custom params
      spinup_years = 100
    )
  ))
)

cat("\n\nCustom parameters timing:\n")
print(timed_custom)


# Test 3: Compare results to verify parameter passing is working
cat("\n\n=== VERIFICATION: Parameter passing is working ===\n")
cat("Q-model SOC in year 2020 (PLOT_001):\n")
cat(sprintf("  Default (u01=0.0157): %.2f t C/ha\n", 
            results_default$summary$q_model[results_default$summary$plot_id == "PLOT_001" & 
                                              results_default$summary$year == 2020]))
cat(sprintf("  Custom  (u01=0.020):  %.2f t C/ha\n", 
            results_custom$summary$q_model[results_custom$summary$plot_id == "PLOT_001" & 
                                             results_custom$summary$year == 2020]))
cat("\nIf these numbers are DIFFERENT, parameter passing works! ✓\n")


# Test 4: Verify fix for years extraction with edge case
cat("\n\n=== VERIFICATION: Years extraction fix ===\n")
tryCatch({
  # This would fail in old version if running only monthly model
  results_edge <- run_soc_models_oneplot_calibration(
    input_df = input_data[input_data$plot_id == "PLOT_001", ],
    site_row = site_data[site_data$plot_id == "PLOT_001", ],
    models = "rothc",  # Only monthly model
    spinup_years = 100
  )
  cat("Years extraction: SUCCESS ✓\n")
  cat("Years found:", paste(head(results_edge$year), collapse=", "), "...\n")
}, error = function(e) {
  cat("Years extraction: FAILED ✗\n")
  cat("Error:", e$message, "\n")
})






# =============================================================================
# COMPREHENSIVE PARAMETER PASSING TEST - ALL MODELS
# =============================================================================
# This tests that custom parameters actually affect model outputs for each model
# =============================================================================

# Assumes you've already sourced:
# - wrapper_calibration_ready.R
# - All model functions (Yasso07, Yasso15, Yasso20, RothC, Q-model)
# - Input_matching_functions.R
# And have input_data and site_data loaded

# Use just one plot for speed
test_input <- input_data[input_data$plot_id == "PLOT_001", ]
test_site <- site_data[site_data$plot_id == "PLOT_001", ]

# =============================================================================
# TEST 1: YASSO07
# =============================================================================
cat("--- TEST 1: YASSO07 ---\n")

# Default parameters
result_y07_default <- run_soc_models_oneplot_calibration(
  input_df = test_input,
  site_row = test_site,
  models = "yasso07",
  params_list = NULL,  # Use defaults
  spinup_years = 100
)

# Custom parameters (modify decomposition rate of A pool)
custom_y07 <- YASSO07_DEFAULT_PARAMS
custom_y07$alpha_A <- custom_y07$alpha_A * 1.5  # Modify existing parameter

result_y07_custom <- run_soc_models_oneplot_calibration(
  input_df = test_input,
  site_row = test_site,
  models = "yasso07",
  params_list = list(yasso07 = custom_y07),
  spinup_years = 100
)

# Compare final year
soc_default <- result_y07_default$yasso07[result_y07_default$year == 2020]
soc_custom <- result_y07_custom$yasso07[result_y07_custom$year == 2020]
difference <- abs(soc_default - soc_custom)

cat(sprintf("  Default params: %.2f t C/ha\n", soc_default))
cat(sprintf("  Custom params:  %.2f t C/ha\n", soc_custom))
cat(sprintf("  Difference:     %.2f t C/ha\n", difference))
cat(sprintf("  Status: %s\n\n", if(difference > 1) "✓ PASS" else "✗ FAIL - No difference detected"))


# =============================================================================
# TEST 2: YASSO15
# =============================================================================
cat("--- TEST 2: YASSO15 ---\n")

result_y15_default <- run_soc_models_oneplot_calibration(
  input_df = test_input,
  site_row = test_site,
  models = "yasso15",
  params_list = NULL,
  spinup_years = 100
)

# Custom parameters (modify temperature sensitivity)
custom_y15 <- YASSO15_DEFAULT_PARAMS
custom_y15$beta1 <- custom_y15$beta1 * 1.3  # Increase temperature sensitivity

result_y15_custom <- run_soc_models_oneplot_calibration(
  input_df = test_input,
  site_row = test_site,
  models = "yasso15",
  params_list = list(yasso15 = custom_y15),
  spinup_years = 100
)

soc_default <- result_y15_default$yasso15[result_y15_default$year == 2020]
soc_custom <- result_y15_custom$yasso15[result_y15_custom$year == 2020]
difference <- abs(soc_default - soc_custom)

cat(sprintf("  Default params: %.2f t C/ha\n", soc_default))
cat(sprintf("  Custom params:  %.2f t C/ha\n", soc_custom))
cat(sprintf("  Difference:     %.2f t C/ha\n", difference))
cat(sprintf("  Status: %s\n\n", if(difference > 1) "✓ PASS" else "✗ FAIL - No difference detected"))


# =============================================================================
# TEST 3: YASSO20
# =============================================================================
cat("--- TEST 3: YASSO20 ---\n")

result_y20_default <- run_soc_models_oneplot_calibration(
  input_df = test_input,
  site_row = test_site,
  models = "yasso20",
  params_list = NULL,
  spinup_years = 100
)

# Custom parameters (modify woody decomposition)
custom_y20 <- YASSO20_DEFAULT_PARAMS
custom_y20$alpha_W <- custom_y20$alpha_W * 3

result_y20_custom <- run_soc_models_oneplot_calibration(
  input_df = test_input,
  site_row = test_site,
  models = "yasso20",
  params_list = list(yasso20 = custom_y20),
  spinup_years = 100
)

soc_default <- result_y20_default$yasso20[result_y20_default$year == 2020]
soc_custom <- result_y20_custom$yasso20[result_y20_custom$year == 2020]
difference <- abs(soc_default - soc_custom)

cat(sprintf("  Default params: %.2f t C/ha\n", soc_default))
cat(sprintf("  Custom params:  %.2f t C/ha\n", soc_custom))
cat(sprintf("  Difference:     %.2f t C/ha\n", difference))
cat(sprintf("  Status: %s\n\n", if(difference > 1) "✓ PASS" else "✗ FAIL - No difference detected"))


# =============================================================================
# TEST 4: ROTHC
# =============================================================================
cat("--- TEST 4: ROTHC ---\n")

result_rothc_default <- run_soc_models_oneplot_calibration(
  input_df = test_input,
  site_row = test_site,
  models = "rothc",
  params_list = NULL,
  spinup_years = 100
)

# Custom parameters (modify DPM decomposition rate)
custom_rothc <- ROTHC_DEFAULT_PARAMS
custom_rothc$k_HUM <- custom_rothc$k_HUM * 3.0

result_rothc_custom <- run_soc_models_oneplot_calibration(
  input_df = test_input,
  site_row = test_site,
  models = "rothc",
  params_list = list(rothc = custom_rothc),
  spinup_years = 100
)

soc_default <- result_rothc_default$rothc[result_rothc_default$year == 2020]
soc_custom <- result_rothc_custom$rothc[result_rothc_custom$year == 2020]
difference <- abs(soc_default - soc_custom)

cat(sprintf("  Default params: %.2f t C/ha\n", soc_default))
cat(sprintf("  Custom params:  %.2f t C/ha\n", soc_custom))
cat(sprintf("  Difference:     %.2f t C/ha\n", difference))
cat(sprintf("  Status: %s\n\n", if(difference > 1) "✓ PASS" else "✗ FAIL - No difference detected"))


# =============================================================================
# TEST 5: Q-MODEL
# =============================================================================
cat("--- TEST 5: Q-MODEL ---\n")

result_q_default <- run_soc_models_oneplot_calibration(
  input_df = test_input,
  site_row = test_site,
  models = "q_model",
  params_list = NULL,
  spinup_years = 100
)

# Custom parameters (modify temperature sensitivity)
custom_q <- Q_MODEL_DEFAULT_PARAMS
custom_q$u01 <- 0.025  # Increase from default 0.0157

result_q_custom <- run_soc_models_oneplot_calibration(
  input_df = test_input,
  site_row = test_site,
  models = "q_model",
  params_list = list(q_model = custom_q),
  spinup_years = 100
)

soc_default <- result_q_default$q_model[result_q_default$year == 2020]
soc_custom <- result_q_custom$q_model[result_q_custom$year == 2020]
difference <- abs(soc_default - soc_custom)

cat(sprintf("  Default params: %.2f t C/ha\n", soc_default))
cat(sprintf("  Custom params:  %.2f t C/ha\n", soc_custom))
cat(sprintf("  Difference:     %.2f t C/ha\n", difference))
cat(sprintf("  Status: %s\n\n", if(difference > 1) "✓ PASS" else "✗ FAIL - No difference detected"))


# =============================================================================
# TEST 6: MULTIPLE MODELS SIMULTANEOUSLY
# =============================================================================
cat("--- TEST 6: MULTIPLE MODELS WITH CUSTOM PARAMS ---\n")

# Custom parameters for multiple models at once
custom_params_all <- list(
  yasso07 = custom_y07,
  yasso15 = custom_y15,
  q_model = custom_q
)

result_multi <- run_soc_models_oneplot_calibration(
  input_df = test_input,
  site_row = test_site,
  models = c("yasso07", "yasso15", "q_model"),
  params_list = custom_params_all,
  spinup_years = 100
)

cat("  Custom parameters for 3 models simultaneously:\n")
cat(sprintf("    Yasso07: %.2f t C/ha\n", result_multi$yasso07[result_multi$year == 2020]))
cat(sprintf("    Yasso15: %.2f t C/ha\n", result_multi$yasso15[result_multi$year == 2020]))
cat(sprintf("    Q-model: %.2f t C/ha\n", result_multi$q_model[result_multi$year == 2020]))
cat("  Status: ✓ PASS (no crashes)\n\n")

