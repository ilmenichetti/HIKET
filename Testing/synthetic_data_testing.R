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


dir.create("./Testing", showWarnings = FALSE)


# =============================================================================
# STEP 1: Generate Synthetic Input Data
# =============================================================================

generate_input_template <- function(plot_ids = c("PLOT_001", "PLOT_002", "PLOT_003", "PLOT_004"),
                                    years = 2000:2020) {
  
  grid <- expand.grid(
    plot_id = plot_ids,
    year    = years,
    month   = 1:12,
    stringsAsFactors = FALSE
  )
  grid <- grid[order(grid$plot_id, grid$year, grid$month), ]
  n <- nrow(grid)
  
  # AWEN proportions by litter class.
  # NWL weights are derived from foliage (1.2) + fine_roots (0.8) + understorey (0.3).
  # CWL weights are derived from stems (0.1) + coarse_roots (0.3).
  # FWL = branches only.
  awen_props <- list(
    nwl = c(A = 0.49, W = 0.09, E = 0.02, N = 0.40),
    fwl = c(A = 0.50, W = 0.05, E = 0.03, N = 0.42),
    cwl = c(A = 0.45, W = 0.03, E = 0.01, N = 0.51)
  )
  
  # Annual litter totals per class (Mg C/ha/yr)
  litter_totals <- list(nwl = 2.3, fwl = 0.4, cwl = 0.4)
  
  # Per-plot base temperature offsets (simulate site variability)
  base_temps <- c(3.5, 2.0, 4.5, 1.5)
  names(base_temps) <- plot_ids
  
  # Per-plot precip scaling
  precip_scales <- c(1.00, 1.10, 0.90, 1.15)
  names(precip_scales) <- plot_ids
  
  # Per-plot litter scaling (simulate different productivity)
  litter_scales <- c(1.00, 1.10, 0.85, 1.20)
  names(litter_scales) <- plot_ids
  
  # Generate litter columns (monthly = annual/12)
  for (cls in names(litter_totals)) {
    total_monthly <- litter_totals[[cls]] / 12
    props <- awen_props[[cls]]
    
    for (frac in c("A", "W", "E", "N")) {
      col_name   <- paste0("C_", cls, "_", frac)
      base_value <- total_monthly * props[frac]
      
      plot_effect <- litter_scales[grid$plot_id]
      year_effect <- 1 + 0.1 * sin(2 * pi * (grid$year - min(years)) / 10)
      
      grid[[col_name]] <- round(base_value * plot_effect * year_effect, 6)
    }
  }
  
  # Climate
  base_temp_vec <- base_temps[grid$plot_id]
  grid$temp_air <- round(base_temp_vec + 12 * sin(2 * pi * (grid$month - 4) / 12) +
                           rnorm(n, 0, 1), 1)
  
  monthly_precip <- c(40, 35, 35, 40, 45, 55, 70, 75, 60, 55, 50, 45)
  grid$precip <- round(monthly_precip[grid$month] *
                         precip_scales[grid$plot_id] *
                         (1 + rnorm(n, 0, 0.2)), 1)
  
  grid$evap <- round(pmax(5, 20 + 40 * sin(2 * pi * (grid$month - 4) / 12) +
                            rnorm(n, 0, 5)), 1)
  
  # SOC observations: 2 measurements per plot (beginning and end of time series),
  # recorded in July (month 7). Values are realistic Finnish forest SOC
  # (30-80 tC/ha), with 5% plot-level noise. All other rows are NA.
  soc_base <- c(55, 48, 62, 70)
  names(soc_base) <- plot_ids
  obs_years <- c(min(years), max(years))
  
  grid$soc_obs_tCha <- NA_real_
  for (pid in plot_ids) {
    for (oy in obs_years) {
      row_idx <- which(grid$plot_id == pid & grid$year == oy & grid$month == 7L)
      if (length(row_idx) == 1L)
        grid$soc_obs_tCha[row_idx] <- round(soc_base[pid] * (1 + rnorm(1, 0, 0.05)), 1)
    }
  }
  
  grid
}

# Generate site data
generate_site_template <- function(plot_ids = c("PLOT_001", "PLOT_002", "PLOT_003", "PLOT_004")) {
  data.frame(
    plot_id     = plot_ids,
    pine_prop   = c(0.30, 0.50, 0.20, 0.70),
    spruce_prop = c(0.60, 0.40, 0.70, 0.20),
    birch_prop  = c(0.10, 0.10, 0.10, 0.10),
    clay        = c(18,   22,   15,   25),
    soil_depth  = c(23,   20,   25,   18),
    stringsAsFactors = FALSE
  )
}


# Generate data
input_data <- generate_input_template(
  plot_ids = c("PLOT_001", "PLOT_002", "PLOT_003", "PLOT_004"),
  years = 2000:2020
)
site_data <- generate_site_template(
  plot_ids = c("PLOT_001", "PLOT_002", "PLOT_003", "PLOT_004")
)

cat(sprintf("Generated input data: %d rows (months)\n", nrow(input_data)))
cat(sprintf("Time span: %d-%d\n", min(input_data$year), max(input_data$year)))
cat(sprintf("SOC observations: %d rows\n", sum(!is.na(input_data$soc_obs_tCha))))


# -----------------------------------------------------
# Climate of the scenario
# -----------------------------------------------------

png("./Testing/synthetic_climate.png", width = 1000, height = 350)
par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))

boxplot(input_data$temp_air ~ input_data$month,
        main = "Monthly temperature (°C)", xlab = "Month", ylab = "Temp (°C)")

boxplot(input_data$precip ~ input_data$month,
        main = "Monthly precipitation (mm)", xlab = "Month", ylab = "Precipitation (mm)")

boxplot(input_data$evap ~ input_data$month,
        main = "Monthly evaporation (mm)", xlab = "Month", ylab = "ET (mm)")

par(mfrow = c(1, 1))
dev.off()


# -----------------------------------------------------
# Monthly litter totals by class (A + W + E + N)
# -----------------------------------------------------

total_litter_class <- function(df, cls) {
  rowSums(df[, paste0("C_", cls, "_", c("A","W","E","N")), drop = FALSE],
          na.rm = TRUE)
}

classes      <- c("nwl", "fwl", "cwl")
class_labels <- c("Non-woody litter (NWL)", "Fine woody litter (FWL)",
                  "Coarse woody litter (CWL)")

png("./Testing/synthetic_litter_by_month.png", width = 1000, height = 380)
par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))

for (i in seq_along(classes)) {
  boxplot(
    total_litter_class(input_data, classes[i]) ~ input_data$month,
    main = paste("Monthly litter -", class_labels[i]),
    xlab = "Month",
    ylab = "Litter input (Mg C/ha/month)"
  )
}

par(mfrow = c(1, 1))
dev.off()


png("./Testing/synthetic_litter_by_class.png", width = 600, height = 500)
cols <- c("forestgreen", "sienna", "grey40")

boxplot(
  lapply(classes, function(cls) total_litter_class(input_data, cls)),
  names = class_labels,
  las   = 2,
  col   = cols,
  ylab  = "Mean monthly litter (Mg C/ha/month)",
  main  = "Relative magnitude of litter inputs by class"
)
dev.off()


# =============================================================================
# WRITE THE TEMPLATE WITH THE INPUT DATA
# =============================================================================
write.csv(input_data, "synthetic_input_data_template.csv", row.names = FALSE)
write.csv(site_data,  "synthetic_site_data_template.csv",  row.names = FALSE)

cat("Saved: synthetic_input_data_template.csv\n")
cat("Saved: synthetic_site_data_template.csv\n")
cat("Saved figures to ./Testing/\n")