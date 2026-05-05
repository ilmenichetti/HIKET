# =============================================================================
# generate_synthetic_data.R
#
# Generates synthetic_input_data_template.csv and synthetic_site_data_template.csv
# for testing the SOC model calibration framework.
#
# Designed to be realistic for Finnish boreal forests:
#   - 4 plots, 1990-2020, monthly resolution
#   - Climate based on Finnish temperature/precipitation gradients
#   - Litter inputs consistent with boreal forest productivity
#   - SOC observations generated from Yasso07 at default parameters
#
# Sources of structure beyond random noise:
#   1. Plot-level disturbance offset (deviation from steady-state at t0)
#      Mimics harvest history -- some plots are below/above SS at sim start
#   2. Multiplicative measurement bias per plot
#      Mimics bulk density estimation errors in field SOC measurement
#   3. Proportional Gaussian observation noise
#
# TRUE parameter values stored as attributes for calibration validation.
# =============================================================================

library(dplyr)

source("./Model_functions/input_compatibility_layer.R")
source("./Model_functions/Decomposition_functions/Yasso/yasso07_wrapper.R")
dyn.load("./Model_functions/Decomposition_functions/Yasso/yasso07.so")

set.seed(42)

# =============================================================================
# TRUE noise / bias parameters
# (what the calibration should approximately recover)
# =============================================================================

TRUE_SIGMA_OBS  <- 0.05   # 5%  proportional Gaussian observation noise
TRUE_SIGMA_INIT <- 0.10   # 10% proportional initialization uncertainty

# Plot-level multiplicative measurement bias
# (systematic over- or under-estimation of SOC stock at field measurement)
# Centred on 1, modest spread -- plausible for bulk density uncertainty
TRUE_MEAS_BIAS <- c(
  PLOT_001 = 1.00,   # unbiased reference plot
  PLOT_002 = 1.08,   # slight overestimation (low bulk density assumed)
  PLOT_003 = 0.93,   # slight underestimation (stony soil, BD overestimated)
  PLOT_004 = 1.04
)

# Plot-level disturbance multiplier on initial SOC
# Values < 1: plot was recently disturbed (harvest), SOC below SS
# Values > 1: plot is accumulating (post-disturbance recovery), SOC above SS
TRUE_DISTURBANCE <- c(
  PLOT_001 = 0.85,   # harvested ~15 years before sim start
  PLOT_002 = 1.00,   # close to steady state
  PLOT_003 = 0.70,   # heavily disturbed, still recovering
  PLOT_004 = 1.10    # old growth, slight accumulation
)


# =============================================================================
# 1. Site data
# =============================================================================

site_data <- data.frame(
  plot_id     = c("PLOT_001","PLOT_002","PLOT_003","PLOT_004"),
  pine_prop   = c(0.30, 0.20, 0.40, 0.25),
  spruce_prop = c(0.60, 0.70, 0.50, 0.65),
  birch_prop  = c(0.10, 0.10, 0.10, 0.10),
  clay        = c(18,   12,   25,   15),
  soil_depth  = c(23,   18,   30,   20),
  stringsAsFactors = FALSE
)

# Normalise species proportions
s <- site_data$pine_prop + site_data$spruce_prop + site_data$birch_prop
site_data$pine_prop   <- site_data$pine_prop   / s
site_data$spruce_prop <- site_data$spruce_prop / s
site_data$birch_prop  <- site_data$birch_prop  / s


# =============================================================================
# 2. Monthly input data (climate + litter)
# =============================================================================

plot_ids <- c("PLOT_001","PLOT_002","PLOT_003","PLOT_004")
years    <- 1990:2020

# Climate baseline by plot (Finnish gradient: colder/wetter north to south)
plot_climate <- data.frame(
  plot_id    = plot_ids,
  base_temp  = c( 3.5,  2.5,  2.0,  1.5),   # mean annual temp (C)
  precip_mult= c( 1.00, 1.05, 1.10, 0.95)   # precipitation multiplier
)

# Litter totals (Mg C/ha/yr) -- boreal forest range
# Slightly different productivity per plot
litter_totals_base <- list(
  foliage      = 1.20,
  branches     = 0.40,
  stems        = 0.10,
  fine_roots   = 0.80,
  coarse_roots = 0.30,
  understorey  = 0.30
)
plot_litter_mult <- c(PLOT_001=1.00, PLOT_002=1.05, PLOT_003=0.90, PLOT_004=1.10)

# AWEN fractions by cohort (Finnish boreal, Liski et al. style)
awen_props <- list(
  foliage      = c(A=0.52, W=0.10, E=0.02, N=0.36),
  branches     = c(A=0.50, W=0.05, E=0.03, N=0.42),
  stems        = c(A=0.45, W=0.02, E=0.01, N=0.52),
  fine_roots   = c(A=0.44, W=0.05, E=0.01, N=0.50),
  coarse_roots = c(A=0.45, W=0.03, E=0.01, N=0.51),
  understorey  = c(A=0.50, W=0.15, E=0.03, N=0.32)
)

# Build grid
grid <- expand.grid(
  plot_id = plot_ids,
  year    = years,
  month   = 1:12,
  stringsAsFactors = FALSE
)
grid <- grid[order(grid$plot_id, grid$year, grid$month), ]
n    <- nrow(grid)

# -- Climate --
clim_base <- plot_climate[match(grid$plot_id, plot_climate$plot_id), ]

# Temperature: seasonal cycle + long-term warming trend + interannual noise
# Warming trend: +0.04 C/year (consistent with observed Finnish trends)
warming_trend <- 0.04 * (grid$year - min(years))
seasonal_temp <- 12 * sin(2 * pi * (grid$month - 4) / 12)

grid$temp_air <- round(
  clim_base$base_temp + seasonal_temp + warming_trend + rnorm(n, 0, 1.2),
  1
)

# Precipitation: seasonal cycle + plot multiplier + interannual noise
monthly_precip_base <- c(40, 35, 35, 40, 45, 55, 70, 75, 60, 55, 50, 45)
grid$precip <- round(
  monthly_precip_base[grid$month] *
    clim_base$precip_mult *
    (1 + rnorm(n, 0, 0.20)),
  1
)
grid$precip <- pmax(grid$precip, 0)

# Evapotranspiration: temperature-driven with floor
grid$evap <- round(
  pmax(5, 20 + 40 * sin(2 * pi * (grid$month - 4) / 12) + rnorm(n, 0, 5)),
  1
)

# -- Litter inputs --
# Annual litter varies with a decadal cycle (stand age dynamics) + small noise
year_effect <- 1 + 0.08 * sin(2 * pi * (grid$year - min(years)) / 12)

for (coh in names(litter_totals_base)) {
  base_annual <- litter_totals_base[[coh]] / 12   # monthly
  props       <- awen_props[[coh]]
  pmult       <- plot_litter_mult[grid$plot_id]
  
  for (frac in c("A","W","E","N")) {
    col <- paste0("C_", coh, "_", frac)
    grid[[col]] <- round(base_annual * props[frac] * pmult * year_effect, 6)
  }
}

# -- SOC observations placeholder --
grid$soc_obs_tCha <- NA_real_


# =============================================================================
# 3. Generate SOC observations from Yasso07 with realistic perturbations
# =============================================================================

soc_meas_years <- c(1997, 2006)
soc_meas_month <- 7

Yasso07_climate <- map_climate_yasso07(grid)
Yasso07_inputs  <- map_inputs_yasso07(grid)
p               <- YASSO07_DEFAULT_PARAMS

for (pid in plot_ids) {
  
  clim   <- Yasso07_climate[Yasso07_climate$plot_id == pid, ]
  inputs <- Yasso07_inputs[Yasso07_inputs$plot_id   == pid, ]
  
  # Stage 1: xi
  xi_array <- compute_xi_yasso07(
    temp_mean = clim$temp_mean,
    temp_amp  = clim$temp_amplitude,
    precip    = clim$precip,
    beta1     = p[["beta1"]],
    beta2     = p[["beta2"]],
    gamma     = p[["gamma"]]
  )
  xi_mean <- mean(xi_array)
  
  # Stage 2: steady-state pools
  C_ss <- yasso07_steady_state(
    params   = p,
    nwl_mean = colMeans(inputs[, c("nwl_A","nwl_W","nwl_E","nwl_N")]),
    fwl_mean = colMeans(inputs[, c("fwl_A","fwl_W","fwl_E","fwl_N")]),
    cwl_mean = colMeans(inputs[, c("cwl_A","cwl_W","cwl_E","cwl_N")]),
    xi_mean  = xi_mean
  )
  
  # Stage 3: transient simulation
  # Apply disturbance multiplier to initial pools -- this is the key perturbation.
  # The model runs from a disturbed initial state, not steady-state.
  # This is what sigma_init in the likelihood is meant to absorb.
  C_init <- C_ss * TRUE_DISTURBANCE[pid]
  
  run_out <- yasso07_run(
    input_df = inputs,
    params   = p,
    C_init   = C_init,
    xi_array = xi_array
  )
  
  sim_years <- clim$year
  
  for (k in seq_along(soc_meas_years)) {
    
    yy      <- soc_meas_years[k]
    idx_grid <- which(grid$plot_id == pid &
                        grid$year  == yy &
                        grid$month == soc_meas_month)
    idx_sim  <- which(sim_years == yy)
    
    if (length(idx_grid) != 1 || length(idx_sim) != 1) next
    
    soc_true <- run_out$total_soc[idx_sim]
    
    # Apply multiplicative measurement bias (systematic, same both years)
    soc_measured <- soc_true * TRUE_MEAS_BIAS[pid]
    
    # Add proportional Gaussian noise
    # First observation: obs noise + residual init uncertainty
    # Subsequent observations: obs noise only
    sigma <- if (k == 1L) {
      soc_measured * sqrt(TRUE_SIGMA_OBS^2 + TRUE_SIGMA_INIT^2)
    } else {
      soc_measured * TRUE_SIGMA_OBS
    }
    
    grid$soc_obs_tCha[idx_grid] <- round(rnorm(1, mean = soc_measured, sd = sigma), 2)
  }
}


# =============================================================================
# 4. Summary and validation
# =============================================================================

cat("=== Synthetic data summary ===\n\n")
cat(sprintf("Plots: %d | Years: %d-%d | Months: %d total rows\n\n",
            length(plot_ids), min(years), max(years), nrow(grid)))

cat("SOC observations:\n")
print(grid %>%
        filter(!is.na(soc_obs_tCha)) %>%
        select(plot_id, year, month, soc_obs_tCha) %>%
        arrange(plot_id, year))

cat("\nTrue perturbation parameters:\n")
cat(sprintf("  sigma_obs  = %.2f (%.0f%% proportional noise)\n",
            TRUE_SIGMA_OBS, TRUE_SIGMA_OBS * 100))
cat(sprintf("  sigma_init = %.2f (%.0f%% initialization uncertainty)\n",
            TRUE_SIGMA_INIT, TRUE_SIGMA_INIT * 100))
cat("\nPlot-level disturbance multipliers (deviation from steady-state):\n")
for (pid in plot_ids) {
  cat(sprintf("  %s: %.2f (%+.0f%% from SS)\n",
              pid, TRUE_DISTURBANCE[pid],
              (TRUE_DISTURBANCE[pid] - 1) * 100))
}
cat("\nPlot-level measurement bias:\n")
for (pid in plot_ids) {
  cat(sprintf("  %s: %.2f (%+.0f%%)\n",
              pid, TRUE_MEAS_BIAS[pid],
              (TRUE_MEAS_BIAS[pid] - 1) * 100))
}

cat("\nNote: measurement bias is NOT captured by the current likelihood --\n")
cat("it will appear as residual misfit. This is intentional and realistic.\n")


# =============================================================================
# 5. Write output
# =============================================================================

input_data <- grid

write.csv(input_data, "synthetic_input_data_template.csv", row.names = FALSE)
write.csv(site_data,  "synthetic_site_data_template.csv",  row.names = FALSE)

cat("\nWritten:\n")
cat("  synthetic_input_data_template.csv\n")
cat("  synthetic_site_data_template.csv\n")