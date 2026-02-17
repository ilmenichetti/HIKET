library(dplyr)
library(ggplot2)

# Load input compatibility layer
source("./Model_functions/input_compatibility.R")

# Read input data
input_raw <- read.csv("synthetic_input_data_template.csv")
site_raw  <- read.csv("synthetic_site_data_template.csv")

str(input_raw)

# -----------------------------------------------------------------------------
# Map all inputs and climates
# -----------------------------------------------------------------------------

monthly_climate <- map_climate_monthly(input_raw)
Yasso07_climate <- map_climate_yasso07(input_raw)
Yasso07_inputs  <- map_inputs_yasso07(input_raw)
Yasso20_inputs  <- map_inputs_yasso20(input_raw)

plots <- unique(Yasso07_climate$plot_id)


# =============================================================================
# Yasso07
# =============================================================================

dyn.load("./Model_functions/Decomposition_functions/Yasso/yasso07.so")
source("./Model_functions/Decomposition_functions/Yasso/yasso07_wrapper.R")

params <- YASSO07_DEFAULT_PARAMS

result_Yasso07 <- do.call(rbind, lapply(plots, function(pid) {
  
  # Subset climate and litter inputs to this plot
  clim   <- Yasso07_climate[Yasso07_climate$plot_id == pid, ]
  inputs <- Yasso07_inputs[Yasso07_inputs$plot_id   == pid, ]
  
  # Stage 1: Annual climate modifier (scalar per year)
  xi_array <- compute_xi_yasso07(
    temp_mean = clim$temp_mean,
    temp_amp  = clim$temp_amplitude,
    precip    = clim$precip,
    beta1     = params["beta1"],
    beta2     = params["beta2"],
    gamma     = params["gamma"]
  )
  
  # Stage 2: Steady-state initial pools using plot-specific mean climate and inputs
  C_init <- yasso07_steady_state(
    params   = params,
    nwl_mean = colMeans(inputs[, c("nwl_A","nwl_W","nwl_E","nwl_N")]),
    fwl_mean = colMeans(inputs[, c("fwl_A","fwl_W","fwl_E","fwl_N")]),
    cwl_mean = colMeans(inputs[, c("cwl_A","cwl_W","cwl_E","cwl_N")]),
    xi_mean  = mean(xi_array)   # mean over all years of this plot
  )
  
  # Stage 3: Forward transient simulation year by year
  out <- yasso07_run(
    input_df = inputs,
    params   = params,
    C_init   = C_init,
    xi_array = xi_array
  )
  
  cbind(plot_id = pid, out)
}))

print(result_Yasso07)


# =============================================================================
# Yasso15
# =============================================================================

dyn.load("./Model_functions/Decomposition_functions/Yasso/yasso15.so")
source("./Model_functions/Decomposition_functions/Yasso/yasso15_wrapper.R")

params <- YASSO15_DEFAULT_PARAMS
leac   <- 0.0

result_Yasso15 <- do.call(rbind, lapply(plots, function(pid) {
  
  # Subset climate and litter inputs to this plot
  # Yasso15 uses the same annual climate/inputs format as Yasso07
  clim   <- Yasso07_climate[Yasso07_climate$plot_id == pid, ]
  inputs <- Yasso07_inputs[Yasso07_inputs$plot_id   == pid, ]
  
  # Stage 1: Annual climate modifier, separately for AWE, N and H pool groups
  xi_arrays <- compute_xi_yasso15(
    temp_mean = clim$temp_mean,
    temp_amp  = clim$temp_amplitude,
    precip    = clim$precip,
    params    = params
  )
  cat("precip_mean:", mean(clim$precip), "\n")
  
  # Stage 2: Steady-state initial pools using plot-specific mean climate and inputs
  C_init <- yasso15_steady_state(
    params      = params,
    nwl_mean    = colMeans(inputs[, c("nwl_A","nwl_W","nwl_E","nwl_N")]),
    fwl_mean    = colMeans(inputs[, c("fwl_A","fwl_W","fwl_E","fwl_N")]),
    cwl_mean    = colMeans(inputs[, c("cwl_A","cwl_W","cwl_E","cwl_N")]),
    xi_means    = list(
      xi_awe = mean(xi_arrays$xi_awe),
      xi_n   = mean(xi_arrays$xi_n),
      xi_h   = mean(xi_arrays$xi_h)
    ),
    leac        = leac,
    precip_mean = mean(clim$precip)   # mean annual precip for leaching term
  )
  
  cat("plot:", pid, "\n")
  cat("xi_means awe/n/h:", 
      mean(xi_arrays$xi_awe), 
      mean(xi_arrays$xi_n), 
      mean(xi_arrays$xi_h), "\n")
  cat("C_init sum:", sum(C_init), "\n")
  
  cat("nwl_mean:", colMeans(inputs[, c("nwl_A","nwl_W","nwl_E","nwl_N")]), "\n")
  cat("nrow(inputs):", nrow(inputs), "\n")
  # Stage 3: Forward transient simulation year by year
  out <- yasso15_run(
    input_df  = inputs,
    params    = params,
    C_init    = C_init,
    xi_arrays = xi_arrays,
    precip    = clim$precip,   # annual precip vector for leaching
    leac      = leac
  )
  
  cbind(plot_id = pid, out)
}))

print(result_Yasso15)


# =============================================================================
# Yasso20
# =============================================================================

# yasso15.so already loaded above — Yasso20 reuses the same Fortran entry points
source("./Model_functions/Decomposition_functions/Yasso/yasso20_wrapper.R")

params <- YASSO20_DEFAULT_PARAMS
leac   <- 0.0

result_Yasso20 <- do.call(rbind, lapply(plots, function(pid) {
  
  # Subset monthly climate and annual litter inputs to this plot
  clim   <- monthly_climate[monthly_climate$plot_id == pid, ]
  inputs <- Yasso20_inputs[Yasso20_inputs$plot_id   == pid, ]
  
  # Stage 1: Annual climate modifier from all 12 monthly temperatures
  # (replaces the four-point Gaussian approximation used in Yasso07/15)
  xi_arrays <- compute_xi_yasso20(
    climate_df = clim,
    params     = params
  )
  
  # Stage 2: Steady-state initial pools
  # Build a synthetic single "mean year" (12 rows) from plot-specific
  # monthly temperature and precipitation averages across all years
  mean_clim <- data.frame(
    year     = 9999L,
    month    = 1:12,
    temp_air = tapply(clim$temp_air, clim$month, mean, na.rm = TRUE),
    precip   = mean(clim$precip)   # monthly mean as-is, NOT divided by 12
    # yasso20_steady_state sums over 12 months -> ~600mm
  )
  
  cat("precip_mean:", mean(tapply(mean_clim$precip, mean_clim$year, sum, na.rm=TRUE)), "\n")
  
  C_init <- yasso20_steady_state(
    climate_df = mean_clim,
    params     = params,
    nwl_mean   = colMeans(inputs[, c("nwl_A","nwl_W","nwl_E","nwl_N")]),
    fwl_mean   = colMeans(inputs[, c("fwl_A","fwl_W","fwl_E","fwl_N")]),
    cwl_mean   = colMeans(inputs[, c("cwl_A","cwl_W","cwl_E","cwl_N")]),
    leac       = leac
  )
  
  cat("plot:", pid, "\n")
  cat("xi_means awe/n/h:", 
      mean(xi_arrays$xi_awe), 
      mean(xi_arrays$xi_n), 
      mean(xi_arrays$xi_h), "\n")
  cat("C_init sum:", sum(C_init), "\n")
  
  cat("nwl_mean:", colMeans(inputs[, c("nwl_A","nwl_W","nwl_E","nwl_N")]), "\n")
  cat("nrow(inputs):", nrow(inputs), "\n")
  
  # Stage 3: Forward transient simulation year by year
  out <- yasso20_run(
    climate_df = clim,
    inputs_df  = inputs,
    params     = params,
    C_init     = C_init,
    leac       = leac
  )
  
  cbind(plot_id = pid, out)
}))

print(result_Yasso20)







# =============================================================================
# Yasso07 official R reference implementation
# =============================================================================

source("./Model_functions/Decomposition_functions/Yasso_original/y07_subroutine.r")

# Build 44-element parameter vector matching the official indexing
# [1-4]: alpha_A/W/E/N, [5-16]: transfers, [17-18]: beta1/2,
# [26]: gamma, [35]: alpha_H, [36]: p_H, [39-41]: delta1/delta2/r
PA <- numeric(44)
PA[1]  <- YASSO07_DEFAULT_PARAMS["alpha_A"]
PA[2]  <- YASSO07_DEFAULT_PARAMS["alpha_W"]
PA[3]  <- YASSO07_DEFAULT_PARAMS["alpha_E"]
PA[4]  <- YASSO07_DEFAULT_PARAMS["alpha_N"]
PA[5]  <- YASSO07_DEFAULT_PARAMS["p_WA"]
PA[6]  <- YASSO07_DEFAULT_PARAMS["p_EA"]
PA[7]  <- YASSO07_DEFAULT_PARAMS["p_NA"]
PA[8]  <- YASSO07_DEFAULT_PARAMS["p_AW"]
PA[9]  <- YASSO07_DEFAULT_PARAMS["p_EW"]
PA[10] <- YASSO07_DEFAULT_PARAMS["p_NW"]
PA[11] <- YASSO07_DEFAULT_PARAMS["p_AE"]
PA[12] <- YASSO07_DEFAULT_PARAMS["p_WE"]
PA[13] <- YASSO07_DEFAULT_PARAMS["p_NE"]
PA[14] <- YASSO07_DEFAULT_PARAMS["p_AN"]
PA[15] <- YASSO07_DEFAULT_PARAMS["p_WN"]
PA[16] <- YASSO07_DEFAULT_PARAMS["p_EN"]
PA[17] <- YASSO07_DEFAULT_PARAMS["beta1"]
PA[18] <- YASSO07_DEFAULT_PARAMS["beta2"]
PA[26] <- YASSO07_DEFAULT_PARAMS["gamma"]
PA[35] <- YASSO07_DEFAULT_PARAMS["alpha_H"]
PA[36] <- YASSO07_DEFAULT_PARAMS["p_H"]
PA[39] <- YASSO07_DEFAULT_PARAMS["delta1"]
PA[40] <- YASSO07_DEFAULT_PARAMS["delta2"]
PA[41] <- YASSO07_DEFAULT_PARAMS["r"]

result_Yasso07_ref <- do.call(rbind, lapply(plots, function(pid) {
  
  clim   <- Yasso07_climate[Yasso07_climate$plot_id == pid, ]
  inputs <- Yasso07_inputs[Yasso07_inputs$plot_id   == pid, ]
  years  <- clim$year
  
  # Reuse same steady-state C_init as Fortran Yasso07
  xi_ref <- compute_xi_yasso07(
    temp_mean = clim$temp_mean,
    temp_amp  = clim$temp_amplitude,
    precip    = clim$precip,
    beta1     = YASSO07_DEFAULT_PARAMS["beta1"],
    beta2     = YASSO07_DEFAULT_PARAMS["beta2"],
    gamma     = YASSO07_DEFAULT_PARAMS["gamma"]
  )
  C_init_ref <- yasso07_steady_state(
    params   = YASSO07_DEFAULT_PARAMS,
    nwl_mean = colMeans(inputs[, c("nwl_A","nwl_W","nwl_E","nwl_N")]),
    fwl_mean = colMeans(inputs[, c("fwl_A","fwl_W","fwl_E","fwl_N")]),
    cwl_mean = colMeans(inputs[, c("cwl_A","cwl_W","cwl_E","cwl_N")]),
    xi_mean  = mean(xi_ref)
  )
  
  # Transient simulation: step through years, accumulating pools
  # Three cohorts (nwl/fwl/cwl) run separately then summed
  C_nwl <- C_init_ref[1:5]
  C_fwl <- C_init_ref[6:10]
  C_cwl <- C_init_ref[11:15]
  
  out <- do.call(rbind, lapply(seq_along(years), function(i) {
    MT <- clim$temp_mean[i]
    TA <- clim$temp_amplitude[i]
    PR <- clim$precip[i]
    
    C_nwl <<- yasso07.light(MT, TA, PR, C_nwl,
                            c(as.numeric(inputs[i, c("nwl_A","nwl_W","nwl_E","nwl_N")]), 0),
                            0, PA, 1)
    C_fwl <<- yasso07.light(MT, TA, PR, C_fwl,
                            c(as.numeric(inputs[i, c("fwl_A","fwl_W","fwl_E","fwl_N")]), 0),
                            2, PA, 1)
    C_cwl <<- yasso07.light(MT, TA, PR, C_cwl,
                            c(as.numeric(inputs[i, c("cwl_A","cwl_W","cwl_E","cwl_N")]), 0),
                            15, PA, 1)
    
    pools <- C_nwl + C_fwl + C_cwl
    data.frame(year = years[i], total_soc = sum(pools))
  }))
  
  cbind(plot_id = pid, out)
}))



# =============================================================================
# Model comparison plot: total SOC over time, one panel per plot
# =============================================================================

# Combine results from all three models into a single long data frame
comparison <- rbind(
  result_Yasso07[,     c("plot_id", "year", "total_soc")] |> cbind(model = "Yasso07 (Fortran)"),
  result_Yasso15[,     c("plot_id", "year", "total_soc")] |> cbind(model = "Yasso15"),
  result_Yasso20[,     c("plot_id", "year", "total_soc")] |> cbind(model = "Yasso20"),
  result_Yasso07_ref[, c("plot_id", "year", "total_soc")] |> cbind(model = "Yasso07 (ref R)")
)

model_colors <- c(
  "Yasso07 (Fortran)" = "#2166ac",
  "Yasso15"           = "#d6604d",
  "Yasso20"           = "#1a9641",
  "Yasso07 (ref R)"   = "#756bb1"
)

p <- ggplot(comparison, aes(x = year, y = total_soc, color = model)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ plot_id, ncol = 2, scales = "free_y") +
  scale_color_manual(values = model_colors) +
  labs(
    title   = "Yasso model comparison — total SOC over time",
    x       = "Year",
    y       = expression("Total SOC (kg C m"^{-2}*")"),
    color   = "Model"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey92"),
    strip.text       = element_text(face = "bold"),
    legend.position  = "bottom",
    panel.grid.minor = element_blank()
  )

print(p)
