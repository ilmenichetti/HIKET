# =============================================================================
# RothC Wrapper Functions
# =============================================================================
# Architecture mirrors Yasso07/15/20 three-stage pattern:
#
#   Stage 1 (outside MCMC, once per plot):
#     compute_xi_rothc() -- reads temp_air, precip, evap, soil_cover, clay,
#                           depth directly from climate_monthly data frame
#                           (as produced by map_climate_monthly(site_df=))
#
#   Stage 2 (inside MCMC, parameter-dependent):
#     rothc_steady_state() -- analytical 4-pool steady state + IOM fixed-point
#
#   Stage 3 (inside MCMC, forward simulation):
#     rothc_run() -- calls Fortran step function month by month with xi array
#
# Standard options: opt_RMmoist=1, opt_SMDbare=1, minRM_Moist=0.2
#
# Rate constants (fixed, not calibrated):
#   DPM_k = 10.0,  RPM_k = 0.3,  Bio_k = 0.66,  Hum_k = 0.02  [yr^-1]
#
# Clay partitioning (Section 1.7 of manual):
#   x = 1.67 * (1.85 + 1.60 * exp(-0.0786 * clay))
#   f_co2 = x / (x + 1)
#   f_bio = 0.46 / (x + 1)
#   f_hum = 0.54 / (x + 1)
# =============================================================================

# Fixed rate constants [yr^-1]
.ROTHC_K <- c(DPM = 10.0, RPM = 0.3, Bio = 0.66, Hum = 0.02)


# =============================================================================
# Stage 1: compute_xi_rothc
# =============================================================================
# Precomputes the monthly combined rate modifier xi = RM_Tmp * RM_Moist * RM_PC
# from the climate data frame alone (no model parameters involved).
#
# Reads all required columns from climate_monthly, which must be produced by
# map_climate_monthly(site_df=) so that clay, depth and soil_cover are present.
#
# The TSMD accumulation is path-dependent (carries month to month), so this
# must iterate sequentially. It is parameter-free, so it runs ONCE outside
# the MCMC loop -- identical to how compute_xi_yasso07() is used.
#
# Arguments:
#   climate_monthly -- data.frame, one row per month, columns:
#                        temp_air   : mean air temperature [deg C]
#                        precip     : rainfall [mm]
#                        evap       : open-pan evaporation [mm]
#                        soil_cover : 1 = vegetated, 0 = bare
#                        clay       : clay content [%]  (constant per plot)
#                        depth      : soil depth [cm]   (constant per plot)
#                      All present when map_climate_monthly(site_df=) is used.
#   opt_RMmoist     -- moisture option (1 = standard v1.0.0; default 1)
#   minRM_Moist     -- minimum RM_Moist floor (default 0.2)
#
# Returns:
#   numeric vector of length nrow(climate_monthly) with monthly xi values
# =============================================================================
compute_xi_rothc <- function(climate_monthly,
                             opt_RMmoist = 1L,
                             minRM_Moist = 0.2) {
  
  n_months <- nrow(climate_monthly)
  xi       <- numeric(n_months)
  
  # Extract columns -- clay and depth constant per plot, take first value
  clay  <- climate_monthly$clay[1]
  depth <- climate_monthly$depth[1]
  
  # Soil water capacity parameters (standard v1.0.0)
  # Max TSMD for 0-23 cm layer, scaled to actual depth
  smd_15bar_base <- -(20.0 + 1.3 * clay - 0.01 * clay^2)
  smd_15bar      <- smd_15bar_base * depth / 23.0
  smd_1bar       <- 0.444 * smd_15bar
  smd_bare       <- 0.556 * smd_15bar
  
  smd_acc <- 0.0  # running state; starts at field capacity
  
  temp  <- climate_monthly$temp_air
  rain  <- climate_monthly$precip
  evap  <- climate_monthly$evap
  cover <- as.integer(climate_monthly$soil_cover)
  
  for (i in seq_len(n_months)) {
    
    # RM_Tmp: Arrhenius-type (eq. 1 in manual)
    T <- temp[i]
    rm_tmp <- if (T < -5.0) 0.0 else 47.91 / (exp(106.06 / (T + 18.27)) + 1.0)
    
    # RM_Moist: accumulate TSMD then evaluate
    df <- rain[i] - 0.75 * evap[i]
    if (cover[i] == 1L) {
      smd_acc <- max(smd_15bar, min(0.0, smd_acc + df))
    } else {
      smd_acc <- max(min(smd_bare, smd_acc), min(0.0, smd_acc + df))
    }
    
    rm_moist <- if (smd_acc > smd_1bar) {
      1.0
    } else {
      minRM_Moist + (1.0 - minRM_Moist) *
        (smd_15bar - smd_acc) / (smd_15bar - smd_1bar)
    }
    
    # RM_PC: plant cover factor (Section 1.6.3)
    rm_pc <- if (cover[i] == 1L) 0.6 else 1.0
    
    xi[i] <- rm_tmp * rm_moist * rm_pc
  }
  
  xi
}


# =============================================================================
# Stage 2: rothc_steady_state
# =============================================================================
# Analytical steady-state for the 4 active pools (DPM, RPM, Bio, Hum) under a
# constant effective rate modifier xi_mean, plus IOM via Falloon fixed-point.
#
# At steady state for DPM and RPM (receive only external input):
#   C* = I / (xi * k)
#
# For Bio and Hum, solve the 2x2 linear system:
#   k_Bio*(1-f_bio)*C_Bio - f_bio*k_Hum*C_Hum = f_bio * S_ext
#  -f_hum*k_Bio*C_Bio + k_Hum*(1-f_hum)*C_Hum = f_hum * S_ext
# where S_ext = k_DPM*C_DPM + k_RPM*C_RPM
#
# IOM: Falloon (1998) fixed-point  IOM = 0.049 * (C_active + IOM)^1.139
#
# Arguments:
#   xi_mean      -- mean of the monthly xi array (from compute_xi_rothc)
#   litter_input -- named vector: DPM and RPM annual C input [t C ha^-1 yr^-1]
#   clay         -- clay content [%]
#
# Returns:
#   named numeric vector: DPM, RPM, Bio, Hum, IOM  [t C ha^-1]
# =============================================================================
rothc_steady_state <- function(xi_mean, litter_input, clay) {
  
  k <- .ROTHC_K
  
  # Clay partitioning fractions (Section 1.7)
  x     <- 1.67 * (1.85 + 1.60 * exp(-0.0786 * clay))
  f_bio <- 0.46 / (x + 1)
  f_hum <- 0.54 / (x + 1)
  
  C_DPM_ss <- litter_input[["DPM"]] / (xi_mean * k["DPM"])
  C_RPM_ss <- litter_input[["RPM"]] / (xi_mean * k["RPM"])
  
  S_ext <- k["DPM"] * C_DPM_ss + k["RPM"] * C_RPM_ss
  
  A <- matrix(c(
    k["Bio"] * (1 - f_bio), -f_bio * k["Hum"],
    -f_hum * k["Bio"],        k["Hum"] * (1 - f_hum)
  ), nrow = 2, byrow = TRUE)
  
  sol      <- solve(A, c(f_bio * S_ext, f_hum * S_ext))
  C_Bio_ss <- sol[1]
  C_Hum_ss <- sol[2]
  
  C_active <- unname(C_DPM_ss + C_RPM_ss + C_Bio_ss + C_Hum_ss)
  
  iom_fun   <- function(iom) iom - 0.049 * (C_active + iom)^1.139
  iom_upper <- max(0.049 * (C_active * 10)^1.139, 10)
  IOM_ss <- tryCatch(
    uniroot(iom_fun, lower = 0, upper = iom_upper, tol = 1e-10)$root,
    error = function(e) 0.049 * C_active^1.139
  )
  
  c(DPM = unname(C_DPM_ss),
    RPM = unname(C_RPM_ss),
    Bio = unname(C_Bio_ss),
    Hum = unname(C_Hum_ss),
    IOM = IOM_ss)
}


# =============================================================================
# Stage 3: rothc_run
# =============================================================================
# Forward transient simulation. Iterates the step function over all months,
# using the precomputed xi array. Falls back to pure R if Fortran not loaded.
#
# Arguments:
#   pools0         -- named vector: initial DPM, RPM, Bio, Hum, IOM
#   xi             -- numeric vector of monthly xi values (from compute_xi_rothc)
#   litter_monthly -- data.frame: C_DPM, C_RPM [t C ha^-1 month^-1],
#                     one row per month (from map_inputs_rothc, expanded monthly)
#   clay           -- clay content [%]
#   obs_times      -- integer vector of months at which to record SOC
#
# Returns:
#   data.frame with columns: DPM, RPM, Bio, Hum, IOM, total_soc
#   one row per obs_times entry  [t C ha^-1]
# =============================================================================
rothc_run <- function(pools0, xi, litter_monthly, clay, obs_times) {
  
  n_months <- length(xi)
  if (nrow(litter_monthly) != n_months)
    stop("xi length and litter_monthly rows must match")
  
  DPM <- pools0[["DPM"]]
  RPM <- pools0[["RPM"]]
  Bio <- pools0[["Bio"]]
  Hum <- pools0[["Hum"]]
  IOM <- pools0[["IOM"]]
  
  x     <- 1.67 * (1.85 + 1.60 * exp(-0.0786 * clay))
  f_bio <- 0.46 / (x + 1)
  f_hum <- 0.54 / (x + 1)
  
  k  <- .ROTHC_K
  dt <- 1.0 / 12.0
  
  use_fortran <- is.loaded("rothc_step_f")
  obs_set     <- as.integer(obs_times)
  n_obs       <- length(obs_set)
  results     <- matrix(0, nrow = n_obs, ncol = 5,
                        dimnames = list(NULL, c("DPM","RPM","Bio","Hum","IOM")))
  ri          <- 1L
  
  for (i in seq_len(n_months)) {
    
    xi_i     <- xi[i]
    C_in_DPM <- litter_monthly$C_DPM[i]
    C_in_RPM <- litter_monthly$C_RPM[i]
    
    if (use_fortran) {
      res <- .Fortran("rothc_step_f",
                      DPM      = as.double(DPM),
                      RPM      = as.double(RPM),
                      Bio      = as.double(Bio),
                      Hum      = as.double(Hum),
                      IOM      = as.double(IOM),
                      xi       = as.double(xi_i),
                      clay     = as.double(clay),
                      C_in_DPM = as.double(C_in_DPM),
                      C_in_RPM = as.double(C_in_RPM))
      DPM <- res$DPM; RPM <- res$RPM; Bio <- res$Bio; Hum <- res$Hum
      
    } else {
      DPM1 <- DPM * exp(-xi_i * k["DPM"] * dt)
      RPM1 <- RPM * exp(-xi_i * k["RPM"] * dt)
      Bio1 <- Bio * exp(-xi_i * k["Bio"] * dt)
      Hum1 <- Hum * exp(-xi_i * k["Hum"] * dt)
      dsum <- (DPM - DPM1) + (RPM - RPM1) + (Bio - Bio1) + (Hum - Hum1)
      DPM  <- DPM1 + C_in_DPM
      RPM  <- RPM1 + C_in_RPM
      Bio  <- Bio1 + f_bio * dsum
      Hum  <- Hum1 + f_hum * dsum
    }
    
    if (ri <= length(obs_set) && i == obs_set[ri]) {
      results[ri, ] <- c(DPM, RPM, Bio, Hum, IOM)
      ri <- ri + 1L
    }
  }
  
  df           <- as.data.frame(results)
  df$total_soc <- df$DPM + df$RPM + df$Bio + df$Hum + df$IOM
  df
}