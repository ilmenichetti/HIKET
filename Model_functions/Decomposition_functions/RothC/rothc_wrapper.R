# =============================================================================
# RothC Wrapper Functions
# =============================================================================
#
# Compile:  system("R CMD SHLIB rothc_step_f.f90")
# Load:     dyn.load("rothc_step_f.so")
#
# Three-stage call pattern:
#   1. compute_xi_rothc()    -- climate modifier time series (R)
#   2. rothc_steady_state()  -- initial pool state C* (R, analytical)
#   3. rothc_run()           -- transient forward simulation (Fortran + R fallback)
#
# IMPORTANT: Unlike the previous version, compute_xi_rothc() is now
# parameter-dependent (a_T, b_T, c_T, r_SMD are calibration parameters).
# It must therefore be called INSIDE the MCMC loop at each proposal,
# mirroring how compute_xi_yasso07() is used.
#
# Clay partitioning (Section 1.7 of manual):
#   x         = 1.67 * (1.85 + 1.60 * exp(-0.0786 * clay))   [hardcoded]
#   f_bio_eff = f_bio / (x + 1)                               [clay-scaled]
#   f_hum_eff = f_hum / (x + 1)                               [clay-scaled]
#   f_co2_eff = 1 - f_bio_eff - f_hum_eff
#
# f_bio and f_hum are free calibration parameters (numerators of the
# clay partitioning fractions). The clay-dependent denominator (x + 1)
# is retained to preserve soil-texture sensitivity.
#
# Hardcoded (not calibrated):
#   evap fraction : 0.75   (open-pan evaporation coefficient)
#   r_bare        : 0.556  (bare soil TSMD limit as fraction of TSMD_15bar)
#   xi_min        : 0.20   (moisture modifier floor)
#   clay pedotransfer coefficients: 20.0, 1.3, 0.01, reference depth 23.0 cm
#   clay texture factor coefficients: 1.67, 1.85, 1.60, 0.0786
# =============================================================================


ROTHC_PARAM_NAMES <- c(
  "k_DPM", "k_RPM", "k_Bio", "k_Hum",  # decomposition rates [yr^-1]
  "a_T", "b_T", "c_T",                  # temperature response (Arrhenius)
  "r_SMD",                              # moisture response shape (TSMD breakpoint)
  "f_bio", "f_hum"                      # humification / inter-pool transfer numerators
)

ROTHC_DEFAULT_PARAMS <- c(
  10.00, 0.30, 0.66, 0.02,             # k_DPM, k_RPM, k_Bio, k_Hum
  47.91, 106.06, 18.27,                # a_T, b_T, c_T
  0.444,                               # r_SMD
  0.46, 0.54                           # f_bio, f_hum
)
names(ROTHC_DEFAULT_PARAMS) <- ROTHC_PARAM_NAMES


# =============================================================================
# Stage 1: compute_xi_rothc
# =============================================================================
# Computes the monthly combined rate modifier xi = RM_Tmp * RM_Moist * RM_PC
# from the climate data frame and the current parameter proposal.
#
# NOTE: Now parameter-dependent -- must be called inside the MCMC loop.
# The TSMD accumulation is path-dependent (sequential month-to-month state)
# but depends only on climate, not on SOC pools.
#
# Arguments:
#   climate_monthly -- data.frame, one row per month, columns:
#                        temp_air   : mean air temperature [deg C]
#                        precip     : rainfall [mm]
#                        evap       : open-pan evaporation [mm]
#                        soil_cover : 1 = vegetated, 0 = bare
#                        clay       : clay content [%]  (constant per plot)
#                        depth      : soil depth [cm]   (constant per plot)
#   a_T, b_T, c_T  -- Arrhenius temperature response parameters (calibrated)
#                     xi_T = a_T / (exp(b_T / (T + c_T)) + 1)
#   r_SMD          -- TSMD breakpoint as fraction of TSMD_15bar (calibrated)
#                     xi_M = 1.0  when TSMD > r_SMD * TSMD_15bar  (wet)
#                     xi_M declines linearly to xi_min as TSMD -> TSMD_15bar (dry)
#
# Returns:
#   numeric vector of length nrow(climate_monthly) with monthly xi values
# =============================================================================
compute_xi_rothc <- function(climate_monthly,
                             a_T   = 47.91,
                             b_T   = 106.06,
                             c_T   = 18.27,
                             r_SMD = 0.444) {

  n_months <- nrow(climate_monthly)
  xi       <- numeric(n_months)

  # Hardcoded constants
  xi_min <- 0.20    # moisture modifier floor
  r_bare <- 0.556   # bare soil TSMD limit as fraction of TSMD_15bar
  f_evap <- 0.75    # open-pan evaporation coefficient

  # Clay and depth are constant per plot
  clay  <- climate_monthly$clay[1]
  depth <- climate_monthly$depth[1]

  # Soil water capacity (standard RothC v1.0.0, Section 1.6.2)
  # Max TSMD for 0-23 cm layer, scaled to actual depth
  smd_15bar_base <- -(20.0 + 1.3 * clay - 0.01 * clay^2)
  smd_15bar      <- smd_15bar_base * depth / 23.0
  smd_1bar       <- r_SMD * smd_15bar   # calibratable breakpoint
  smd_bare       <- r_bare * smd_15bar  # hardcoded bare soil limit

  smd_acc <- 0.0  # running TSMD state; initialised at field capacity

  temp  <- climate_monthly$temp_air
  rain  <- climate_monthly$precip
  evap  <- climate_monthly$evap
  cover <- as.integer(climate_monthly$soil_cover)

  for (i in seq_len(n_months)) {

    # RM_Tmp: Arrhenius-type temperature modifier (Section 1.6.1)
    T      <- temp[i]
    rm_tmp <- if (T < -5.0) 0.0 else a_T / (exp(b_T / (T + c_T)) + 1.0)

    # RM_Moist: accumulate TSMD then evaluate piecewise linear response
    df <- rain[i] - f_evap * evap[i]
    if (cover[i] == 1L) {
      smd_acc <- max(smd_15bar, min(0.0, smd_acc + df))
    } else {
      smd_acc <- max(min(smd_bare, smd_acc), min(0.0, smd_acc + df))
    }

    rm_moist <- if (smd_acc > smd_1bar) {
      1.0
    } else {
      xi_min + (1.0 - xi_min) *
        (smd_15bar - smd_acc) / (smd_15bar - smd_1bar)
    }

    # RM_PC: plant cover modifier (Section 1.6.3) -- hardcoded, not calibrated
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
#   C* = I / (xi_mean * k)
#
# For Bio and Hum, solve the 2x2 linear system arising from mutual recycling:
#   k_Bio*(1-f_bio_eff)*C_Bio - f_bio_eff*k_Hum*C_Hum = f_bio_eff * S_ext
#  -f_hum_eff*k_Bio*C_Bio + k_Hum*(1-f_hum_eff)*C_Hum = f_hum_eff * S_ext
# where S_ext = k_DPM*C_DPM + k_RPM*C_RPM
#
# IOM: Falloon (1998) fixed-point  IOM = 0.049 * (C_active + IOM)^1.139
#
# Arguments:
#   params       -- named numeric vector (ROTHC_PARAM_NAMES)
#   xi_mean      -- mean of the monthly xi array (from compute_xi_rothc)
#   litter_input -- named vector: DPM and RPM annual C input [t C ha^-1 yr^-1]
#   clay         -- clay content [%]
#
# Returns:
#   named numeric vector: DPM, RPM, Bio, Hum, IOM  [t C ha^-1]
# =============================================================================
rothc_steady_state <- function(params, xi_mean, litter_input, clay) {

  k_DPM <- params[["k_DPM"]]
  k_RPM <- params[["k_RPM"]]
  k_Bio <- params[["k_Bio"]]
  k_Hum <- params[["k_Hum"]]
  f_bio <- params[["f_bio"]]
  f_hum <- params[["f_hum"]]

  # Clay-scaled effective partitioning fractions (Section 1.7)
  # Numerators f_bio, f_hum are free parameters; denominator retains clay effect
  x         <- 1.67 * (1.85 + 1.60 * exp(-0.0786 * clay))
  f_bio_eff <- f_bio / (x + 1)
  f_hum_eff <- f_hum / (x + 1)

  # Steady state for pools receiving only external input
  C_DPM_ss <- litter_input[["DPM"]] / (xi_mean * k_DPM)
  C_RPM_ss <- litter_input[["RPM"]] / (xi_mean * k_RPM)

  # External decomposition flux driving Bio and Hum
  S_ext <- k_DPM * C_DPM_ss + k_RPM * C_RPM_ss

  # 2x2 linear system for Bio and Hum at steady state
  A <- matrix(c(
    k_Bio * (1 - f_bio_eff),  -f_bio_eff * k_Hum,
    -f_hum_eff * k_Bio,        k_Hum * (1 - f_hum_eff)
  ), nrow = 2, byrow = TRUE)

  sol      <- solve(A, c(f_bio_eff * S_ext, f_hum_eff * S_ext))
  C_Bio_ss <- sol[1]
  C_Hum_ss <- sol[2]

  # IOM: Falloon (1998) fixed-point iteration
  C_active  <- unname(C_DPM_ss + C_RPM_ss + C_Bio_ss + C_Hum_ss)
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
# Forward transient simulation. Iterates the monthly step function using the
# precomputed xi array. Calls Fortran if loaded; falls back to pure R otherwise.
#
# Arguments:
#   pools0         -- named vector: initial DPM, RPM, Bio, Hum, IOM [t C ha^-1]
#   params         -- named parameter vector (ROTHC_PARAM_NAMES)
#   xi             -- numeric vector of monthly xi values (from compute_xi_rothc)
#   litter_monthly -- data.frame: C_DPM, C_RPM [t C ha^-1 month^-1],
#                     one row per month
#   clay           -- clay content [%]
#   obs_times      -- integer vector of month indices at which to record SOC
#
# Returns:
#   data.frame with columns: DPM, RPM, Bio, Hum, IOM, total_soc
#   one row per obs_times entry  [t C ha^-1]
# =============================================================================
rothc_run <- function(pools0, params, xi, litter_monthly, clay, obs_times) {

  n_months <- length(xi)
  if (nrow(litter_monthly) != n_months)
    stop("xi length and litter_monthly rows must match")

  k_DPM <- params[["k_DPM"]]
  k_RPM <- params[["k_RPM"]]
  k_Bio <- params[["k_Bio"]]
  k_Hum <- params[["k_Hum"]]
  f_bio <- params[["f_bio"]]
  f_hum <- params[["f_hum"]]

  # Clay-scaled effective partitioning fractions
  x         <- 1.67 * (1.85 + 1.60 * exp(-0.0786 * clay))
  f_bio_eff <- f_bio / (x + 1)
  f_hum_eff <- f_hum / (x + 1)

  dt <- 1.0 / 12.0

  DPM <- pools0[["DPM"]]
  RPM <- pools0[["RPM"]]
  Bio <- pools0[["Bio"]]
  Hum <- pools0[["Hum"]]
  IOM <- pools0[["IOM"]]

  use_fortran <- is.loaded("rothc_step_f")
  obs_set     <- as.integer(obs_times)
  n_obs       <- length(obs_set)
  results     <- matrix(0.0, nrow = n_obs, ncol = 5,
                        dimnames = list(NULL, c("DPM","RPM","Bio","Hum","IOM")))
  ri <- 1L

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
                      k_DPM    = as.double(k_DPM),
                      k_RPM    = as.double(k_RPM),
                      k_Bio    = as.double(k_Bio),
                      k_Hum    = as.double(k_Hum),
                      f_bio    = as.double(f_bio_eff),
                      f_hum    = as.double(f_hum_eff),
                      C_in_DPM = as.double(C_in_DPM),
                      C_in_RPM = as.double(C_in_RPM))
      DPM <- res$DPM; RPM <- res$RPM; Bio <- res$Bio; Hum <- res$Hum

    } else {
      DPM1 <- DPM * exp(-xi_i * k_DPM * dt)
      RPM1 <- RPM * exp(-xi_i * k_RPM * dt)
      Bio1 <- Bio * exp(-xi_i * k_Bio * dt)
      Hum1 <- Hum * exp(-xi_i * k_Hum * dt)
      dsum <- (DPM - DPM1) + (RPM - RPM1) + (Bio - Bio1) + (Hum - Hum1)
      DPM  <- DPM1 + C_in_DPM
      RPM  <- RPM1 + C_in_RPM
      Bio  <- Bio1 + f_bio_eff * dsum
      Hum  <- Hum1 + f_hum_eff * dsum
    }

    if (ri <= n_obs && i == obs_set[ri]) {
      results[ri, ] <- c(DPM, RPM, Bio, Hum, IOM)
      ri <- ri + 1L
    }
  }

  df           <- as.data.frame(results)
  df$total_soc <- df$DPM + df$RPM + df$Bio + df$Hum + df$IOM
  df
}
