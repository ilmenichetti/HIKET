library(dplyr)
library(ggplot2)

source("./Model_functions/input_compatibility_layer.R")
source("./Model_functions/Decomposition_functions/RothC/rothc_wrapper.R")

input_raw <- read.csv("synthetic_input_data_template.csv")
site_raw  <- read.csv("synthetic_site_data_template.csv")

str(input_raw)

# -----------------------------------------------------------------------------
# Map all inputs and climates
# -----------------------------------------------------------------------------

# site_df supplies clay and soil_depth into the monthly climate data frame
monthly_climate <- map_climate_monthly(input_raw, site_df = site_raw)
RothC_inputs    <- map_inputs_rothc(input_raw)

plots <- unique(monthly_climate$plot_id)


# =============================================================================
# RothC
# =============================================================================

dyn.load("./Model_functions/Decomposition_functions/RothC/rothc_step_f.so")

params <- list(PL_DPM_f = 0.20)
result_RothC <- do.call(rbind, lapply(plots, function(pid) {
  
  clim <- monthly_climate[monthly_climate$plot_id == pid, ]
  inp  <- RothC_inputs[RothC_inputs$plot_id       == pid, ]
  
  clay  <- clim$clay[1]
  depth <- clim$depth[1]
  
  xi <- compute_xi_rothc(climate_monthly = clim)
  
  litter_ss <- c(
    DPM = mean(inp$C_DPM) * 12,
    RPM = mean(inp$C_RPM) * 12
  )
  xi_mean <- mean(xi[seq_len(min(12L, length(xi)))])
  
  C_init <- rothc_steady_state(
    xi_mean      = xi_mean,
    litter_input = litter_ss,
    clay         = clay
  )
  
  cat("plot:", pid, "\n")
  cat("xi_mean:", round(xi_mean, 4), "\n")
  cat("C_init (DPM/RPM/Bio/Hum/IOM):",
      paste(round(C_init, 4), collapse = " / "), "\n")
  cat("C_init total:", round(sum(C_init), 3), "\n")
  
  litter_monthly <- data.frame(
    C_DPM = inp$C_DPM,
    C_RPM = inp$C_RPM
  )
  
  obs_months <- seq(12L, nrow(litter_monthly), by = 12L)
  
  cbind(plot_id = pid,
        data.frame(year = unique(inp$year)),
        rothc_run(
          pools0         = C_init,
          xi             = xi,
          litter_monthly = litter_monthly,
          clay           = clay,
          obs_times      = obs_months
        ))
}))

print(result_RothC)







# =============================================================================
# RothC reference validation
# =============================================================================
# Calls the original Rothamsted RothC.for Fortran subroutine directly via
# .Fortran() and compares pool-by-pool output against our rothc_run() wrapper.
#
# The official subroutine signature (from RothC.for):
#
#   Subroutine RothC(timeFact, DPM, RPM, Bio, Hum, IOM, SOC, total_CO2,
#                    DPM_Rage, RPM_Rage, Bio_Rage, Hum_Rage, Total_Rage,
#                    modernC, clay, depth, TEMP, RAIN, PEVAP, PC,
#                    PL_DPM_f, PL_RPM_f, OA_DPM_f, OA_RPM_f, OA_Bio_f, OA_Hum_f,
#                    C_Inp, OA_Inp, SMD, RM_TMP, RM_Moist, RM_PC,
#                    opt_RMmoist, opt_SMDbare, silt, BD, OC, minRM_Moist)
#
# Key design notes:
#   - Called one month at a time; SMD is passed in/out and carried between calls
#   - Takes total C_inp + PL_DPM_f/PL_RPM_f fractions (not pre-split)
#   - We use opt_RMmoist=1, opt_SMDbare=1 (standard v1.0.0, matching our wrapper)
#   - Radiocarbon arguments (modernC, *_Rage) are dummied out: modernC=1.0,
#     Rage values initialised to 0.0 and ignored
#   - OA_Inp=0 throughout (no organic amendment in our dataset)
# =============================================================================

# Load the compiled original Fortran subroutine
# Compile once with:
#   gfortran -shared -fPIC -o ./Model_functions/Decomposition_functions/RothC_original/RothC_ref.so 
#            ./Model_functions/Decomposition_functions/RothC_original/RothC.for
dyn.load("./Model_functions/Decomposition_functions/RothC_original/RothC_ref.so")

# -----------------------------------------------------------------------------
# Helper: call RothC.for for one month, carrying SMD as state
# -----------------------------------------------------------------------------
# Returns a list with updated pool values and the output SMD.
# The caller is responsible for passing smd_prev from the previous call.
.rothc_ref_step <- function(DPM, RPM, Bio, Hum, IOM,
                            temp, rain, evap, pc,
                            C_inp, PL_DPM_f, PL_RPM_f,
                            clay, depth, smd_prev,
                            modernC     = 1.0,
                            total_CO2   = 0.0,
                            opt_RMmoist = 1L,
                            opt_SMDbare = 1L,
                            silt        = 50.0,
                            BD          = 1.2,
                            OC          = 1.0,
                            minRM_Moist = 0.2) {
  
  res <- .Fortran("rothc",
                  timeFact   = as.integer(12),
                  DPM        = as.double(DPM),
                  RPM        = as.double(RPM),
                  Bio        = as.double(Bio),
                  Hum        = as.double(Hum),
                  IOM        = as.double(IOM),
                  SOC        = as.double(0),
                  total_CO2  = as.double(total_CO2),
                  DPM_Rage   = as.double(0),
                  RPM_Rage   = as.double(0),
                  Bio_Rage   = as.double(0),
                  Hum_Rage   = as.double(0),
                  Total_Rage = as.double(0),
                  modernC    = as.double(modernC),
                  clay       = as.double(clay),
                  depth      = as.double(depth),
                  TEMP       = as.double(temp),
                  RAIN       = as.double(rain),
                  PEVAP      = as.double(evap),
                  PC         = as.integer(pc),
                  PL_DPM_f   = as.double(PL_DPM_f),
                  PL_RPM_f   = as.double(PL_RPM_f),
                  OA_DPM_f   = as.double(0),
                  OA_RPM_f   = as.double(0),
                  OA_Bio_f   = as.double(0),
                  OA_Hum_f   = as.double(1),   # OA_Inp=0 so fractions irrelevant; sum must = 1
                  C_Inp      = as.double(C_inp),
                  OA_Inp     = as.double(0),
                  SMD        = as.double(smd_prev),
                  RM_TMP     = as.double(0),
                  RM_Moist   = as.double(0),
                  RM_PC      = as.double(0),
                  opt_RMmoist = as.integer(opt_RMmoist),
                  opt_SMDbare = as.integer(opt_SMDbare),
                  silt       = as.double(silt),
                  BD         = as.double(BD),
                  OC         = as.double(OC),
                  minRM_Moist = as.double(minRM_Moist)
  )
  
  list(DPM = res$DPM, RPM = res$RPM, Bio = res$Bio, Hum = res$Hum,
       IOM = res$IOM, SOC = res$SOC, total_CO2 = res$total_CO2,
       SMD = res$SMD, RM_TMP = res$RM_TMP, RM_Moist = res$RM_Moist,
       RM_PC = res$RM_PC)
}


# -----------------------------------------------------------------------------
# rothc_ref_run: run the original Fortran over a full time series
# -----------------------------------------------------------------------------
# Arguments:
#   pools0         -- named vector: DPM, RPM, Bio, Hum, IOM
#   climate_monthly -- data.frame from map_climate_monthly(site_df=):
#                      temp_air, precip, evap, soil_cover, clay, depth
#   inputs_monthly  -- data.frame from map_inputs_rothc():
#                      C_DPM, C_RPM (used to back-calculate C_inp + PL_DPM_f)
#   obs_times      -- integer vector of months at which to record output
#
# Returns:
#   data.frame with columns: DPM, RPM, Bio, Hum, IOM, total_soc
rothc_ref_run <- function(pools0, climate_monthly, inputs_monthly, obs_times) {
  
  n_months <- nrow(climate_monthly)
  if (nrow(inputs_monthly) != n_months)
    stop("climate_monthly and inputs_monthly row counts must match")
  
  clay  <- climate_monthly$clay[1]
  depth <- climate_monthly$depth[1]
  
  DPM <- pools0[["DPM"]]
  RPM <- pools0[["RPM"]]
  Bio <- pools0[["Bio"]]
  Hum <- pools0[["Hum"]]
  IOM <- pools0[["IOM"]]
  smd <- 0.0   # initialise at field capacity, as in Shell.for
  
  obs_set <- as.integer(obs_times)
  n_obs   <- length(obs_set)
  results <- matrix(0, nrow = n_obs, ncol = 5,
                    dimnames = list(NULL, c("DPM","RPM","Bio","Hum","IOM")))
  ri <- 1L
  
  for (i in seq_len(n_months)) {
    
    C_DPM  <- inputs_monthly$C_DPM[i]
    C_RPM  <- inputs_monthly$C_RPM[i]
    C_inp  <- C_DPM + C_RPM
    
    # Recover PL fractions from pre-split inputs; guard against zero input
    if (C_inp > 1e-10) {
      PL_DPM_f <- C_DPM / C_inp
      PL_RPM_f <- C_RPM / C_inp
    } else {
      PL_DPM_f <- 0.5
      PL_RPM_f <- 0.5
    }
    
    step <- .rothc_ref_step(
      DPM      = DPM,
      RPM      = RPM,
      Bio      = Bio,
      Hum      = Hum,
      IOM      = IOM,
      temp     = climate_monthly$temp_air[i],
      rain     = climate_monthly$precip[i],
      evap     = climate_monthly$evap[i],
      pc       = as.integer(climate_monthly$soil_cover[i]),
      C_inp    = C_inp,
      PL_DPM_f = PL_DPM_f,
      PL_RPM_f = PL_RPM_f,
      clay     = clay,
      depth    = depth,
      smd_prev = smd
    )
    
    DPM <- step$DPM
    RPM <- step$RPM
    Bio <- step$Bio
    Hum <- step$Hum
    smd <- step$SMD   # carry SMD to next month
    
    if (ri <= n_obs && i == obs_set[ri]) {
      results[ri, ] <- c(DPM, RPM, Bio, Hum, IOM)
      ri <- ri + 1L
    }
  }
  
  df           <- as.data.frame(results)
  df$total_soc <- df$DPM + df$RPM + df$Bio + df$Hum + df$IOM
  df
}


# =============================================================================
# Full time series comparison: wrapper vs original Fortran
# =============================================================================

result_RothC_ref <- do.call(rbind, lapply(plots, function(pid) {
  
  clim <- monthly_climate[monthly_climate$plot_id == pid, ]
  inp  <- RothC_inputs[RothC_inputs$plot_id       == pid, ]
  
  clay  <- clim$clay[1]
  depth <- clim$depth[1]
  
  
  # -----------------------------------------------------------------------------
  # KEY POINT: for this test we use the same initialization for both the wrapper and the original Fortran,
  # so that any differences in output reflect only structural/numerical differences, not initialization 
  # -----------------------------------------------------------------------------
  
  # --- INITIALIZATION PART 1: compute mean xi for steady-state, as new implementation ---
  xi      <- compute_xi_rothc(clim)
  xi_mean <- mean(xi[seq_len(min(12L, length(xi)))])
  
  # --- INITIALIZATION PART 2: compute steady-state pool sizes, as new implementation ---
   litter_ss <- c(DPM = mean(inp$C_DPM) * 12,
                 RPM = mean(inp$C_RPM) * 12)
  
  # rothc_steady_state() solves analytically for the pool sizes where
  # decomposition exactly balances litter input at the mean xi.
  # IOM is computed internally via the Falloon equation:
  #   IOM = 0.049 * TOC^1.139
  # where TOC = C_active + IOM, solved self-consistently.
  # These pools are passed to BOTH the wrapper and the original Fortran,
  # so both runs start from identical initial conditions -- which is why
  # the comparison reflects only structural/numerical differences, not
  # initialization differences.
  C_init <- rothc_steady_state(xi_mean, litter_ss, clay)
  
  cat("plot:", pid, " (ref)\n")
  
  # --- SIMULATION: run original Fortran from C_init ---
  # SMD is initialized to 0.0 inside rothc_ref_run() (field capacity),
  # matching Shell.for. It is then carried forward month by month as
  # internal state -- this is the second piece of initialization.
  run_ref <- rothc_ref_run(
    pools0          = C_init,
    climate_monthly = clim,
    inputs_monthly  = inp,
    obs_times       = seq(12L, nrow(inp), by = 12L)
  )
  
  cbind(plot_id = pid,
        data.frame(year = unique(inp$year)),
        run_ref)
}))

# =============================================================================
# Comparison plot: wrapper vs original Fortran
# =============================================================================

comparison <- rbind(
  result_RothC[,     c("plot_id","year","total_soc")] |> cbind(model = "RothC (wrapper)"),
  result_RothC_ref[, c("plot_id","year","total_soc")] |> cbind(model = "RothC (ref. Fortran)")
)

model_colors <- c(
  "RothC (wrapper)"          = "#2166ac",
  "RothC (ref. Fortran)" = "#6baed6"
)
model_linetypes <- c(
  "RothC (wrapper)"          = "solid",
  "RothC (ref. Fortran)" = "dashed"
)

p <- ggplot(comparison, aes(x = year, y = total_soc, color = model, linewidth = model)) +
  geom_line(aes(linetype = model)) +
  scale_color_manual(values = c(
    "RothC (wrapper)"          = "#d73027",
    "RothC (ref. Fortran)" = "#4575b4"
  )) +
  scale_linetype_manual(values = c(
    "RothC (wrapper)"          = "solid",
    "RothC (ref. Fortran)" = "solid"
  )) +
  scale_linewidth_manual(values = c(
    "RothC (wrapper)"          = 0.8,
    "RothC (ref. Fortran)" = 2.5
  )) +
  facet_wrap(~ plot_id, ncol = 2, scales = "free_y") +
  labs(
    title    = "RothC — total SOC over time",
    x        = "Year",
    y        = expression("Total SOC (t C ha"^{-1}*")"),
    color    = "Implementation",
    linetype = "Implementation",
    linewidth = "Implementation"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey92"),
    strip.text       = element_text(face = "bold"),
    legend.position  = "bottom",
    panel.grid.minor = element_blank()
  )
print(p)

ggsave("./Testing/rothc_comparison_and_validation.png",
       plot = p, width = 10, height = 8, dpi = 300)
