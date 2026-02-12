# =============================================================================
# YASSO20 — explicit-parameter version (Yasso07-style), with validations
# =============================================================================
# Goals of this refactor:
#  1) Replace theta[1:35] indexing with named parameters (params$alpha_A, etc.)
#  2) Add strong input + parameter validation to fail fast with clear errors
#  3) Make the code readable for debugging (explanatory comments)
#
# Model structure:
#  Pools: A, W, E, N, H  (AWEN + Humus)
#  Dynamics (continuous within year): dC/dt = A * C + b
#  Annual step uses exact matrix exponential solution.
#
# Litter size classes:
#  nwl (diameter 0), fwl (2), cwl (15)
#  Each is simulated separately and summed at the end.
# =============================================================================

# -----------------------------------------------------------------------------
# Default parameters (explicit names)
# -----------------------------------------------------------------------------
YASSO20_DEFAULT_PARAMS <- list(
  # Baseline decomposition coefficients (positive rates; climate scales them)
  alpha_A = 0.51,
  alpha_W = 5.19,
  alpha_E = 0.13,
  alpha_N = 0.10,
  alpha_H = 0.0015,
  
  # Transfer coefficients among AWEN (pattern matrix A_p)
  # Convention: p_XY = fraction from source Y transferred into destination X
  # (Rows = destination, columns = source)
  p_WA = 0.50,  # W -> A
  p_EA = 0.00,  # E -> A
  p_NA = 1.00,  # N -> A
  
  p_AW = 1.00,  # A -> W
  p_EW = 0.99,  # E -> W
  p_NW = 0.00,  # N -> W
  
  p_AE = 0.00,  # A -> E
  p_WE = 0.00,  # W -> E
  p_NE = 0.00,  # N -> E
  
  p_AN = 0.00,  # A -> N
  p_WN = 0.163,  # W -> N
  p_EN = 0.00,  # E -> N
  
  # Temperature response (quadratic in T, averaged across months)
  beta1_AWE = 0.158,
  beta2_AWE = -0.002,
  beta1_N   = 0.17,
  beta2_N   = -0.005,
  beta1_H   = 0.067,
  beta2_H   = 0.00,
  
  # Precipitation response (damps/enables decomposition)
  gamma_AWE = -1.44,
  gamma_N   = -2.00,
  gamma_H   = -6.90,
  
  # Humus formation: transfer from A/W/E/N into H
  p_H = 0.0042,
  
  # Wood size modifier parameters
  phi1 = -2.55,
  phi2 = 1.24,
  r    = 0.25
)

# =============================================================================
# VALIDATION HELPERS
# =============================================================================

# -----------------------------------------------------------------------------
# Validate parameter list
# -----------------------------------------------------------------------------
# What we check (fail fast with informative errors):
# - all required names present
# - all values are finite numerics
# - rates (alpha_*) are non-negative
# - transfer coefficients are non-negative
# - humus transfer p_H is non-negative
# - precipitation modifiers are safe for typical P (won't create NA/Inf)
# - size modifier base stays positive for the diameters you will use
# -----------------------------------------------------------------------------
validate_yasso20_params <- function(params, diameters_to_check = c(0, 2, 15)) {
  
  if (!is.list(params)) stop("params must be a list", call. = FALSE)
  
  required <- c(
    # rates
    "alpha_A","alpha_W","alpha_E","alpha_N","alpha_H",
    # transfers (AWEN)
    "p_WA","p_EA","p_NA","p_AW","p_EW","p_NW","p_AE","p_WE","p_NE","p_AN","p_WN","p_EN",
    # temp response
    "beta1_AWE","beta2_AWE","beta1_N","beta2_N","beta1_H","beta2_H",
    # precip response
    "gamma_AWE","gamma_N","gamma_H",
    # humus transfer + size
    "p_H","phi1","phi2","r"
  )
  
  missing <- setdiff(required, names(params))
  if (length(missing) > 0) {
    stop("Yasso20 params missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  
  # Ensure all required are numeric & finite
  bad_type <- required[!vapply(params[required], function(x) is.numeric(x) && length(x) == 1, logical(1))]
  if (length(bad_type) > 0) {
    stop("Yasso20 params must be numeric scalars: ", paste(bad_type, collapse = ", "), call. = FALSE)
  }
  
  bad_finite <- required[!vapply(params[required], function(x) is.finite(x), logical(1))]
  if (length(bad_finite) > 0) {
    stop("Yasso20 params contain NA/Inf: ", paste(bad_finite, collapse = ", "), call. = FALSE)
  }
  
  # Rates should be >= 0 (we also abs() them later, but explicit is better)
  rate_names <- c("alpha_A","alpha_W","alpha_E","alpha_N","alpha_H")
  if (any(vapply(params[rate_names], function(x) x < 0, logical(1)))) {
    stop("Yasso20 alpha_* rates should be >= 0 (got negative).", call. = FALSE)
  }
  
  # Transfer coefficients should be >= 0
  transfer_names <- c("p_WA","p_EA","p_NA","p_AW","p_EW","p_NW","p_AE","p_WE","p_NE","p_AN","p_WN","p_EN","p_H")
  if (any(vapply(params[transfer_names], function(x) x < 0, logical(1)))) {
    stop("Yasso20 transfer coefficients p_* should be >= 0.", call. = FALSE)
  }
  
  # Optional sanity: precipitation modifier expression should be finite for typical P
  # (We don't force sign of gamma; your defaults are negative and that is OK.)
  P_test <- c(0, 300, 600, 1200, 2000)
  gnames <- c("gamma_AWE","gamma_N","gamma_H")
  for (gn in gnames) {
    vals <- 1 - exp(params[[gn]] * P_test / 1000)
    if (any(!is.finite(vals))) {
      stop("Precip modifier produced non-finite values for ", gn, ". Check gamma.", call. = FALSE)
    }
  }
  
  # Size modifier base must be positive for diameters used (if diameter > 0)
  for (d in diameters_to_check) {
    if (d > 0) {
      base <- 1 + params$phi1*d + params$phi2*d^2
      if (!is.finite(base) || base <= 0) {
        stop(sprintf("Size modifier base <= 0 for diameter %.3f (base=%.3f). Check phi1/phi2.", d, base),
             call. = FALSE)
      }
    }
  }
  
  invisible(TRUE)
}

# -----------------------------------------------------------------------------
# Validate input_df for yasso20_run()
# -----------------------------------------------------------------------------
# input_df is ANNUAL rows, but must contain:
# - precip (annual)
# - T_1..T_12 (monthly mean temps for each year)
# - litter columns: nwl_*, fwl_*, cwl_* (annual AWEN inputs)
# -----------------------------------------------------------------------------
validate_yasso20_input_df <- function(input_df) {
  if (!is.data.frame(input_df)) stop("input_df must be a data.frame", call. = FALSE)
  
  req <- c("precip", paste0("T_", 1:12))
  missing <- setdiff(req, names(input_df))
  if (length(missing) > 0) {
    stop("Yasso20 input_df missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  
  if (any(!is.finite(input_df$precip))) stop("input_df$precip contains NA/Inf.", call. = FALSE)
  
  Tmat <- as.matrix(input_df[, paste0("T_", 1:12)])
  if (any(!is.finite(Tmat))) stop("input_df contains NA/Inf in T_1..T_12.", call. = FALSE)
  
  # Check litter columns exist for at least one litter type
  litter_present <- any(grepl("^(nwl|fwl|cwl)_[AWEN]$", names(input_df)))
  if (!litter_present) {
    stop("input_df has no litter columns like nwl_A..nwl_N / fwl_* / cwl_*.", call. = FALSE)
  }
  
  # Optional: warn if any negative litter inputs
  lit_cols <- grep("^(nwl|fwl|cwl)_[AWEN]$", names(input_df), value = TRUE)
  if (any(input_df[, lit_cols, drop = FALSE] < 0, na.rm = TRUE)) {
    warning("Some litter inputs are negative. This is unusual; check mapping/units.")
  }
  
  invisible(TRUE)
}

# =============================================================================
# CLIMATE MODIFIERS
# =============================================================================

# -----------------------------------------------------------------------------
# Temperature modifier
# -----------------------------------------------------------------------------
# Computes mean over months of exp(beta1*T + beta2*T^2)
# This is nonlinear, so temperature seasonality matters.
yasso20_temp_mod <- function(monthly_temp, params) {
  tem  <- mean(exp(params$beta1_AWE*monthly_temp + params$beta2_AWE*monthly_temp^2))
  temN <- mean(exp(params$beta1_N  *monthly_temp + params$beta2_N  *monthly_temp^2))
  temH <- mean(exp(params$beta1_H  *monthly_temp + params$beta2_H  *monthly_temp^2))
  list(tem = tem, temN = temN, temH = temH)
}

# -----------------------------------------------------------------------------
# Precipitation modifier
# -----------------------------------------------------------------------------
# Applies factor (1 - exp(gamma * P/1000)) to each temperature modifier.
yasso20_precip_mod <- function(P, params, tem, temN, temH) {
  tem  <- tem  * (1 - exp(params$gamma_AWE * P/1000))
  temN <- temN * (1 - exp(params$gamma_N   * P/1000))
  temH <- temH * (1 - exp(params$gamma_H   * P/1000))
  list(tem = tem, temN = temN, temH = temH)
}

# -----------------------------------------------------------------------------
# Wood size modifier
# -----------------------------------------------------------------------------
# Returns scalar in (0,1], reducing rates for larger diameter.
yasso20_size_mod <- function(d, params) {
  if (is.null(d) || d <= 0) return(1)
  base <- 1 + params$phi1*d + params$phi2*d^2
  if (!is.finite(base) || base <= 0) {
    stop(sprintf("Size modifier base <= 0 (base=%.3f, d=%.3f). Check phi1/phi2/r.", base, d),
         call. = FALSE)
  }
  min(1, base^(-abs(params$r)))
}

# =============================================================================
# BUILD SYSTEM MATRIX
# =============================================================================

# -----------------------------------------------------------------------------
# Build A matrix for dC/dt = A*C + b
# -----------------------------------------------------------------------------
# Structure:
# - A_p sets the transfer pattern among pools (unitless)
# - k sets the actual decomposition rates (1/time)
# - A = A_p %*% k
yasso20_build_matrix <- function(params, tem, temN, temH, size_dep, leac = 0, precip = 0) {
  
  # Baseline decomposition rates (positive; climate scales them)
  alpha <- abs(c(params$alpha_A, params$alpha_W, params$alpha_E, params$alpha_N, params$alpha_H))
  
  # Pattern matrix (rows=destination, cols=source; pools ordered A,W,E,N,H)
  A_p <- matrix(c(
    -1,          params$p_WA, params$p_EA, params$p_NA, 0,
    params$p_AW, -1,          params$p_EW, params$p_NW, 0,
    params$p_AE, params$p_WE, -1,          params$p_NE, 0,
    params$p_AN, params$p_WN, params$p_EN, -1,          0,
    params$p_H,  params$p_H,  params$p_H,  params$p_H, -1
  ), 5, 5, byrow = TRUE)
  
  # Rate diagonal: climate + size scaling
  k <- diag(c(
    tem  * alpha[1] * size_dep,
    tem  * alpha[2] * size_dep,
    tem  * alpha[3] * size_dep,
    temN * alpha[4] * size_dep,
    temH * alpha[5]           # humus not size-scaled
  ))
  
  A <- A_p %*% k
  
  # Leaching (not used in your current calls, but left here)
  # NOTE: check sign convention if you activate leaching.
  for (i in 1:4) {
    A[i, i] <- A[i, i] + leac * precip / 1000
  }
  
  A
}

# =============================================================================
# ANNUAL STEP (exact integration)
# =============================================================================

yasso20_annual_step <- function(pools, input_awen, params, monthly_temp, precip,
                                diameter = 0, leac = 0, ss = FALSE) {
  
  # ---- Validate inputs for this step (fail fast) ----
  if (length(pools) != 5 || any(!is.finite(pools))) stop("`pools` must be length-5 finite numeric.", call. = FALSE)
  if (length(input_awen) != 5 || any(!is.finite(input_awen))) stop("`input_awen` must be length-5 finite numeric.", call. = FALSE)
  if (length(monthly_temp) != 12 || any(!is.finite(monthly_temp))) stop("`monthly_temp` must be length-12 finite numeric.", call. = FALSE)
  if (!is.finite(precip)) stop("`precip` must be finite numeric.", call. = FALSE)
  
  # ---- Climate modifiers ----
  tmods <- yasso20_temp_mod(monthly_temp, params)
  pmods <- yasso20_precip_mod(precip, params, tmods$tem, tmods$temN, tmods$temH)
  
  # Guard against NA/Inf modifiers (prevents cryptic if() errors)
  if (!is.finite(pmods$tem) || !is.finite(pmods$temN) || !is.finite(pmods$temH)) {
    stop("Non-finite climate modifier(s) (tem/temN/temH). Check precip + temps.", call. = FALSE)
  }
  
  # ---- Size modifier ----
  size_dep <- yasso20_size_mod(diameter, params)
  
  # ---- Build annual A matrix ----
  A <- yasso20_build_matrix(params, pmods$tem, pmods$temN, pmods$temH, size_dep, leac, precip)
  
  b <- input_awen
  
  # If decomposition is effectively off, just add input (numerical safeguard)
  if (pmods$tem <= 1e-16) {
    new_pools <- pools + b
    resp <- 0
    return(list(pools = new_pools, respiration = resp))
  }
  
  if (ss) {
    # Steady state for this year's forcing: solve(-A * C = b)
    new_pools <- as.vector(solve(-A, b))
  } else {
    # Exact solution over one year for linear ODE with constant A and b
    expA <- as.matrix(Matrix::expm(A))
    Ainv <- solve(A)
    new_pools <- as.vector(expA %*% pools + (expA - diag(5)) %*% Ainv %*% b)
  }
  
  # Respiration by mass balance: old + input - new
  resp <- sum(pools) + sum(b) - sum(new_pools)
  
  list(
    pools = pmax(0, new_pools),
    respiration = max(0, resp)
  )
}

# =============================================================================
# STEADY STATE (analytical)
# =============================================================================
yasso20_steady_state <- function(params, input_awen, monthly_temp, precip, diameter = 0, leac = 0) {
  
  if (length(input_awen) != 5 || any(!is.finite(input_awen))) stop("`input_awen` must be length-5 finite numeric.", call. = FALSE)
  if (length(monthly_temp) != 12 || any(!is.finite(monthly_temp))) stop("`monthly_temp` must be length-12 finite numeric.", call. = FALSE)
  if (!is.finite(precip)) stop("`precip` must be finite numeric.", call. = FALSE)
  
  tmods <- yasso20_temp_mod(monthly_temp, params)
  pmods <- yasso20_precip_mod(precip, params, tmods$tem, tmods$temN, tmods$temH)
  size_dep <- yasso20_size_mod(diameter, params)
  
  if (!is.finite(pmods$tem) || !is.finite(pmods$temN) || !is.finite(pmods$temH)) {
    stop("Non-finite climate modifier(s) in steady state.", call. = FALSE)
  }
  
  A <- yasso20_build_matrix(params, pmods$tem, pmods$temN, pmods$temH, size_dep, leac, precip)
  as.vector(solve(-A, input_awen))
}

# =============================================================================
# MAIN RUN (multi-litter)
# =============================================================================
yasso20_run <- function(params = YASSO20_DEFAULT_PARAMS, C0 = NULL, input_df, spinup = TRUE) {
  
  if (!is.list(params)) params <- as.list(params)
  
  # Validate params + input at the top (so errors are near the call site)
  validate_yasso20_params(params)
  validate_yasso20_input_df(input_df)
  
  n_years <- nrow(input_df)
  
  # Fixed diameters for the 3 litter types (same choices as Yasso07)
  wood_sizes <- list(nwl = 0, fwl = 2, cwl = 15)
  litter_types <- c("nwl", "fwl", "cwl")
  
  # Storage for total pools across litter types
  pools_mat <- matrix(0, n_years, 5, dimnames = list(NULL, c("A","W","E","N","H")))
  respiration <- numeric(n_years)
  
  # Extract climate series once
  monthly_temp <- as.matrix(input_df[, paste0("T_", 1:12)])
  precip <- input_df$precip
  
  for (lt in litter_types) {
    
    # Skip litter types not present in input_df
    if (!(paste0(lt, "_A") %in% names(input_df))) next
    
    # Annual inputs for this litter type (vectors length n_years)
    input_A <- input_df[[paste0(lt, "_A")]]
    input_W <- input_df[[paste0(lt, "_W")]]
    input_E <- input_df[[paste0(lt, "_E")]]
    input_N <- input_df[[paste0(lt, "_N")]]
    
    # ---- Initialize pools for this litter type ----
    # If C0 not provided, initialize at steady state for mean forcing.
    if (is.null(C0) || !lt %in% names(C0)) {
      
      avg_input <- c(mean(input_A), mean(input_W), mean(input_E), mean(input_N), 0)
      
      C0_lt <- yasso20_steady_state(
        params,
        avg_input,
        colMeans(monthly_temp),
        mean(precip),
        diameter = wood_sizes[[lt]],
        leac = 0
      )
      
    } else {
      C0_lt <- C0[[lt]]
      if (length(C0_lt) != 5 || any(!is.finite(C0_lt))) {
        stop("C0 for litter type ", lt, " must be a finite length-5 vector.", call. = FALSE)
      }
    }
    
    pools_lt <- C0_lt
    
    # Store pools for this litter type
    lt_pools_mat <- matrix(0, n_years, 5)
    lt_pools_mat[1, ] <- pools_lt
    
    # ---- Time loop over years ----
    # Note: updates t -> t+1, so we run 1:(n_years-1)
    for (t in 1:(n_years - 1)) {
      
      b <- c(input_A[t], input_W[t], input_E[t], input_N[t], 0)
      
      res <- yasso20_annual_step(
        pools = pools_lt,
        input_awen = b,
        params = params,
        monthly_temp = monthly_temp[t, ],
        precip = precip[t],
        diameter = wood_sizes[[lt]],
        leac = 0,
        ss = FALSE
      )
      
      pools_lt <- res$pools
      lt_pools_mat[t + 1, ] <- pools_lt
      
      # Add respiration for this litter type into total
      respiration[t + 1] <- respiration[t + 1] + res$respiration
    }
    
    # Add this litter type pools into overall pools
    pools_mat <- pools_mat + lt_pools_mat
  }
  
  total_soc <- rowSums(pools_mat)
  
  list(
    pools = pools_mat,
    total_soc = total_soc,
    respiration = respiration,
    pool_names = c("A","W","E","N","H")
  )
}
