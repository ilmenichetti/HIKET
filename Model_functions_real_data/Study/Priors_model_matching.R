# =============================================================================
# Priors_model_matching.R
#
# PURPOSE:
#   Study the original Yasso15/20 posterior distributions from the model
#   repositories and derive prior specifications for the HIKET calibration.
#   Prints ready-to-paste sigma_ppm, free_defaults, and param_spec blocks
#   for each model's calibration script. No files written to disk.
#
# INPUTS (all under ./Yasso_original/):
#   Yasso15.dat   -- 10001 rows x 35 cols (row 1 = MAP; rows 2-10001 = MCMC)
#   Yasso20.dat   -- 10000 rows x 35 cols (MCMC samples)
#   y15par.csv    -- Yasso15 MAP values with original parameter names
#   ParY20.dat    -- Yasso20 MAP values, 1 row x 35 cols (original ordering)
#
# PRIOR DESIGN PHILOSOPHY:
#   Climate/size parameters: informative priors from global calibration.
#     Posterior SDs from repository give sigma_ppm in the correct space.
#     Prevents gamma/sigma_input compensation loop and delta2 NaN pathology.
#   Transfer fractions: weakly informative (sigma_ppm = 1.0).
#     Finnish data drives pool flow structure; structural comparison axis.
#   Nuisance (sigma_init, sigma_input): sigma_ppm = 0.5 (log scale).
#
# COLUMN REMAPPING:
#   Original Yasso15/20 ordering vs YASSO15_PARAM_NAMES differ at positions
#   17-35. ORIG_TO_LORENZO maps .dat columns to YASSO15_PARAM_NAMES positions.
#   Original:  17-21=w1-w5, 22-23=beta1/2, 24-25=betaN1/2, 26-27=betaH1/2,
#              28-30=gamma/gammaN/gammaH, 31-32=p_H/alpha_H, 33-35=delta1/2/r
#   Lorenzo's: 17-18=beta1/2, 19=gamma, 20-21=betaN1/2, 22=gammaN,
#              23-24=betaH1/2, 25=gammaH, 26-27=p_H/alpha_H, 28-30=delta1/2/r,
#              31-35=w1-w5
#
# TRANSFORM CHANGES REQUIRED IN CALIBRATION SCRIPTS:
#   delta2: unconstrained -> log    (enforce > 0; NaN bug for d > ~0.4cm)
#   r:      unconstrained -> log    (only -ABS(r) used in Fortran)
#   beta1:  unconstrained -> log    (enforce > 0; T increases decomp in Finland)
#   betaN1: unconstrained -> log    (same)
#   betaH1: unconstrained -> log    (same)
#   delta1 and gamma: remain unconstrained; sign enforced via tight prior centre
#
# YASSO20 STRUCTURAL ZEROS:
#   p_NW, p_AE, p_WE, p_NE, p_AN, p_EN, w1-w5 are fixed at 0 in the original
#   Yasso20 calibration (SD = 0). These must move from param_spec (free) to
#   FIXED_RATE_NAMES in run_Yasso20_calibration.R.
# =============================================================================

DAT_DIR <- "./Model_functions_real_data/Decomposition_functions/Yasso_original"

source("./Model_functions_real_data/Decomposition_functions/Yasso/yasso07_wrapper.R")
source("./Model_functions_real_data/Decomposition_functions/Yasso/yasso15_wrapper.R")
source("./Model_functions_real_data/Decomposition_functions/Yasso/yasso20_wrapper.R")


# =============================================================================
# 1.  Column remapping vector
# =============================================================================
# For each position i in YASSO15_PARAM_NAMES, ORIG_TO_LORENZO[i] is the
# column to read from the .dat file (original ordering).

ORIG_TO_LORENZO <- c(
  1:16,               # alphas + fractions: unchanged
  22, 23,             # beta1,  beta2        (orig 22-23)
  28,                 # gamma                (orig 28)
  24, 25,             # betaN1, betaN2       (orig 24-25)
  29,                 # gammaN               (orig 29)
  26, 27,             # betaH1, betaH2       (orig 26-27)
  30,                 # gammaH               (orig 30)
  31, 32,             # p_H, alpha_H         (orig 31-32)
  33, 34, 35,         # delta1, delta2, r    (orig 33-35)
  17, 18, 19, 20, 21  # w1-w5               (orig 17-21)
)
stopifnot(length(ORIG_TO_LORENZO) == 35L)
stopifnot(sort(ORIG_TO_LORENZO) == 1:35)


# =============================================================================
# 2.  Read and remap posterior samples
# =============================================================================

y15_raw     <- read.table(file.path(DAT_DIR, "Yasso15.dat"), header = FALSE)
y15_raw     <- y15_raw[-1L, ]
y15_samples <- as.data.frame(y15_raw[, ORIG_TO_LORENZO])
colnames(y15_samples) <- YASSO15_PARAM_NAMES

y20_raw     <- read.table(file.path(DAT_DIR, "Yasso20.dat"), header = FALSE)
y20_samples <- as.data.frame(y20_raw[, ORIG_TO_LORENZO])
colnames(y20_samples) <- YASSO20_PARAM_NAMES

cat(sprintf("Yasso15 posterior: %d samples x %d parameters\n",
            nrow(y15_samples), ncol(y15_samples)))
cat(sprintf("Yasso20 posterior: %d samples x %d parameters\n",
            nrow(y20_samples), ncol(y20_samples)))


# =============================================================================
# 3.  Read MAP values
# =============================================================================

# Yasso15: named CSV -- direct name lookup (avoids positional confusion)
y15_map_raw  <- read.csv(file.path(DAT_DIR, "y15par.csv"))
y15_map_orig <- setNames(y15_map_raw$value, y15_map_raw$name)

NAME_MAP_Y15 <- c(
  aA="alpha_A", aW="alpha_W", aE="alpha_E", aN="alpha_N",
  pWA="p_WA",   pEA="p_EA",   pNA="p_NA",   pAW="p_AW",
  pEW="p_EW",   pNW="p_NW",   pAE="p_AE",   pWE="p_WE",
  pNE="p_NE",   pAN="p_AN",   pWN="p_WN",   pEN="p_EN",
  w1="w1",   w2="w2",   w3="w3",   w4="w4",   w5="w5",
  b1="beta1",   b2="beta2",
  bN1="betaN1", bN2="betaN2",
  bH1="betaH1", bH2="betaH2",
  g="gamma",    gN="gammaN",  gH="gammaH",
  pH="p_H",     aH="alpha_H",
  th1="delta1", th2="delta2", r="r"
)
y15_map <- setNames(y15_map_orig[names(NAME_MAP_Y15)], NAME_MAP_Y15)

# Yasso20: single-row .dat, original ordering
y20_map_raw <- suppressWarnings(
  read.table(file.path(DAT_DIR, "ParY20.dat"), header = FALSE))
y20_map <- setNames(as.numeric(y20_map_raw[1L, ORIG_TO_LORENZO]),
                    YASSO20_PARAM_NAMES)

# Yasso07: published MAP from y07par_gui.csv (no MCMC posterior available).
# Validation: confirm wrapper defaults match published values.
# y07par_gui uses different names and contains NA rows for parameters not
# present in Yasso07 (shared file format with later model versions).
y07_pub_raw  <- read.csv(file.path(DAT_DIR, "y07par_gui.csv"))
y07_pub_orig <- setNames(y07_pub_raw$value, y07_pub_raw$name)

NAME_MAP_Y07 <- c(
  aA="alpha_A",  aW="alpha_W",  aE="alpha_E",  aN="alpha_N",
  p1WA="p_WA",   p2EA="p_EA",   p3NA="p_NA",   p4AW="p_AW",
  p5EW="p_EW",   p6NW="p_NW",   p7AE="p_AE",   p8WE="p_WE",
  p9NE="p_NE",   p10AN="p_AN",  p11WN="p_WN",  p12EN="p_EN",
  b1="beta1",    b2="beta2",    gamma="gamma",
  aH="alpha_H",  pH="p_H",
  phi1="delta1", phi2="delta2", r="r"
)
y07_map <- setNames(y07_pub_orig[names(NAME_MAP_Y07)], NAME_MAP_Y07)

# Compare against wrapper defaults (should be identical)
y07_wrapper_vals <- YASSO07_DEFAULT_PARAMS[names(y07_map)]
y07_discrepancy  <- abs(y07_wrapper_vals - y07_map) / (abs(y07_map) + 1e-12)
y07_mismatch     <- y07_discrepancy > 1e-4

cat("\n=== Yasso07 published vs wrapper defaults ===\n")
df07_check <- data.frame(
  parameter   = names(y07_map),
  published   = round(y07_map, 7),
  wrapper     = round(y07_wrapper_vals, 7),
  rel_diff_pct = round(y07_discrepancy * 100, 4),
  match       = !y07_mismatch,
  row.names   = NULL
)
print(df07_check, row.names = FALSE)
if (any(y07_mismatch))
  warning(sprintf("Yasso07: %d parameters differ between y07par_gui.csv and wrapper defaults:\n  %s",
                  sum(y07_mismatch), paste(names(y07_map)[y07_mismatch], collapse = ", ")))
if (!any(y07_mismatch))
  cat("  All values match -- wrapper defaults consistent with published parameters.\n")


# =============================================================================
# 4.  Sanity check: posterior means vs MAP
# =============================================================================

key_check <- c("beta1","gamma","gammaN","gammaH","delta1","delta2","r","w1","w3")

cat("\n=== SANITY CHECK: posterior means vs MAP ===\n")
cat("Yasso15  (Mean should match y15par.csv; MAP and Mean should agree):\n")
print(round(rbind(MAP  = y15_map[key_check],
                  Mean = colMeans(y15_samples)[key_check]), 4))

cat("\nYasso20  (MAP and Mean should be close):\n")
print(round(rbind(MAP  = y20_map[key_check],
                  Mean = colMeans(y20_samples)[key_check]), 4))

cat("\nYasso07  (no posterior; published vs wrapper defaults):\n")
key07 <- intersect(key_check, names(y07_map))
print(round(rbind(published = y07_map[key07],
                  wrapper   = y07_wrapper_vals[key07]), 4))


# =============================================================================
# 5.  Posterior summary
# =============================================================================

post_summary <- function(samples, map_vals, model_name) {
  mn    <- colMeans(samples)
  sds   <- apply(samples, 2L, sd)
  fixed <- sds < 1e-8
  
  df <- data.frame(
    parameter = names(mn),
    map       = round(map_vals[names(mn)], 5),
    mean      = round(mn, 5),
    sd        = round(sds, 5),
    fixed     = fixed,
    row.names = NULL
  )
  cat(sprintf("\n=== %s posterior summary ===\n", model_name))
  print(df, row.names = FALSE)
  
  if (any(fixed))
    cat(sprintf(
      "\n  %s structurally fixed (SD~0):  %s\n",
      model_name, paste(names(mn)[fixed], collapse = ", ")))
  
  invisible(df)
}

s15 <- post_summary(y15_samples, y15_map, "Yasso15")
s20 <- post_summary(y20_samples, y20_map, "Yasso20")


# =============================================================================
# 6.  Sigma_ppm computation
# =============================================================================
# Log-transformed parameters (beta1, betaN1, betaH1, delta2, r):
#   sigma_ppm = SD(log(samples)) -- correct for unconstrained log space
# Unconstrained parameters (gamma, gammaN, gammaH, beta2, etc., delta1):
#   sigma_ppm = SD(samples) in physical space
# gammaH cap at 1.5: parameter weakly identified at Finnish P (~600mm) because
#   [1-exp(gammaH*0.6)] ~ 1 for any gammaH < -5. Large SD is an artefact of
#   near-unidentifiability, not genuine uncertainty. Yasso20 better constrained
#   (SD=1.39); use that as the cap for both models.

GAMMATH_SD_CAP <- 1.5
FRACTION_SIGMA <- 1.0
LOG_PARAMS     <- c("beta1","betaN1","betaH1","delta2","r")

CLIMATE_SIZE <- c("beta1","beta2","gamma",
                  "betaN1","betaN2","gammaN",
                  "betaH1","betaH2","gammaH",
                  "delta1","delta2","r")

compute_sigma_ppm <- function(samples, map_vals, model_name) {
  sds   <- apply(samples, 2L, sd)
  sigma <- setNames(rep(FRACTION_SIGMA, ncol(samples)), colnames(samples))
  
  for (p in intersect(CLIMATE_SIZE, colnames(samples))) {
    if (sds[p] < 1e-8) {
      sigma[p] <- 0.10
      next
    }
    if (p %in% LOG_PARAMS) {
      x <- samples[[p]]
      x <- x[x > 0 & is.finite(x)]
      sigma[p] <- if (length(x) > 1L) sd(log(x)) else 0.10
    } else {
      sigma[p] <- sds[p]
    }
  }
  
  if ("gammaH" %in% names(sigma))
    sigma["gammaH"] <- min(sigma["gammaH"], GAMMATH_SD_CAP)
  
  sigma["sigma_init"]  <- 0.5
  sigma["sigma_input"] <- 0.5
  
  cat(sprintf("\n=== %s sigma_ppm ===\n", model_name))
  df <- data.frame(
    parameter = names(sigma),
    post_sd   = round(sds[names(sigma)], 5),
    sigma_ppm = round(sigma, 5),
    space     = ifelse(names(sigma) %in% LOG_PARAMS,      "log",
                       ifelse(names(sigma) %in% CLIMATE_SIZE,    "physical",
                              ifelse(names(sigma) %in% c("sigma_init",
                                                         "sigma_input"), "log", "—"))),
    note      = ifelse(names(sigma) %in% CLIMATE_SIZE,    "from posterior",
                       ifelse(names(sigma) %in% c("sigma_init",
                                                  "sigma_input"), "fixed",
                              "wide (fractions)")),
    row.names = NULL
  )
  print(df, row.names = FALSE)
  invisible(sigma)
}

sigma15 <- compute_sigma_ppm(y15_samples, y15_map, "Yasso15")
sigma20 <- compute_sigma_ppm(y20_samples, y20_map, "Yasso20")


# =============================================================================
# 7.  Free defaults from posterior means
# =============================================================================
# Prior centres for climate/size params set to posterior mean; fractions
# keep model defaults. Values in physical space; to_unconstrained() in the
# calibration script converts to unconstrained space for best_x.

print_defaults <- function(samples, map_vals, model_name,
                           default_params, param_names) {
  mn       <- colMeans(samples)
  defaults <- default_params[param_names]
  
  for (p in intersect(CLIMATE_SIZE, param_names))
    defaults[p] <- mn[p]
  
  if ("sigma_init"  %in% names(defaults)) defaults["sigma_init"]  <- 1.0
  if ("sigma_input" %in% names(defaults)) defaults["sigma_input"] <- 1.0
  
  cat(sprintf("\n=== %s free_defaults (physical space, update in calibration script) ===\n",
              model_name))
  df <- data.frame(
    parameter     = names(defaults),
    model_default = round(default_params[names(defaults)], 5),
    prior_centre  = round(defaults, 5),
    source        = ifelse(names(defaults) %in% CLIMATE_SIZE,
                           "posterior mean", "model default"),
    row.names = NULL
  )
  print(df, row.names = FALSE)
  invisible(defaults)
}

defaults15 <- print_defaults(y15_samples, y15_map, "Yasso15",
                             YASSO15_DEFAULT_PARAMS, YASSO15_PARAM_NAMES)
defaults20 <- print_defaults(y20_samples, y20_map, "Yasso20",
                             YASSO20_DEFAULT_PARAMS, YASSO20_PARAM_NAMES)

# Yasso07: no posterior; free_defaults = published values (wrapper defaults)
cat("
=== Yasso07 free_defaults (physical space) ===
")
y07_free_names <- c(
  "p_WA","p_EA","p_NA","p_AW","p_EW","p_NW",
  "p_AE","p_WE","p_NE","p_AN","p_WN","p_EN",
  "beta1","beta2","gamma","delta1","delta2","r"
)
df07_def <- data.frame(
  parameter    = y07_free_names,
  prior_centre = round(YASSO07_DEFAULT_PARAMS[y07_free_names], 6),
  source       = ifelse(y07_free_names %in% CLIMATE_SIZE,
                        "published MAP (Tuomi 2008-2010)", "published MAP"),
  row.names = NULL
)
print(df07_def, row.names = FALSE)


# =============================================================================
# 8.  Yasso20 fixed parameters
# =============================================================================

y20_fixed <- names(which(apply(y20_samples, 2L, sd) < 1e-8))

cat("\n=== YASSO20 STRUCTURAL ZEROS ===\n")
cat("Move these from param_spec to FIXED_RATE_NAMES in run_Yasso20_calibration.R:\n")
cat(sprintf("  c(%s)\n", paste(sprintf('"%s"', y20_fixed), collapse = ", ")))
cat("Fixed values (all 0.0):\n")
print(round(y20_map[y20_fixed], 6))


# =============================================================================
# 9.  Transform changes required
# =============================================================================

cat("
=== PARAM_SPEC CHANGES (all three calibration scripts) ===

Replace the current grouped entries for delta/size and climate with
individual entries as shown below.

  # --- Size modifier (all models) ---
  list(names = 'delta1', type = 'unconstrained')  # can be negative
  list(names = 'delta2', type = 'log')             # enforce > 0; NaN fix
  list(names = 'r',      type = 'log')             # only -ABS(r) used

  # --- Climate: Yasso07 ---
  list(names = 'beta1',  type = 'log')
  list(names = 'beta2',  type = 'unconstrained')
  list(names = 'gamma',  type = 'unconstrained')   # < 0 enforced by tight prior

  # --- Climate: Yasso15/20 (replace the three existing group entries) ---
  list(names = 'beta1',   type = 'log')
  list(names = 'beta2',   type = 'unconstrained')
  list(names = 'gamma',   type = 'unconstrained')
  list(names = 'betaN1',  type = 'log')
  list(names = 'betaN2',  type = 'unconstrained')
  list(names = 'gammaN',  type = 'unconstrained')
  list(names = 'betaH1',  type = 'log')
  list(names = 'betaH2',  type = 'unconstrained')
  list(names = 'gammaH',  type = 'unconstrained')
")


# =============================================================================
# 10.  Ready-to-paste sigma_ppm blocks
# =============================================================================

cat_sigma_block <- function(sigma, model_name) {
  params_in_model <- intersect(CLIMATE_SIZE, names(sigma))
  cat(sprintf(
    "\n# ---- sigma_ppm for %s (paste into calibration script) ----\n",
    model_name))
  cat("sigma_ppm <- setNames(rep(1.0, N_FREE), FREE_NAMES)\n")
  for (p in params_in_model)
    cat(sprintf('sigma_ppm["%s"] <- %.5f\n', p, sigma[p]))
  cat('sigma_ppm["sigma_init"]  <- 0.5\n')
  cat('sigma_ppm["sigma_input"] <- 0.5\n')
  cat("# Transfer fractions remain at 1.0\n")
}

cat("\n=== READY-TO-PASTE SIGMA_PPM BLOCKS ===\n")
cat_sigma_block(sigma15, "Yasso15")
cat_sigma_block(sigma20, "Yasso20")

cat('
# ---- sigma_ppm for Yasso07 (paste into calibration script) ----
# Climate: Tuomi et al. 2008/2009; size: Tuomi et al. 2010 Table 4
sigma_ppm <- setNames(rep(1.0, N_FREE), FREE_NAMES)
sigma_ppm["beta1"]  <- 0.20    # log space; exp(0.2*T) at T~10C gives 2x range
sigma_ppm["beta2"]  <- 0.05    # physical space; quadratic, tight
sigma_ppm["gamma"]  <- 0.30    # physical space; keeps modifier in [0.3,0.6] at Finnish P
sigma_ppm["delta1"] <- 0.15    # physical space; Tuomi 2010 CI +/-0.16
sigma_ppm["delta2"] <- 0.10    # log space;      Tuomi 2010 CI +/-0.10
sigma_ppm["r"]      <- 0.015   # log space;      Tuomi 2010 CI +/-0.013
sigma_ppm["sigma_init"]  <- 0.5
sigma_ppm["sigma_input"] <- 0.5
# Transfer fractions remain at 1.0
')

message("\nPriors_model_matching.R complete -- no files written to disk.")