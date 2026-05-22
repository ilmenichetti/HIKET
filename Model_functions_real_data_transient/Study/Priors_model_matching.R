# =============================================================================
# Priors_model_matching.R
#
# PURPOSE:
#   Derive and validate prior specifications for all three Yasso models in the
#   HIKET calibration. Prints ready-to-paste sigma_ppm, free_defaults, and
#   param_spec blocks for each calibration script. No files written to disk.
#
# PRIMARY MAP SOURCES:
#   Yasso20: Yasso20_sample_parameters.rda (Ryassofortran / FMI)
#            -- preferred over ParY20.dat; FMI-validated point estimate.
#   Yasso15: y15par.csv (named CSV from model repository)
#   Yasso07: y07par_gui.csv (published MAP; Tuomi et al. 2008-2010)
#            Wrapper defaults cross-checked against this file.
#
# SIGMA_PPM SOURCE:
#   Yasso15/20: full MCMC posterior from Yasso15.dat / Yasso20.dat
#               (SD in log/physical space depending on transform type)
#   Yasso07: literature CIs (Tuomi 2008-2010); no posterior available
#
# COLUMN REMAPPING:
#   Original Yasso15/20 .dat column ordering differs from YASSO*_PARAM_NAMES.
#   ORIG_TO_LORENZO maps .dat columns to our parameter name ordering.
#   Original:  17-21=w1-w5, 22-23=beta1/2, 24-25=betaN1/2, 26-27=betaH1/2,
#              28-30=gamma/gammaN/gammaH, 31-32=p_H/alpha_H, 33-35=delta1/2/r
#   HIKET:     17-18=beta1/2, 19=gamma, 20-21=betaN1/2, 22=gammaN,
#              23-24=betaH1/2, 25=gammaH, 26-27=p_H/alpha_H, 28-30=delta1/2/r,
#              31-35=w1-w5
#
# KEY DESIGN DECISIONS:
#   * Climate/size params: informative priors from posterior SDs (Yasso15/20)
#     or literature CIs (Yasso07). Prevents gamma/sigma_input compensation loop
#     and delta2 NaN pathology for d > ~0.4cm.
#   * Transfer fractions: sigma_ppm = 1.0 (weakly informative). Finnish data
#     should drive pool flow structure -- the primary structural comparison axis.
#   * Log transforms enforced for: beta1, betaN1, betaH1, delta2, r.
#     Yasso07 r stored negative in wrapper -- abs() fix required.
#   * gammaH SD capped at 1.5: near-unidentifiable at Finnish P (~600mm),
#     large SD is artefact of flat likelihood surface, not genuine uncertainty.
#   * Yasso20 structural zeros: p_NW, p_AE, p_WE, p_NE, p_AN, p_EN, w1-w5
#     have SD=0 in posterior -- must be fixed, not calibrated.
# =============================================================================

# When sourced from a calibration script to set DEFAULT_PARAMS only, set this
# flag first to skip the slow posterior-reading and sigma_ppm sections:
#   DEFAULTS_ONLY <- TRUE; source("path/to/Priors_model_matching.R")
# Running the script standalone (Rscript or interactive) runs everything.
DEFAULTS_ONLY <- exists("DEFAULTS_ONLY") && isTRUE(DEFAULTS_ONLY)

DAT_DIR  <- "./Model_functions_real_data/Decomposition_functions/Yasso_original"
RDA_FILE <- file.path(DAT_DIR, "Yasso20_sample_parameters.rda")

if (!DEFAULTS_ONLY) {
  source("./Model_functions_real_data/Decomposition_functions/Yasso/yasso07_wrapper.R")
  source("./Model_functions_real_data/Decomposition_functions/Yasso/yasso15_wrapper.R")
  source("./Model_functions_real_data/Decomposition_functions/Yasso/yasso20_wrapper.R")
}


# =============================================================================
# 1.  Column remapping vector (.dat files -> HIKET parameter ordering)
# =============================================================================

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

# Short-name -> full HIKET name mapping (shared across all Yasso models)
NAME_MAP <- c(
  aA="alpha_A",  aW="alpha_W",  aE="alpha_E",  aN="alpha_N",
  pWA="p_WA",    pEA="p_EA",    pNA="p_NA",    pAW="p_AW",
  pEW="p_EW",    pNW="p_NW",    pAE="p_AE",    pWE="p_WE",
  pNE="p_NE",    pAN="p_AN",    pWN="p_WN",    pEN="p_EN",
  w1="w1",       w2="w2",       w3="w3",       w4="w4",   w5="w5",
  b1="beta1",    b2="beta2",
  bN1="betaN1",  bN2="betaN2",
  bH1="betaH1",  bH2="betaH2",
  g="gamma",     gN="gammaN",   gH="gammaH",
  pH="p_H",      aH="alpha_H",
  th1="delta1",  th2="delta2",  r="r"
)


# =============================================================================
# 2.  Read posterior samples (.dat files -- needed for sigma_ppm computation)
# =============================================================================
# Skipped when DEFAULTS_ONLY = TRUE (sigma_ppm not needed, only MAP values).

if (!DEFAULTS_ONLY) {
  y15_raw     <- read.table(file.path(DAT_DIR, "Yasso15.dat"), header = FALSE)
  y15_raw     <- y15_raw[-1L, ]           # row 1 = MAP; rows 2-10001 = MCMC
  y15_samples <- as.data.frame(y15_raw[, ORIG_TO_LORENZO])
  colnames(y15_samples) <- YASSO15_PARAM_NAMES

  y20_raw     <- read.table(file.path(DAT_DIR, "Yasso20.dat"), header = FALSE)
  y20_samples <- as.data.frame(y20_raw[, ORIG_TO_LORENZO])
  colnames(y20_samples) <- YASSO20_PARAM_NAMES

  cat(sprintf("Yasso15 posterior: %d samples x %d parameters\n",
              nrow(y15_samples), ncol(y15_samples)))
  cat(sprintf("Yasso20 posterior: %d samples x %d parameters\n",
              nrow(y20_samples), ncol(y20_samples)))
}


# =============================================================================
# 3.  MAP values
# =============================================================================

# --- Yasso20: primary source = Ryassofortran .rda (FMI-validated) ---
local({
  e <- new.env()
  load(RDA_FILE, envir = e)
  sp <<- e$sample_parameters
})
y20_map_rda <- setNames(sp[names(NAME_MAP)], NAME_MAP)[YASSO20_PARAM_NAMES]

# Secondary: ParY20.dat (cross-validation only)
y20_map_dat <- tryCatch({
  raw <- suppressWarnings(read.table(file.path(DAT_DIR, "ParY20.dat"),
                                     header = FALSE))
  setNames(as.numeric(raw[1L, ORIG_TO_LORENZO]), YASSO20_PARAM_NAMES)
}, error = function(e) {
  message("ParY20.dat not found -- skipping cross-validation")
  NULL
})

# Cross-validate rda vs dat for key climate/size params
if (!is.null(y20_map_dat)) {
  key_params <- c("beta1","beta2","gamma","betaN1","betaN2","gammaN",
                  "betaH1","betaH2","gammaH","delta1","delta2","r",
                  "p_H","alpha_H")
  cat("\n=== Yasso20 MAP: rda vs ParY20.dat ===\n")
  cmp <- data.frame(
    parameter = key_params,
    rda       = round(y20_map_rda[key_params], 5),
    dat       = round(y20_map_dat[key_params], 5),
    rel_diff  = round(abs(y20_map_rda[key_params] - y20_map_dat[key_params]) /
                        (abs(y20_map_dat[key_params]) + 1e-12) * 100, 2),
    row.names = NULL
  )
  print(cmp, row.names = FALSE)
  cat("  (rel_diff in %; rda used as primary source)\n")
}

# Use rda as the authoritative MAP for Yasso20
y20_map <- y20_map_rda

# --- Yasso15: named CSV ---
y15_map_raw  <- read.csv(file.path(DAT_DIR, "y15par.csv"))
y15_map_orig <- setNames(y15_map_raw$value, y15_map_raw$name)
y15_map      <- setNames(y15_map_orig[names(NAME_MAP)], NAME_MAP)[YASSO15_PARAM_NAMES]

# --- Yasso07: published MAP (y07par_gui.csv) ---
y07_pub_raw  <- read.csv(file.path(DAT_DIR, "y07par_gui.csv"))
y07_pub_orig <- setNames(y07_pub_raw$value, y07_pub_raw$name)

NAME_MAP_Y07 <- c(
  aA="alpha_A",   aW="alpha_W",   aE="alpha_E",   aN="alpha_N",
  p1WA="p_WA",    p2EA="p_EA",    p3NA="p_NA",    p4AW="p_AW",
  p5EW="p_EW",    p6NW="p_NW",    p7AE="p_AE",    p8WE="p_WE",
  p9NE="p_NE",    p10AN="p_AN",   p11WN="p_WN",   p12EN="p_EN",
  b1="beta1",     b2="beta2",     gamma="gamma",
  aH="alpha_H",   pH="p_H",
  phi1="delta1",  phi2="delta2",  r="r"
)
y07_map <- setNames(y07_pub_orig[names(NAME_MAP_Y07)], NAME_MAP_Y07)

# Validate wrapper defaults vs y07par_gui.csv
# NOTE: r is stored negative in YASSO07_DEFAULT_PARAMS (-0.307) but published
#       as positive (|r|). The log transform in param_spec requires positive
#       free_defaults["r"] -- apply abs() in the calibration script.
y07_wrapper_vals <- YASSO07_DEFAULT_PARAMS[names(y07_map)]
y07_discrepancy  <- abs(abs(y07_wrapper_vals) - abs(y07_map)) /
  (abs(y07_map) + 1e-12)
y07_mismatch     <- y07_discrepancy > 1e-3

cat("\n=== Yasso07 published vs wrapper defaults ===\n")
cat("  (|r| compared for both; r stored negative in wrapper -- abs() fix required)\n")
df07_check <- data.frame(
  parameter    = names(y07_map),
  published    = round(y07_map, 7),
  wrapper      = round(y07_wrapper_vals, 7),
  wrapper_abs  = round(abs(y07_wrapper_vals), 7),
  rel_diff_pct = round(y07_discrepancy * 100, 3),
  match        = !y07_mismatch,
  row.names    = NULL
)
print(df07_check, row.names = FALSE)
if (any(y07_mismatch))
  warning(sprintf("Yasso07: %d parameters differ >0.1%%: %s",
                  sum(y07_mismatch),
                  paste(names(y07_map)[y07_mismatch], collapse = ", ")))



# =============================================================================
# 4-9.  Posterior analysis (skipped when DEFAULTS_ONLY = TRUE)
# =============================================================================

if (!DEFAULTS_ONLY) {

# =============================================================================
# 4.  Sanity check: posterior means vs MAP
# =============================================================================

key_check <- c("beta1","gamma","gammaN","gammaH","delta1","delta2","r",
               "p_H","alpha_H")

cat("\n=== SANITY CHECK: posterior means vs MAP ===\n")
cat("Yasso15:\n")
kp15 <- intersect(key_check, colnames(y15_samples))
print(round(rbind(MAP  = y15_map[kp15],
                  Mean = colMeans(y15_samples)[kp15]), 5))

cat("\nYasso20 (rda = primary MAP; dat = secondary):\n")
kp20 <- intersect(key_check, colnames(y20_samples))
out20 <- rbind(MAP_rda = y20_map_rda[kp20],
               Mean    = colMeans(y20_samples)[kp20])
if (!is.null(y20_map_dat)) out20 <- rbind(out20, MAP_dat = y20_map_dat[kp20])
print(round(out20, 5))

cat("\nYasso07 (published MAP vs wrapper |value|):\n")
kp07 <- intersect(key_check, names(y07_map))
print(round(rbind(published   = y07_map[kp07],
                  wrapper_abs = abs(y07_wrapper_vals[kp07])), 5))

cat("\nShared climate/size parameters across models (rda/published MAP):\n")
shared <- c("beta1","beta2","gamma","delta1","delta2","r","p_H","alpha_H")
shared15 <- intersect(shared, names(y15_map))
shared20 <- intersect(shared, names(y20_map))
shared07 <- intersect(shared, names(y07_map))
all_shared <- unique(c(shared15, shared20, shared07))
cmp_shared <- data.frame(
  parameter = all_shared,
  Yasso07   = round(abs(y07_map)[all_shared], 5),
  Yasso15   = round(y15_map[all_shared], 5),
  Yasso20   = round(y20_map[all_shared], 5),
  row.names = NULL
)
print(cmp_shared, row.names = FALSE)
cat("  (Yasso07 r shown as |r|; stored negative in wrapper)\n")


# =============================================================================
# 5.  Posterior summary (SD = structural zero detection)
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
    cat(sprintf("\n  STRUCTURAL ZEROS (fixed in calibration): %s\n",
                paste(names(mn)[fixed], collapse = ", ")))
  invisible(df)
}

s15 <- post_summary(y15_samples, y15_map, "Yasso15")
s20 <- post_summary(y20_samples, y20_map, "Yasso20")


# =============================================================================
# 6.  Sigma_ppm computation from posterior
# =============================================================================
# Log-transformed in param_spec (enforce positivity): beta1, betaN1, betaH1,
#   delta2, r -> sigma_ppm = SD(log(samples))
# Unconstrained: gamma, gammaN, gammaH, beta2, betaN2, betaH2, delta1
#   -> sigma_ppm = SD(samples) in physical space
# gammaH capped at 1.5: near-unidentifiable at Finnish P (~600mm); large SD
#   is artefact of flat likelihood surface at extreme negative values.

GAMMATH_SD_CAP <- 1.5
LOG_PARAMS     <- c("beta1","betaN1","betaH1","delta2","r")
CLIMATE_SIZE   <- c("beta1","beta2","gamma",
                    "betaN1","betaN2","gammaN",
                    "betaH1","betaH2","gammaH",
                    "delta1","delta2","r")

compute_sigma_ppm <- function(samples, model_name) {
  sds   <- apply(samples, 2L, sd)
  sigma <- setNames(rep(1.0, ncol(samples)), colnames(samples))
  
  for (p in intersect(CLIMATE_SIZE, colnames(samples))) {
    if (sds[p] < 1e-8) { sigma[p] <- 0.05; next }   # structural zero: tight
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
    space     = ifelse(names(sigma) %in% LOG_PARAMS,                  "log",
                       ifelse(names(sigma) %in% CLIMATE_SIZE,                "physical",
                              ifelse(names(sigma) %in% c("sigma_init","sigma_input"),"log",
                                     "—"))),
    source    = ifelse(names(sigma) %in% CLIMATE_SIZE,                "posterior SD",
                       ifelse(names(sigma) %in% c("sigma_init","sigma_input"),"fixed 0.5",
                              "wide (fractions)")),
    row.names = NULL
  )
  print(df, row.names = FALSE)
  invisible(sigma)
}

sigma15 <- compute_sigma_ppm(y15_samples, "Yasso15")
sigma20 <- compute_sigma_ppm(y20_samples, "Yasso20")


# =============================================================================
# 7.  Free defaults from MAP (prior centres in physical space)
# =============================================================================
# Climate/size params: MAP from primary source (rda for Y20, y15par for Y15,
#   published for Y07). Transfer fractions: keep model defaults.
# sigma_init / sigma_input: start at 1.0 (no a priori scaling).
# Yasso07 r: abs() applied here and flagged for calibration script.

print_defaults <- function(map_vals, model_name, default_params, param_names,
                           apply_abs_r = FALSE) {
  defaults <- default_params[param_names]
  for (p in intersect(CLIMATE_SIZE, param_names))
    defaults[p] <- map_vals[p]
  if (apply_abs_r && "r" %in% names(defaults))
    defaults["r"] <- abs(defaults["r"])
  if ("sigma_init"  %in% names(defaults)) defaults["sigma_init"]  <- 1.0
  if ("sigma_input" %in% names(defaults)) defaults["sigma_input"] <- 1.0
  
  cat(sprintf("\n=== %s free_defaults (physical space) ===\n", model_name))
  df <- data.frame(
    parameter     = names(defaults),
    model_default = round(default_params[names(defaults)], 6),
    prior_centre  = round(defaults, 6),
    source        = ifelse(names(defaults) %in% CLIMATE_SIZE, "MAP (primary)",
                           ifelse(names(defaults) %in% c("sigma_init","sigma_input"),
                                  "fixed 1.0", "model default")),
    row.names = NULL
  )
  if (apply_abs_r && "r" %in% rownames(df))
    df$source[df$parameter == "r"] <- "MAP (abs applied; wrapper stores -|r|)"
  print(df, row.names = FALSE)
  invisible(defaults)
}

defaults15 <- print_defaults(y15_map, "Yasso15",
                             YASSO15_DEFAULT_PARAMS, YASSO15_PARAM_NAMES)
defaults20 <- print_defaults(y20_map, "Yasso20",
                             YASSO20_DEFAULT_PARAMS, YASSO20_PARAM_NAMES)

# Yasso07: MAP from y07par_gui.csv; abs() on r
y07_free_names <- c(
  "p_WA","p_EA","p_NA","p_AW","p_EW","p_NW",
  "p_AE","p_WE","p_NE","p_AN","p_WN","p_EN",
  "beta1","beta2","gamma","delta1","delta2","r"
)
defaults07 <- YASSO07_DEFAULT_PARAMS[y07_free_names]
for (p in intersect(CLIMATE_SIZE, y07_free_names))
  defaults07[p] <- y07_map[p]
defaults07["r"] <- abs(defaults07["r"])

cat("\n=== Yasso07 free_defaults (physical space) ===\n")
df07_def <- data.frame(
  parameter    = names(defaults07),
  prior_centre = round(defaults07, 6),
  source       = ifelse(names(defaults07) %in% CLIMATE_SIZE,
                        "published MAP (Tuomi 2008-2010)",
                        "model default"),
  row.names = NULL
)
df07_def$source[df07_def$parameter == "r"] <-
  "published MAP, abs() applied (wrapper stores -|r|)"
print(df07_def, row.names = FALSE)


# =============================================================================
# 8.  Yasso20 structural zeros
# =============================================================================

y20_fixed <- names(which(apply(y20_samples, 2L, sd) < 1e-8))
cat("\n=== YASSO20 STRUCTURAL ZEROS ===\n")
cat("Move from param_spec to FIXED_RATE_NAMES in run_Yasso20_calibration.R:\n")
cat(sprintf("  c(%s)\n", paste(sprintf('"%s"', y20_fixed), collapse = ", ")))
cat("Values (should all be 0.0 from rda):\n")
print(round(y20_map[y20_fixed], 6))


# =============================================================================
# 9.  Ready-to-paste sigma_ppm blocks
# =============================================================================

cat_sigma_block <- function(sigma, model_name) {
  params <- intersect(CLIMATE_SIZE, names(sigma))
  cat(sprintf(
    "\n# ---- sigma_ppm: %s (paste into calibration script) ----\n", model_name))
  cat("sigma_ppm <- setNames(rep(1.0, N_FREE), FREE_NAMES)\n")
  for (p in params)
    cat(sprintf('sigma_ppm["%s"] <- %.5f\n', p, sigma[p]))
  cat('sigma_ppm["sigma_init"]  <- 0.5\n')
  cat('sigma_ppm["sigma_input"] <- 0.5\n')
  cat("# Transfer fractions remain at 1.0\n")
}

cat("\n=== READY-TO-PASTE SIGMA_PPM BLOCKS ===\n")
cat_sigma_block(sigma15, "Yasso15")
cat_sigma_block(sigma20, "Yasso20")

cat('
# ---- sigma_ppm: Yasso07 (Tuomi et al. 2008/2009/2010 CIs) ----
sigma_ppm <- setNames(rep(1.0, N_FREE), FREE_NAMES)
sigma_ppm["beta1"]  <- 0.20    # log space; exp(0.2*T) at T~10C: 2x range
sigma_ppm["beta2"]  <- 0.05    # physical space; quadratic, tight
sigma_ppm["gamma"]  <- 0.30    # physical space; keeps modifier in [0.3,0.6] at Finnish P
sigma_ppm["delta1"] <- 0.15    # physical space; Tuomi 2010 Table 4: +/-0.16 cm^-1
sigma_ppm["delta2"] <- 0.10    # log space;      Tuomi 2010 Table 4: +/-0.10 cm^-2
sigma_ppm["r"]      <- 0.015   # log space;      Tuomi 2010 Table 4: +/-0.013
sigma_ppm["sigma_init"]  <- 0.5
sigma_ppm["sigma_input"] <- 0.5
# Transfer fractions remain at 1.0
')


# =============================================================================
# 10.  Param_spec changes required in calibration scripts
# =============================================================================

cat('
=== PARAM_SPEC CHANGES (all Yasso calibration scripts) ===

Replace current grouped entries for size/climate with individual entries:

  # --- Size modifier (all models) ---
  list(names = "delta1", type = "unconstrained")  # can be negative
  list(names = "delta2", type = "log")             # enforce > 0 (NaN fix for d>0.4cm)
  list(names = "r",      type = "log")             # only -ABS(r) used in Fortran

  # --- Climate: Yasso07 ---
  list(names = "beta1",  type = "log")             # T response; enforce > 0
  list(names = "beta2",  type = "unconstrained")
  list(names = "gamma",  type = "unconstrained")   # sign enforced by tight prior

  # --- Climate: Yasso15/20 ---
  list(names = "beta1",   type = "log")
  list(names = "beta2",   type = "unconstrained")
  list(names = "gamma",   type = "unconstrained")
  list(names = "betaN1",  type = "log")
  list(names = "betaN2",  type = "unconstrained")
  list(names = "gammaN",  type = "unconstrained")
  list(names = "betaH1",  type = "log")
  list(names = "betaH2",  type = "unconstrained")
  list(names = "gammaH",  type = "unconstrained")

  # --- Yasso07 calibration script only: one-line fix after free_defaults ---
  free_defaults["r"] <- abs(free_defaults["r"])
  # (YASSO07_DEFAULT_PARAMS stores r as -0.307; log transform requires positive)
')

} # end if (!DEFAULTS_ONLY)


# =============================================================================
# 11.  Yasso15 wrapper validation (analogous to Yasso07 check, Section 3)
# =============================================================================
# Flags any discrepancy between YASSO15_DEFAULT_PARAMS (hardcoded in wrapper)
# and the authoritative y15par.csv values. Run standalone to verify.

if (!DEFAULTS_ONLY) {
  y15_wrapper_vals <- YASSO15_DEFAULT_PARAMS[names(y15_map)]
  y15_discrepancy  <- abs(y15_wrapper_vals - y15_map) / (abs(y15_map) + 1e-12)
  y15_mismatch     <- y15_discrepancy > 1e-3
  cat("\n=== Yasso15 y15par.csv vs wrapper defaults ===\n")
  print(data.frame(
    parameter    = names(y15_map),
    published    = round(y15_map, 7),
    wrapper      = round(y15_wrapper_vals, 7),
    rel_diff_pct = round(y15_discrepancy * 100, 3),
    match        = !y15_mismatch,
    row.names    = NULL
  ), row.names = FALSE)
  if (any(y15_mismatch))
    warning(sprintf("Yasso15: %d parameters differ >0.1%%: %s",
                    sum(y15_mismatch),
                    paste(names(y15_map)[y15_mismatch], collapse = ", ")))
}


# =============================================================================
# 12.  Assign authoritative DEFAULT_PARAMS from source files
# =============================================================================
# These two lines are the reason calibration scripts source this file.
# They overwrite the aliased wrapper values with MAP values read in Section 3:
#   YASSO15_DEFAULT_PARAMS  <- y15_map      (source: y15par.csv)
#   YASSO20_DEFAULT_PARAMS  <- y20_map_rda  (source: Yasso20_sample_parameters.rda)
# Yasso07 wrapper already matches y07par_gui.csv (validated in Section 3).
#
# USAGE in calibration scripts (source AFTER wrappers):
#   DEFAULTS_ONLY <- TRUE
#   source("<path>/Priors_model_matching.R")

YASSO15_DEFAULT_PARAMS <- y15_map
YASSO20_DEFAULT_PARAMS <- y20_map_rda

message("\nPriors_model_matching.R complete -- no files written to disk.")