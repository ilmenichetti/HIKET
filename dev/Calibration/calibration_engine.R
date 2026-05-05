# =============================================================================
# calibration_engine.R
#
# Model-agnostic calibration engine for the HIKET multi-model SOC framework.
#
# PURPOSE:
#   Provides all machinery that is shared across the five decomposition models
#   (Yasso07, Yasso15, Yasso20, RothC, Q-model). Model-specific scripts source
#   this file and supply only what genuinely differs between models: a parameter
#   specification table, default parameter values, and three forward-run
#   functions (compute_xi, steady_state, run_model).
#
# DESIGN PRINCIPLE:
#   The engine knows nothing about soil carbon dynamics. It only knows about
#   parameter transformations, MCMC sampling, and statistical diagnostics.
#   All physical model logic lives in the model-specific wrappers, which are
#   the same functions used for forward (non-calibration) runs. This means
#   a function validated in forward mode is exactly the function evaluated
#   inside the likelihood -- there is no separate "calibration version".
#
# EXPORTED FUNCTIONS (called by model scripts):
#   build_transforms()             Section 2 -- generate transform functions
#                                               from a parameter spec table
#   run_prior_predictive_match()   Section 3 -- choose prior width sigma_ppm
#   make_likelihood()              Section 4 -- generic likelihood factory
#   run_mcmc_chains()              Section 5 -- MCMC orchestration
#   run_diagnostics()              Section 6 -- R-hat, ESS, trace/marginal PDFs
#   save_results()                 Section 7 -- serialise outputs
#
# DEPENDENCIES:
#   BayesianTools  -- DEzs sampler, createPrior, createBayesianSetup
#   parallel       -- mclapply for the inner plot loop
#   compiler       -- cmpfun for byte-compiling the likelihood
#   coda           -- gelman.diag, effectiveSize, mcmc.list
#
# =============================================================================

library(BayesianTools)
library(parallel)
library(compiler)
library(coda)


# =============================================================================
# 1.  Primitive transform helpers
#
#     These three families of functions implement the building blocks that
#     build_transforms() (Section 2) assembles into full parameter transforms.
#
#     WHY TRANSFORMATIONS ARE NEEDED:
#       DEzs and other gradient-free samplers work best in unconstrained
#       real space (R^n). Our physical parameters have constraints:
#         - flow fractions must be positive and sum to <= budget B
#         - rate constants must be positive
#         - single proportions must lie in (0,1)
#       We transform to unconstrained space for sampling, then back to physical
#       space for model evaluation. The Jacobian of this transformation must be
#       included in the likelihood to preserve the correct posterior geometry.
#
#     THE THREE PRIMITIVE TYPES:
#       stick_break / stick_break_inv :  compositional constraints (sum <= B)
#       exp / log                      :  positivity constraints
#       inv_logit / logit              :  unit interval constraints
# =============================================================================

# --- Logit / inverse-logit ---
# Standard sigmoid pair. Used directly for single proportions (logit type)
# and as the building block inside stick-breaking.
logit     <- function(p) log(p / (1 - p))
inv_logit <- function(x) 1 / (1 + exp(-x))


# --- Stick-breaking transform (unconstrained -> physical) ---
#
# Converts a vector of unconstrained reals v into a vector of positive
# fractions p that automatically satisfies sum(p) <= budget.
#
# The metaphor: imagine a stick of length `budget`. At each step i,
# we break off a fraction of the remaining stick, where the fraction
# is determined by inv_logit(v[i]). The remaining stick after each
# break is budget - sum(p[1..i-1]).
#
# This guarantees:
#   p[i] >= 0 for all i        (inv_logit maps to (0,1), times rem > 0)
#   sum(p) <= budget            (we never break off more than what remains)
#
# Used for Yasso AWEN inter-pool transfer fractions, where each column
# (e.g. p_AW, p_AE, p_AN) must sum to <= B = 1 - p_H.
stick_break <- function(v, budget) {
  p   <- numeric(length(v))
  rem <- budget                        # remaining stick length
  for (i in seq_along(v)) {
    p[i] <- inv_logit(v[i]) * rem      # break off a piece
    rem  <- rem - p[i]                 # shrink remaining
  }
  p
}

# --- Stick-breaking inverse (physical -> unconstrained) ---
#
# Reverses stick_break(). Used to initialise the sampler at published
# default parameters (best_x = to_unconstrained(free_defaults)).
#
# At each step i, the logit of (p[i] / remaining) gives the unconstrained
# value that would have produced p[i] via stick_break().
stick_break_inv <- function(p, budget) {
  v   <- numeric(length(p))
  rem <- budget
  for (i in seq_along(p)) {
    v[i] <- logit(p[i] / rem)          # invert the fraction step
    rem  <- rem - p[i]
  }
  v
}

# --- Log Jacobian of stick-breaking transform ---
#
# When we change variables from unconstrained v to physical p, we must
# multiply the prior/likelihood by the absolute Jacobian determinant
# |d_p / d_v| to account for the volume distortion.
#
# For stick-breaking, the Jacobian factorises as a product over steps:
#   |d_p[i] / d_v[i]| = p[i]/rem[i] * (1 - p[i]/rem[i]) * rem[i]
#                      = frac * (1 - frac) * rem
#
# We sum log terms here; these are added to the log-likelihood in the
# likelihood factory (Section 4). The correction is included in the
# likelihood rather than the prior to keep the prior object clean and
# allow BayesianTools to use it without modification.
log_jac_stick <- function(p, budget) {
  rem <- budget
  jac <- 0
  for (i in seq_along(p)) {
    frac <- p[i] / rem
    jac  <- jac + log(frac * (1 - frac) * rem)
    rem  <- rem - p[i]
  }
  jac
}


# =============================================================================
# 2.  Transform builder
#
#     Accepts a param_spec list and returns a named list containing:
#       $to_original(x)       : unconstrained vector -> physical vector
#       $to_unconstrained(p)  : physical vector -> unconstrained vector
#       $log_jacobian(x, p)   : log |d_physical / d_unconstrained|
#       $param_names          : character vector of all free parameter names
#       $n_params             : total number of free parameters
#
#     HOW param_spec WORKS:
#       param_spec is a list of groups, each describing a set of parameters
#       that share the same transformation type and constraint. The engine
#       loops over these groups to assemble the full transforms. This means
#       model scripts declare what each parameter IS (its type and constraint)
#       rather than writing custom transform functions.
#
#     SUPPORTED TYPES:
#       "stick_break"   : list(names = c(...), type = "stick_break", budget = B)
#                         For compositional groups: fractions summing to <= B.
#                         Used for Yasso AWEN inter-pool transfer columns.
#
#       "log"           : list(names = c(...), type = "log")
#                         For strictly positive parameters.
#                         Used for rate constants (alpha_i, k_i) and
#                         observation error terms (sigma_obs, sigma_init).
#
#       "unconstrained" : list(names = c(...), type = "unconstrained")
#                         For parameters with no constraint (identity transform).
#                         Used for temperature/precipitation response parameters
#                         (beta1, beta2, gamma) which can take any real value.
#                         Contributes zero to the Jacobian.
#
#       "logit"         : list(names = c(...), type = "logit")
#                         For single proportions in (0,1).
#                         Used for RothC DPM:RPM ratio, f_bio, f_hum.
#
#     EXTENSIBILITY:
#       Adding a new parameter type requires adding one case to the three
#       switch() statements below and one entry to the Jacobian switch.
#       All other engine code is unaffected.
# =============================================================================

build_transforms <- function(param_spec) {
  
  # Flatten all parameter names in declaration order.
  # This defines the layout of the unconstrained vector x passed to the sampler.
  all_names <- unlist(lapply(param_spec, `[[`, "names"))
  n_params  <- length(all_names)
  
  # --- to_original: unconstrained R^n -> physical (constrained) space ---
  #
  # Called on every likelihood evaluation (millions of times during MCMC).
  # Byte-compiled via cmpfun() in make_likelihood(); not compiled here because
  # it is also called during diagnostics and prior predictive matching where
  # compilation overhead is not worth it.
  to_original <- function(x) {
    if (is.null(names(x))) names(x) <- all_names
    p <- numeric(n_params)
    names(p) <- all_names
    for (grp in param_spec) {
      nms <- grp$names
      switch(grp$type,                     # dispatch by group type
             stick_break   = { p[nms] <- stick_break(x[nms], grp$budget) },
             log           = { p[nms] <- exp(x[nms]) },
             unconstrained = { p[nms] <- x[nms] },
             logit         = { p[nms] <- inv_logit(x[nms]) }
      )
    }
    p
  }
  
  # --- to_unconstrained: physical -> unconstrained R^n ---
  #
  # Called once before sampling to initialise best_x from default parameters,
  # and once per prior predictive matching sample. Not performance-critical.
  to_unconstrained <- function(p) {
    if (is.null(names(p))) names(p) <- all_names
    x <- numeric(n_params)
    names(x) <- all_names
    for (grp in param_spec) {
      nms <- grp$names
      switch(grp$type,                     # inverse transform per group
             stick_break   = { x[nms] <- stick_break_inv(p[nms], grp$budget) },
             log           = { x[nms] <- log(p[nms]) },
             unconstrained = { x[nms] <- p[nms] },
             logit         = { x[nms] <- logit(p[nms]) }
      )
    }
    x
  }
  
  # --- log_jacobian: log |d_physical / d_unconstrained| ---
  #
  # Summed across all parameter groups. The sign convention follows the
  # change-of-variables formula: we need the log determinant of the Jacobian
  # of the forward map (unconstrained -> physical). For diagonal transforms
  # (all types here) this is a sum of log absolute partial derivatives.
  #
  # By transform type:
  #   stick_break   : log_jac_stick() -- see Section 1
  #   log           : d(exp(x))/dx = exp(x) = p, so log|J| = sum(log p)
  #   unconstrained : identity, log|J| = 0
  #   logit         : d(inv_logit(x))/dx = p(1-p), so log|J| = sum(log p(1-p))
  log_jacobian <- function(x, p) {
    if (is.null(names(x))) names(x) <- all_names
    if (is.null(names(p))) names(p) <- all_names
    jac <- 0
    for (grp in param_spec) {
      nms <- grp$names
      switch(grp$type,                     # accumulate log |d p / d x|
             stick_break   = { jac <- jac + log_jac_stick(p[nms], grp$budget) },
             log           = { jac <- jac + sum(log(p[nms])) },
             unconstrained = { },                                # identity contributes 0
             logit         = { jac <- jac + sum(log(p[nms] * (1 - p[nms]))) }
      )
    }
    jac
  }
  
  # Return all four outputs as a named list.
  # Model scripts typically unpack these immediately after calling build_transforms().
  list(
    to_original      = to_original,
    to_unconstrained = to_unconstrained,
    log_jacobian     = log_jacobian,
    param_names      = all_names,
    n_params         = n_params
  )
}


# =============================================================================
# 3.  Prior predictive matching
#
#     PURPOSE:
#       The Gaussian prior is centred on best_x (the published defaults in
#       unconstrained space) with width sigma_ppm. We want this width to be
#       wide enough that the prior is weakly informative -- i.e. the prior
#       predictive distribution of SOC should have a coefficient of variation
#       (CV) at least as large as the observed CV of SOC across plots.
#       This ensures the prior does not artificially constrain the posterior.
#
#     PROCEDURE:
#       For each candidate sigma in sigma_grid:
#         1. Draw n_samples random parameter vectors from N(best_x, sigma^2 I).
#         2. Transform each to physical space via to_original().
#         3. Run the forward model for one representative plot.
#         4. Compute CV of the resulting SOC predictions.
#       Find the smallest sigma such that the prior predictive CV >= observed CV.
#
#     NOTE ON SINGLE-PLOT EVALUATION:
#       We evaluate the prior predictive on a single plot (plots_real[1] in the
#       model scripts) for computational efficiency. This is sufficient because
#       we are matching the marginal prior predictive variance, not the full
#       joint distribution.
#
#     ARGUMENTS:
#       best_x       : starting point / prior centre (unconstrained space)
#       to_original  : transform function from build_transforms()
#       run_one_plot : function(p_free) -> scalar SOC prediction
#                      Defined in each model script; wraps the model's
#                      compute_xi / steady_state / run_model pipeline
#                      for a single fixed plot.
#       obs_cv       : observed CV = sd(SOC_obs) / mean(SOC_obs)
#       sigma_grid   : candidate prior widths to evaluate (default 0.1 to 3.0)
#       n_samples    : Monte Carlo samples per sigma value (default 200)
#       min_valid    : minimum number of finite predictions required to
#                      compute a CV estimate; if fewer, that sigma is skipped
# =============================================================================

run_prior_predictive_match <- function(best_x, to_original, run_one_plot,
                                       obs_cv,
                                       sigma_grid = seq(0.1, 3.0, by = 0.1),
                                       n_samples  = 200L,
                                       min_valid  = 20L) {
  
  message("Running prior predictive matching...")
  message(sprintf("  Target CV (observed): %.3f", obs_cv))
  
  # Evaluate predicted CV for each candidate sigma
  ppm_cv <- sapply(seq_along(sigma_grid), function(j) {
    sv <- sigma_grid[j]
    cat(sprintf("\r  sigma = %.1f (%d/%d)...", sv, j, length(sigma_grid)))
    flush.console()
    
    # Draw n_samples perturbations around best_x and run forward model
    preds <- vapply(seq_len(n_samples), function(i) {
      x_s <- best_x + rnorm(length(best_x), 0, sv)    # perturb in unconstrained space
      # tryCatch: invalid parameter combinations (e.g. negative xi) return NA
      tryCatch(run_one_plot(to_original(x_s)), error = function(e) NA_real_)
    }, numeric(1))
    
    # Require a minimum number of valid predictions for a reliable CV estimate
    valid <- preds[is.finite(preds) & preds > 0]
    if (length(valid) < min_valid) NA_real_ else sd(valid) / mean(valid)
  })
  cat("\n")
  
  # Select the smallest sigma whose prior predictive CV meets the target
  idx       <- which(ppm_cv >= obs_cv)[1]
  sigma_ppm <- if (is.na(idx)) {
    warning("Prior predictive matching did not reach target CV -- defaulting to 1.0")
    1.0
  } else {
    sigma_grid[idx]
  }
  
  message(sprintf("  Selected sigma_ppm = %.1f\n", sigma_ppm))
  sigma_ppm
}


# =============================================================================
# 3.5  Save processed inputs
#
#     PURPOSE:
#       Serialises all pre-processed data structures that entered the likelihood
#       to disk. Serves two purposes:
#         1. Enables post-hoc analysis (posterior predictive checks) without
#            re-sourcing the calibration script or raw data files.
#         2. Provides a validation record: the saved structures reflect exactly
#            what the model saw, independently verifiable against the raw inputs.
#
#     WHAT IS SAVED:
#       All per-plot lookup structures built during data preparation, plus
#       the prior centre and starting point for reference. Raw CSV data is
#       NOT saved -- only the processed structures that entered the likelihood.
#
#     ARGUMENTS:
#       inputs_list : named list containing:
#         $climate_by_plot : named list of per-plot climate data frames
#         $inputs_by_plot  : named list of per-plot litter input data frames
#         $litter_means    : named list of per-plot mean litter fractions
#         $obs_meta        : named list of per-plot observation metadata
#         $plots           : character vector of plot IDs used in likelihood
#         $plots_real      : character vector of unique real plot IDs
#         $obs_cv          : observed SOC coefficient of variation
#         $free_defaults   : default free parameters in physical space
#         $best_x          : default free parameters in unconstrained space
#       run_config  : list(MODEL_NAME, RUN_ID, DIR_INPUTS)
#
#     RETURNS:
#       Path to saved file (invisible)
# =============================================================================

save_inputs <- function(inputs_list, run_config) {
  
  MODEL_NAME <- run_config$MODEL_NAME
  RUN_ID     <- run_config$RUN_ID
  DIR_INPUTS <- run_config$DIR_INPUTS
  
  # Attach provenance metadata so the file is self-documenting
  inputs_list$provenance <- list(
    model       = MODEL_NAME,
    run_id      = RUN_ID,
    timestamp   = Sys.time(),
    n_plots     = length(inputs_list$plots),
    n_real      = length(inputs_list$plots_real),
    sigma_ppm   = inputs_list$sigma_ppm,   # add this
    r_version   = R.version$version.string
  )
  
  out_path <- file.path(DIR_INPUTS,
                        sprintf("%s_inputs_%s.rds", MODEL_NAME, RUN_ID))
  saveRDS(inputs_list, out_path)
  
  message(sprintf("Inputs saved: [%s]", out_path))
  invisible(out_path)
}

# =============================================================================
# 4.  Generic likelihood factory
#
#     PURPOSE:
#       Returns a byte-compiled likelihood function with all model callables
#       and pre-computed data structures captured in its closure. The factory
#       pattern (function that returns a function) serves two purposes:
#         1. Each chain gets its own independent closure with no shared state,
#            which is required for parallel chain execution.
#         2. n_cores is baked in at construction time; cmpfun() can have
#            problems resolving free variables, so we force() everything
#            explicitly before compiling.
#
#     THE THREE MODEL INTERFACE FUNCTIONS:
#       compute_xi(clim, model_params)
#           Takes plot-level climate data and the current parameter proposal.
#           Returns the environmental modifier time series xi(t) -- a vector
#           of length n_years (or n_months for RothC).
#           This is the SAME function used for forward runs, not a wrapper.
#
#       steady_state(model_params, litter_means, xi_mean)
#           Computes the analytical steady-state pool distribution C* from
#           mean climate and mean litter inputs. Used to initialise each plot
#           at t=0, avoiding a costly numerical spinup.
#           Called once per plot per likelihood evaluation.
#
#       run_model(inputs, model_params, C_init, xi_array)
#           Runs the full transient simulation for one plot given C* and xi(t).
#           Returns a data frame with columns year | pool1 | ... | total_soc.
#           Single Fortran call per plot; avoids R<->Fortran round-trips per
#           timestep.
#
#     OBSERVATION ERROR MODEL:
#       SOC_obs ~ Normal(SOC_hat, sd_vec^2)
#       where sd_vec = SOC_hat * sigma_obs                 (subsequent obs)
#             sd_vec = SOC_hat * sqrt(sigma_obs^2 +        (first observation
#                                     sigma_init^2)         per plot)
#
#       sigma_obs is FIXED across all models (passed as sigma_obs_fixed),
#       derived from the observed SOC coefficient of variation (obs_cv).
#       Fixing sigma_obs ensures structural model error is not absorbed into
#       the likelihood and shows up in residuals instead, enabling fair
#       comparison across model structures.
#
#       sigma_init remains a free parameter: it captures initialisation
#       uncertainty which genuinely differs by model (number of pools,
#       steady-state assumptions) and should be inferred from the data.
#
#       sigma_obs_fixed is set once in the model script as obs_cv, which is
#       computed before the prior predictive matching and reused here.
#
#     JACOBIAN:
#       The log Jacobian from build_transforms() is added to the summed
#       log-likelihood here. It compensates for the volume distortion of the
#       unconstrained -> physical transformation, ensuring that the sampler
#       explores the correct posterior rather than a distorted version of it.
#       Including it in the likelihood (rather than the prior density) keeps
#       the prior object usable by BayesianTools without modification.
#
#     -Inf GUARDS:
#       Any numerical failure (non-finite xi, negative C_init, non-finite SOC)
#       returns -Inf for that plot. This prevents invalid parameter combinations
#       from crashing the sampler; DEzs will simply reject such proposals.
#       The global check at the end catches the case where a small number of
#       plots fail due to numerical edge cases.
#
#     PERFORMANCE:
#       The inner mclapply() loop over plots is the only meaningful source of
#       parallelism. With n_cores = 90 (Roihu) and 3700 plots, each likelihood
#       evaluation distributes ~41 plot evaluations per core. The Fortran
#       model functions do the heavy lifting; R overhead per plot is minimal.
#
#     ARGUMENTS (all captured in closure):
#       n_cores          : cores for mclapply plot loop
#       to_original      : from build_transforms()
#       log_jacobian     : from build_transforms()
#       assemble_params  : function(p_free) -> full model parameter vector
#                          Defined in each model script. Combines free
#                          (calibrated) and fixed parameters in the order
#                          expected by the Fortran wrappers.
#       compute_xi       : model interface function (see above)
#       steady_state     : model interface function (see above)
#       run_model        : model interface function (see above)
#       sigma_obs_fixed  : fixed proportional observation error (scalar).
#                          Set to obs_cv in each model script. Same value
#                          used for all models to enable fair structural
#                          comparison.
#       plots            : character vector of plot IDs to evaluate
#                          May be a replicated subset (testing) or the full
#                          3700-plot dataset (production).
#       climate_by_plot  : named list of per-plot climate data frames,
#                          pre-split once before MCMC to avoid O(n) scanning
#       inputs_by_plot   : named list of per-plot litter input data frames
#       litter_means     : named list of per-plot mean litter fractions,
#                          used by steady_state()
#       obs_meta         : named list of per-plot observation metadata:
#                          $idx      : integer indices into the simulation
#                                      output matching observation years
#                          $soc_obs  : observed SOC values (tC/ha)
#                          $is_first : logical vector, TRUE for the first
#                                      observation per plot
# =============================================================================

make_likelihood <- function(n_cores,
                            to_original,
                            log_jacobian,
                            assemble_params,
                            compute_xi,
                            steady_state,
                            run_model,
                            sigma_obs_fixed,
                            plots,
                            climate_by_plot,
                            inputs_by_plot,
                            litter_means,
                            obs_meta,
                            steady_state_n = NULL) {
  
  # Force all arguments into the closure explicitly before byte-compiling.
  # cmpfun() freezes the environment at compile time; without force(),
  # lazy evaluation can cause closures to share the wrong bindings,
  # especially when make_likelihood() is called inside a loop (e.g. per chain).
  force(n_cores); force(to_original); force(log_jacobian)
  force(assemble_params); force(compute_xi); force(steady_state)
  force(run_model); force(sigma_obs_fixed); force(plots)
  force(climate_by_plot)
  force(inputs_by_plot); force(litter_means); force(obs_meta)
  force(steady_state_n)                   # capture for closure
  
  # Byte-compile the inner function. This gives a modest speedup (~10-20%) for
  # the R-level bookkeeping (parameter assembly, Jacobian, result aggregation).
  # The dominant cost is the Fortran model calls, which are unaffected.
  cmpfun(function(x) {
    
    # --- Step 1: transform proposal to physical space ---
    p_free       <- to_original(x)
    sigma_init   <- p_free["sigma_init"]
    model_params <- assemble_params(p_free)   # adds fixed parameters
    
    # Jacobian correction for the unconstrained -> physical transformation.
    # Must be added once per likelihood evaluation, not once per plot.
    log_jac <- log_jacobian(x, p_free)
    
    # --- Step 2: parallel plot loop ---
    # Each worker evaluates one plot's contribution to the total log-likelihood.
    # Workers share no state; all data is accessed via pre-built named lists.
    log_liks <- parallel::mclapply(plots, function(pid) {
      
      # Retrieve pre-split data for this plot (O(1) lookup)
      clim   <- climate_by_plot[[pid]]
      inputs <- inputs_by_plot[[pid]]
      lm     <- litter_means[[pid]]
      meta   <- obs_meta[[pid]]
      
      # Guard: if observation year indices are missing, skip this plot
      if (any(is.na(meta$idx))) return(-Inf)
      
      # --- Stage 1: environmental modifier time series ---
      # compute_xi is the model's own forward-run function, unchanged.
      xi_array <- tryCatch(
        compute_xi(clim, model_params),
        error = function(e) NULL)
      if (is.null(xi_array))                    return(-Inf)  # numerical failure
      n_ss <- if (is.null(steady_state_n)) length(xi_array)    # full series
      else min(steady_state_n, length(xi_array))   # or first N years
      xi_mean <- mean(xi_array[seq_len(n_ss)])
      if (!is.finite(xi_mean) || xi_mean <= 0)  return(-Inf)  # degenerate xi
      
      # --- Stage 2: steady-state initialisation ---
      # Analytical solution for C* given current parameters and mean climate.
      # This replaces a costly numerical spinup (1000+ years).
      C_init <- tryCatch(
        steady_state(model_params, lm, xi_mean),
        error = function(e) NULL)
      if (is.null(C_init) ||
          any(!is.finite(C_init)) ||
          any(C_init < 0))                      return(-Inf)  # invalid pools
      
      # --- Stage 3: transient simulation ---
      # Single Fortran call per plot; returns full SOC trajectory.
      run_out <- tryCatch(
        run_model(inputs, model_params, C_init, xi_array),
        error = function(e) NULL)
      if (is.null(run_out) ||
          any(!is.finite(run_out$total_soc)))   return(-Inf)
      
      # Extract model predictions at observation years
      SOC_hat <- run_out$total_soc[meta$idx]
      if (any(!is.finite(SOC_hat)) ||
          any(SOC_hat <= 0))                    return(-Inf)  # non-physical SOC
      
      # --- Observation error model ---
      # Multiplicative normal: sd proportional to predicted SOC.
      # sigma_obs is fixed (= obs_cv, set once in the model script) so that
      # structural error is not absorbed and remains visible in residuals.
      # sigma_init is free: captures initialisation uncertainty at t=0.
      # sd_vec <- ifelse(meta$is_first,
      #                  SOC_hat * sqrt(sigma_obs_fixed^2 + sigma_init^2),
      #                  SOC_hat * sigma_obs_fixed)
      # 
      # # Sum log-likelihood contributions across all observations for this plot
      # sum(dnorm(meta$soc_obs, mean = SOC_hat, sd = sd_vec, log = TRUE))
      
      #attempting a fix, model is overpredicting badly
      sigma_vec <- ifelse(meta$is_first,
                          sqrt(sigma_obs_fixed^2 + sigma_init^2),
                          sigma_obs_fixed)
      sum(dnorm(log(meta$soc_obs), mean = log(SOC_hat), sd = sigma_vec, log = TRUE))
      
      
    }, mc.cores = n_cores)
    
    # --- Step 3: aggregate across plots ---
    log_liks <- unlist(log_liks)
    
    # If any plot returned -Inf, the whole proposal is rejected.
    # NOTE: a small number of -Inf returns from numerically marginal plots
    # might be acceptable in production; consider changing to
    # sum(log_liks[is.finite(log_liks)]) if this causes excessive rejection.
    if (any(!is.finite(log_liks))) return(-Inf)
    
    # Total log-likelihood + Jacobian correction
    sum(log_liks) + log_jac
  })
}


# =============================================================================
# 5.  MCMC orchestration
#
#     PURPOSE:
#       Runs N_CHAINS independent DEzs chains and collects results.
#       Each chain has its own log file written to DIR_LOGS.
#
#     SAMPLER CHOICE:
#       DEzs (Differential Evolution with snooker updates) is a
#       population-based sampler well-suited to correlated, moderately
#       high-dimensional posteriors like ours (20-31 free parameters).
#       It is implemented in BayesianTools and requires no gradient information,
#       which is appropriate for our black-box Fortran model functions.
#
#     PARALLELISATION STRATEGY:
#       Single-nesting (current): chains run sequentially (lapply), each chain
#       uses all available cores for the plot loop (mclapply). This is the
#       only safe strategy on macOS, where nested forking is not supported.
#
#       Hybrid (Roihu / Linux): replace the outer lapply with
#       mclapply(mc.cores = N_CHAINS) and set CORES_PER_CHAIN = 90.
#       This gives N_CHAINS * CORES_PER_CHAIN total cores in use simultaneously.
#       Nested forking works on Linux; 4 chains * 90 cores = 360 of 384 cores.
#
#     LOGGING:
#       Each chain writes a progress line to its own log file every N_LOG
#       evaluations, recording: timestamp, evaluation count, log-likelihood,
#       evaluations per second, and elapsed time. This allows monitoring of
#       long Roihu runs without attaching to the R session.
#
#     ARGUMENTS:
#       ll_fn          : likelihood function from make_likelihood()
#       prior          : BayesianTools prior object (Gaussian in unconstrained
#                        space, centred on best_x with width sigma_ppm)
#       best_x         : chain starting point in unconstrained space
#                        (= to_unconstrained(published_defaults))
#       param_names    : character vector of free parameter names (for labels)
#       mcmc_settings  : list(N_CHAINS, N_ITER, N_BURNIN, N_LOG)
#       run_config     : list(MODEL_NAME, RUN_ID, DIR_LOGS,
#                             CORES_PER_CHAIN, n_plots)
#
#     RETURNS:
#       List of N_CHAINS BayesianTools chain result objects.
# =============================================================================

run_mcmc_chains <- function(ll_fn, prior, best_x, param_names,
                            mcmc_settings, run_config) {
  
  N_CHAINS        <- mcmc_settings$N_CHAINS
  N_ITER          <- mcmc_settings$N_ITER
  N_BURNIN        <- mcmc_settings$N_BURNIN
  N_LOG           <- mcmc_settings$N_LOG
  MODEL_NAME      <- run_config$MODEL_NAME
  RUN_ID          <- run_config$RUN_ID
  DIR_LOGS        <- run_config$DIR_LOGS
  CORES_PER_CHAIN <- run_config$CORES_PER_CHAIN
  
  message("-------------------------------------------------------------")
  message(sprintf("Starting MCMC: %d chains x %d iterations", N_CHAINS, N_ITER))
  message(sprintf("  Cores per chain: %d | Burn-in: %d | Log every: %d",
                  CORES_PER_CHAIN, N_BURNIN, N_LOG))
  message("-------------------------------------------------------------\n")
  
  set.seed(2025)   # reproducibility; fixed across all model scripts
  
  # Outer loop: sequential on Mac, replace with mclapply on Roihu.
  # Each iteration runs one complete chain.
  chain_results <- lapply(seq_len(N_CHAINS), function(chain_id) {
    
    # Per-chain log file. Named with model, chain index, and run ID so
    # multiple model runs can share the same log directory without collision.
    log_file <- file.path(
      DIR_LOGS,
      sprintf("%s_chain_%02d_%s.log", MODEL_NAME, chain_id, RUN_ID))
    
    # Write chain header to log file
    cat(sprintf("==========================================================\n"),
        file = log_file)
    cat(sprintf("  %s  Chain %d / %d\n", MODEL_NAME, chain_id, N_CHAINS),
        file = log_file, append = TRUE)
    cat(sprintf("  Run ID:  %s\n", RUN_ID),    file = log_file, append = TRUE)
    cat(sprintf("  Started: %s\n", Sys.time()), file = log_file, append = TRUE)
    cat(sprintf("  Plots:   %d | Iter: %d | Cores: %d\n",
                run_config$n_plots, N_ITER, CORES_PER_CHAIN),
        file = log_file, append = TRUE)
    cat(sprintf("==========================================================\n\n"),
        file = log_file, append = TRUE)
    
    message(sprintf("  Chain %d / %d  -->  %s", chain_id, N_CHAINS, log_file))
    
    # Per-chain evaluation counter and timer.
    # These are mutable via <<- inside the logging wrapper below.
    call_count    <- 0L                          # eval counter (mutated)
    t_chain_start <- proc.time()[["elapsed"]]    # chain wall-clock start
    
    # Logging wrapper: wraps ll_fn with a counter and periodic file writes.
    # This is the function actually passed to BayesianTools; the sampler sees
    # only a standard likelihood function with no side effects other than
    # the log file writes. The wrapper adds negligible overhead (~microseconds).
    ll_fn_logged <- function(x) {
      call_count <<- call_count + 1L
      result <- ll_fn(x)
      if (call_count %% N_LOG == 0L) {
        elapsed <- proc.time()[["elapsed"]] - t_chain_start
        cat(sprintf("[%s]  eval %5d  |  ll: %10.3f  |  %.2f eval/s  |  %.1f min\n",
                    format(Sys.time(), "%H:%M:%S"),
                    call_count, result,
                    call_count / elapsed, elapsed / 60),
            file = log_file, append = TRUE)
      }
      result
    }
    
    # Create BayesianSetup for this chain.
    # parallel = FALSE: we handle parallelism inside ll_fn via mclapply;
    # BayesianTools' built-in parallel support is not used.
    setup <- createBayesianSetup(
      likelihood = ll_fn_logged,
      prior      = prior,
      names      = param_names,
      parallel   = FALSE)
    
    # Run DEzs. nrChains = 1 because we manage the chain loop ourselves,
    # which gives us per-chain logging, independent closures, and the ability
    # to swap lapply for mclapply for Roihu without touching sampler code.
    result <- runMCMC(
      bayesianSetup = setup,
      sampler       = "DEzs",
      settings      = list(
        iterations = N_ITER,   # total iterations including burn-in
        burnin     = N_BURNIN, # discarded at the start
        nrChains   = 1L,       # one chain per call; see above
        eps        = 0.001))   # DEzs jitter term for ergodicity
    
    elapsed_chain <- proc.time()[["elapsed"]] - t_chain_start
    cat(sprintf("\nChain %d finished: %s  (%.1f min)\n",
                chain_id, Sys.time(), elapsed_chain / 60),
        file = log_file, append = TRUE)
    message(sprintf("  Chain %d / %d complete (%.1f min)",
                    chain_id, N_CHAINS, elapsed_chain / 60))
    result
  })
  
  chain_results
}


# =============================================================================
# 6.  Convergence diagnostics
#
#     PURPOSE:
#       Compute standard MCMC convergence diagnostics and generate output files
#       for all models in a consistent format. All files are written to DIR_DIAG
#       and prefixed with MODEL_NAME and RUN_ID so multiple model runs can share
#       the same directory.
#
#     DIAGNOSTICS COMPUTED:
#       Gelman-Rubin R-hat  : Values < 1.05 indicate good convergence.
#                             Values > 1.10 trigger a warning and suggest
#                             the run should be extended.
#       Effective sample size (ESS): Values < 200 indicate poor mixing.
#       Posterior quantile table: 2.5%, 25%, 50%, 75%, 97.5% in physical space.
#       Trace plots (PDF): Time series of selected parameters for each chain.
#       Marginal plots (PDF): Posterior vs prior density for selected parameters.
#
#     PHYSICAL SPACE REPORTING:
#       All diagnostics are reported in physical (constrained) parameter space
#       via to_original(). R-hat and ESS are computed in unconstrained space
#       (as the sampler sees them); the posterior summary is in physical space
#       for interpretability.
#
#     ARGUMENTS:
#       chain_results    : list of BayesianTools chain objects from run_mcmc_chains()
#       to_original      : transform function from build_transforms()
#       param_names      : character vector of free parameter names
#       best_x           : prior centre (unconstrained), used to sample prior
#                          draws for the marginal comparison plots
#       sigma_ppm        : prior width, used to sample prior draws
#       run_config       : list(MODEL_NAME, RUN_ID, DIR_DIAG, N_CHAINS,
#                               N_ITER, N_BURNIN, n_plots)
#       highlight_params : character vector of parameter names to include in
#                          trace and marginal plots. If NULL, all parameters
#                          are plotted. For models with 31 parameters, plotting
#                          all would produce unwieldy PDFs; recommend passing
#                          climate, modifier, and nuisance parameters only.
#
#     RETURNS:
#       Named list (invisible):
#         $all_raw   : matrix of raw (unconstrained) post-burnin samples
#         $all_phys  : matrix of physical-space post-burnin samples
#         $gr        : gelman.diag() output object, or NULL if it failed
#         $ess_vals  : effectiveSize() output vector, or NULL if it failed
# =============================================================================

run_diagnostics <- function(chain_results, to_original, param_names,
                            best_x, sigma_ppm, run_config,
                            mcmc_settings,
                            highlight_params = NULL) {
  
  MODEL_NAME <- run_config$MODEL_NAME
  RUN_ID     <- run_config$RUN_ID
  DIR_DIAG   <- run_config$DIR_DIAG
  N_CHAINS   <- mcmc_settings$N_CHAINS
  N_ITER     <- mcmc_settings$N_ITER
  N_BURNIN   <- mcmc_settings$N_BURNIN
  
  message("-------------------------------------------------------------")
  message("Computing convergence diagnostics...")
  message(sprintf("  Output: %s/%s_*_%s.*", DIR_DIAG, MODEL_NAME, RUN_ID))
  message("-------------------------------------------------------------\n")
  
  # --- Validate chain completion ---
  # Check that each chain returned a valid BayesianTools object and that the
  # majority of samples are finite. A chain is considered passing if >= 95%
  # of samples are fully finite; below this threshold the chain is flagged
  # but not automatically discarded (the researcher may still want to inspect).
  n_ok <- 0L
  for (i in seq_len(N_CHAINS)) {
    ok_i <- inherits(chain_results[[i]], "mcmcSamplerList") ||
      inherits(chain_results[[i]], "mcmcSampler")
    if (ok_i) {
      samp_i <- getSample(chain_results[[i]], coda = FALSE)
      pct    <- 100 * sum(apply(samp_i, 1, function(r) all(is.finite(r)))) /
        max(nrow(samp_i), 1)                        # fraction finite rows
      message(sprintf("  Chain %d: %d post-burnin samples, %.0f%% finite  [%s]",
                      i, nrow(samp_i), pct,
                      if (pct >= 95) "PASS" else "WARN"))
      n_ok <- n_ok + 1L
    } else {
      message(sprintf("  Chain %d: FAILED (not a valid chain object)", i))
    }
  }
  if (n_ok < N_CHAINS)
    stop(sprintf("%d / %d chains failed. Aborting diagnostics.", N_CHAINS - n_ok, N_CHAINS))
  
  # --- Combine post-burnin samples across chains ---
  # all_raw: unconstrained space (as sampled by DEzs) -- used for R-hat and ESS
  # all_phys: physical space (via to_original) -- used for posterior summaries
  all_raw  <- do.call(rbind,
                      lapply(chain_results, function(r) getSample(r, coda = FALSE)))
  all_phys <- t(apply(all_raw, 1, to_original))
  colnames(all_phys) <- param_names
  
  # coda::mcmc.list: required format for gelman.diag()
  chain_coda <- coda::mcmc.list(lapply(chain_results, function(r)
    coda::mcmc(getSample(r, coda = FALSE))))
  
  # --- Gelman-Rubin R-hat ---
  # multivariate = FALSE: compute univariate R-hat for each parameter separately.
  # The multivariate psrf is also printed in the file for reference.
  # tryCatch: gelman.diag() can fail if chains are degenerate; we report NULL
  # in that case and continue rather than stopping the entire diagnostics run.
  gr <- tryCatch(gelman.diag(chain_coda, multivariate = FALSE),
                 error = function(e) {
                   warning("gelman.diag() failed: ", conditionMessage(e))
                   NULL
                 })
  
  rhat_file <- file.path(DIR_DIAG,
                         sprintf("%s_rhat_%s.txt", MODEL_NAME, RUN_ID))
  sink(rhat_file)
  cat(sprintf("%s -- Gelman-Rubin R-hat\n", MODEL_NAME))
  cat(sprintf("Run: %s  |  Date: %s\n", RUN_ID, Sys.time()))
  cat(sprintf("Plots: %d | Chains: %d | Iter: %d | Burnin: %d\n\n",
              run_config$n_plots, N_CHAINS, N_ITER, N_BURNIN))
  if (!is.null(gr)) {
    rv  <- gr$psrf[, "Point est."]
    ruc <- gr$psrf[, "Upper C.I."]
    n_good <- sum(rv < 1.05, na.rm = TRUE)
    n_ok_r <- sum(rv >= 1.05 & rv < 1.10, na.rm = TRUE)
    n_warn <- sum(rv >= 1.10, na.rm = TRUE)
    cat(sprintf("Good (<1.05): %d | OK (1.05-1.10): %d | Warn (>1.10): %d\n\n",
                n_good, n_ok_r, n_warn))
    cat(sprintf("%-15s  %-10s  %-10s  %-8s\n",
                "Parameter", "Point est.", "Upper C.I.", "Status"))
    cat(strrep("-", 48), "\n")
    for (nm in param_names) {
      sta <- if (is.na(rv[nm])) "N/A" else
        if (rv[nm] < 1.05) "Good" else
          if (rv[nm] < 1.10) "OK" else "Warning"
      cat(sprintf("%-15s  %-10.4f  %-10.4f  %-8s\n", nm, rv[nm], ruc[nm], sta))
    }
    cat(sprintf("\nMultivariate psrf: %.4f\n", gr$mpsrf))
  } else {
    cat("R-hat computation failed -- see warnings.\n")
  }
  sink()
  
  if (!is.null(gr)) {
    n_warn <- sum(gr$psrf[,1] >= 1.10, na.rm = TRUE)
    message(sprintf("R-hat: %.3f -- %.3f  (%d warnings)  [%s]",
                    min(gr$psrf[,1], na.rm=TRUE), max(gr$psrf[,1], na.rm=TRUE),
                    n_warn, rhat_file))
    if (n_warn > 0)
      message(sprintf("  WARNING: %d parameter(s) with R-hat > 1.10 -- consider longer run", n_warn))
  }
  
  # --- Effective Sample Size ---
  # Computed on combined (across-chain) unconstrained samples.
  # ESS < 200 per parameter is a common threshold for reliable posterior summaries.
  ess_vals <- tryCatch(effectiveSize(coda::mcmc(all_raw)),
                       error = function(e) {
                         warning("effectiveSize() failed: ", conditionMessage(e))
                         NULL
                       })
  
  ess_file <- file.path(DIR_DIAG,
                        sprintf("%s_ess_%s.txt", MODEL_NAME, RUN_ID))
  sink(ess_file)
  cat(sprintf("%s -- Effective Sample Size\n", MODEL_NAME))
  cat(sprintf("Run: %s  |  Date: %s\n\n", RUN_ID, Sys.time()))
  if (!is.null(ess_vals)) {
    cat(sprintf("Total post-burnin samples: %d\n", nrow(all_raw)))
    cat(sprintf("ESS range: %.0f -- %.0f\n", min(ess_vals), max(ess_vals)))
    cat(sprintf("Below 200:  %d / %d\n\n", sum(ess_vals < 200), length(ess_vals)))
    cat(sprintf("%-15s  %-10s\n", "Parameter", "ESS"))
    cat(strrep("-", 28), "\n")
    for (nm in param_names)
      cat(sprintf("%-15s  %-10.0f\n", nm, ess_vals[nm]))
  } else {
    cat("ESS computation failed -- see warnings.\n")
  }
  sink()
  
  if (!is.null(ess_vals)) {
    n_low <- sum(ess_vals < 200, na.rm = TRUE)
    message(sprintf("ESS:   %.0f -- %.0f  (%d < 200)  [%s]",
                    min(ess_vals), max(ess_vals), n_low, ess_file))
    if (n_low > 0)
      message("  WARNING: Low ESS -- consider running longer chains")
  }
  
  # --- Posterior quantile summary (physical space) ---
  post_q    <- apply(all_phys, 2, quantile,
                     probs = c(0.025, 0.25, 0.50, 0.75, 0.975), na.rm = TRUE)
  post_file <- file.path(DIR_DIAG,
                         sprintf("%s_posterior_summary_%s.txt", MODEL_NAME, RUN_ID))
  sink(post_file)
  cat(sprintf("%s -- Posterior summary (physical space)\n", MODEL_NAME))
  cat(sprintf("Run: %s  |  Samples: %d\n\n", RUN_ID, nrow(all_phys)))
  cat(sprintf("%-15s  %8s  %8s  %8s  %8s  %8s\n",
              "Parameter", "2.5%", "25%", "50%", "75%", "97.5%"))
  cat(strrep("-", 62), "\n")
  for (nm in param_names)
    cat(sprintf("%-15s  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f\n", nm,
                post_q["2.5%",nm], post_q["25%",nm], post_q["50%",nm],
                post_q["75%",nm], post_q["97.5%",nm]))
  sink()
  message(sprintf("Posterior summary: [%s]", post_file))
  
  # --- Trace and marginal plots ---
  # If highlight_params is NULL, plot all parameters.
  # For models with 31 parameters, plotting all produces large PDFs;
  # passing climate + nuisance parameters keeps the PDFs concise and
  # focuses attention on the most policy-relevant quantities.
  plot_params <- if (is.null(highlight_params)) param_names else highlight_params
  n_pp        <- length(plot_params)
  n_cols      <- 2L
  n_rows      <- ceiling(n_pp / n_cols)
  cols        <- c("steelblue","firebrick","darkgreen","darkorange")[seq_len(N_CHAINS)]
  
  # Trace plots: post-burnin iteration vs parameter value, one line per chain
  trace_file <- file.path(DIR_DIAG,
                          sprintf("%s_traces_%s.pdf", MODEL_NAME, RUN_ID))
  pdf(trace_file, width = 10, height = max(4, n_rows * 2.5))
  par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 3, 1))
  for (nm in plot_params) {
    # Extract physical-space trace for this parameter from each chain
    traces <- lapply(chain_results, function(r) {
      phys <- t(apply(getSample(r, coda = FALSE), 1, to_original))  # back-transform
      phys[, nm]
    })
    ylim_r <- range(unlist(traces), na.rm = TRUE)
    plot(traces[[1]], type = "l", col = cols[1], ylim = ylim_r,
         xlab = "Post-burnin iteration", ylab = nm,
         main = sprintf("%s | %s", MODEL_NAME, nm))
    for (ci in seq_len(N_CHAINS)[-1]) lines(traces[[ci]], col = cols[ci])
    legend("topright", legend = sprintf("Chain %d", seq_len(N_CHAINS)),
           col = cols, lwd = 1, bty = "n", cex = 0.7)
  }
  dev.off()
  message(sprintf("Traces:    [%s]", trace_file))
  
  # Marginal plots: posterior vs prior density for selected parameters.
  # Prior draws are sampled from N(best_x, sigma_ppm^2 I) and transformed
  # to physical space. This gives an honest picture of what the data contributed
  # relative to the prior.
  set.seed(99)   # reproducible prior draws for the marginal plots
  prior_raw  <- matrix(rnorm(3000 * length(best_x), mean = best_x, sd = sigma_ppm),
                       nrow = 3000)
  prior_phys <- t(apply(prior_raw, 1, to_original))
  
  marg_file <- file.path(DIR_DIAG,
                         sprintf("%s_marginals_%s.pdf", MODEL_NAME, RUN_ID))
  pdf(marg_file, width = 10, height = max(4, n_rows * 2.5))
  par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 3, 1))
  for (nm in plot_params) {
    post_v  <- all_phys[is.finite(all_phys[, nm]), nm]
    prior_v <- prior_phys[is.finite(prior_phys[, nm]), nm]
    d_post  <- density(post_v)                   # posterior density
    d_prior <- density(prior_v)                  # prior density
    xlim_r  <- range(c(d_post$x, d_prior$x))
    ylim_r  <- range(c(d_post$y, d_prior$y))
    plot(d_post, col = "steelblue", lwd = 2,
         xlim = xlim_r, ylim = ylim_r,
         main = sprintf("%s | %s", MODEL_NAME, nm),
         xlab = nm, ylab = "Density")
    lines(d_prior, col = "grey60", lwd = 1, lty = 2)
    legend("topright",
           legend = c("Posterior", "Prior"),
           col    = c("steelblue", "grey60"),
           lwd    = c(2, 1), lty = c(1, 2), bty = "n", cex = 0.7)
  }
  dev.off()
  message(sprintf("Marginals: [%s]", marg_file))
  
  # Return samples invisibly; model scripts can use these for further analysis
  # without re-reading the saved files.
  invisible(list(
    all_raw  = all_raw,
    all_phys = all_phys,
    gr       = gr,
    ess_vals = ess_vals
  ))
}


# =============================================================================
# 7.  Save results
#
#     PURPOSE:
#       Serialise chain objects, posterior samples, and run metadata to disk.
#       Three separate files are written to DIR_RUNS:
#         MODEL_chains_RUNID.rds    : raw BayesianTools chain objects
#                                     (needed for chain extension / restart)
#         MODEL_posterior_RUNID.rds : physical-space posterior samples matrix
#                                     (ready for downstream analysis)
#         MODEL_metadata_RUNID.rds  : run configuration and summary statistics
#                                     (for bookkeeping and reproducibility)
#
#     METADATA CONTENTS:
#       All fields needed to reproduce or extend the run: model name, run ID,
#       timestamp, plot count, chain/iteration settings, prior width, wallclock
#       time, R version, platform, and summary convergence statistics.
#
#     ARGUMENTS:
#       chain_results  : from run_mcmc_chains()
#       all_phys       : physical-space samples matrix from run_diagnostics()
#       run_config     : list(MODEL_NAME, RUN_ID, DIR_RUNS, CORES_PER_CHAIN,
#                             n_plots)
#       mcmc_settings  : list(N_CHAINS, N_ITER, N_BURNIN, N_LOG)
#       sigma_ppm      : prior width from run_prior_predictive_match()
#       t_elapsed      : wallclock seconds from system.time()
#       gr             : gelman.diag() output from run_diagnostics(), or NULL
#       ess_vals       : effectiveSize() output from run_diagnostics(), or NULL
#
#     RETURNS:
#       metadata list (invisible)
# =============================================================================

save_results <- function(chain_results, all_phys, run_config,
                         mcmc_settings, sigma_ppm, sigma_obs_fixed,
                         t_elapsed, gr = NULL, ess_vals = NULL) {
  
  MODEL_NAME <- run_config$MODEL_NAME
  RUN_ID     <- run_config$RUN_ID
  DIR_RUNS   <- run_config$DIR_RUNS
  
  metadata <- list(
    model           = MODEL_NAME,
    run_id          = RUN_ID,
    timestamp       = Sys.time(),
    n_plots         = run_config$n_plots,
    n_chains        = mcmc_settings$N_CHAINS,
    n_iter          = mcmc_settings$N_ITER,
    n_burnin        = mcmc_settings$N_BURNIN,
    cores_per_chain = run_config$CORES_PER_CHAIN,
    sigma_ppm       = sigma_ppm,
    sigma_obs_fixed = sigma_obs_fixed,
    wallclock_min   = t_elapsed / 60,
    r_version       = R.version$version.string,
    platform        = .Platform$OS.type,
    # Convergence summaries: stored as ranges for quick inspection
    rhat_range      = if (!is.null(gr))       range(gr$psrf[,1], na.rm=TRUE) else NA,
    ess_range        = if (!is.null(ess_vals)) range(ess_vals, na.rm=TRUE)   else NA,
    n_samples       = nrow(all_phys)
  )
  
  # Save three files; all prefixed with MODEL_NAME and RUN_ID
  saveRDS(chain_results,
          file.path(DIR_RUNS, sprintf("%s_chains_%s.rds",    MODEL_NAME, RUN_ID)))
  saveRDS(all_phys,
          file.path(DIR_RUNS, sprintf("%s_posterior_%s.rds", MODEL_NAME, RUN_ID)))
  saveRDS(metadata,
          file.path(DIR_RUNS, sprintf("%s_metadata_%s.rds",  MODEL_NAME, RUN_ID)))
  
  message(sprintf("\nSaved to %s:", DIR_RUNS))
  message(sprintf("  %s_chains_%s.rds",    MODEL_NAME, RUN_ID))
  message(sprintf("  %s_posterior_%s.rds", MODEL_NAME, RUN_ID))
  message(sprintf("  %s_metadata_%s.rds",  MODEL_NAME, RUN_ID))
  
  invisible(metadata)
}