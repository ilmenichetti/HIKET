# =============================================================================
# calibration_engine_transient.R
#
# Transient pre-initialization extension of calibration_engine.R.
#
# WHAT THIS FILE DOES:
#   Sources the original calibration_engine.R, then redefines make_likelihood()
#   to add the transient_init flag. All other engine functions (build_transforms,
#   run_mcmc_chains, run_diagnostics, save_results, etc.) are inherited unchanged.
#
# NEW ARGUMENT in make_likelihood():
#   transient_init (logical, default FALSE):
#     FALSE -- original behaviour: sigma_init added in quadrature to sd_vec
#              at the first observation per plot (likelihood patch).
#     TRUE  -- transient pre-run behaviour: sigma_init enters the model
#              dynamics via the pre-run in the model-specific init function
#              (see *_wrapper_transient.R). sd_vec uses sigma_obs_fixed
#              for ALL observations; sigma_init no longer patches the likelihood.
#
#   With transient_init = FALSE the engine is bit-for-bit identical to the
#   original. This flag is the revert switch.
#
# REJECTION LOGIC:
#   Only three hard rejections are applied, all strictly necessary for
#   arithmetic validity:
#     (1) run_model() returned NULL (Fortran crash or tryCatch error)
#     (2) SOC_hat is non-finite (Inf/NaN) at observation years -- dnorm undefined
#     (3) SOC_hat <= 0 at observation years -- undefined in multiplicative error model
#
#   No physical-ceiling guard (e.g. SOC > 1000) is applied. Physically
#   impossible predictions receive a very large negative log-likelihood and
#   are rejected by MCMC naturally; a hard -Inf ceiling is unnecessary and
#   would asymmetrically penalise models (e.g. Yasso20) whose recycling
#   architecture explores a wider region of parameter space. Post-calibration-
#   window instability is never evaluated in the likelihood: the trajectory
#   beyond meta$idx is computed but not touched here. It appears in full in
#   the predictive scripts, where structural divergence across models is a
#   primary scientific output.
#
#   This matches Viskari (2022), who applies only simplex checks (outflow
#   fractions sum <= 1) and no physical ceiling.
# =============================================================================

source("./Calibration_real_data/calibration_engine.R")


# Override make_likelihood() with transient_init support.
# Only the sd_vec line and function signature differ from the original.
make_likelihood <- function(n_cores,
                            to_original,
                            log_jacobian,
                            assemble_params,
                            compute_xi,
                            compute_xi_mean,
                            steady_state,
                            run_model,
                            sigma_obs_fixed,
                            plots,
                            climate_by_plot,
                            inputs_by_plot,
                            litter_means,
                            obs_meta,
                            steady_state_n  = NULL,
                            transient_init  = FALSE) {
  
  force(n_cores); force(to_original); force(log_jacobian)
  force(assemble_params); force(compute_xi); force(compute_xi_mean)
  force(steady_state); force(run_model); force(sigma_obs_fixed); force(plots)
  force(climate_by_plot); force(inputs_by_plot); force(litter_means)
  force(obs_meta); force(steady_state_n); force(transient_init)
  
  cmpfun(function(x) {
    
    p_free       <- to_original(x)
    sigma_init   <- p_free["sigma_init"]
    model_params <- assemble_params(p_free)
    log_jac      <- log_jacobian(x, p_free)
    
    log_liks <- parallel::mclapply(plots, function(pid) {
      
      clim   <- climate_by_plot[[pid]]
      inputs <- inputs_by_plot[[pid]]
      lm     <- litter_means[[pid]]
      meta   <- obs_meta[[pid]]
      
      if (any(is.na(meta$idx))) return(-Inf)
      
      xi_array <- tryCatch(compute_xi(clim, model_params),
                           error = function(e) NULL)
      if (is.null(xi_array)) return(-Inf)
      
      n_ss <- if (is.null(steady_state_n)) nrow(clim)
      else min(steady_state_n, nrow(clim))
      
      xi_for_ss <- tryCatch(
        compute_xi_mean(clim[seq_len(n_ss), , drop = FALSE], model_params),
        error = function(e) NULL)
      if (is.null(xi_for_ss) || !is_valid_xi(xi_for_ss)) return(-Inf)
      
      C_init <- tryCatch(
        steady_state(model_params, lm, xi_for_ss),
        error = function(e) NULL)
      if (is.null(C_init) ||
          any(!is.finite(C_init)) ||
          any(C_init < 0))                        return(-Inf)
      
      run_out <- tryCatch(
        run_model(inputs, model_params, C_init, xi_array),
        error = function(e) NULL)
      if (is.null(run_out)) return(-Inf)
      
      # Extract predicted SOC at observation years only.
      # Trajectory values beyond meta$idx are not evaluated here; post-
      # calibration-window behaviour appears in the predictive scripts.
      SOC_hat <- run_out$total_soc[meta$idx]
      if (any(!is.finite(SOC_hat)) ||
          any(SOC_hat <= 0))                      return(-Inf)
      
      # -------------------------------------------------------------------
      # Observation error model
      # -------------------------------------------------------------------
      # transient_init = FALSE (original):
      #   sigma_init added in quadrature at first obs -- likelihood patch for
      #   uncertain analytical steady-state initialisation.
      #
      # transient_init = TRUE (new):
      #   sigma_init already propagated physically through the 68-year pre-run
      #   (see *_wrapper_transient.R). Adding it again in sd_vec would
      #   double-count. All observations use sigma_obs_fixed only.
      # -------------------------------------------------------------------
      sd_vec <- if (transient_init) {
        SOC_hat * sigma_obs_fixed
      } else {
        ifelse(meta$is_first,
               SOC_hat * sqrt(sigma_obs_fixed^2 + sigma_init^2),
               SOC_hat * sigma_obs_fixed)
      }

      # C5: per-observation SD inflation (down-weight the suspect 1985 campaign).
      # meta$sigma_infl is 1 everywhere except 1985 obs (SIGMA_1985_INFL). Absent
      # in older input bundles -> treated as 1 (no effect), so this is backward-safe.
      if (!is.null(meta$sigma_infl)) sd_vec <- sd_vec * meta$sigma_infl

      sum(dnorm(meta$soc_obs, mean = SOC_hat, sd = sd_vec, log = TRUE))
      
    }, mc.cores = n_cores)
    
    log_liks <- unlist(log_liks)
    if (any(!is.finite(log_liks))) return(-Inf)
    sum(log_liks) + log_jac
  })
}

# =============================================================================
# assert_inputs_current()  (2026-08-04)
#
# FAIL-FAST GUARD AGAINST A STALE INPUT BUNDLE.
#
# Data/ is gitignored, so the SOC-baseline swap and the 1985 litter
# reconstruction (both 2026-08-04) do NOT travel with `git pull`. On a cluster
# updated from git alone, the code would be current while
# Data/model_inputs/*.csv were months old -- and the run would complete
# NORMALLY, producing plausible-looking posteriors fitted to the wrong target.
# That is the same failure shape as the stale-.so incident: right code, wrong
# data, silent garbage, discovered only after the compute is spent.
#
# Two cheap invariants separate current from stale inputs unambiguously:
#   SOC target  : homogenized median ~64 tC/ha; the superseded (stoniness-
#                 inflated, 1m-extrapolated) target ran ~100.
#   1985 litter : reconstructed ~1.75 tC/ha/yr; the raw artefactual first year
#                 of the Tupek series is ~0.09.
#
# Called immediately after the input bundle is read in every run_*_calibration.R.
# =============================================================================
assert_inputs_current <- function(input_raw,
                                  soc_median_max = 80,
                                  litter_1985_min = 0.5) {
  soc <- input_raw$soc_obs_tCha
  soc <- soc[!is.na(soc)]
  if (!length(soc)) stop("assert_inputs_current: no soc_obs_tCha in the bundle.")
  soc_med <- median(soc)

  lit_cols <- grep("^C_(nwl|fwl|cwl)_", names(input_raw), value = TRUE)
  j85 <- input_raw[input_raw$year == 1985L, lit_cols, drop = FALSE]
  # rows are MONTHLY (annual/12) -> sum the 12 months to get tC/ha/yr per plot
  j85_annual <- tapply(rowSums(j85), input_raw$plot_id[input_raw$year == 1985L], sum)
  j85_med <- median(j85_annual, na.rm = TRUE)

  message(sprintf("Input currency check: median soc_obs = %.1f tC/ha | 1985 litter = %.2f tC/ha/yr",
                  soc_med, j85_med))

  if (soc_med > soc_median_max)
    stop(sprintf(paste0("STALE INPUT BUNDLE: median soc_obs = %.1f tC/ha (expected ~64, limit %.0f).\n",
                        "  This looks like the SUPERSEDED stoniness-inflated SOC target.\n",
                        "  Re-run Data/Data_work.R and rsync Data/model_inputs/ -- Data/ is gitignored\n",
                        "  and does NOT arrive via git pull."),
                 soc_med, soc_median_max))

  if (!is.finite(j85_med) || j85_med < litter_1985_min)
    stop(sprintf(paste0("STALE INPUT BUNDLE: median 1985 litter = %.3f tC/ha/yr (expected ~1.75).\n",
                        "  The artefactual first year of the Tupek series (~0.09) is still present;\n",
                        "  the 1986-1990 backcast has not been applied. Re-run Data/Data_work.R\n",
                        "  and rsync Data/model_inputs/."),
                 j85_med))

  invisible(TRUE)
}
