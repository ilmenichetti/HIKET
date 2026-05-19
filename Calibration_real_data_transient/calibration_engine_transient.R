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
      if (is.null(run_out) ||
          any(!is.finite(run_out$total_soc)))     return(-Inf)

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
      #   (see *_transient_init() in model wrappers). Adding it again in
      #   sd_vec would double-count. All observations use sigma_obs_fixed only.
      # -------------------------------------------------------------------
      sd_vec <- if (transient_init) {
        SOC_hat * sigma_obs_fixed
      } else {
        ifelse(meta$is_first,
               SOC_hat * sqrt(sigma_obs_fixed^2 + sigma_init^2),
               SOC_hat * sigma_obs_fixed)
      }

      sum(dnorm(meta$soc_obs, mean = SOC_hat, sd = sd_vec, log = TRUE))

    }, mc.cores = n_cores)

    log_liks <- unlist(log_liks)
    if (any(!is.finite(log_liks))) return(-Inf)
    sum(log_liks) + log_jac
  })
}
