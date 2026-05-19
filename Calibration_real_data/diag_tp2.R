# =============================================================================
# diag_tp2.R
#
# Diagnostic script for the TP2 posterior predictive pipeline.
# Run from the HIKET root directory:
#   apptainer_wrapper exec Rscript Calibration_real_data/diag_tp2.R
#
# Tests the full single-plot forward pass outside mclapply to identify
# which step fails and what tp2_run actually returns.
# =============================================================================

source('./Calibration_real_data/calibration_engine.R')
source('./Model_functions_real_data/input_compatibility_layer.R')
source('./Model_functions_real_data/Decomposition_functions/SimpleModels/tp2_wrapper.R')
source('./Model_functions_real_data/Decomposition_functions/Yasso/yasso07_wrapper.R')

RUN_ID     <- "20260512_193632"
inputs_pkg <- readRDS(file.path('./Data/model_inputs',
                                sprintf('TP2_inputs_%s.rds', RUN_ID)))
posterior  <- readRDS(file.path('./Calibration_real_data/runs',
                                sprintf('TP2_posterior_%s.rds', RUN_ID)))

cat(sprintf('Posterior: %d samples x %d params\n',
            nrow(posterior), ncol(posterior)))
cat('Param names:', names(posterior), '\n')

pid    <- inputs_pkg$plots_real[[1]]
clim   <- inputs_pkg$climate_by_plot[[pid]]
inputs <- inputs_pkg$inputs_by_plot[[pid]]
lm     <- inputs_pkg$litter_means[[pid]]
p      <- unlist(posterior[1, ])

cat('\nParameter values (draw 1):\n')
print(round(p, 5))

# --- Step 1: xi_mean for steady state ---
xi_ss <- tryCatch(
  compute_xi_mean_yasso07(clim[1:20, ], p['beta1'], p['beta2'], p['gamma']),
  error = function(e) { cat('xi_ss ERROR:', e$message, '\n'); NULL })
cat(sprintf('\nxi_ss: %s  finite: %s\n', xi_ss,
            if (!is.null(xi_ss)) is.finite(xi_ss) else 'NULL'))

# --- Step 2: steady state ---
C_init <- tryCatch(
  tp2_steady_state(p, lm, xi_ss),
  error = function(e) { cat('C_init ERROR:', e$message, '\n'); NULL })
cat('C_init:', C_init, '\n')
cat('C_init finite:', all(is.finite(C_init)), '  all >= 0:', all(C_init >= 0), '\n')

# --- Step 3: xi time series ---
xi_array <- tryCatch(
  compute_xi_yasso07(clim$temp_mean, clim$temp_amplitude, clim$precip,
                     p['beta1'], p['beta2'], p['gamma']),
  error = function(e) { cat('xi_array ERROR:', e$message, '\n'); NULL })
cat(sprintf('xi_array length: %d  finite: %s  range: [%.4f, %.4f]\n',
            length(xi_array), all(is.finite(xi_array)),
            min(xi_array), max(xi_array)))

# --- Step 4: tp2_run ---
run_out <- tryCatch(
  tp2_run(inputs, p, C_init, xi_array),
  error = function(e) { cat('tp2_run ERROR:', e$message, '\n'); NULL })

cat('\nrun_out is.null:', is.null(run_out), '\n')
if (!is.null(run_out)) {
  cat('run_out names:', names(run_out), '\n')
  cat('run_out nrow:', nrow(run_out), '\n')
  cat('run_out class:', class(run_out), '\n')
  str(run_out)
}

# --- Step 5: data.frame construction as in predictive script ---
df_test <- tryCatch(
  data.frame(
    plot_id   = pid,
    year      = run_out$year,
    draw      = 1L,
    A         = run_out$A,
    H         = run_out$H,
    total_soc = run_out$total_soc
  ),
  error = function(e) { cat('data.frame ERROR:', e$message, '\n'); NULL })

if (!is.null(df_test)) {
  cat(sprintf('\ndata.frame OK: %d rows x %d cols\n', nrow(df_test), ncol(df_test)))
  print(head(df_test, 3))
} else {
  cat('data.frame: FAILED\n')
}


# --- Step 6: test inside mclapply with error capture ---
cat('\nTesting inside mclapply (3 plots)...\n')
library(parallel)

model_params <- unlist(posterior[1, ])

test_results <- mclapply(inputs_pkg$plots_real[1:3], function(pid) {
  tryCatch({
    clim   <- climate_by_plot[[pid]]
    inputs <- inputs_by_plot[[pid]]
    lm     <- litter_means[[pid]]
    
    xi_array <- compute_xi_tp2_engine(clim, model_params)
    n_ss     <- min(20L, nrow(clim))
    xi_ss    <- compute_xi_mean_tp2_engine(clim[seq_len(n_ss), , drop=FALSE],
                                           model_params)
    C_init   <- steady_state_tp2_engine(model_params, lm, xi_ss)
    run_out  <- tp2_run_engine(inputs, model_params, C_init, xi_array)
    
    data.frame(plot_id=pid, year=run_out$year, draw=1L,
               A=run_out$A, H=run_out$H, total_soc=run_out$total_soc)
    
  }, error = function(e) {
    cat(sprintf("  Worker error [%s]: %s\n", pid, e$message))
    NULL
  })
}, mc.cores = 3L)

cat('Results per plot:\n')
for (i in seq_along(test_results))
  cat(sprintf('  Plot %s: %s\n', inputs_pkg$plots_real[[i]],
              if (is.null(test_results[[i]])) 'NULL' else
                sprintf('data.frame %d rows', nrow(test_results[[i]]))))



# Extract so mclapply workers can find them
climate_by_plot <- inputs_pkg$climate_by_plot
inputs_by_plot  <- inputs_pkg$inputs_by_plot
litter_means    <- inputs_pkg$litter_means

# Also define engine wrappers (as in predictive script)
compute_xi_tp2_engine <- function(clim, model_params) {
  compute_xi_yasso07(clim$temp_mean, clim$temp_amplitude, clim$precip,
                     model_params["beta1"], model_params["beta2"], model_params["gamma"])
}
compute_xi_mean_tp2_engine <- function(clim_ss, model_params) {
  compute_xi_mean_yasso07(clim_ss, model_params["beta1"],
                          model_params["beta2"], model_params["gamma"])
}
steady_state_tp2_engine <- function(model_params, lm, xi_mean) {
  tp2_steady_state(model_params, lm, xi_mean)
}
tp2_run_engine <- function(inputs, model_params, C_init, xi_array) {
  tp2_run(inputs, model_params, C_init, xi_array)
}