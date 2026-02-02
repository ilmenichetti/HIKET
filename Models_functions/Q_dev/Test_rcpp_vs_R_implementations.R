# =============================================================================
# TEST RCPP SPEEDUP
# =============================================================================
source("./Models_functions/Q_transient_T_hybrid.R") # load Q-model functions and he wrapper with the hybrid initialization approach for Q-model
source("./Models_functions/Q_transient_T_hybrid_RCPP.R") # same but with compiled bits with rcpp
source("./Models_functions/synthetic_data_template") # load the data (it wastes quite some calculations running also some tests...)


# OPTION 1: Test with R compiler (quick, 2-3x speedup)
# -----------------------------------------------------
library(compiler)
enableJIT(3)

cat("Testing with R compiler enabled...\n")
system.time({
  results_compiled <- run_soc_models_oneplot(
    input_df = input_data,
    site_row = site_data,
    models = c("q_model"),
    spinup_years = 20
  )
})

enableJIT(0)

cat("Testing with R compiler disabled\n")
system.time({
  results_compiled_0 <- run_soc_models_oneplot(
    input_df = input_data,
    site_row = site_data,
    models = c("q_model"),
    spinup_years = 20
  )
})

# OPTION 2: Test with Rcpp version (5-15x speedup)
# -------------------------------------------------
source("./Models_functions/Q_transient_T_hybrid_RCPP.R")

cat("\nTesting with Rcpp...\n")
system.time({
  results_rcpp <- q_model_run_hybrid_rcpp(
    input_df = map_input_to_q_model(input_data),
    site_row = site_data,
    spinup_years = 20
  )
})

# Compare results (should be nearly identical)
cat("\nFinal SOC comparison:\n")
cat(sprintf("  Compiled R: %.2f gC/m²\n", results_compiled$summary$q_model[1]))
cat(sprintf("  Compiled R: %.2f gC/m²\n", results_compiled_0$summary$q_model[1]))
cat(sprintf("  Rcpp:       %.2f gC/m²\n", results_rcpp$total_soc[1]))

