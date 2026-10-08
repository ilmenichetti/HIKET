# =============================================================================
# eqinit_nesting_test.R -- the equilibrium arm must NEST in production (2026-10-08)
# -----------------------------------------------------------------------------
# GATE before launching the equilibrium-init arm (NEXT_RUN_equilibrium_init.md).
# The arm is production with sigma_init fixed at 1, so at ANY parameter vector
# with sigma_init = 1 the two arms must give the SAME data log-likelihood. The
# likelihood functions differ only in the Jacobian (one free parameter fewer), so
# the comparison is ll_fn(x) - log_jacobian(x): the data term.
#
# Faithful: sources each REAL calibration script up to the MCMC launch (as
# preflight_prior_pushforward.R does), once per arm, in separate R processes
# because HIKET_EQUILIBRIUM_INIT is read at source time. Production run config
# (correlated likelihood, sigma_tot 0.800, 1985 inflation 1) is set here.
#
# Usage (project root):
#   Rscript doublechecks/eqinit_nesting_test.R [MODEL ...]   # default: all six
# Side effects removed at the end: the input bundle + sanity PNG each setup writes.
# =============================================================================

args   <- commandArgs(trailingOnly = TRUE)
WORKER <- Sys.getenv("EQTEST_WORKER", "")
OUTDIR <- Sys.getenv("EQTEST_OUTDIR", file.path(tempdir(), "eqinit_nesting"))
dir.create(OUTDIR, showWarnings = FALSE)
N_THETA <- 4L   # defaults + 3 prior draws, all with sigma_init = 1

# === worker: source one script up to the MCMC launch, evaluate the data ll ===
if (nzchar(WORKER)) {
  MODEL <- WORKER
  arm   <- if (Sys.getenv("HIKET_EQUILIBRIUM_INIT") == "1") "eq" else "prod"
  src   <- readLines(sprintf("Calibration_real_data_transient/run_%s_transient_calibration.R", MODEL))
  cutix <- grep("^t_run <- system.time\\(\\{", src)[1]
  e <- new.env(parent = globalenv())
  sink(file.path(OUTDIR, sprintf("%s_%s.log", MODEL, arm)))
  source(textConnection(paste(src[seq_len(cutix - 1L)], collapse = "\n")), local = e)
  sink()
  th_file <- file.path(OUTDIR, sprintf("%s_theta.rds", MODEL))
  if (arm == "prod") {
    # thetas in PHYSICAL space, drawn from the production prior, sigma_init := 1
    set.seed(2025)
    X <- rbind(e$best_x, t(replicate(N_THETA - 1L,
               rnorm(length(e$best_x), e$best_x, e$sigma_ppm * 0.3))))
    TH <- t(apply(X, 1, function(x) { p <- e$to_original(x); p["sigma_init"] <- 1; p }))
    saveRDS(TH, th_file)
  } else TH <- readRDS(th_file)
  ll_data <- apply(TH, 1, function(th) {
    th <- th[e$FREE_NAMES]
    x  <- e$to_unconstrained(th)
    e$ll_fn(x) - e$log_jacobian(x, e$to_original(x))
  })
  saveRDS(list(ll = ll_data, n_free = length(e$FREE_NAMES), run_id = e$RUN_ID,
               model_name = e$MODEL_NAME),
          file.path(OUTDIR, sprintf("%s_%s.rds", MODEL, arm)))
  quit(save = "no")
}

# === driver ===================================================================
MODELS <- if (length(args)) args else c("SP1", "TP2", "TP3", "Yasso07", "Yasso15", "Yasso20")
base_env <- c("HIKET_CORRELATED_LIK=1", "HIKET_SIGMA_TOT=0.800", "HIKET_SIGMA_1985_INFL=1",
              "SLURM_CPUS_PER_TASK=4", paste0("EQTEST_OUTDIR=", OUTDIR))
res <- list()
for (m in MODELS) {
  for (arm in c("prod", "eq")) {
    env <- c(base_env, paste0("EQTEST_WORKER=", m),
             paste0("HIKET_EQUILIBRIUM_INIT=", if (arm == "eq") "1" else "0"))
    rc <- system2("Rscript", c("--no-save", "doublechecks/eqinit_nesting_test.R"), env = env)
    if (rc != 0) stop(sprintf("%s %s worker failed (log in %s)", m, arm, OUTDIR))
  }
  p <- readRDS(file.path(OUTDIR, sprintf("%s_prod.rds", m)))
  q <- readRDS(file.path(OUTDIR, sprintf("%s_eq.rds",   m)))
  # remove the setup side effects (bundle + sanity PNG) of both arms
  for (r in list(p, q))
    unlink(c(Sys.glob(sprintf("Data/model_inputs/%s_inputs_%s.rds", r$model_name, r$run_id)),
             Sys.glob(sprintf("Calibration_real_data_transient/diagnostics/%s/*_%s.png",
                              r$model_name, r$run_id))))
  unlink(sprintf("Calibration_real_data_transient/diagnostics/%s_eqinit", m), recursive = TRUE)
  d <- max(abs(p$ll - q$ll))
  res[[m]] <- data.frame(model = m, n_free_prod = p$n_free, n_free_eq = q$n_free,
                         ll_prod_1 = p$ll[1], max_abs_diff = d,
                         pass = all(is.finite(p$ll)) && d < 1e-6 && q$n_free == p$n_free - 1L)
  print(res[[m]])
}
out <- do.call(rbind, res)
cat("\n=== NESTING TEST ===\n"); print(out, row.names = FALSE)
if (!all(out$pass)) quit(status = 1L)
