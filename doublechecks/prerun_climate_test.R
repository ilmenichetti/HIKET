# =============================================================================
# prerun_climate_test.R -- gate for the pre-run climate window (2026-10-12)
# -----------------------------------------------------------------------------
# Three arms per model, each sourcing the REAL calibration script up to the MCMC
# launch in its own R process (HIKET_PRERUN_CLIM_YEARS is read at source time):
#   off   : unset                      -> production
#   ident : HIKET_PRERUN_CLIM_YEARS=1985-2004 -> the new pathway (climate_monthly.csv,
#           mapper, registry, precip_mean) on the production window: MUST equal off
#   on    : HIKET_PRERUN_CLIM_YEARS=1961-1984 -> the proposed change
# Reported: the data log-likelihood at defaults + 3 prior draws (off == ident to
# 1e-9 is the PASS criterion), the pre-run climate modifier and the 1985 starting
# stock at the CURRENT production posterior median (how far the change moves it).
# Production config set here (correlated likelihood, sigma_tot 0.800, 1985 infl 1).
#
# Usage (project root):  Rscript doublechecks/prerun_climate_test.R [MODEL ...]
# Side effects removed at the end: the input bundle + sanity PNG each setup writes.
# =============================================================================

args   <- commandArgs(trailingOnly = TRUE)
WORKER <- Sys.getenv("PCTEST_WORKER", "")
OUTDIR <- Sys.getenv("PCTEST_OUTDIR", file.path(tempdir(), "prerun_climate"))
dir.create(OUTDIR, showWarnings = FALSE)
N_THETA <- 4L

# === worker ===================================================================
if (nzchar(WORKER)) {
  MODEL <- WORKER; ARM <- Sys.getenv("PCTEST_ARM")
  src   <- readLines(sprintf("Calibration_real_data_transient/run_%s_transient_calibration.R", MODEL))
  cutix <- grep("^t_run <- system.time\\(\\{", src)[1]
  e <- new.env(parent = globalenv())
  sink(file.path(OUTDIR, sprintf("%s_%s.log", MODEL, ARM)))
  source(textConnection(paste(src[seq_len(cutix - 1L)], collapse = "\n")), local = e)
  sink()
  # the engine functions handed to make_likelihood (as attribution_decomposition.R)
  st <- grep("ll_fn <- make_likelihood", src)[1]; op <- 0L
  for (i in seq(st, length(src))) {
    ch <- strsplit(src[i], "")[[1]]; op <- op + sum(ch == "(") - sum(ch == ")")
    if (op == 0L) { en <- i; break }
  }
  ml <- as.list(str2lang(sub("^\\s*ll_fn\\s*<-\\s*", "", paste(src[seq(st, en)], collapse = "\n"))))[-1]
  f  <- lapply(setNames(nm = c("assemble_params", "compute_xi_mean", "steady_state")),
               function(k) eval(ml[[k]], envir = e))
  n_ss <- eval(ml[["steady_state_n"]], envir = e)

  set.seed(2025)
  X  <- rbind(e$best_x, t(replicate(N_THETA - 1L,
              rnorm(length(e$best_x), e$best_x, e$sigma_ppm * 0.3))))
  ll <- apply(X, 1, function(x) e$ll_fn(x))

  # production posterior median -> pre-run xi and 1985 starting stock per plot
  source("manuscript/figures/run_ids.R", local = TRUE)
  post <- readRDS(sprintf("Calibration_real_data_transient/runs/%s_posterior_%s.rds", MODEL, RID[[MODEL]]))
  pm   <- apply(post, 2, median)
  mp   <- f$assemble_params(setNames(as.numeric(pm), colnames(post)))
  pl   <- names(e$litter_means)
  st   <- t(vapply(pl, function(pid) {
    clim <- e$climate_by_plot[[pid]]
    xi   <- f$compute_xi_mean(hiket_ss_clim(pid, clim, min(n_ss, nrow(clim))), mp)
    C0   <- f$steady_state(mp, e$litter_means[[pid]], xi)
    c(xi = unname(unlist(xi))[1], C1985 = sum(C0))
  }, numeric(2)))
  saveRDS(list(ll = ll, st = st, label = hiket_prerun_label(), run_id = e$RUN_ID,
               model_name = e$MODEL_NAME, prod_rid = RID[[MODEL]]),
          file.path(OUTDIR, sprintf("%s_%s.rds", MODEL, ARM)))
  quit(save = "no")
}

# === driver ===================================================================
MODELS <- if (length(args)) args else c("SP1", "TP2", "TP3", "Yasso07", "Yasso15", "Yasso20")
ARMS   <- c(off = "", ident = "1985-2004", on = "1961-1984")
base_env <- c("HIKET_CORRELATED_LIK=1", "HIKET_SIGMA_TOT=0.800", "HIKET_SIGMA_1985_INFL=1",
              "HIKET_EQUILIBRIUM_INIT=0", "SLURM_CPUS_PER_TASK=4", paste0("PCTEST_OUTDIR=", OUTDIR))
res <- list()
for (m in MODELS) {
  R <- list()
  for (a in names(ARMS)) {
    env <- c(base_env, paste0("PCTEST_WORKER=", m), paste0("PCTEST_ARM=", a),
             paste0("HIKET_PRERUN_CLIM_YEARS=", ARMS[[a]]))
    rc <- system2("Rscript", c("--no-save", "doublechecks/prerun_climate_test.R"), env = env)
    if (rc != 0) stop(sprintf("%s %s worker failed (log in %s)", m, a, OUTDIR))
    R[[a]] <- readRDS(file.path(OUTDIR, sprintf("%s_%s.rds", m, a)))
    unlink(c(Sys.glob(sprintf("Data/model_inputs/%s_inputs_%s.rds", R[[a]]$model_name, R[[a]]$run_id)),
             Sys.glob(sprintf("Calibration_real_data_transient/diagnostics/%s/*_%s.png",
                              R[[a]]$model_name, R[[a]]$run_id))))
  }
  d_id <- max(abs(R$off$ll - R$ident$ll))
  ratio <- R$on$st[, "C1985"] / R$off$st[, "C1985"]
  res[[m]] <- data.frame(
    model = m, labels = paste(R$off$label, R$ident$label, R$on$label, sep = " | "),
    ll_off = R$off$ll[1], ident_maxdiff = d_id, dll_on = R$on$ll[1] - R$off$ll[1],
    xi_off = median(R$off$st[, "xi"]), xi_on = median(R$on$st[, "xi"]),
    C1985_off = mean(R$off$st[, "C1985"]), C1985_on = mean(R$on$st[, "C1985"]),
    C_ratio_med = median(ratio), C_ratio_p05 = quantile(ratio, .05), C_ratio_p95 = quantile(ratio, .95),
    pass = all(is.finite(c(R$off$ll, R$on$ll))) && d_id < 1e-9, row.names = NULL)
  print(res[[m]])
}
out <- do.call(rbind, res)
cat("\n=== PRE-RUN CLIMATE TEST (stock at the production posterior median,", R$off$prod_rid, ") ===\n")
print(out, row.names = FALSE, digits = 4)
saveRDS(out, "doublechecks/figures/prerun_climate_test.rds")
if (!all(out$pass)) quit(status = 1L)
