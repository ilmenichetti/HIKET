# =============================================================================
# run_ablation.R   (2026-08-04)
#
# Short-chain MCMC ablations for the LIKELIHOOD-SIDE factors that the forward-only
# test (c3_preinit_shape_ablation.R) cannot reach.
#
# Why this is needed. Six things changed at once since the last posteriors
# (20260710_*): the SOC target (2026-08-04), and C1/C2/C3/C4b/C5 (2026-07-16).
# c3_preinit_shape_ablation.R attributes C3 cleanly because C3 changes the FORWARD
# run, so it can be isolated at fixed parameters. C5 cannot: it changes only the
# observation SD, so it moves nothing until the sampler responds to it. Same for
# the SOC target. Those need actual (short) chains.
#
# Design. Each config re-runs the REAL calibration script with one factor varied
# through an environment variable; defaults are unchanged when unset, so a config
# with no overrides reproduces production behaviour exactly.
#
#   HIKET_SIGMA_1985_INFL : C5 strength (production 2.0; 1.0 = C5 off)
#   HIKET_PREINIT_LINEAR  : 1 = C3 off (old linear ramp)
#   HIKET_N_CHAINS / _N_ITER / _N_BURNIN : short-run settings
#
# These are ABLATIONS, not production posteriors: fewer chains and iterations,
# enough to compare posterior locations of the well-identified parameters
# (sigma_init, sigma_input) and the campaign-level fit, not to publish.
#
# Usage:
#   Rscript doublechecks/run_ablation.R [MODEL] [N_ITER] [N_CHAINS]
#   e.g.  Rscript doublechecks/run_ablation.R TP2 6000 3
# =============================================================================

args    <- commandArgs(trailingOnly = TRUE)
MODEL   <- if (length(args) >= 1) args[[1]] else "TP2"
N_ITER  <- if (length(args) >= 2) args[[2]] else "6000"
N_CHAIN <- if (length(args) >= 3) args[[3]] else "3"
# Optional 4th arg: comma-separated config names to run (default: all). Lets a
# single config be added to a model whose suite already ran, without repeating it.
ONLY    <- if (length(args) >= 4) strsplit(args[[4]], ",")[[1]] else NULL
stopifnot(MODEL %in% c("SP1","TP2","TP3","Yasso07","Yasso15","Yasso20"))

script <- sprintf("Calibration_real_data_transient/run_%s_transient_calibration.R", MODEL)
if (!file.exists(script)) stop("Run from the project root.")

OUT <- file.path("doublechecks", "ablation_logs")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# --- the configs -------------------------------------------------------------
# A0 is the reference (production settings, short chains). Every other config
# differs from A0 in exactly ONE factor, so differences are attributable.
configs <- list(
  A0_reference   = c(),                                   # C5 = 2.0, C3 on
  A1_C5_off      = c(HIKET_SIGMA_1985_INFL = "1.0"),      # trust 1985 fully
  A2_C5_1.5      = c(HIKET_SIGMA_1985_INFL = "1.5"),
  A3_C5_3.0      = c(HIKET_SIGMA_1985_INFL = "3.0"),
  A4_C3_off      = c(HIKET_PREINIT_LINEAR  = "1"),        # old linear ramp
  # A5 is not a weighting variant: it WITHHOLDS the 1985 campaign entirely, so the
  # model's 1985 prediction becomes out-of-sample. This is the only config that can
  # test C5's PREMISE (are the VMI8 stocks really low?) rather than its consequences.
  # Score it with ablation_fit_by_campaign.R, which re-derives the full obs set and
  # therefore scores 1985 as a held-out campaign.
  A5_LOCO_no1985 = c(HIKET_DROP_CAMPAIGN = "1985")
)

if (!is.null(ONLY)) {
  miss <- setdiff(ONLY, names(configs))
  if (length(miss)) stop("Unknown config(s): ", paste(miss, collapse = ", "))
  configs <- configs[ONLY]
}

base_env <- c(HIKET_N_ITER   = N_ITER,
              HIKET_N_CHAINS = N_CHAIN,
              HIKET_N_BURNIN = as.character(max(200, floor(as.numeric(N_ITER) / 5))))

cat(sprintf("\n=== ABLATION SUITE :: %s  (%s iter x %s chains, %d configs) ===\n",
            MODEL, N_ITER, N_CHAIN, length(configs)))

results <- list()
for (nm in names(configs)) {
  envv <- c(base_env, configs[[nm]])
  logf <- file.path(OUT, sprintf("%s_%s.log", MODEL, nm))
  cat(sprintf("\n--- %-16s %s\n", nm,
              if (length(configs[[nm]]))
                paste(names(configs[[nm]]), unlist(configs[[nm]]), sep = "=", collapse = " ")
              else "(production defaults)"))
  t0 <- Sys.time()
  st <- system2("Rscript", script, env = paste0(names(envv), "=", envv),
                stdout = logf, stderr = logf)
  mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  # RUN_ID is echoed in the script banner: "Run: YYYYmmdd_HHMMSS"
  lg  <- readLines(logf, warn = FALSE)
  rid <- sub(".*Run:\\s*([0-9_]+).*", "\\1", grep("\\|\\s*Run:", lg, value = TRUE)[1])
  cat(sprintf("    exit %d | %.1f min | RUN_ID %s\n", st, mins, rid))
  results[[nm]] <- data.frame(config = nm, exit = st, minutes = round(mins, 1),
                              run_id = rid, log = logf, stringsAsFactors = FALSE)
}

res    <- do.call(rbind, results)
idx_f  <- file.path(OUT, sprintf("%s_ablation_index.csv", MODEL))
# A partial run (ONLY set) must EXTEND the index, not replace it -- otherwise
# adding one config to a finished suite would erase the other four.
if (file.exists(idx_f)) {
  prev <- read.csv(idx_f, stringsAsFactors = FALSE)
  prev <- prev[!(prev$config %in% res$config), , drop = FALSE]
  common <- union(names(prev), names(res))
  for (nm in setdiff(common, names(prev))) prev[[nm]] <- NA
  for (nm in setdiff(common, names(res)))  res[[nm]]  <- NA
  res <- rbind(prev[, common, drop = FALSE], res[, common, drop = FALSE])
  res <- res[order(res$config), , drop = FALSE]
}
write.csv(res, idx_f, row.names = FALSE)
cat("\n=== SUITE DONE ===\n"); print(res, row.names = FALSE)
cat(sprintf("\nIndex: %s\n", file.path(OUT, sprintf("%s_ablation_index.csv", MODEL))))
cat("Compare posteriors with doublechecks/summarise_ablation.R\n")
