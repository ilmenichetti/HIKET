# =============================================================================
# summarise_ablation.R   (2026-08-04)
#
# Side-by-side comparison of the ablation posteriors produced by run_ablation.R.
#
# Reads doublechecks/ablation_logs/<MODEL>_ablation_index.csv, loads each config's
# posterior matrix from Calibration_real_data_transient/runs/, and prints one
# table per parameter: median [2.5%, 97.5%] for every config, plus the shift
# relative to the A0 reference.
#
# What to look at first:
#   sigma_init  -- the pre-run scaler. The forward test showed it out-leverages
#                  C3 by ~an order of magnitude on the 1985 stock, so this is
#                  where C5 should show up: down-weighting 1985 frees sigma_init.
#   sigma_input -- the input multiplier; the flux_pair bound should keep it
#                  physical (< ~3.6) in every config.
#
# Usage:  Rscript doublechecks/summarise_ablation.R [MODEL]
# =============================================================================

args  <- commandArgs(trailingOnly = TRUE)
MODEL <- if (length(args) >= 1) args[[1]] else "TP2"

idx_f <- file.path("doublechecks", "ablation_logs",
                   sprintf("%s_ablation_index.csv", MODEL))
if (!file.exists(idx_f)) stop("No ablation index for ", MODEL, " -- run run_ablation.R first.")
idx <- read.csv(idx_f, stringsAsFactors = FALSE)

load_post <- function(run_id) {
  # quarantine first (see quarantine_ablation_runs.R), then the production dir
  cand <- c(file.path("doublechecks", "ablation_runs",
                      sprintf("%s_posterior_%s.rds", MODEL, run_id)),
            file.path("Calibration_real_data_transient", "runs",
                      sprintf("%s_posterior_%s.rds", MODEL, run_id)))
  f <- cand[file.exists(cand)][1]
  if (is.na(f)) return(NULL)
  readRDS(f)
}

posts <- lapply(idx$run_id, load_post)
names(posts) <- idx$config
ok <- !vapply(posts, is.null, logical(1))
if (any(!ok)) {
  cat("Missing posteriors for:", paste(idx$config[!ok], collapse = ", "), "\n")
  cat("(a config whose chains failed will have exit != 0 in the index)\n\n")
}
posts <- posts[ok]
if (!length(posts)) stop("No posteriors could be loaded.")

pars <- colnames(posts[[1]])
ref  <- if ("A0_reference" %in% names(posts)) "A0_reference" else names(posts)[1]

cat(sprintf("\n=====================================================================\n"))
cat(sprintf("ABLATION POSTERIOR COMPARISON  ::  %s\n", MODEL))
cat(sprintf("reference config = %s ; shifts are (config - reference) medians\n", ref))
cat(sprintf("SHORT chains -- for comparing locations, NOT production posteriors\n"))
cat(sprintf("=====================================================================\n"))

for (p in pars) {
  cat(sprintf("\n%s\n", p))
  cat(sprintf("  %-16s %10s %10s %10s %10s\n",
              "config", "median", "2.5%", "97.5%", "vs ref"))
  med_ref <- median(posts[[ref]][, p])
  for (nm in names(posts)) {
    v <- posts[[nm]][, p]
    q <- quantile(v, c(.025, .975))
    cat(sprintf("  %-16s %10.4f %10.4f %10.4f %+10.4f\n",
                nm, median(v), q[1], q[2], median(v) - med_ref))
  }
}

# --- headline: the two parameters the ablations are really about -------------
cat(sprintf("\n=====================================================================\n"))
cat("HEADLINE\n")
for (p in intersect(c("sigma_init", "sigma_input"), pars)) {
  cat(sprintf("\n  %s (median by config):\n", p))
  for (nm in names(posts))
    cat(sprintf("    %-16s %8.4f\n", nm, median(posts[[nm]][, p])))
}
cat(sprintf("=====================================================================\n"))
