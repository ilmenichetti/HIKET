# =============================================================================
# test_correlated_likelihood.R   (2026-08-19)
#
# THE NON-NEGOTIABLE GATE. With every tau set to 0, Sigma collapses to
# sigma_e^2 I and the correlated likelihood MUST reproduce the independent one
# EXACTLY -- not approximately, not up to a constant. If it does not, the
# implementation is wrong and no cluster time should be spent.
#
# LOCAL UNIT TEST, no cluster. It runs the REAL ll_fn from the REAL calibration
# script at the REAL data, so it exercises the production code path rather than
# a reimplementation of it.
#
# ⚠ EACH CONFIGURATION RUNS IN ITS OWN SUBPROCESS, and that is not fussiness.
#   The engine's `source()` of correlated_likelihood.R is nested inside a
#   source(local = e), and a nested source() defaults to globalenv() -- so two
#   configurations in one session silently overwrite each other's switches. The
#   first version of this test "failed" for exactly that reason. Subprocesses
#   also reproduce production, where every run is a fresh R process.
#
# Checks:
#   1. tau = 0 identity at five parameter vectors            (THE GATE)
#   2. the incompatibility guards actually fire
#   3. the real tau's give a finite, different value
#
# Usage:  Rscript doublechecks/test_correlated_likelihood.R [MODEL]   (default SP1)
# =============================================================================

M <- local({ a <- commandArgs(trailingOnly = TRUE); if (length(a)) a[[1]] else "SP1" })
WORKER <- "doublechecks/_corr_lik_worker.R"

run <- function(env) {
  out <- suppressWarnings(system2("Rscript", c(WORKER, M), env = env,
                                  stdout = TRUE, stderr = TRUE))
  st  <- attr(out, "status")
  v   <- grep("^RESULT ", out, value = TRUE)
  list(ok = is.null(st) || st == 0, vals = as.numeric(sub("^RESULT \\d+ ", "", v)),
       log = out)
}

ZERO <- c("HIKET_CORRELATED_LIK=1", "HIKET_SIGMA_1985_INFL=1",
          "HIKET_TAU_R=0", "HIKET_TAU_P=0", "HIKET_TAU_C_1985=0",
          "HIKET_TAU_C_BASE=0", "HIKET_SIGMA_TOT=0.800")
INDEP <- c("HIKET_SIGMA_1985_INFL=1", "HIKET_SIGMA_TOTAL=0.800")
REAL  <- c("HIKET_CORRELATED_LIK=1", "HIKET_SIGMA_1985_INFL=1")

cat(sprintf("\n=== correlated likelihood: unit tests on %s ===\n\n", M))

# ---- 1. THE GATE -----------------------------------------------------------
cat("1. tau = 0 must reproduce the independent likelihood EXACTLY\n\n")
a <- run(INDEP); b <- run(ZERO)
if (!a$ok || !length(a$vals)) { cat(tail(a$log, 25), sep = "\n"); stop("independent build failed") }
if (!b$ok || !length(b$vals)) { cat(tail(b$log, 25), sep = "\n"); stop("tau=0 build failed") }

# A rejected draw gives -Inf in BOTH paths; that is agreement, not a difference.
# abs(-Inf - -Inf) is NaN, so the two cases are separated explicitly -- one path
# rejecting while the other does not is a HARD failure, not a tolerance question.
both_rej <- !is.finite(a$vals) & !is.finite(b$vals) & a$vals == b$vals
one_rej  <- xor(is.finite(a$vals), is.finite(b$vals))
d <- ifelse(both_rej, 0, abs(a$vals - b$vals))

cat(sprintf("%-6s %20s %20s %14s\n", "draw", "independent", "correlated (tau=0)", "difference"))
for (k in seq_along(d))
  cat(sprintf("%-6d %20.10f %20.10f %14s\n", k, a$vals[k], b$vals[k],
              if (both_rej[k]) "both -Inf" else sprintf("%.2e", d[k])))
if (any(one_rej))
  stop("draw(s) ", paste(which(one_rej), collapse = ", "),
       ": one path rejected and the other did not -- the guards have diverged.")
TOL <- 1e-8
cat(sprintf("\n   %s worst |difference| = %.3e  (tolerance %.0e)\n\n",
            if (max(d) < TOL) "PASS --" else "*** FAIL ***", max(d), TOL))
if (max(d) >= TOL)
  stop("tau = 0 does NOT reproduce the independent likelihood. Do not run anything.")

# ---- 2. THE GUARDS ---------------------------------------------------------
cat("2. incompatibility guards -- each MUST refuse\n\n")
guard <- function(label, extra, expect) {
  r  <- run(c(ZERO, extra))
  hit <- any(grepl(expect, r$log, fixed = TRUE))
  cat(sprintf("   %-24s %s\n", label, if (!r$ok && hit) "PASS" else "*** FAIL ***"))
  !r$ok && hit
}
g <- c(guard("Student-t stacking",    "HIKET_LIK_DF=6",          "do NOT stack"),
       guard("SIGMA_1985_INFL != 1",  "HIKET_SIGMA_1985_INFL=2", "requires SIGMA_1985_INFL = 1"),
       guard("multiplicative normal", "HIKET_LOGNORMAL_LIK=0",   "requires the LOG-NORMAL"))
if (!all(g)) stop("a guard failed to fire.")

# ---- 3. THE REAL SETTINGS --------------------------------------------------
cat("\n3. the adopted tau's (R 0.117 | P 0.396 | C 0.06/0.03 | sigma_e 0.685)\n\n")
r <- run(REAL)
if (!r$ok || !length(r$vals)) { cat(tail(r$log, 25), sep = "\n"); stop("real-tau build failed") }
for (l in grep("^SETUP ", r$log, value = TRUE))
  cat("   ", sub("^SETUP\\s*", "", l), "\n", sep = "")
cat(sprintf("   independent  %16.4f\n   correlated   %16.4f    (difference %+.1f)\n",
            a$vals[1], r$vals[1], r$vals[1] - a$vals[1]))
if (!is.finite(r$vals[1])) stop("the correlated likelihood is not finite at best_x.")
tm <- function(z) as.numeric(sub("^TIMING ", "", grep("^TIMING ", z$log, value = TRUE)[1]))
cat(sprintf("\n   cost per evaluation: %.3f s independent -> %.3f s correlated (%+.0f%%)\n",
            tm(a), tm(r), 100 * (tm(r) / tm(a) - 1)))

cat("\n   ⚠ That difference is NOT a model comparison -- the normalising constant\n")
cat("     moved. Compare runs on RMSE distributions, never on log-likelihood.\n")
cat("\nALL TESTS PASSED\n")
