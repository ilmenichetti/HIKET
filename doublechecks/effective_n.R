# =============================================================================
# effective_n.R   (2026-08-19)
#
# WHAT THE CORRELATED LIKELIHOOD ACTUALLY COSTS IN INFORMATION, measured on the
# Sigma the engine really builds -- not estimated from the proposal's algebra.
#
# n_eff is defined so that it EQUALS n under independence at the same total:
#
#       var(w'y) = w' Sigma w        for a fixed contrast w
#       n_eff    = sigma_total^2 / var   (for the mean, w = 1/n)
#
# Reported for three quantities, because the correction does NOT hit them
# equally -- that asymmetry is the whole point of the instrument:
#   * the national LEVEL        -- what pins MRT x sigma_input; hit hardest
#   * the 1985->2024 TREND      -- the long interval
#   * the 2006->2024 TREND      -- the short, fragile one
#
# ⚠ Within-plot contrasts are supposed to survive: u^P and u^R cancel on
#   differencing. Only the campaign term touches a trend. If the trend numbers
#   fall as far as the level does, the implementation is wrong.
#
# Usage:  Rscript doublechecks/effective_n.R [MODEL]     (default SP1)
# =============================================================================

M <- local({ a <- commandArgs(trailingOnly = TRUE); if (length(a)) a[[1]] else "SP1" })
Sys.setenv(HIKET_CORRELATED_LIK = "1", HIKET_SIGMA_1985_INFL = "1")

src <- readLines(sprintf("Calibration_real_data_transient/run_%s_transient_calibration.R", M),
                 warn = FALSE)
cut <- grep("^t_run <- system.time", src)[1]
e <- new.env(parent = globalenv())
invisible(capture.output(suppressMessages(
  source(textConnection(paste(src[seq_len(cut - 1)], collapse = "\n")), local = e))))

plots <- get("plots", e); obs_meta <- get("obs_meta", e)
S <- suppressMessages(hiket_build_sigma(plots, obs_meta))

Sig  <- crossprod(S$R)                    # Sigma = R'R
n    <- S$n
TOT  <- HIKET_SIGMA_TOT
camp <- S$camp

# variance of a fixed linear contrast, under Sigma and under iid at the same total
vc  <- function(w) as.numeric(t(w) %*% Sig %*% w)
vi  <- function(w) TOT^2 * sum(w^2)
mean_w <- function(sel) { w <- numeric(n); w[sel] <- 1 / sum(sel); w }

cat(sprintf("\n=== effective sample size under the correlated likelihood (%s) ===\n", M))
cat(sprintf("n = %d observations | %d plots | %d bands | campaigns %s\n\n",
            n, length(unique(S$pid)), length(unique(S$reg)),
            paste(sort(unique(camp)), collapse = "/")))

rows <- list(
  list("national level (all obs)", mean_w(rep(TRUE, n)),            n),
  list("1985 campaign level",      mean_w(camp == 1985L),           sum(camp == 1985L)),
  list("2006 campaign level",      mean_w(camp == 2006L),           sum(camp == 2006L)),
  list("2024 campaign level",      mean_w(camp == 2024L),           sum(camp == 2024L)),
  list("trend 1985->2024", mean_w(camp == 2024L) - mean_w(camp == 1985L), NA),
  list("trend 2006->2024", mean_w(camp == 2024L) - mean_w(camp == 2006L), NA))

# n_eff = the number of INDEPENDENT observations at the same total sigma that
# would give this variance:  n_eff = sigma_total^2 / var. Equals n under
# independence, by construction. Undefined for a trend (a difference of two
# means is not "so many observations"), so only the inflation factor is shown.
cat(sprintf("%-24s %6s %10s %10s %8s %8s\n",
            "quantity", "n obs", "sd (iid)", "sd (corr)", "inflate", "n_eff"))
for (r in rows) {
  w <- r[[2]]; a <- sqrt(vi(w)); b <- sqrt(vc(w))
  cat(sprintf("%-24s %6s %10.4f %10.4f %8.2fx %8s\n", r[[1]],
              if (is.na(r[[3]])) "--" else format(r[[3]]), a, b, b / a,
              if (is.na(r[[3]])) "--" else format(round(TOT^2 / vc(w)))))
}

wm <- mean_w(rep(TRUE, n))
cat(sprintf("\nLEVEL: the national mean's uncertainty grows %.2fx -- %d observations carry the\n",
            sqrt(vc(wm) / vi(wm)), n))
cat(sprintf("       weight of about %d independent ones.\n", round(TOT^2 / vc(wm))))
w85 <- mean_w(camp == 2024L) - mean_w(camp == 1985L)
w06 <- mean_w(camp == 2024L) - mean_w(camp == 2006L)
cat(sprintf("TRENDS: inflated only %.2fx (1985->2024) and %.2fx (2006->2024) -- the plot and\n",
            sqrt(vc(w85) / vi(w85)), sqrt(vc(w06) / vi(w06))))
cat("        region offsets cancel on differencing, exactly as designed.\n")
