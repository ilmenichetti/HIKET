# =============================================================================
# test_flux_pair_transform.R
#   Verifies the new `flux_pair` bounded transform in calibration_engine.R:
#     (1) round-trip  to_original(to_unconstrained(p)) == p
#     (2) fluxes stay inside the window [a,b] for ALL unconstrained inputs
#     (3) analytical log-Jacobian == numerical log|det d(si,sinit)/d(x1,x2)|
#   Run:  Rscript doublechecks/test_flux_pair_transform.R
# =============================================================================

source("Calibration_real_data/calibration_engine.R")

a <- 0.05; b <- 8.7; Jb <- 2.8            # window (tC/ha/yr) + mean-litter bridge

param_spec <- list(
  list(names = c("sigma_input", "sigma_init"),
       type = "flux_pair", window = c(a, b), J_bar = Jb)
)
tr <- build_transforms(param_spec)

ok <- TRUE
report <- function(tag, pass) {
  cat(sprintf("  [%s] %s\n", if (pass) "PASS" else "FAIL", tag))
  if (!pass) ok <<- FALSE
}

# --- (1) round-trip over a grid of physical (sigma_input, sigma_init) ---------
# choose sigma pairs whose implied fluxes are safely inside (a,b)
si_grid    <- c(0.30, 1.00, 2.50, 3.00)          # F_now = si*Jb in (0.84, 8.4)
sinit_grid <- c(0.20, 0.70, 1.00, 1.50)          # F_1917 = sinit*F_now
max_rt <- 0
for (si in si_grid) for (sinit in sinit_grid) {
  F_now  <- si * Jb
  F_1917 <- sinit * F_now
  if (F_now <= a || F_now >= b || F_1917 <= a || F_1917 >= b) next  # skip out-of-window
  p  <- c(sigma_input = si, sigma_init = sinit)
  x  <- tr$to_unconstrained(p)
  p2 <- tr$to_original(x)
  max_rt <- max(max_rt, abs(p2["sigma_input"] - si), abs(p2["sigma_init"] - sinit))
}
report(sprintf("round-trip max error = %.2e", max_rt), max_rt < 1e-10)

# --- (2) fluxes in-window for arbitrary unconstrained draws -------------------
set.seed(1)
in_window <- TRUE
for (i in 1:2000) {
  x <- c(sigma_input = rnorm(1, 0, 6), sigma_init = rnorm(1, 0, 6))  # extreme draws
  p <- tr$to_original(x)
  F_now  <- p["sigma_input"] * Jb
  F_1917 <- p["sigma_init"] * F_now
  if (F_now  <= a || F_now  >= b) in_window <- FALSE
  if (F_1917 <= a || F_1917 >= b) in_window <- FALSE
}
report("both fluxes inside (a,b) for 2000 extreme draws", in_window)

# --- (3) analytical vs numerical log-Jacobian --------------------------------
# numerical |det| of the map (x1,x2) -> (sigma_input, sigma_init) via central diff
num_logjac <- function(x, eps = 1e-6) {
  base <- tr$to_original(x)
  J <- matrix(0, 2, 2)
  for (j in 1:2) {
    xp <- x; xm <- x
    xp[j] <- xp[j] + eps; xm[j] <- xm[j] - eps
    dp <- (tr$to_original(xp) - tr$to_original(xm)) / (2 * eps)
    J[, j] <- c(dp["sigma_input"], dp["sigma_init"])
  }
  log(abs(det(J)))
}
max_jac_err <- 0
for (x1 in c(-3, -1, 0, 1, 2)) for (x2 in c(-2, 0, 1.5)) {
  x <- c(sigma_input = x1, sigma_init = x2)
  p <- tr$to_original(x)
  ana <- tr$log_jacobian(x, p)
  num <- num_logjac(x)
  max_jac_err <- max(max_jac_err, abs(ana - num))
}
report(sprintf("analytical vs numerical log-Jac max error = %.2e", max_jac_err),
       max_jac_err < 1e-5)

cat(if (ok) "\nALL FLUX_PAIR TESTS PASSED\n" else "\nFLUX_PAIR TESTS FAILED\n")
quit(status = if (ok) 0 else 1)
