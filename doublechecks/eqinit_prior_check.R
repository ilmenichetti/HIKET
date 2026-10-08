# =============================================================================
# eqinit_prior_check.R -- same sigma_input prior in both arms (2026-10-08)
# -----------------------------------------------------------------------------
# The prior is Gaussian in SAMPLING space and the engine adds log_jacobian() to
# the log-likelihood, so the effective prior on x is N(x) * |J(x)|. In the
# transient arm (flux_pair) |J| = g'(x1) g'(x2) / (J_bar * F_now); integrating
# x2 out leaves N(x1) g'(x1) / F_now for sigma_input. The equilibrium arm's
# flux_now must reproduce exactly that marginal -- including the 1/F_now.
#
# Checks, per model, the effective sigma_input density of the two arms on a grid
# (transient: x2 integrated numerically), plus the naive alternative WITHOUT the
# -log(F_now) term, to show what the fix prevents. No data, no model runs.
# Usage: Rscript doublechecks/eqinit_prior_check.R   (project root)
# =============================================================================

source("./Calibration_real_data/calibration_engine.R")

J_BAR <- 2.511   # cross-plot mean litter (units bridge); any positive value works
specs <- c(SP1 = "SP1", TP2 = "TP2", TP3 = "TP3",
           Yasso07 = "YASSO07", Yasso15 = "YASSO15", Yasso20 = "YASSO20")

tv <- function(f, g, dx) 0.5 * sum(abs(f - g)) * dx   # total-variation distance

out <- lapply(names(specs), function(m) {
  source(sprintf("./Prior_specs/%s_priors.R", m))
  P   <- specs[[m]]
  win <- get(paste0(P, "_INPUT_FLUX_WINDOW"))
  fd  <- get(paste0(P, "_FREE_DEFAULTS"))[c("sigma_input", "sigma_init")]
  sd  <- get(paste0(P, "_SIGMA_PPM"))[c("sigma_input", "sigma_init")]

  tr_pair <- build_transforms(list(list(names = c("sigma_input", "sigma_init"),
                                        type = "flux_pair", window = win, J_bar = J_BAR)))
  tr_now  <- build_transforms(list(list(names = "sigma_input",
                                        type = "flux_now", window = win, J_bar = J_BAR)))
  mu <- tr_pair$to_unconstrained(fd)
  stopifnot(abs(tr_now$to_unconstrained(fd["sigma_input"]) - mu[1]) < 1e-12)   # same centre

  x1 <- seq(mu[1] - 6 * sd[1], mu[1] + 6 * sd[1], length.out = 801); dx1 <- diff(x1[1:2])
  x2 <- seq(mu[2] - 6 * sd[2], mu[2] + 6 * sd[2], length.out = 801); dx2 <- diff(x2[1:2])

  # transient arm: effective density on (x1, x2), x2 integrated out
  f_pair <- vapply(x1, function(a) {
    sum(vapply(x2, function(b) {
      x <- c(sigma_input = a, sigma_init = b)
      exp(dnorm(a, mu[1], sd[1], log = TRUE) + dnorm(b, mu[2], sd[2], log = TRUE) +
          tr_pair$log_jacobian(x, tr_pair$to_original(x)))
    }, numeric(1))) * dx2
  }, numeric(1))
  # equilibrium arm (flux_now) and the naive version without -log(F_now)
  f_now <- vapply(x1, function(a) {
    x <- c(sigma_input = a); p <- tr_now$to_original(x)
    exp(dnorm(a, mu[1], sd[1], log = TRUE) + tr_now$log_jacobian(x, p))
  }, numeric(1))
  f_naive <- vapply(x1, function(a) {
    F <- bounded_fwd(a, win[1], win[2])
    exp(dnorm(a, mu[1], sd[1], log = TRUE) + log_jac_bounded(F, win[1], win[2]))
  }, numeric(1))
  nrm <- function(f) f / (sum(f) * dx1)
  f_pair <- nrm(f_pair); f_now <- nrm(f_now); f_naive <- nrm(f_naive)
  med <- function(f) tr_now$to_original(c(sigma_input = x1[which(cumsum(f) * dx1 >= 0.5)[1]]))
  data.frame(model = m,
             tv_flux_now = tv(f_pair, f_now, dx1),
             tv_naive    = tv(f_pair, f_naive, dx1),
             median_si_transient = med(f_pair),
             median_si_flux_now  = med(f_now),
             median_si_naive     = med(f_naive))
})
out <- do.call(rbind, out)
cat("\n=== sigma_input effective prior: transient vs equilibrium arm ===\n")
print(out, row.names = FALSE, digits = 4)
pass <- all(out$tv_flux_now < 1e-3)
cat(sprintf("\nflux_now matches the transient marginal: %s\n", if (pass) "PASS" else "FAIL"))
if (!pass) quit(status = 1L)
