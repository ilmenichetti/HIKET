# =============================================================================
# correlated_likelihood.R   (2026-08-19)
#
# THE CORRELATED-ERROR LIKELIHOOD. Design settled 2026-08-18; full write-up in
# manuscript/HIKET_correlated_likelihood_proposal.tex (8 pp). Read that first --
# this file implements it, it does not argue for it.
#
#     log y_ij = log f_ij(theta) + u^R_r(i) + u^P_i + u^C_j + e_ij
#
# WHY. The current likelihood treats 1205 plot-years as independent. They are
# not: plot ICC 0.57-0.71, latitude-band means varying 5-6x more than iid allows,
# campaign levels that are not comparable. The national LEVEL is therefore
# measured far more precisely than the data support -- and the level is what pins
# the MRT x sigma_input ridge. 1205 observations carry the weight of about 260
# for the national mean.
#
# NOTHING IS ESTIMATED. All four variances are FIXED and the three offsets are
# MARGINALISED -- integrated out, never given a value -- leaving only their
# footprint in the covariance:
#
#     Sigma = sigma_e^2 I + tau_P^2 Z_P Z_P' + tau_R^2 Z_R Z_R' + tau_C^2 Z_C Z_C'
#
# CONSEQUENCE FOR THE IMPLEMENTATION: Sigma does not depend on theta, so it is
# CONSTANT. Factorise once (~12 MB dense at n = 1205), then each evaluation is a
# single triangular solve, ~1 ms. Woodbury / Sherman-Morrison are NOT needed.
# But the campaign term couples every observation to every other one in its
# campaign, so the per-plot sum becomes ONE GLOBAL SOLVE -- that is the only
# structural change to the engine.
#
# ⚠ THE NON-NEGOTIABLE CHECK: with all tau = 0, this must reproduce the current
#   log-likelihood EXACTLY. doublechecks/test_correlated_likelihood.R does that,
#   and it is a LOCAL UNIT TEST -- run it before any cluster time is spent.
#
# ⚠ THE VARIANCE IS SPLIT, NEVER ADDED. sigma_total = 0.800 is UNCHANGED;
#   sigma_e = 0.685 is what remains after the three shared terms are taken out.
#   Every observation's independent noise therefore FALLS (0.800 -> 0.685), which
#   is why the scheme is ~12% TIGHTER on the trend than the SIGMA_1985_INFL = 2
#   it replaces, while being 2.4x MORE sceptical of the 1985 LEVEL. Those are two
#   different questions with opposite answers -- see doublechecks/campaign_tau_bracket.R.
#
# ⚠ THE ICC RATIO TRAP. The observations' own ICC is 0.799. Transplanting that
#   RATIO onto a total of 0.800 gives sigma_e = 0.359, BELOW the models' within-plot
#   error -- recreating the overconfidence in the TREND instead of the level.
#   Carry ABSOLUTE standard deviations across, never ratios.
#
# ⚠ DOES NOT STACK WITH: Student-t (HIKET_LIK_DF), because a t does not decompose
#   into shared + independent Gaussian; and SIGMA_1985_INFL, because tau_C REPLACES
#   it (stacking gives an effective tau_C of 0.085, the broad scheme by accident).
#   Both are refused loudly below rather than silently combined.
#
# ⚠ LOG-LIKELIHOODS ARE NOT COMPARABLE across this change -- the normalising
#   constant moves. Compare runs on RMSE distributions (S14_rmse_posterior).
#
# SWITCH: HIKET_CORRELATED_LIK=1 enables. Default OFF, so production behaviour is
# unchanged until deliberately flipped. On Roihu it must be passed as
# SINGULARITYENV_HIKET_CORRELATED_LIK -- a bare export never reaches R inside the
# r-env container.
# =============================================================================

.envnum2 <- function(k, d) {
  v <- suppressWarnings(as.numeric(Sys.getenv(k, NA)))
  if (is.na(v)) d else v
}

.hiket_correlated_lik <- identical(Sys.getenv("HIKET_CORRELATED_LIK"), "1")

# --- the fixed variance components ------------------------------------------
# tau_R  latitude-band means, campaign-centred            MEASURED
# tau_P  2006-2024 pair covariance -- the only pair       MEASURED
#        untouched by the 1985 problems
# tau_C  prescribed; 1985 carries twice the others        ASSUMED (see below)
# sigma_e  the REMAINDER of the unchanged total           DERIVED
HIKET_TAU_R      <- .envnum2("HIKET_TAU_R",      0.117)
HIKET_TAU_P      <- .envnum2("HIKET_TAU_P",      0.396)
HIKET_TAU_C_1985 <- .envnum2("HIKET_TAU_C_1985", 0.060)
HIKET_TAU_C_BASE <- .envnum2("HIKET_TAU_C_BASE", 0.030)
HIKET_SIGMA_TOT  <- .envnum2("HIKET_SIGMA_TOT",  0.800)

# sigma_e is derived from the BASE campaign term, so 1985's extra scepticism sits
# on top rather than being paid for by the other campaigns. 1985's implied total
# is then 0.802 -- deliberate, and negligible.
.hiket_sigma_e <- local({
  v <- HIKET_SIGMA_TOT^2 - HIKET_TAU_R^2 - HIKET_TAU_P^2 - HIKET_TAU_C_BASE^2
  if (v <= 0)
    stop("correlated_likelihood: the shared terms exceed the total variance; ",
         "sigma_e would be imaginary. Split the total, never add to it.")
  sqrt(v)
})

# =============================================================================
# hiket_latitude_bands()
#
# Eight equal-count bands on ETRS northing, the grouping tau_R = 0.117 was
# measured on. Equal COUNT, not equal width: Finland's plots are not uniformly
# distributed in latitude and equal-width bands would put a handful of plots in
# the northernmost one, making its mean noise rather than signal.
# =============================================================================
hiket_latitude_bands <- function(y_etrs, n_band = 8L) {
  if (any(!is.finite(y_etrs)))
    stop("hiket_latitude_bands: non-finite northing for ",
         sum(!is.finite(y_etrs)), " plot(s); the regional term cannot be built.")
  br <- quantile(y_etrs, probs = seq(0, 1, length.out = n_band + 1L), names = FALSE)
  br[1] <- -Inf; br[length(br)] <- Inf
  setNames(as.integer(cut(y_etrs, breaks = br, labels = FALSE)), names(y_etrs))
}

# =============================================================================
# hiket_build_sigma()
#
# Builds the observation index and the Cholesky factor of Sigma, ONCE.
#
# ⚠ ROW ORDER IS LOAD-BEARING. The residual vector the likelihood assembles is
#   unlist(lapply(plots, ...)) -- plots in order, and within each plot its
#   obs_meta$soc_obs in order. This function walks the SAME two loops, so the
#   two orders cannot drift. Never reorder one without the other.
# =============================================================================
hiket_build_sigma <- function(plots, obs_meta) {
  pid  <- unlist(lapply(plots, function(p) rep(p, length(obs_meta[[p]]$soc_obs))))
  camp <- unlist(lapply(plots, function(p) as.integer(obs_meta[[p]]$campaign)))
  reg  <- unlist(lapply(plots, function(p)
                        rep(as.integer(obs_meta[[p]]$region),
                            length(obs_meta[[p]]$soc_obs))))
  if (length(camp) != length(pid) || any(is.na(camp)))
    stop("correlated_likelihood: obs_meta$campaign missing or wrong length. ",
         "Add `campaign = obs_plot$year` to obs_meta in the run script.")
  if (any(is.na(reg)))
    stop("correlated_likelihood: obs_meta$region missing. ",
         "Add `region = <band>` to obs_meta in the run script.")

  n <- length(pid)
  # tau_C is campaign-specific: 1985 carries twice the others.
  tc <- ifelse(camp == 1985L, HIKET_TAU_C_1985, HIKET_TAU_C_BASE)

  # Sigma = sigma_e^2 I + tau_P^2 [same plot] + tau_R^2 [same band] + tc_i tc_j [same campaign]
  S <- diag(.hiket_sigma_e^2, n, n)
  S <- S + HIKET_TAU_P^2 * outer(pid,  pid,  "==")
  S <- S + HIKET_TAU_R^2 * outer(reg,  reg,  "==")
  S <- S + outer(tc, tc, "*") * outer(camp, camp, "==")

  R <- tryCatch(chol(S), error = function(e)
    stop("correlated_likelihood: Sigma is not positive definite (", conditionMessage(e),
         "). Check that the tau's are non-negative and sigma_e > 0."))

  message(sprintf(paste0("[ERROR MODEL] CORRELATED likelihood: n = %d obs, %d plots, %d bands, %d campaigns\n",
                         "              tau_R %.3f | tau_P %.3f | tau_C %.3f (1985) / %.3f | sigma_e %.3f",
                         " | total %.3f"),
                  n, length(unique(pid)), length(unique(reg)), length(unique(camp)),
                  HIKET_TAU_R, HIKET_TAU_P, HIKET_TAU_C_1985, HIKET_TAU_C_BASE,
                  .hiket_sigma_e, HIKET_SIGMA_TOT))

  list(R = R, n = n, logdet = 2 * sum(log(diag(R))),
       pid = pid, camp = camp, reg = reg)
}

# =============================================================================
# hiket_corr_ll()
#
# Gaussian log-density of the log-residual vector under the constant Sigma.
#
# ⚠ MATCHES THE INDEPENDENT PATH'S CONVENTION: the existing log-normal branch is
#   dnorm() on log residuals and omits the -sum(log y) Jacobian (a constant in
#   theta). This does the same, so the tau = 0 identity is EXACT rather than
#   exact-up-to-a-constant.
# =============================================================================
hiket_corr_ll <- function(r, S) {
  z <- backsolve(S$R, r, transpose = TRUE)     # z = R'^{-1} r  =>  r' Sigma^{-1} r = z'z
  -0.5 * (sum(z * z) + S$logdet + S$n * log(2 * pi))
}

# =============================================================================
# hiket_total_sigma()
#
# The TOTAL observation+model error on the log scale, resolved the same way the
# likelihood resolves it, for use by the PREDICTIVE stage.
#
# ⚠ THE KEY FACT, and the reason the predictive stage needs almost no change:
#   the correlated likelihood SPLITS a fixed total, it does not add to it. So the
#   MARGINAL variance of a SINGLE observation is
#       tau_R^2 + tau_P^2 + tau_C^2 + sigma_e^2 = sigma_total^2
#   -- identical to the independent model. A posterior-predictive interval for
#   one plot-year is therefore the SAME width either way, and per-observation
#   coverage cannot distinguish the two error models. What the shared offsets
#   change is JOINT statements over many observations (the national level, a
#   campaign mean, a trend) -- see doublechecks/effective_n.R, where the level
#   inflates 2.18x while a single observation does not move at all.
#
# ⚠ This ALSO closes a long-standing defect that predates the correlated work:
#   the predictive scripts quoted "95% coverage" from the spread of the model
#   MEAN across draws, with NO observation error injected, giving ~0.05. That
#   column was relabelled "Param cov" as a stopgap; a true posterior-predictive
#   interval needs this sigma. Both are now reported side by side.
# =============================================================================
hiket_total_sigma <- function(sigma_obs_fixed = NA_real_) {
  if (.hiket_correlated_lik) return(HIKET_SIGMA_TOT)
  s <- suppressWarnings(as.numeric(Sys.getenv("HIKET_SIGMA_TOTAL", NA)))
  if (is.finite(s)) return(s)
  sigma_obs_fixed
}
