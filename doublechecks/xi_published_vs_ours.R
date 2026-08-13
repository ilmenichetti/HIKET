# =============================================================================
# xi_published_vs_ours.R — WHY THE CLIMATE RESPONSE OWNS YASSO07's MRT GAP
#                          AND CANNOT OWN THE OTHER TWO. (2026-08-13)
#
# THE PUZZLE. The MRT variance decomposition (F14, lower row) and the swap test
# (mrt_attribution.R) agree that Yasso07's whole MRT gap is the CLIMATE response
# (~97-109%), while Yasso15's is mostly transfer FRACTIONS (88%) and Yasso20's is
# mixed. That asymmetry looked like a prior artefact. It is not: it is STRUCTURAL.
#
# THE MECHANISM. Yasso07 applies ONE climate modifier xi to every pool, so xi is a
# pure rescaling of time and
#           MRT = MRT_reference / xi     (exactly)
# Any climate misfit therefore converts one-for-one into MRT, and the model has
# nowhere else to put it. Yasso15/Yasso20 carry THREE pool-specific modifiers
# (xi_awe, xi_n, xi_h, each with its own beta/gamma). MRT is dominated by the slow
# HUMUS pool, which has its own betaH1 -- so climate can no longer rescale the
# system uniformly, the components can even move in OPPOSITE directions, and the
# leverage that Yasso07 hands to climate is structurally unavailable.
#
# WHAT IT PRINTS (run 20260812_0809*, sigma 0.80 + Student-t nu=6):
#   Yasso07   xi 0.855 -> 1.822  (x2.13)   MRT 33.47 -> 15.20   predicted-from-xi 15.7
#   Yasso15   awe x0.98  n x1.21  h x1.04  MRT 30.38 -> 22.07
#   Yasso20   awe x0.76  n x1.16  h x1.06  MRT 19.03 -> 17.52
#
# TWO READINGS THAT MATTER FOR THE PAPER.
#   1. Yasso07's PUBLISHED xi is 0.855 < 1 at the Finnish mean climate -- its
#      published parameterisation says Finland decomposes SLOWER than its reference
#      condition, while Yasso15/20 already say faster (1.16-1.66). Yasso07 starts
#      furthest from what the Finnish data want AND is the only one able to travel
#      the whole distance in a single parameter.
#   2. Calibration INVERTS the family ordering. Published 33.5 / 30.4 / 19.0
#      (07 > 15 > 20); ours 15.2 / 22.1 / 17.5 (15 > 20 > 07). Yasso07 goes from
#      the slowest published model to the fastest calibrated one.
#
# ⚠ OPEN: is xi = 1.82 defensible? It needs beta1 0.0987 -> 0.1578 (+60%), and S12
# puts that ~4.4 prior sigma off centre against the corrected (Tuomi-95%-as-2sigma)
# width. Either the Finnish data genuinely demand it, or the single-xi structure is
# absorbing a misfit that belongs elsewhere. Pair with the FMI warming-rate check.
#
# Run from repo root:  Rscript doublechecks/xi_published_vs_ours.R
# =============================================================================

suppressMessages(library(BayesianTools))
source("doublechecks/intrinsic_mrt_lib.R")     # setup(), mrt_fun(), ref
options(width = 110)

RID <- c(Yasso07 = "20260812_080941", Yasso15 = "20260812_080940", Yasso20 = "20260812_080940")

cat("\n", strrep("=", 92),
    "\nCLIMATE MODIFIER AT THE REFERENCE CLIMATE: published -> our posterior median\n",
    strrep("=", 92), "\n", sep = "")
cat(sprintf("reference: T = %.2f C | amplitude = %.2f C | precip = %.0f mm\n\n",
            ref$clim$temp_mean, ref$clim$temp_amplitude, ref$clim$precip))

OUT <- list()
for (M in names(RID)) {
  e <- setup(M); f <- mrt_fun(M, e)
  p_pub <- e$.to_original(get("best_x", e))
  smp <- getSample(readRDS(sprintf("Calibration_real_data_transient/runs/%s_posterior_%s.rds",
                                   M, RID[[M]])))
  sh <- intersect(names(p_pub), colnames(smp))
  p_our <- p_pub; p_our[sh] <- apply(smp[, sh, drop = FALSE], 2, median)

  # Yasso07: a single scalar xi. Yasso15/20: a named list of three.
  if (M == "Yasso07") {
    cxm <- get("compute_xi_mean_yasso07", e)
    g <- function(p) { mp <- e$.assemble(p)
      c(all = unname(cxm(ref$clim, mp[["beta1"]], mp[["beta2"]], mp[["gamma"]]))) }
  } else {
    cxm <- get("compute_xi_mean_yasso15", e)
    YP  <- get(sprintf("%s_PARAM_NAMES", toupper(M)), e)
    g <- function(p) { mp <- e$.assemble(p)
      unlist(lapply(cxm(clim_ss = ref$clim, params = mp[YP]), unname)) }
  }

  a <- g(p_pub); b <- g(p_our)
  cat(sprintf("--- %-8s  MRT %6.2f -> %6.2f  (x%.2f)   beta1 %.4f -> %.4f\n",
              M, f(p_pub), f(p_our), f(p_our)/f(p_pub), p_pub[["beta1"]], p_our[["beta1"]]))
  for (k in names(a))
    cat(sprintf("      xi[%-6s] %6.3f -> %6.3f   x%.2f\n", k, a[k], b[k], b[k]/a[k]))
  if (M == "Yasso07")
    cat(sprintf("      => MRT predicted from xi ALONE: %.1f yr (actual %.1f) -- xi owns the whole gap\n",
                f(p_pub) * a[["all"]] / b[["all"]], f(p_our)))
  OUT[[M]] <- list(xi_pub = a, xi_our = b, mrt_pub = f(p_pub), mrt_our = f(p_our),
                   beta1 = c(pub = p_pub[["beta1"]], our = p_our[["beta1"]]))
}

cat("\n", strrep("=", 92), "\nREAD-OUT\n", strrep("=", 92), "\n", sep = "")
cat("A single xi (Yasso07) rescales time, so MRT = MRT_ref / xi and climate has TOTAL leverage.\n")
cat("Three pool-specific xi (Yasso15/20) cannot: the slow humus pool sets MRT and carries its own\n")
cat("betaH1, so climate misfit is split and can partly cancel. 'MRT too short' is therefore NOT one\n")
cat("phenomenon shared by the family and must not be discussed as if it were.\n")
saveRDS(OUT, "doublechecks/xi_published_vs_ours.rds")
