# =============================================================================
# mrt_attribution.R — WHICH parameters make our MRT shorter than published?
#
# THE PUZZLE. Intrinsic MRT (unit input, fixed reference climate, pure
# steady-state routine) comes out at:
#     Yasso07  15.2 vs published point 33.5   (46% of it)
#     Yasso15  22.1 vs 30.4                   (73%)
#     Yasso20  17.7 vs 19.0                   (93%)
# The gap is NOT uniform, so it is not one shared mechanism such as "sigma_input
# is the loosest parameter and absorbs everything". And the transfer fractions are
# prior-pinned (KL < 1 nat, most within 0.5 sigma of their centres), so they
# should not be able to move MRT far -- yet Yasso07 is at less than half.
#
# THE TEST. MRT here is a function of the free parameters only. Split them into
# groups and swap one group at a time between OUR posterior median and the
# PUBLISHED defaults, holding everything else at the other source:
#
#     published   : all published            (= the "published POINT" number)
#     +fractions  : published, but OUR p_*   -> isolates the transfer fractions
#     +climate    : published, but OUR betas -> isolates the climate response
#     +other      : published, but OUR rest  -> woody size, p_H, alpha_H, ...
#     ours        : all ours                 (= the "OUR posterior" number)
#
# Whichever single swap collapses the gap identifies the culprit. If none does but
# the total does, the cause is an interaction and no single parameter group owns it.
#
# WHY IT MATTERS. The planned next lever is tightening the sigma_input prior,
# predicted to bring Yasso15 to ~34 yr. That presumes the fast MRT is bought by
# sigma_input's looseness. If instead the climate response is driving it, that
# experiment is aimed at the wrong parameter and would only degrade the fit.
#
# Run from repo root:  Rscript doublechecks/mrt_attribution.R
# =============================================================================

suppressMessages(library(BayesianTools))
options(width = 120)
source("doublechecks/intrinsic_mrt_lib.R")   # setup(), mrt_fun(), ref  (shared)

RID <- c(Yasso07 = "20260812_080941", Yasso15 = "20260812_080940", Yasso20 = "20260812_080940")

# --- parameter groups --------------------------------------------------------
is_fraction <- function(n) grepl("^p_", n)                     # transfer fractions incl. p_H
is_climate  <- function(n) grepl("^(beta|gamma)", n)           # beta1/2, gamma, N and H variants
is_size     <- function(n) n %in% c("delta1","delta2","r","w1","w2","w3","w4","w5")

cat("\n", strrep("=", 104), "\nMRT ATTRIBUTION: which parameter group carries the gap?\n",
    strrep("=", 104), "\n", sep = "")

res <- list()
for (M in names(RID)) {
  e <- setup(M); f <- mrt_fun(M, e)
  p_def <- e$.to_original(get("best_x", e))                    # published defaults

  smp <- getSample(readRDS(sprintf("Calibration_real_data_transient/runs/%s_posterior_%s.rds",
                                   M, RID[[M]])))
  shared <- intersect(names(p_def), colnames(smp))             # the FREE parameters
  p_our  <- p_def
  p_our[shared] <- apply(smp[, shared, drop = FALSE], 2, median)

  swap <- function(sel) { p <- p_def; nm <- shared[sel(shared)]
                          if (length(nm)) p[nm] <- p_our[nm]; f(p) }
  v <- c(published  = f(p_def),
         `+fractions` = swap(is_fraction),
         `+climate`   = swap(is_climate),
         `+size`      = swap(is_size),
         `+rest`      = swap(function(n) !is_fraction(n) & !is_climate(n) & !is_size(n)),
         ours        = f(p_our))

  gap <- v["ours"] - v["published"]
  cat(sprintf("\n=== %s ===   published %.2f -> ours %.2f  (gap %+.2f yr)\n", M,
              v["published"], v["ours"], gap))
  cat(sprintf("  %-12s %8s %10s   %s\n", "swap", "MRT", "shift", "share of gap"))
  for (nm in names(v)[2:5]) {
    sh <- v[[nm]] - v["published"]
    cat(sprintf("  %-12s %8.2f %+10.2f   %s\n", nm, v[[nm]], sh,
                if (abs(gap) > 1e-9) sprintf("%+5.0f%%", 100*sh/gap) else "--"))
  }
  cat(sprintf("  %-12s %8.2f %+10.2f   %s\n", "free params", NA, NA,
              paste(shared, collapse = " ")))
  res[[M]] <- v
}

cat("\n", strrep("=", 104), "\nREAD-OUT\n", strrep("=", 104), "\n", sep = "")
cat("A single group near 100% owns the gap. Groups summing to far from 100% mean the\n")
cat("effect is INTERACTIVE -- swapping one at a time understates it, and no single\n")
cat("prior can be blamed. Either way, if 'climate' rather than the fractions carries\n")
cat("the gap, the planned sigma_input tightening is aimed at the wrong parameter.\n")
saveRDS(res, "doublechecks/mrt_attribution.rds")
