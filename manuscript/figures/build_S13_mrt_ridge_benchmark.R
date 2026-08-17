setwd("/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling")
# =============================================================================
# S13 -- F14's RIDGE, FOR THE BENCHMARK TRIO (SP1 / TP2 / TP3).
#
# The appendix companion to F14. Identical construction, identical reference
# condition, identical estimator and colour scale, so the two figures can be laid
# side by side and read as one comparison. Everything in F14's header applies here
# and is not repeated; only what DIFFERS for the simple models is documented below.
#
# WHY IT IS WORTH A PANEL. F14 establishes that MRT and sigma_input trade off along
# a ridge in the Yasso family (r = -0.79/-0.73/-0.57) and that the simple models do
# not do this (r ~ -0.26/-0.31, from doublechecks/ridge_test.R section 1). That
# contrast is currently asserted from a CORRELATION ALONE. This figure shows the
# geometry the correlation summarises -- a round blob against Yasso's diagonal
# smear -- and prices the same "what would a longer MRT cost" question for the
# trio, on the same axes and the same colour scale.
#
# --- FOUR THINGS DIFFER FROM F14 -------------------------------------------
#
# 1. MRT IS COMPUTED THROUGH THE MODEL'S OWN *_steady_state ROUTINE, NOT THE
#    CLOSED FORM IN ridge_test.R -- and that is not cosmetic. ridge_test computes
#      SP1: 1/alpha | TP2: 1/aA + pH/aH | TP3: 1/aA + pS/aS + pS*pH/aH
#    which is MRT at xi = 1, i.e. WITH THE CLIMATE RESPONSE SWITCHED OFF. All three
#    models put xi on every pool rate (TP3 since C2, 2026-07-16), so exactly as in
#    Yasso07, xi is a pure rescaling of time and MRT = MRT_ref / xi EXACTLY. Leaving
#    xi out therefore drops a whole per-draw factor from the x coordinate -- and
#    since the calibration moves xi, it drops a factor that MOVES. Here xi is
#    included, at the same dataset-mean reference climate F14 uses, so the x axis
#    is the same physical quantity in both figures.
#    ⚠ CONSEQUENCE: the correlations printed at the end are NOT the -0.26/-0.31 of
#    ridge_test.R. Both are correct; they are different quantities. If they diverge
#    materially, the climate response is doing work in the split and the storyline
#    note's "the simple models show no such ridge" needs restating on this basis.
#
# 2. sigma_input IS FORCED TO 1 INSIDE THE MRT EVALUATION. The simple models'
#    steady states are written as J = lm$J_total_mean * sigma_input, so an unmodified
#    call would put sigma_input on BOTH axes and manufacture the ridge the figure is
#    testing for. Unit input (J_total_mean = 1) with sigma_input = 1 gives the pure
#    generator property, matching Yasso's intrinsic MRT.
#    ⚠ And NOT the engine's steady_state_*_engine binding: that wraps *_transient_init
#    (a 1917 equilibrium plus a 68-yr ramp), which is contaminated by sigma_init. Same
#    trap documented for Yasso in CLAUDE.md; the pure routines are used here.
#
# 3. THE DASHED REFERENCE IS THE ICBM ANCHOR, NOT A PUBLISHED PARAMETERISATION.
#    The trio has no published calibration to compare against -- by design, since C1
#    anchors their kinetics externally on ICBM (Andren & Katterer 1997: k1 = 0.8,
#    k2 = 0.00605, h = 0.13, bulk MRT 1/k1 + h/k2 = 22.74 yr at the Ultuna xi = 1).
#    So the line is drawn at the anchor evaluated EXACTLY as F14 draws the published
#    point: MRT of e$.to_original(best_x), i.e. of the prior centre, pushed through
#    the same reference condition. It is the same kind of object as F14's dashed line
#    -- the external, non-SOC-fitted reference the calibration is allowed to leave --
#    and it is what the "transferability diagnostic" of C1 is measured against.
#    NB it is NOT 22.74: the reference climate is Finland, not Ultuna, so the anchor
#    lands at 22.74 * xi_Ultuna / xi_Finland. The printout reports both.
#
# 4. THE VARIANCE PANEL GAINS A "rates" GROUP AND LOSES "woody size". In Yasso the
#    a-vector is FIXED, so MRT variance can only come from fractions, climate and the
#    woody-size submodel. The trio has no size submodel, and its slow rates are free
#    under a very-informative prior -- so "rates" is a live contributor here and is
#    the group that carries C1's transferability question. SP1 has no fractions at
#    all (single pool); its bar chart correctly shows only two groups.
#
# --- ✅ RUN IDs NOW COME FROM ONE JOB (fixed 2026-08-17) ----------------------
# HISTORICAL: SP1 used to sit on 20260810_152914 (sigma 0.72, Gaussian) while TP2/TP3
# used 20260812_0809* (sigma 0.80, Student-t nu = 6), because SP1 OOMed at a chain
# boundary in job 597028. Log-likelihoods are NOT comparable across different sigma
# -- the normalising constant changes -- so SP1's deficits were on a different
# likelihood from TP2/TP3's and could not be compared in absolute units across panels.
# That caveat is now RETIRED: all three take the corrected-target run via run_ids.R,
# one error model throughout, and F14 does the same, so S13 and F14 are again strictly
# comparable. Each panel is still stamped with its RUN_ID.
#
# Run from repo root:  Rscript manuscript/figures/build_S13_mrt_ridge_benchmark.R
# =============================================================================

suppressMessages(library(BayesianTools))
suppressMessages(library(fields))   # interp.surface.grid
suppressMessages(library(pals))     # kovesi diverging rainbow
# Sourced for setup() and -- the point -- for `ref`, the SAME dataset-mean climate
# F14's MRT is evaluated at. Rebuilding it locally from a simple model's own
# environment would give a near-identical number and no guarantee of it.
source("doublechecks/intrinsic_mrt_lib.R")
source("manuscript/figures/model_palette.R")

source("manuscript/figures/run_ids.R")   # auto-selects current RUN_IDs
RID <- RID[c("SP1", "TP2", "TP3")]       # this figure is the benchmark trio only
ICBM <- list(k1 = 0.8, k2 = 0.00605, h = 0.13, xi_Ultuna = 0.9397)
NB    <- 60L      # ESTIMATION grid cells per axis      (as F14)
NFINE <- 420L     # DISPLAY grid, bilinear              (as F14)
RMAX  <- 3L       # max adaptive window radius          (as F14)
DLL   <- 25       # colour span, ll units below the best (as F14)
# ⚠ MINN = 20, NOT F14's 5 -- the ONE deliberate departure, and it is an estimator
# choice, not a change to the quantity plotted. The trio's likelihood surface is far
# STEEPER in these coordinates than Yasso's: it falls 42-44 ll units across the
# visited region (SP1 18) where the Yasso ridge is nearly flat over its whole length.
# A max over ~5 draws carries several units of noise, invisible against a flat surface
# but rendered as red-on-green speckle against a steep one. Requiring 20 draws behind
# every plotted cell -- at the SAME grid resolution and the SAME 25-unit colour scale --
# removes it. The change is strictly CONSERVATIVE: it masks more cells, never fewer,
# and a max over more draws is a tighter lower bound on the profile, so the surface
# shown is if anything closer to the truth.
# TODO when both figures are rebuilt on the corrected-target run: settle on one MINN
# for F14 and S13, and update F14's caption ("fewer than five times") to match.
MINN  <- 20L
set.seed(2025)

# --- intrinsic MRT for a simple model ----------------------------------------
# Unit litter, dataset-mean climate, sigma_input neutralised, PURE steady state.
mrt_fun_simple <- function(M, e) {
  ss  <- get(sprintf("%s_steady_state", tolower(M)), envir = e)
  cxm <- get(sprintf("compute_xi_mean_%s_engine", tolower(M)), envir = e)
  LM1 <- list(J_total_mean = 1)                       # unit input
  function(p) {
    mp <- e$.assemble(p)
    mp["sigma_input"] <- 1                            # see header note 2
    xi <- cxm(clim_ss = ref$clim, model_params = mp)
    sum(ss(mp, LM1, xi))
  }
}
# xi at the reference climate, for the anchor bookkeeping in the printout
xi_fun_simple <- function(M, e) {
  cxm <- get(sprintf("compute_xi_mean_%s_engine", tolower(M)), envir = e)
  function(p) unname(cxm(clim_ss = ref$clim, model_params = e$.assemble(p)))
}

MRT_ICBM <- 1/ICBM$k1 + ICBM$h/ICBM$k2               # 22.74 yr, at xi = 1 (Ultuna)
cat(sprintf("\nICBM anchor: 1/k1 + h/k2 = %.2f yr at xi = 1 (Ultuna r = 1)\n\n", MRT_ICBM))

D <- list(); PUB <- setNames(numeric(0), character(0))
for (M in names(RID)) {
  e <- setup(M); f <- mrt_fun_simple(M, e); xf <- xi_fun_simple(M, e)
  bx   <- get("best_x", e)
  free <- names(bx)                                   # order matters for to_original
  p_anchor <- e$.to_original(bx)                      # prior centre == the ICBM anchor
  PUB[M]   <- f(p_anchor)

  # ⚠ start = 2 -- drop the first retained iteration of every internal DEzs chain.
  # Same start-up-state artefact documented at length in build_F14_mrt_ridge.R: 15
  # rows per model sitting 100-200 ll below the bulk. Extracting per sampler also
  # avoids the 1-in-3 thinning getSample applies to a sampler LIST.
  .o <- readRDS(sprintf("Calibration_real_data_transient/runs/%s_chains_%s.rds", M, RID[[M]]))
  s  <- do.call(rbind, lapply(.o, function(z) getSample(z, parametersOnly = FALSE, start = 2)))
  # CHAINS are in SAMPLING space; to_original resolves the transforms and the
  # coupled sigma_input/sigma_init flux_pair construction.
  ph  <- lapply(seq_len(nrow(s)), function(k)
                 tryCatch(e$.to_original(s[k, free]), error = function(z) NULL))
  mrt <- vapply(ph, function(q) if (is.null(q)) NA_real_ else
                  tryCatch(f(q), error = function(z) NA_real_), numeric(1))
  si  <- vapply(ph, function(q) if (is.null(q)) NA_real_ else q[["sigma_input"]], numeric(1))
  xiv <- vapply(ph, function(q) if (is.null(q)) NA_real_ else
                  tryCatch(xf(q), error = function(z) NA_real_), numeric(1))
  ll  <- s[, "Llikelihood"]
  ok  <- is.finite(mrt) & is.finite(si) & is.finite(ll)
  PH  <- do.call(rbind, ph[ok])                       # physical draws, for the variance panel
  D[[M]] <- list(d = data.frame(mrt = mrt[ok], si = si[ok], ll = ll[ok], xi = xiv[ok]),
                 PH = PH, f = f, anchor = p_anchor, xi_anchor = xf(p_anchor))
  cat(sprintf("%-4s %6d draws | MRT %.1f-%.1f (median %.1f) | anchor %.1f | ll range %.0f\n",
              M, sum(ok), min(D[[M]]$d$mrt), max(D[[M]]$d$mrt), median(D[[M]]$d$mrt),
              PUB[[M]], diff(range(D[[M]]$d$ll))))
  cat(sprintf("       xi at reference climate: anchor %.3f -> posterior median %.3f (x%.2f)",
              D[[M]]$xi_anchor, median(D[[M]]$d$xi), median(D[[M]]$d$xi)/D[[M]]$xi_anchor))
  cat(sprintf(" | anchor MRT check %.2f x %.4f/%.4f = %.2f\n",
              MRT_ICBM, ICBM$xi_Ultuna, D[[M]]$xi_anchor,
              MRT_ICBM * ICBM$xi_Ultuna / D[[M]]$xi_anchor))
}

# --- first-order variance indices of MRT, per model --------------------------
# "rates" replaces Yasso's "woody size": the trio has no size submodel, and its slow
# rates ARE free (very-informative prior), so they can contribute. Fractions exist
# only for TP2 (p_H) and TP3 (p_S, p_H); SP1 returns NA there and the bar is dropped.
grp_of <- function(n)
  ifelse(grepl("^p_", n), "fractions",
  ifelse(grepl("^(beta|gamma)", n), "climate",
  ifelse(grepl("^alpha", n), "rates", NA)))
GC <- c(climate = "#1b9e77", fractions = "#d95f02", rates = "#7570b3")

VAR <- list()
for (M in names(RID)) {
  PH <- D[[M]]$PH; f <- D[[M]]$f
  med <- apply(PH, 2, median)
  kk  <- if (nrow(PH) > 4000L) round(seq(1, nrow(PH), length.out = 4000L)) else seq_len(nrow(PH))
  V   <- var(D[[M]]$d$mrt[kk])
  VAR[[M]] <- vapply(names(GC), function(g) {
    nm <- colnames(PH)[which(grp_of(colnames(PH)) == g)]
    if (!length(nm)) return(NA_real_)
    v <- vapply(kk, function(k) { q <- med; q[nm] <- PH[k, nm]
      tryCatch(f(q), error = function(z) NA_real_) }, numeric(1))
    var(v, na.rm = TRUE) / V }, numeric(1))
}

# =============================================================================
png("manuscript/figures/S13_mrt_ridge_benchmark.png",
    width = 13.2, height = 9.0, units = "in", res = 200)
layout(matrix(1:6, 2, 3, byrow = TRUE), heights = c(1.35, 1))
par(mar = c(4.4, 4.8, 5.0, 1.2), mgp = c(2.8, 0.7, 0), las = 1)
# Same reversed Kovesi bgymr rainbow as F14 -- the two figures MUST share a colour
# map or the side-by-side comparison they exist for is not readable.
pal <- rev(pals::kovesi.diverging_rainbow_bgymr_45_85_c67(256))

for (M in names(RID)) {                       # --- row 1: the ridge ---
  d <- D[[M]]$d
  xr <- range(c(d$mrt, PUB[[M]])) + c(-0.6, 0.6); yr <- range(d$si)
  nb <- NB
  bx <- cut(d$mrt, seq(xr[1], xr[2], length.out = nb+1), labels = FALSE)
  by <- cut(d$si,  seq(yr[1], yr[2], length.out = nb+1), labels = FALSE)
  # MAX per cell (a sampled lower bound on the profile likelihood), with the
  # adaptive window that keeps full resolution where draws are dense and only
  # smooths the sparse tails. Rationale in full in build_F14_mrt_ridge.R.
  Zc <- matrix(NA_real_, nb, nb); Nc <- matrix(0L, nb, nb)
  Z <- tapply(d$ll, list(bx, by), max); N <- tapply(d$ll, list(bx, by), length)
  Zc[cbind(as.integer(rownames(Z))[row(Z)], as.integer(colnames(Z))[col(Z)])] <- Z
  Nc[cbind(as.integer(rownames(N))[row(N)], as.integer(colnames(N))[col(N)])] <-
    ifelse(is.na(N), 0L, N)
  shiftm <- function(M, di, dj, fill) {          # M shifted by (di, dj), edges filled
    out <- matrix(fill, nrow(M), ncol(M))
    si <- max(1L, 1L - di):min(nrow(M), nrow(M) - di)
    sj <- max(1L, 1L - dj):min(ncol(M), ncol(M) - dj)
    out[si + di, sj + dj] <- M[si, sj]; out
  }
  Zi <- Zc; Zi[is.na(Zi)] <- -Inf
  Zcum <- Zi; Ncum <- Nc                          # radius 0
  Zf <- matrix(NA_real_, nb, nb); RUSED <- matrix(NA_integer_, nb, nb)
  take <- function(r) { hit <- is.na(Zf) & Ncum >= MINN & is.finite(Zcum)
                        Zf[hit] <<- Zcum[hit]; RUSED[hit] <<- r }
  take(0L)
  for (r in seq_len(RMAX)) {                      # add only the new ring at each radius
    for (di in -r:r) for (dj in -r:r) if (max(abs(di), abs(dj)) == r) {
      Zcum <- pmax(Zcum, shiftm(Zi, di, dj, -Inf)); Ncum <- Ncum + shiftm(Nc, di, dj, 0L)
    }
    take(r)
  }
  cat(sprintf("  %-4s window radius used: median %d, max %d cells (cell = %.2f%% of axis)\n",
              M, median(RUSED, na.rm = TRUE), suppressWarnings(max(RUSED, na.rm = TRUE)),
              100/nb))
  best <- max(Zf, na.rm = TRUE)
  Zp   <- pmax(Zf, best - DLL)                    # deficit, clamped

  gx <- seq(xr[1], xr[2], length.out = nb); gy <- seq(yr[1], yr[2], length.out = nb)
  fx <- seq(xr[1], xr[2], length.out = NFINE); fy <- seq(yr[1], yr[2], length.out = NFINE)
  Zi2 <- fields::interp.surface.grid(list(x = gx, y = gy, z = Zp),
                                     grid.list = list(x = fx, y = fy))$z

  image(fx, fy, Zi2, col = pal, zlim = c(best - DLL, best),
        xlab = "intrinsic MRT (yr)", ylab = expression(sigma[input]),
        main = "", cex.main = 1.0)
  box()
  contour(MASS::kde2d(d$mrt, d$si, n = 60), add = TRUE, drawlabels = FALSE,
          col = adjustcolor("white", 0.75), lwd = 1.3, nlevels = 5)
  points(median(d$mrt), median(d$si), pch = 21, bg = "white", col = "black", cex = 1.5, lwd = 2)
  abline(v = PUB[[M]], lty = 2, lwd = 2, col = "grey15")
  title(main = sprintf("%s  --  best attainable fit over (MRT, sigma_input)", M),
        line = 2.65, cex.main = 1.05)
  mtext(sprintf("posterior reaches %.1f-%.1f yr;  ICBM anchor %.1f (dashed);  run %s",
                min(d$mrt), max(d$mrt), PUB[[M]], RID[[M]]),
        side = 3, line = 1.45, cex = 0.60, col = "grey35")
  mtext(sprintf("colour = log-likelihood deficit: 0 (blue = best fit) to %d+ (red);  ll spans %.0f here",
                DLL, diff(range(Zf, na.rm = TRUE))),
        side = 3, line = 0.45, cex = 0.58, col = "grey45")
}

par(mar = c(4.4, 4.8, 3.4, 1.2))
for (M in names(RID)) {                       # --- row 2: variance of MRT ---
  v <- 100 * VAR[[M]]; v <- v[is.finite(v)]
  bp <- barplot(v, col = GC[names(v)], border = NA, ylim = c(0, max(115, max(v)*1.15)),
                ylab = "share of posterior MRT variance (%)",
                main = sprintf("%s  --  what makes MRT vary", M), cex.main = 1.0, cex.names = 0.95)
  abline(h = 100, lty = 3, col = "grey50")
  text(bp, v, sprintf("%.0f%%", v), pos = 3, cex = 0.85, font = 2)
  mtext(sprintf("first-order indices, total %.0f%%  (>100%% = the groups interact)", sum(v)),
        side = 3, line = 0.3, cex = 0.62, col = "grey35")
}
dev.off()
cat("\nWrote manuscript/figures/S13_mrt_ridge_benchmark.png\n")

cat("\ncorrelation of log MRT with log sigma_input (a ridge would be strongly negative):\n")
cat("⚠ NOT comparable to ridge_test.R section 1 -- that one omits xi (see header note 1).\n")
for (M in names(RID)) { d <- D[[M]]$d
  cat(sprintf("  %-4s r = %+0.3f | sd(log product)/sd(log MRT) = %.2f", M,
              cor(log(d$mrt), log(d$si)), sd(log(d$mrt)+log(d$si))/sd(log(d$mrt))))
  cat(sprintf("   [xi omitted: r = %+0.3f]\n",
              cor(log(d$mrt) + log(d$xi), log(d$si)))) }
saveRDS(D, "manuscript/figures/S13_mrt_ridge_benchmark.rds")
