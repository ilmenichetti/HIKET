setwd("/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling")
# =============================================================================
# F14 -- THE RIDGE: fitness over (MRT, sigma_input), read straight off the MCMC.
#
# Built from the posterior draws themselves rather than from a grid slice. Each
# draw already carries a full physical parameter set AND its log-likelihood
# (chains keep Lposterior / Llikelihood / Lprior), so:
#     x = intrinsic MRT recomputed per draw   (collapses kinetics + fractions +
#                                              climate response into one scalar)
#     y = sigma_input of that draw
#     colour = MAX Llikelihood in the cell (a sampled lower bound on the PROFILE
#              likelihood: the best fit achievable there). Cells with < 5 draws are
#              masked -- unvisited is not the same as unfittable.
#
# WHY THIS BEATS A GRID SLICE. A slice must hold the other parameters somewhere,
# and the Yasso transforms make that awkward -- fractions are stick-breaking (not
# invertible coordinate-wise) and sigma_input/sigma_init are COUPLED in the
# transform (moving one raw coordinate moves both physical values, via the
# flux_pair construction). The draws sidestep all of it: everything varies jointly
# and in physical space, so this is a genuine projection of the posterior, not a
# conditional slice.
#
# HOW TO READ IT. A diagonal elongation means MRT and sigma_input trade off -- the
# data constrain their product, not the split, and tightening sigma_input would
# drag MRT with it. A round blob means they are separately determined and
# constraining one will NOT move the other. Colour tells whether fit actually
# degrades along the ridge or is flat over it.
#
# THE HONEST LIMIT. The posterior only visits where it has mass, so this cannot
# say what MRT = 30 would COST if no draw ever reached 30. The printed reach of
# each posterior states how far it actually got; anything beyond that needs a
# profile likelihood, not this figure.
#
# Run from repo root:  Rscript manuscript/figures/build_F14_mrt_ridge.R
# =============================================================================

suppressMessages(library(BayesianTools))
suppressMessages(library(fields))   # interp.surface.grid
suppressMessages(library(pals))     # kovesi diverging rainbow
source("doublechecks/intrinsic_mrt_lib.R")     # setup(), mrt_fun(), ref
source("manuscript/figures/model_palette.R")

RID <- c(Yasso07 = "20260812_080941", Yasso15 = "20260812_080940", Yasso20 = "20260812_080940")
PUB <- c(Yasso07 = 33.47, Yasso15 = 30.38, Yasso20 = 19.03)   # published POINT MRT
NS  <- Inf                                                     # use EVERY draw (see note)
NB  <- 60L                                                     # ESTIMATION grid cells per axis
NFINE <- 420L                                                  # DISPLAY grid (bilinear interpolation)
RMAX <- 3L                                                     # max adaptive window radius
DLL  <- 25                                                     # colour span, in ll units below the best
set.seed(2025)

# WHY ALL THE DRAWS. One MRT evaluation costs ~1e-4 s, so the whole chain (~75k
# draws) is ~8 s -- there is no reason to subsample. And it is not cosmetic: the
# max-per-cell estimator is biased LOW, and the bias shrinks with the number of
# draws IN THAT CELL. At 6000 draws the dense cells near the mode were nearly
# converged while the sparse tail cells were badly underestimated, so the ridge
# appeared to decay faster than it does -- a bias pointing the convenient way.
# Using every draw shrinks it and, more importantly, evens it out across the plane.
# It does NOT extend the posterior's REACH: the white region stays white, because
# that is where the sampler never went, and only a profile likelihood can price it.

D <- list()
for (M in names(RID)) {
  e <- setup(M); f <- mrt_fun(M, e)
  p_pub <- e$.to_original(get("best_x", e))
  # ⚠ DROP THE FIRST RETAINED ITERATION OF EVERY INTERNAL CHAIN (start = 2).
  # The run stores 5 samplers, each a DEzs with 3 internal chains. Their first retained
  # row is a START-UP STATE, not a posterior sample: 15 rows per model (3 x 5) sitting
  # 100-200 log-likelihood units below everything else, with a clean gap to the bulk
  # (Yasso07 min ll -1346 -> -1170 once dropped). Left in, they appear as isolated
  # "puddles" far from the posterior cloud, they set the reported MRT reach, and -- worst
  # -- they were the ONLY draws near the published MRT for Yasso07/15, so any cost read
  # off them is an artefact of initialisation rather than a statement about the posterior.
  # Extracting per sampler also avoids the 1-in-3 thinning getSample applies to the list.
  .o <- readRDS(sprintf("Calibration_real_data_transient/runs/%s_chains_%s.rds", M, RID[[M]]))
  s  <- do.call(rbind, lapply(.o, function(z) getSample(z, parametersOnly = FALSE, start = 2)))
  # ⚠ The CHAINS are in SAMPLING space while the posterior .rds is physical
  # (exp(-1.8477) = 0.1576 for Yasso07's beta1). Every draw must go through
  # to_original -- which also resolves the stick-breaking fractions and the
  # coupled sigma_input/sigma_init correctly, so the physical values are exact.
  bx   <- get("best_x", e)
  free <- names(bx)                              # order matters for to_original
  i    <- if (is.infinite(NS)) seq_len(nrow(s)) else sample(nrow(s), min(NS, nrow(s)))
  ph   <- lapply(i, function(k) tryCatch(e$.to_original(s[k, free]), error = function(z) NULL))
  mrt  <- vapply(ph, function(q) if (is.null(q)) NA_real_ else
                   tryCatch(f(q), error = function(z) NA_real_), numeric(1))
  si   <- vapply(ph, function(q) if (is.null(q)) NA_real_ else q[["sigma_input"]], numeric(1))
  ll   <- s[i, "Llikelihood"]
  ok <- is.finite(mrt) & is.finite(si) & is.finite(ll)
  PH <- do.call(rbind, ph[ok])                       # physical draws, for panel (b)
  D[[M]] <- list(d = data.frame(mrt = mrt[ok], si = si[ok], ll = ll[ok]),
                 PH = PH, f = f, pub = e$.to_original(bx))
  cat(sprintf("%-9s %5d draws | MRT %.1f-%.1f (median %.1f) | published %.1f | ll range %.0f\n",
              M, sum(ok), min(D[[M]]$d$mrt), max(D[[M]]$d$mrt), median(D[[M]]$d$mrt), PUB[[M]],
              diff(range(D[[M]]$d$ll))))
}

# --- first-order variance indices of MRT, per model --------------------------
grp_of <- function(n)
  ifelse(grepl("^p_", n), "fractions",
  ifelse(grepl("^(beta|gamma)", n), "climate",
  ifelse(n %in% c("delta1","delta2","r","w1","w2","w3","w4","w5"), "woody size", NA)))
GC <- c(climate = "#1b9e77", fractions = "#d95f02", `woody size` = "#66c2a5")

VAR <- list()
for (M in names(RID)) {
  PH <- D[[M]]$PH; f <- D[[M]]$f
  med <- apply(PH, 2, median)
  # variance indices converge in a few thousand draws -- no need for the full chain
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
png("manuscript/figures/F14_mrt_ridge.png", width = 13.2, height = 9.0, units = "in", res = 200)
layout(matrix(1:6, 2, 3, byrow = TRUE), heights = c(1.35, 1))
par(mar = c(4.4, 4.8, 5.0, 1.2), mgp = c(2.8, 0.7, 0), las = 1)
# Kovesi bgymr diverging rainbow, REVERSED: BLUE = best fit, through green/yellow, to
# RED = worst. Chosen over viridis for CONTRAST -- viridis is monotone in lightness, and
# most of this surface sits within a few ll units of the best, so a lightness-only ramp
# renders the plane nearly flat. A rainbow varies in hue as well, which resolves the
# small deficits along the ridge that are the whole point of the panel.
# ⚠ The usual objection to rainbows -- false banding on a continuous field -- is exactly
# what is wanted here: the bands ARE iso-likelihood contours, and Kovesi's construction
# keeps the steps perceptually even, so no band is an artefact of the colour map.
pal <- rev(pals::kovesi.diverging_rainbow_bgymr_45_85_c67(256))

for (M in names(RID)) {                       # --- row 1: the ridge ---
  d <- D[[M]]$d
  xr <- range(c(d$mrt, PUB[[M]])) + c(-0.6, 0.6); yr <- range(d$si)
  nb <- NB
  bx <- cut(d$mrt, seq(xr[1], xr[2], length.out = nb+1), labels = FALSE)
  by <- cut(d$si,  seq(yr[1], yr[2], length.out = nb+1), labels = FALSE)
  # MAX per cell, not mean: this approximates the PROFILE likelihood in these two
  # coordinates -- "the best fit achievable here if everything else rearranges" --
  # whereas the mean answers "what the posterior typically achieves here", which is
  # partly self-referential (the draws were placed by L and p) and is dragged down
  # by the bad corners of the other 18 dimensions.
  # ⚠ It is a LOWER BOUND on the profile, and it degrades where draws are sparse --
  # i.e. in the far tail at high MRT, exactly the region we most want to judge. So
  # cells below MINN draws are masked rather than plotted: an unvisited cell means
  # "the sampler did not go there", NOT "the model cannot fit there".
  # RESOLUTION, AND THE TRADE-OFF IT SITS ON. A hard bin is limited by draws-per-cell:
  # at nb = 60 that is ~21 on average, and simply raising nb punches holes and makes
  # each cell's max a worse lower bound. But widening the window uniformly is no free
  # lunch either -- it buys draws by SPENDING resolution, and a fine display grid then
  # only makes a smoother picture of a blurrier estimate.
  #
  # So the window is ADAPTIVE: at each node it grows from radius 0 until it holds MINN
  # draws (capped at RMAX). Dense regions -- the ridge itself -- keep the full 1/NB
  # resolution; only the sparse tails get smoothed, and exactly as much as they need.
  # The per-model print reports the radius actually used, which IS the honest resolution
  # of each panel. Rationale for taking a max over a neighbourhood at all: the profile
  # likelihood is smooth (an envelope of optima), so nearby draws are informative about
  # this node, and a bin that discards them is simply a worse estimator.
  # Bleed is bounded by the radius used, so a blank cell still means UNVISITED.
  MINN <- 5L
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
  cat(sprintf("  %-9s window radius used: median %d, max %d cells (cell = %.2f%% of axis)\n",
              M, median(RUSED, na.rm = TRUE), suppressWarnings(max(RUSED, na.rm = TRUE)),
              100/nb))
  # COLOUR SCALE. Plot the DEFICIT from the best fit found, clamped at DLL units.
  # A raw linear scale over the full ll range does not work here: the range is 100-200
  # units, owned by a handful of terrible far-tail cells, which flattens the 5-10 units
  # of structure along the ridge -- the only part anyone reads -- into one flat colour.
  # Clamping is also the honest cut: a deficit past ~25 ll units is decisively rejected,
  # so resolving 30 from 200 conveys nothing. Everything worse saturates at the dark end.
  best <- max(Zf, na.rm = TRUE)
  Zp   <- pmax(Zf, best - DLL)

  # INTERPOLATION. The estimate lives on the coarse grid (NB cells, where each cell has
  # enough draws to support a max); the DISPLAY is bilinearly interpolated onto NFINE.
  # This is a display choice, not new information -- it removes the staircase without
  # claiming resolution the draws cannot support, and it is the same continuous surface
  # the eventual profile-likelihood version will produce directly. NA corners propagate,
  # so the support boundary is preserved (eroded by at most one coarse cell), and a blank
  # region still means UNVISITED.
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
  mtext(sprintf("posterior reaches %.1f-%.1f yr;  published %.1f (dashed)",
                min(d$mrt), max(d$mrt), PUB[[M]]), side = 3, line = 1.45, cex = 0.60, col = "grey35")
  mtext(sprintf("colour = log-likelihood deficit: 0 (blue = best fit) to %d+ (red)", DLL),
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
cat("\nWrote manuscript/figures/F14_mrt_ridge.png\n")

cat("\ncorrelation of log MRT with log sigma_input (a ridge would be strongly negative):\n")
for (M in names(RID)) { d <- D[[M]]$d
  cat(sprintf("  %-9s r = %+0.3f | sd(log product)/sd(log MRT) = %.2f\n", M,
              cor(log(d$mrt), log(d$si)), sd(log(d$mrt)+log(d$si))/sd(log(d$mrt)))) }
saveRDS(D, "manuscript/figures/F14_mrt_ridge.rds")
