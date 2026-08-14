# =============================================================================
# ridge_test.R  (2026-08-12)
#
# THE OVERCONFIDENCE DIAGNOSTIC. Establishes, from a finished run, that the
# posterior reports a precision the data cannot support -- and quantifies the
# correlated-error structure responsible.
#
# THE ARGUMENT IT SUPPORTS. The likelihood mainly pins the PRODUCT
# MRT x sigma_input x J_raw (= the observed stock), so the split between MRT and
# sigma_input is decided by whatever else constrains them. If the data genuinely
# could not separate them, the posterior would ride that ridge: strong negative
# corr(log MRT, log sigma_input) and a product much tighter than either factor.
# It does not. Yet ACROSS runs the product is constant to ~4% while the split
# slides by ~2 posterior SD. Between-run sensitivity exceeding within-run
# uncertainty is the operational definition of an overconfident posterior.
#
# WHY: 1269 plot-years are treated as independent and are not (sections 2-4).
#
# Sections
#   1  ridge      -- corr(log MRT, log sigma_input); sd(log product)/sd(log MRT)
#   2  plot       -- ICC and design effect (the plot-persistent component)
#   3  spatial    -- latitude-band mean residuals vs the iid expectation
#   4  scale/tail -- in-sample vs holdout spread; per-campaign spread; kurtosis
#
# MRT is closed-form for SP1/TP2/TP3 (section 1 only); the Yasso trio needs
# doublechecks/intrinsic_mrt.R, whose per-draw values are not paired to
# sigma_input draws, so they are reported there instead.
#
# ⚠ 2026-08-13 -- SECTION 1 IS SIMPLE-MODELS-ONLY AND ITS RESULT DOES NOT GENERALISE.
# Its -0.26/-0.37 was written into CLAUDE.md and memory as "the posterior does not
# explore that ridge" in a paragraph evidenced by Yasso15 -- a model it never touched.
# manuscript/figures/build_F14_mrt_ridge.R closes the gap (per-draw intrinsic MRT
# paired to per-draw sigma_input) and finds the OPPOSITE for Yasso: corr -0.79/-0.73/
# -0.57, ratio 0.62/0.73/0.92 = a real ridge. Not a run effect. Quote section 1 for
# SP1/TP2/TP3 ONLY; for the Yasso family cite F14.
#
# ⚠⚠ 2026-08-13, LATER THE SAME DAY -- AND SECTION 1 IS WRONG FOR THE SIMPLE MODELS TOO.
# Its closed forms are MRT AT xi = 1, i.e. with the climate response switched OFF. All
# three models put xi on every pool rate (TP3 since C2), so MRT = MRT_ref/xi EXACTLY --
# the same structure as Yasso07 -- and the calibration MOVES xi (TP2/TP3 medians 2.92/
# 3.21 at the reference climate, against an anchor of 1.00). Put xi back and the trio
# rides a ridge STRONGER than any Yasso model: corr -0.92/-0.87/-0.86, ratio 0.41/0.49/
# 0.51. Second, smaller cause: -0.263/-0.313 come from the posterior .rds, which
# getSample() thins to 1001 draws; the same closed form over the 225015 chain draws
# gives -0.350/-0.362. THE RIDGE IS UNIVERSAL, NOT A YASSO PROPERTY -- do not quote
# section 1's correlations as evidence that the simple models lack a ridge.
# Correct basis: manuscript/figures/build_S13_mrt_ridge_benchmark.R.
#
# USAGE: Rscript --no-save doublechecks/ridge_test.R
# =============================================================================
suppressMessages(library(BayesianTools))
setwd(Sys.getenv("HIKET_ROOT", "."))

RID <- c(SP1     = "20260810_152914", TP2     = "20260810_152914",
         TP3     = "20260810_152914", Yasso07 = "20260810_152915",
         Yasso15 = "20260810_152916", Yasso20 = "20260810_152917")

resid <- function(M) {
  d <- read.csv(sprintf("Calibration_real_data_transient/diagnostics/%s/%s_residuals_%s.csv",
                        M, M, RID[[M]]))
  d[is.finite(d$residual_log), ]
}
post <- function(M) getSample(readRDS(sprintf(
  "Calibration_real_data_transient/runs/%s_posterior_%s.rds", M, RID[[M]])))

# --- 1. the ridge -------------------------------------------------------------
# ratio -> 0 means a tight ridge IS being explored (product pinned, split free);
# ratio -> ~1 means MRT and sigma_input are pinned essentially independently.
source("Prior_specs/TP2_priors.R"); source("Prior_specs/TP3_priors.R")
cat("\n=== 1. RIDGE (does the posterior explore MRT x sigma_input?) ===\n")
cat(sprintf("%-5s %8s %8s %8s | %9s %8s %11s\n",
            "model", "sd_lMRT", "sd_lSin", "cor", "sd_lProd", "ratio", "prior/post"))
one <- function(lab, MRT, SIN) {
  lM <- log(MRT); lS <- log(SIN)
  cat(sprintf("%-5s %8.4f %8.4f %+8.3f | %9.4f %8.3f %10.1fx\n", lab,
              sd(lM), sd(lS), cor(lM, lS), sd(lM + lS),
              sd(lM + lS)/sd(lM), 0.50/sd(lS)))     # 0.50 = sigma_input prior log-SD
}
s <- post("SP1"); one("SP1", 1/s[, "alpha"], s[, "sigma_input"])
s <- post("TP2"); one("TP2", 1/TP2_ALPHA_A_FIXED + s[, "p_H"]/s[, "alpha_H"], s[, "sigma_input"])
s <- post("TP3"); one("TP3", 1/TP3_ALPHA_A_FIXED + s[, "p_S"]/s[, "alpha_S"] +
                        s[, "p_S"]*s[, "p_H"]/s[, "alpha_H"], s[, "sigma_input"])

# --- 2. plot-persistent component --------------------------------------------
# A plot's residual recurs across ITS campaigns because the cause is a static
# site property no model contains. Design effect = 1+(m-1)*ICC deflates n.
cat("\n=== 2. PLOT-PERSISTENT (campaign means removed) ===\n")
cat(sprintf("%-8s %6s %9s %9s %9s | %7s %9s\n",
            "model", "ICC", "r(85,06)", "r(06,24)", "r(85,24)", "deff", "eff_n"))
for (M in names(RID)) {
  d <- resid(M); d$r <- d$residual_log - ave(d$residual_log, d$year)
  ag <- tapply(d$r, d$plot_id, mean); nn <- tapply(d$r, d$plot_id, length)
  gm <- mean(d$r); k <- length(nn); N <- nrow(d)
  msb <- sum(nn*(ag - gm)^2)/(k - 1)
  msw <- sum((d$r - ag[as.character(d$plot_id)])^2)/(N - k)
  m0  <- (N - sum(nn^2)/N)/(k - 1)
  su2 <- max(0, (msb - msw)/m0); ICC <- su2/(su2 + msw)
  ys <- sort(unique(d$year))
  pr <- function(y1, y2) {
    a <- d[d$year == y1, c("plot_id","r")]; b <- d[d$year == y2, c("plot_id","r")]
    m <- merge(a, b, by = "plot_id"); if (nrow(m) < 20) NA else cor(m$r.x, m$r.y)
  }
  deff <- 1 + (N/k - 1)*ICC
  cat(sprintf("%-8s %6.3f %9.3f %9.3f %9.3f | %7.2f %9.0f\n", M, ICC,
              pr(ys[1],ys[2]), pr(ys[2],ys[3]), pr(ys[1],ys[3]), deff, N/deff))
}

# --- 3. spatial component -----------------------------------------------------
# Latitude bands stand in for NFI regions (assign_nfi_regions.R not yet run).
# The implied regional sd barely moves sigma (0.13 is small next to 0.79) but
# inflates the SE of the NATIONAL MEAN ~2.3x -- which is where the ridge lives.
cat("\n=== 3. SPATIAL (latitude bands as region proxy) ===\n")
for (M in names(RID)) {
  d <- resid(M); d <- d[is.finite(d$lat_WGS84), ]
  d$r <- d$residual_log - ave(d$residual_log, d$year)
  d$band <- cut(d$lat_WGS84, quantile(d$lat_WGS84, seq(0, 1, length = 9)),
                include.lowest = TRUE)
  mb <- tapply(d$r, d$band, mean); nb <- tapply(d$r, d$band, length)
  obs <- var(mb); exp <- mean(var(d$r)/nb)
  se_iid <- sd(d$residual_log)/sqrt(nrow(d))
  se_cor <- sqrt(max(0, obs - exp)/length(mb) + var(d$r)/nrow(d))
  cat(sprintf("%-8s excess var %4.1fx | regional sd %.3f | SE(national mean) %.4f -> %.4f (%.1fx)\n",
              M, obs/exp, sqrt(max(0, obs - exp)), se_iid, se_cor, se_cor/se_iid))
}

# --- 4. scale and tails -------------------------------------------------------
# in-sample spread is optimistic; per-campaign spread shows a single sigma is a
# compromise; kurtosis ~7 shows the Gaussian tails are far too thin.
cat("\n=== 4. SCALE AND TAILS ===\n")
cat(sprintf("%-8s %7s %8s %7s | %7s %7s %7s | %8s %9s\n", "model",
            "calib", "holdout", "ratio", "1985", "2006", "2024", "kurtosis", "sd_d0624"))
for (M in names(RID)) {
  d <- resid(M); ys <- sort(unique(d$year)); r <- d$residual_log
  sy <- sapply(ys, function(y) sd(r[d$year == y]))
  a <- d[d$year == ys[2], c("plot_id","residual_log")]
  b <- d[d$year == ys[3], c("plot_id","residual_log")]
  m <- merge(a, b, by = "plot_id")
  cat(sprintf("%-8s %7.3f %8.3f %7.2f | %7.3f %7.3f %7.3f | %8.2f %9.3f\n", M,
              sd(r[!d$is_holdout]), sd(r[d$is_holdout]),
              sd(r[d$is_holdout])/sd(r[!d$is_holdout]), sy[1], sy[2], sy[3],
              mean((r - mean(r))^4)/sd(r)^4,
              sd(m$residual_log.y - m$residual_log.x)))
}
cat("\nsd_d0624 = spread of the 2006->2024 CHANGE: shared errors cancel in a\n",
    "difference, which is why it is ~2x tighter than the raw residuals.\n", sep = "")
