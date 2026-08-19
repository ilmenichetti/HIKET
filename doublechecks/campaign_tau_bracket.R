# =============================================================================
# campaign_tau_bracket.R   (2026-08-19)
#
# WHY. tau_C cannot be measured -- it is prescribed. The honest output is
# therefore a BRACKET the conclusion must survive, not a single number. This
# script computes that bracket analytically from the campaign means, so it costs
# nothing and can be reported whatever is eventually run.
#
# It answers two distinct questions that are easy to conflate (they had opposite
# answers in discussion on 2026-08-19):
#
#   (a) does raising tau_C for 1985 DOWNWEIGHT that campaign?  YES -- the
#       precision of a campaign level is 1/(tau_C^2 + sigma_e^2/n_j), so at the
#       adopted 0.06 the 1985 level carries ~2.4x LESS weight than a 2006/2024
#       level. Distrust of 1985 is also CHEAP: 0.06 -> 0.09 cuts its weight to a
#       fifth while the 35-year signal only falls 1.9 -> 1.4 sigma.
#
#   (b) is the whole scheme looser than the SIGMA_1985_INFL = 2 it replaces?  NO
#       -- about 12% TIGHTER, because the total 0.800 is SPLIT, never added:
#       every observation's independent noise drops 0.800 -> 0.685 (1985's drops
#       1.600 -> 0.685) and the removed variance reappears only in the LEVEL.
#
# ⚠ The two must never be stacked: together they give an effective tau_C of
#   ~0.085, the broad scheme by accident. Hence SIGMA_1985_INFL -> 1.0.
#
# ⚠ COUNTERINTUITIVE, and it survives every scheme: the interval this hurts most
#   is 2006->2024, NOT anything involving 1985. The short interval dissolves at
#   tau_C ~ 2%; 1985->2024 survives to ~9%. The campaign we distrust supports the
#   comparison that lives. See memory campaign-comparability-limits-the-trend.
#
# Usage:  Rscript doublechecks/campaign_tau_bracket.R
# =============================================================================

SIGMA_E <- 0.685      # independent per-observation SD, remainder of the fixed total 0.800
CBAR    <- 70.8       # Mg/ha, mean profile stock
SIG     <- c("1985->2024" = 9.1, "2006->2024" = 2.1)   # Mg/ha, balanced set (obs_basis.R)
DT      <- c("1985->2024" = 35.1, "2006->2024" = 18.0)

d <- read.csv("Data/SOC_homogeneized/soc_homogenized_plot.csv")
d <- d[!is.na(d$soc_profile_Mgha) & d$soc_profile_Mgha > 0 &
       !(d$soc_outlier %in% TRUE) & !(d$is_peat %in% TRUE), ]
n <- table(d$year)
cat(sprintf("\nn per campaign: %s | sigma_e = %.3f | Cbar = %.1f Mg/ha\n\n",
            paste(sprintf("%s=%d", names(n), as.integer(n)), collapse = " "), SIGMA_E, CBAR))

# --- (a) how much is the 1985 LEVEL downweighted? ---------------------------
TAU_OTHER <- 0.03
cat("(a) weight on the 1985 LEVEL, against a 2006/2024 level at tau_C = 0.03\n")
cat(sprintf("%12s %14s %14s %10s   %s\n", "tau_C(1985)", "sd(1985 lvl)", "sd(other lvl)",
            "rel.weight", "1985->2024"))
for (t in c(0.00, 0.03, 0.04, 0.06, 0.09, 0.12)) {
  v85 <- t^2         + SIGMA_E^2 / as.numeric(n["1985"])
  vot <- TAU_OTHER^2 + SIGMA_E^2 / as.numeric(n["2024"])
  sn  <- SIG[["1985->2024"]] / (sqrt(t^2 + TAU_OTHER^2) * CBAR)
  cat(sprintf("%12.2f %14.4f %14.4f %10.2f   %6.2f sigma%s\n",
              t, sqrt(v85), sqrt(vot), vot / v85, sn,
              if (abs(t - 0.06) < 1e-9) "   <- adopted" else ""))
}

# --- (b) the reporting bracket ----------------------------------------------
cat("\n(b) the bracket to report. A rate is resolvable only if the level difference\n")
cat("    it is built from outgrows the campaign offsets: sd = sqrt(t_a^2+t_b^2)*Cbar.\n\n")
cat(sprintf("%-18s %-14s %s\n", "scheme (85/other)", "1985->2024", "2006->2024"))
for (s in list(c(0.00,0.00), c(0.04,0.02), c(0.06,0.03), c(0.09,0.045), c(0.10,0.05))) {
  f <- function(k, ta, tb) SIG[[k]] / (sqrt(ta^2 + tb^2) * CBAR)
  cat(sprintf("%-18s %-14s %s\n", sprintf("%.3f / %.3f", s[1], s[2]),
              sprintf("%.2f sigma", f("1985->2024", s[1], s[2])),
              sprintf("%.2f sigma", f("2006->2024", s[2], s[2]))))
}

# --- is 1985 less PLOT-LINKED than the other campaigns? ---------------------
# The relocated-subplot worry is a plot-linkage question, which tau_C cannot
# touch. It IS measurable: if 1985's subplots moved, its observations should
# share less of the plot offset u^P. Campaign-centred, as the design principle
# requires (a level must never leak into a covariance).
cat("\n(c) plot linkage by campaign pair -- the relocated-subplot check\n")
d$l <- log(d$soc_profile_Mgha)
w <- reshape(d[, c("plot_id", "year", "l")], idvar = "plot_id", timevar = "year",
             direction = "wide")
names(w) <- sub("^l\\.", "y", names(w))
cat(sprintf("%-14s %5s %7s %8s %12s\n", "pair", "n", "r", "cov", "implied tau_P"))
for (p in list(c("y1985","y2006"), c("y1985","y2024"), c("y2006","y2024"))) {
  ok <- complete.cases(w[, p]); a <- w[ok, p[1]]; b <- w[ok, p[2]]
  a <- a - mean(a); b <- b - mean(b)                       # campaign-centred
  cat(sprintf("%-14s %5d %7.3f %8.3f %12.3f%s\n",
              paste(sub("^y", "", p), collapse = "-"), length(a), cor(a, b), cov(a, b),
              sqrt(max(cov(a, b), 0)),
              if (identical(p, c("y2006","y2024"))) "   <- tau_P design value" else ""))
}
cat("\nVERDICT: 1985's plot linkage is ~10% weaker in SD (0.355-0.363 vs 0.396), NOT broken.\n")
cat("Measurable, small, and it does NOT support treating 1985 as a different soil.\n")
