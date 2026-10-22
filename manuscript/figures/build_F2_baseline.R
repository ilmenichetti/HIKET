# =============================================================================
# build_F2_baseline.R   (rebuilt 2026-10-22)
#
# F2 -- the published model before any calibration. YASSO07 (the model of the
# Finnish inventory; Lorenzo, 2026-10-22) at its PUBLISHED parameter values, started in balance, with the litter of the product and no
# input correction, against the three campaigns. Two failures in one glance:
#   (1) FLAT -- a soil started in balance has no deficit to recover from;
#   (2) OFF LEVEL -- the level is the product's litter times the published transit
#               time: Yasso07 (MTT 33.5 yr) sits ABOVE the soil (+22 -> +12 tC/ha),
#               Yasso20 (MTT 22 yr) BELOW it (-11 -> -21). The gap closes or widens only
#               because the soil gains and the model does not.
#
# WHAT CHANGED (2026-10-22). The old version plotted a June prediction file
# (diagnostics/Yasso20_baseline/, 447 plots, the pre-August SOC target and litter,
# a parameter vector assembled in Prior_specs). This version runs on the current
# production setup, through the same machinery as the published arm of F15
# (doublechecks/published_arm_dat.R: setup()):
#   * model: HIKET_F2_MODEL (default Yasso07; Yasso20 also supported);
#   * parameters: Yasso07 -- its default vector, which IS the published one (Tuomi et
#     al. 2009, 2011; no .dat exists), as in F12 and the published arm of F15.
#     Yasso20 -- ParY20.dat, the FMI MAP (Viskari et al. 2022), with full override and
#     the published arm's budget clamp -- NOT to_original(best_x), which for Yasso20 is
#     a hybrid with Yasso15's fractions;
#   * start: sigma_init = 1, i.e. the balanced start (steady state under the 1985-89
#     litter and the run's own 1961-1984 pre-run climate);
#   * input: sigma_input = 1, the product's litter as it is;
#   * plots: the plots measured in all three campaigns, unweighted -- the basis of
#     the observations it is compared with (obs_basis.R).
# Usage: Rscript manuscript/figures/build_F2_baseline.R   (project root)
# =============================================================================

setwd("/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling")
suppressWarnings(suppressMessages(library(BayesianTools)))
src <- readLines("doublechecks/published_arm_dat.R")
cut <- grep("^MR <- readRDS", src)[1]
eval(parse(text = paste(src[seq_len(cut - 1)], collapse = "\n")), envir = globalenv())
source("manuscript/figures/obs_basis.R")

M  <- Sys.getenv("HIKET_F2_MODEL", "Yasso07")
stopifnot(M %in% c("Yasso07", "Yasso20"))
e  <- setup(M)
p_def <- e$.to_original(get("best_x", e))
if (M == "Yasso20") {
  pub <- setNames(as.numeric(suppressWarnings(read.table(file.path(DAT, "ParY20.dat")))[1, ]), DCOL)
  p   <- p_def; nm <- intersect(names(p), DCOL); p[nm] <- pub[nm]
  mp  <- e$.assemble_params(p)
  for (v in intersect(names(mp), DCOL)) mp[v] <- pub[[v]]   # full override, incl. fixed rates
  mp  <- clamp_budget(mp)
} else {
  mp  <- e$.assemble_params(p_def)                         # Yasso07 defaults = published
}
mp["sigma_init"] <- 1; mp["sigma_input"] <- 1

om   <- get("obs_meta", e)
BAL  <- balanced_plots(om)
ibp <- get("inputs_by_plot", e); cbp <- get("climate_by_plot", e); lms <- get("litter_means", e)
pl  <- intersect(names(ibp), BAL)            # all 310, held-out plots included, as the observations
tr <- parallel::mclapply(pl, function(pid) {
  cl <- cbp[[pid]]; inp <- ibp[[pid]]
  xa <- tryCatch(e$.compute_xi(cl, mp), error = function(z) NULL); if (is.null(xa)) return(NULL)
  xs <- tryCatch(e$.compute_xi_mean(hiket_ss_clim(pid, cl, min(e$.n_ss, nrow(cl))), mp),
                 error = function(z) NULL); if (is.null(xs)) return(NULL)
  C0 <- tryCatch(e$.steady_state(mp, lms[[pid]], xs), error = function(z) NULL)
  if (is.null(C0) || any(!is.finite(C0))) return(NULL)
  ro <- tryCatch(e$.run_model(inp, mp, C0, xa), error = function(z) NULL); if (is.null(ro)) return(NULL)
  t <- ro$total_soc; if (any(!is.finite(t))) return(NULL)
  data.frame(year = inp$year, soc = t)
}, mc.cores = NCORE)
ok <- !vapply(tr, is.null, logical(1))
d  <- do.call(rbind, tr[ok])
tj <- do.call(rbind, lapply(split(d$soc, d$year), function(x) c(m = mean(x), se = sd(x)/sqrt(length(x)))))
tj <- data.frame(year = as.integer(rownames(tj)), m = tj[, "m"], se = tj[, "se"])
cm <- obs_campaigns(om, BAL)
message("F2 ", basis_note(BAL), "; model run on ", sum(ok), " of ", length(pl), " plots")

png("manuscript/figures/F2_baseline.png", width = 8.4, height = 5.2, units = "in", res = 200)
par(mar = c(4.0, 4.6, 3.0, 1.2), mgp = c(2.6, 0.7, 0), las = 1, cex.axis = 1.05, cex.lab = 1.2, cex.main = 1.1)
YL <- range(tj$m - tj$se, tj$m + tj$se, cm$lo, cm$hi); YL <- YL + c(-0.08, 0.42) * diff(YL)   # headroom for the legend
plot(NA, xlim = c(min(tj$year), max(tj$year) + 2), ylim = YL, xlab = "Year",
     ylab = expression("Mean SOC across plots (tC ha"^-1*")"),
     main = sprintf("%s before calibration: published parameters, balanced start", M))
mb <- approx(tj$year, tj$m, cm$year)$y
lab_side <- ifelse(cm$year == max(cm$year), 2, 4)
segments(cm$year, mb, cm$year, cm$m, col = "grey65", lwd = 8, lend = 1)
text(cm$year, (mb + cm$m)/2, sprintf("%+.0f", mb - cm$m), pos = lab_side, offset = 0.5,
     cex = 0.8, col = "grey35", font = 2)
polygon(c(tj$year, rev(tj$year)), c(tj$m - tj$se, rev(tj$m + tj$se)),
        col = adjustcolor("steelblue", 0.30), border = NA)
lines(tj$year, tj$m, col = "steelblue", lwd = 2.6)
arrows(cm$year, cm$lo, cm$year, cm$hi, angle = 90, code = 3, length = 0.05, col = "firebrick", lwd = 1.8)
points(cm$year, cm$m, pch = 19, col = "firebrick", cex = 1.5)
legend("topleft", inset = c(0.01, 0.02), bty = "n", cex = 0.85,
       legend = c(sprintf("%s, published parameters, balanced start, product litter (mean +/- SE, %d plots)", M, sum(ok)),
                  "Observed campaign mean +/- 95% CI",
                  "Model-observation gap"),
       pch = c(NA, 19, NA), lwd = c(2.6, NA, 7), col = c("steelblue", "firebrick", "grey65"))
mtext(basis_note(BAL), side = 1, line = 2.9, adj = 0, cex = 0.6, col = "grey40")
dev.off()
cat("Wrote manuscript/figures/F2_baseline.png\n")
cat(sprintf("Model: %.1f (1985) .. %.1f (2006) .. %.1f (2024);  obs %.1f -> %.1f -> %.1f;  gaps %s\n",
            tj$m[tj$year == 1985], tj$m[tj$year == 2006], tj$m[tj$year == max(tj$year)],
            cm$m[1], cm$m[2], cm$m[3], paste(sprintf("%+.1f", mb - cm$m), collapse = " / ")))
