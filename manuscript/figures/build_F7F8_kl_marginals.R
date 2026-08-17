source("manuscript/figures/run_ids.R")   # auto-selects current RUN_IDs
setwd("/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling")
# F7 + F8 rebuild -- faithful KL + marginal recompute.
# The saved posterior (runs/<MODEL>_posterior_<RUNID>.rds) is already PHYSICAL-space
# (all_phys). We reproduce the prior EXACTLY as the engine does by sourcing each model's
# real calibration setup (preflight-style) to obtain to_original / best_x / sigma_ppm,
# then draw the same seed-99 prior and push it through to_original. We reuse the engine's
# own plot_one_marginal_honest / classify_param / class_cols so styling matches the pipeline.
#   F7 -> merged 3-panel KL barplot, the three Yasso models (shared y-axis).
#   F8 -> readable fraction marginals for BOTH TP3 and Yasso20 (user picks which to feature).

RUNID <- unlist(RID[c("TP3","Yasso07","Yasso15","Yasso20")])

# --- source a model's calibration setup (up to MCMC launch) into a fresh env -------
setup_model <- function(MODEL) {
  sp  <- file.path("Calibration_real_data_transient",
                   sprintf("run_%s_transient_calibration.R", MODEL))
  src <- readLines(sp, warn = FALSE)
  cut <- grep("^t_run <- system.time\\(\\{", src)[1]
  e   <- new.env(parent = globalenv())
  suppressWarnings(suppressMessages(
    source(textConnection(paste(src[seq_len(cut - 1L)], collapse = "\n")), local = e)))
  to_original <- get("to_original", e); best_x <- get("best_x", e)
  sigma_ppm   <- get("sigma_ppm", e);   FREE   <- get("FREE_NAMES", e)
  # prior draws in physical space -- identical recipe to run_diagnostics()
  set.seed(99); n_prior <- 3000L
  pr_raw <- sapply(seq_along(best_x), function(j) rnorm(n_prior, best_x[j], sigma_ppm[j]))
  colnames(pr_raw) <- names(best_x)
  prior_phys <- t(apply(pr_raw, 1, to_original))
  # posterior (already physical)
  post <- readRDS(sprintf("Calibration_real_data_transient/runs/%s_posterior_%s.rds",
                          MODEL, RUNID[[MODEL]]))
  list(e = e, FREE = FREE, prior = prior_phys, post = post,
       classify = get("classify_param", e), ccols = get("class_cols", e),
       clev = get("class_levels", e), marg = get("plot_one_marginal_honest", e))
}

# KL(posterior || prior) per parameter -- copied verbatim from calibration_engine.R
kl_from_kde <- function(post_v, prior_v, n_grid = 512L) {
  post_v <- post_v[is.finite(post_v)]; prior_v <- prior_v[is.finite(prior_v)]
  if (length(post_v) < 10L || length(prior_v) < 10L) return(NA_real_)
  xlp <- quantile(prior_v, c(0.005, 0.995)); xlo <- range(post_v)
  xl  <- c(min(xlp[1], xlo[1]), max(xlp[2], xlo[2]))
  dp  <- density(post_v,  from = xl[1], to = xl[2], n = n_grid)
  dq  <- density(prior_v, from = xl[1], to = xl[2], n = n_grid)
  p <- pmax(dp$y, 1e-10); q <- pmax(dq$y, 1e-10)
  sum(p * log(p / q)) * diff(dp$x[1:2])
}

cat("Setting up models (sourcing calibration setups)...\n")
M <- lapply(c("Yasso07","Yasso15","Yasso20","TP3"), setup_model)
names(M) <- c("Yasso07","Yasso15","Yasso20","TP3")

# ============================ F7: three-Yasso KL ============================
kl_of <- function(mm) {
  vapply(mm$FREE, function(nm)
    kl_from_kde(mm$post[, nm], mm$prior[, nm]), numeric(1))
}
kls   <- lapply(M[c("Yasso07","Yasso15","Yasso20")], kl_of)

# Per-panel y-axis (each scaled to its own max) so short-bar panels aren't mostly
# empty; the dashed 1-nat line anchors the reading. Larger text, tight margins.
png("manuscript/figures/F7_kl_three_yasso.png", width = 11, height = 8, units = "in", res = 200)
par(mfrow = c(3, 1), mar = c(5.0, 4.8, 2.4, 0.8), mgp = c(2.9, 0.7, 0),
    oma = c(2.4, 0, 0, 0), cex.axis = 1.0, cex.lab = 1.2)
for (mn in c("Yasso07","Yasso15","Yasso20")) {
  mm <- M[[mn]]; kv <- kls[[mn]]
  cl <- vapply(names(kv), mm$classify, character(1)); bcol <- mm$ccols[cl]
  barplot(kv, col = bcol, border = NA, las = 2, cex.names = 0.95,
          ylim = c(0, max(kv, na.rm = TRUE) * 1.15),
          ylab = "KL (nats)", main = mn, font.main = 1, cex.main = 1.35)
  abline(h = 1, lty = 2, col = "grey40", lwd = 1.3)
  if (mn == "Yasso07") {
    present <- mm$clev[mm$clev %in% cl]
    legend("topright", legend = present, fill = mm$ccols[present], border = NA,
           bty = "n", cex = 1.05, title = "Parameter class")
  }
}
mtext("KL(posterior || prior) per parameter; note y-axes differ",
      side = 1, line = 0.8, outer = TRUE, cex = 0.85, col = "grey30")
dev.off()
cat("wrote F7_kl_three_yasso.png\n")

# ============================ F8: readable marginals ============================
make_marginals <- function(mn, params, file, ncol, w, h) {
  mm <- M[[mn]]; params <- intersect(params, mm$FREE)
  nr <- ceiling(length(params) / ncol)
  png(file, width = w, height = h, units = "in", res = 200)
  par(mfrow = c(nr, ncol), mar = c(4, 4, 2.6, 1))
  for (nm in params) {
    pv <- mm$post[is.finite(mm$post[, nm]), nm]
    qv <- mm$prior[is.finite(mm$prior[, nm]), nm]
    mm$marg(pv, qv, nm, mn)
  }
  dev.off(); cat("wrote", basename(file), "\n")
}
yfrac <- grep("^p_[A-Z][A-Z]$", M$Yasso20$FREE, value = TRUE)          # 12 transfer fractions
tp3par <- intersect(c("alpha_A","alpha_S","alpha_H","p_S","p_H"), M$TP3$FREE)  # TP3 kinetics
make_marginals("Yasso20", yfrac, "manuscript/figures/F8_marginals_Yasso20.png", ncol = 4, w = 11, h = 8)
make_marginals("TP3",     tp3par, "manuscript/figures/F8_marginals_TP3.png",     ncol = 3, w = 9,  h = 6)
cat("F8 params: TP3 =", paste(tp3par, collapse=","), " | Yasso20 =", paste(yfrac, collapse=","), "\n")
