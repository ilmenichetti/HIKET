source("manuscript/figures/run_ids.R")   # auto-selects current RUN_IDs
setwd("/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling")
# =============================================================================
# S11 + S12 -- WHAT THE DATA ACTUALLY MOVED, IN ALL SIX MODELS AT A GLANCE.
#
# S11  KL(posterior || prior) per free parameter, one panel per model.
# S12  Prior (grey) vs posterior (coloured) boxplots, same layout.
#
# WHY BOTH. KL is one number and it CONFLATES two different things: a parameter
# that shifted, and a parameter that merely narrowed on the same centre. Both give
# a large KL, and they mean opposite things -- the first is the data disagreeing
# with the prior, the second is the data confirming it more sharply. S12 separates
# them by eye: displacement of the box shows the shift, its height the narrowing.
#
# These merge and supersede the two half-views we had: F7 covered the three Yassos
# and S5 the three simple models, so no figure ever showed the ensemble together --
# which is what a structural intercomparison actually needs.
#
# STANDARDISATION (S12). Parameters differ by orders of magnitude (beta1 ~ 0.1,
# gamma ~ -1.3, fractions 0-1, sigma_input ~ 2.6), so each is expressed in units of
# ITS OWN PRIOR: z = (x - mean(prior)) / sd(prior), using the physical-space prior
# draws. The grey prior box is then the same reference everywhere and the posterior
# is read directly as "how many prior sigmas did it move, and how much did it tighten".
#
# WHY IT MATTERS HERE (2026-08-13). Intrinsic MRT is a function of the transfer
# fractions, the CLIMATE response and the woody-size terms -- sigma_input is NOT in
# it. mrt_attribution.R shows Yasso07's entire 18.3 yr MRT gap is climate (97%),
# while Yasso15/Yasso20's smaller gaps are mostly fractions (65%, 79%). Yasso07's
# beta1 sits only ~0.45 prior sigma off centre, yet xi doubles and MRT halves --
# exponential leverage. S12 is where that is visible.
#
# Faithful prior reconstruction, identical to run_diagnostics(): source each model's
# real calibration setup for to_original / best_x / sigma_ppm, draw the same seed-99
# prior, push it through to_original. Posteriors are already in physical space.
#
# Run from repo root:  Rscript manuscript/figures/build_S11S12_prior_posterior.R
# =============================================================================

source("manuscript/figures/model_palette.R")
suppressMessages(library(BayesianTools))
RUNID  <- as.list(RID)
MODELS <- c("SP1","TP2","TP3","Yasso07","Yasso15","Yasso20")

setup_model <- function(MODEL) {
  src <- readLines(file.path("Calibration_real_data_transient",
                             sprintf("run_%s_transient_calibration.R", MODEL)), warn = FALSE)
  cut <- grep("^t_run <- system.time\\(\\{", src)[1]
  e   <- new.env(parent = globalenv())
  suppressWarnings(suppressMessages(
    source(textConnection(paste(src[seq_len(cut - 1L)], collapse = "\n")), local = e)))
  to_original <- get("to_original", e); best_x <- get("best_x", e)
  sigma_ppm   <- get("sigma_ppm", e)
  set.seed(99); n_prior <- 3000L
  pr_raw <- sapply(seq_along(best_x), function(j) rnorm(n_prior, best_x[j], sigma_ppm[j]))
  colnames(pr_raw) <- names(best_x)
  list(prior = t(apply(pr_raw, 1, to_original)),
       post  = readRDS(sprintf("Calibration_real_data_transient/runs/%s_posterior_%s.rds",
                               MODEL, RUNID[[MODEL]])))
}

kl_from_kde <- function(post_v, prior_v, n_grid = 512L) {   # verbatim from calibration_engine.R
  post_v <- post_v[is.finite(post_v)]; prior_v <- prior_v[is.finite(prior_v)]
  if (length(post_v) < 10L || length(prior_v) < 10L) return(NA_real_)
  xlp <- quantile(prior_v, c(0.005, 0.995)); xlo <- range(post_v)
  xl  <- c(min(xlp[1], xlo[1]), max(xlp[2], xlo[2]))
  dp <- density(post_v, from = xl[1], to = xl[2], n = n_grid)
  dq <- density(prior_v, from = xl[1], to = xl[2], n = n_grid)
  p <- pmax(dp$y, 1e-10); q <- pmax(dq$y, 1e-10)
  sum(p * log(p / q)) * diff(dp$x[1:2])
}

# Classes and colours follow the engine's own classify_param()/class_cols() so this
# figure matches the rest of the pipeline -- with ONE deliberate refinement: the
# engine lumps "Climate & size", and we split it. The MRT attribution
# (doublechecks/mrt_attribution.R) shows Yasso07's entire MRT gap is the CLIMATE
# response while the woody-size terms contribute 2%, so merging them would hide the
# result this figure exists to show.
pclass <- function(n)
  ifelse(grepl("^sigma", n),              "Auxiliary (input & init)",
  ifelse(grepl("^p_", n),                 "Transfer fraction",
  ifelse(grepl("^alpha", n),              "Decomposition rate",
  ifelse(n %in% c("delta1","delta2","r","w1","w2","w3","w4","w5"), "Woody size",
  ifelse(grepl("^(beta|gamma)", n),       "Climate response", "Other")))))
PCOL <- c("Decomposition rate"       = "#e7298a",   # engine colour
          "Transfer fraction"        = "#d95f02",   # engine colour
          "Climate response"         = "#1b9e77",   # engine "Climate & size"
          "Woody size"               = "#66c2a5",   # lighter, split out
          "Auxiliary (input & init)" = "#7570b3",   # engine colour
          "Other"                    = "grey60")
HILITE <- "Climate response"

cat("Sourcing calibration setups for all six models...\n")
D <- list()
for (M in MODELS) {
  s  <- setup_model(M)
  sm <- getSample(s$post)
  nm <- intersect(colnames(sm), colnames(s$prior))
  kl <- vapply(nm, function(n) kl_from_kde(sm[, n], s$prior[, n]), numeric(1))
  D[[M]] <- list(nm = nm, kl = setNames(kl, nm), post = sm[, nm, drop = FALSE],
                 prior = s$prior[, nm, drop = FALSE])
  cat(sprintf("  %-9s %2d params | max KL %.1f (%s)\n", M, length(nm),
              max(kl, na.rm = TRUE), nm[which.max(kl)]))
}

lab_axis <- function(bp, nm, cl, ylo, cex = 0.72) {
  text(bp, ylo, nm, srt = 90, adj = 1, xpd = NA, cex = cex,
       col = ifelse(cl == HILITE, PCOL[[HILITE]], "grey25"),
       font = ifelse(cl == HILITE, 2, 1))
}

# ---- one figure per model FAMILY: KL on top, prior-vs-posterior below --------
FAM <- list(
  benchmark = list(models = c("SP1","TP2","TP3"),
                   file = "manuscript/figures/S11_params_benchmark.png",
                   title = "S11   Benchmark models (SP1 / TP2 / TP3)"),
  yasso     = list(models = c("Yasso07","Yasso15","Yasso20"),
                   file = "manuscript/figures/S12_params_yasso.png",
                   title = "S12   Yasso family (Yasso07 / Yasso15 / Yasso20)"))

draw_family <- function(fam) {
  MM <- fam$models
  png(fam$file, width = 13.5, height = 9.6, units = "in", res = 200)
  par(mfrow = c(2,3), mar = c(7.6,4.6,3.2,0.9), mgp = c(2.9,0.7,0), las = 1, oma = c(0,0,3.0,0))
  ymax <- max(unlist(lapply(D[MM], `[[`, "kl")), na.rm = TRUE) * 1.10

  for (M in MM) {                                   # --- top row: KL ---
    k <- sort(D[[M]]$kl, decreasing = TRUE); cl <- pclass(names(k))
    bp <- barplot(k, col = PCOL[cl], border = NA, ylim = c(0, ymax),
                  names.arg = rep("", length(k)), ylab = "KL(posterior || prior)  [nats]",
                  main = sprintf("%s  --  information gain (max %.1f nats)", M, max(k, na.rm = TRUE)))
    lab_axis(bp, names(k), cl, -ymax * 0.02)
    abline(h = 1, lty = 2, col = "grey35")
    if (M == MM[1]) legend("topright", bty = "n", cex = 0.80, fill = PCOL[names(PCOL) != "Other"],
                           border = NA, legend = names(PCOL)[names(PCOL) != "Other"])
  }
  for (M in MM) {                                   # --- bottom row: boxes ---
    nm <- D[[M]]$nm[order(D[[M]]$kl, decreasing = TRUE)]
    cl <- pclass(nm); n <- length(nm)
    z <- lapply(nm, function(q) {
      pr <- D[[M]]$prior[, q]; mu <- mean(pr); sdv <- sd(pr)
      if (!is.finite(sdv) || sdv == 0) sdv <- 1
      list(pr = (pr - mu)/sdv, po = (D[[M]]$post[, q] - mu)/sdv) })
    yr <- range(unlist(lapply(z, function(q) quantile(c(q$pr, q$po), c(0.01, 0.99)))))
    yr <- c(min(-3, yr[1]), max(3, yr[2]))
    plot(NA, xlim = c(0.4, n + 0.6), ylim = yr, xaxt = "n", xlab = "",
         ylab = "standardised by its OWN prior   (prior sigma)",
         main = sprintf("%s  --  prior (grey) vs posterior", M))
    abline(h = 0, col = "grey70"); abline(h = c(-2,2), col = "grey88", lty = 3)
    for (i in seq_len(n)) {
      boxplot(z[[i]]$pr, at = i, add = TRUE, boxwex = 0.78, outline = FALSE, axes = FALSE,
              col = adjustcolor("grey55", 0.40), border = "grey45", medlwd = 1)
      boxplot(z[[i]]$po, at = i, add = TRUE, boxwex = 0.44, outline = FALSE, axes = FALSE,
              col = adjustcolor(PCOL[cl[i]], 0.85), border = PCOL[cl[i]], medlwd = 2)
    }
    lab_axis(seq_len(n), nm, cl, yr[1] - diff(yr)*0.035)
    if (M == MM[1]) legend("topleft", bty = "n", cex = 0.78, fill = adjustcolor("grey55", 0.40),
                           border = "grey45", legend = "prior")
  }
  mtext(paste0(fam$title,
        "   |   top: information gain (dashed = 1 nat)   |   bottom: box displacement = shift, box height = tightening"),
        outer = TRUE, side = 3, line = 0.6, cex = 0.84, font = 2)
  dev.off(); cat("Wrote ", fam$file, "\n", sep = "")
}

for (f in FAM) draw_family(f)

saveRDS(D, "manuscript/figures/S11S12_prior_posterior.rds")

# --- console: shift and tightening for the two MRT-relevant classes ----------
cat("\nMedian |shift| in prior sigma, and posterior/prior SD ratio:\n")
cat(sprintf("%-9s %-26s %-26s\n", "model", "climate response", "transfer fraction"))
for (M in MODELS) {
  nm <- D[[M]]$nm; cl <- pclass(nm)
  f <- function(c0) {
    ii <- which(cl == c0); if (!length(ii)) return("-")
    sh <- vapply(ii, function(i) { pr <- D[[M]]$prior[,nm[i]]
      (median(D[[M]]$post[,nm[i]]) - mean(pr))/sd(pr) }, numeric(1))
    ra <- vapply(ii, function(i) sd(D[[M]]$post[,nm[i]])/sd(D[[M]]$prior[,nm[i]]), numeric(1))
    sprintf("shift %+.2f  width x%.2f", median(sh), median(ra)) }
  cat(sprintf("%-9s %-26s %-26s\n", M, f("climate response"), f("transfer fraction")))
}
