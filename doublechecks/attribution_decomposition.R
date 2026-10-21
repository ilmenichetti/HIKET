# =============================================================================
# attribution_decomposition.R   (2026-10-08, target revised the same day)
#
# ATTRIBUTION OF THE MODELLED SINK TO HISTORICAL CARBON, LITTER INPUTS AND CLIMATE,
# 1917-2084, for each of the six calibrated models (production run).
# Methodology for review: manuscript/appendices/appendix_attribution.tex.
#
# THE QUANTITY: the annual sink F_t = C_t - C_{t-1} itself (not its change), split
# EXACTLY into three parts that add up to it every year:
#
#   historical C  the INITIALISATION story.
#                 1917-1984: the whole pre-run sink -- the soil starts at equilibrium
#                   with the 1917 litter and the climate is held at its 1985-2004
#                   mean, so every bit of sink is the soil catching up with the rise
#                   of litter from the 1917 to the 1985 level.
#                 1985-2084: the sink the soil owes to the deficit it carries into
#                   1985 = sink of the production start MINUS sink of an equilibrium
#                   start (sigma_init = 1, same parameters; a reference state) under the SAME
#                   litter and climate. Exact because the models are linear in
#                   carbon: the difference is the free relaxation of the deficit.
#   litter and climate  SEQUENTIAL split (Lorenzo, 2026-10-09; replaces the Shapley
#                 average). Three runs carry the partition:
#                   run 1  calibrated start, observed litter, observed climate (the model)
#                   run 2  balanced start,   observed litter, observed climate
#                   run 3  balanced start,   observed litter, REFERENCE climate
#                 history = 1 - 2 ; climate = 2 - 3 ; litter = 3 (- run 5, which is 0)
#                 Order: litter first, then climate, so the litter x climate joint effect
#                 (the extra litter decomposing faster in warmer years) goes to CLIMATE --
#                 climate acting on the litter that actually arrived.
#                 Run 5 (balanced, both references) is kept as the check that the balanced
#                 state is balanced; run 4 (reference litter, observed climate) is kept ONLY
#                 to report the size of the joint effect (column 'interaction').
# Climate = all of a model's climate modifiers together (ONE band, Yasso15/20 too).
#
# Every run is the model's OWN multi-year run (calibration-script engine functions).
# CHECKS per model: (1) the pre-run replayed with the run function reproduces the
# model's initial state; (2) the equilibrium start stays at equilibrium under
# reference forcing; (3) the three parts add up to the sink.
#
# SCOPE: balanced plot set (all three campaigns), unweighted national mean;
# projection 2025-2084 exactly as the predictive stage (last 20 climate years
# recycled, litter held at 2024). N_DRAW posterior draws per model.
#
# Usage:  Rscript doublechecks/attribution_decomposition.R [N_DRAW] [MODELS...]
#         HIKET_ATTR_CORES (default 4)
#         HIKET_ATTR_PARAMS=eq  -> PARAMETER-SOURCE TEST (2026-10-13): same runs, same
#           definitions, but the draws come from the BALANCED-START calibration (eqinit
#           posteriors, sigma_init = 1). Inherited part = 0 by construction; the climate and
#           litter parts then show how much they depend on the calibration's initial-state
#           assumption (it enters ONLY through the parameters -- no run is reparameterised).
#         HIKET_ATTR_MODE=recalibrated -> THE RECALIBRATED DECOMPOSITION (2026-10-13): every
#           scenario carries its OWN calibration (Calibration_real_data_transient/
#           scenario_switches.R), so the parts are a comparison between assumptions only:
#             S1  production                  (estimated start, observed litter and climate)
#             S2  <model>_eqinit              (balanced start)
#             S3  <model>_eqinit_climheld     (balanced start, reference climate)
#             S4  <model>_eqinit_litheld      (balanced start, litter held)     [five terms]
#             S1' <model>_climheld            (estimated start, reference climate) [five terms]
#           history = S1 - S2, climate = S2 - S3, litter = S3; joint effects when S4 / S1'
#           exist (litter x climate = S2 - S3 - S4; history x climate = (S1-S2) - (S1'-S3)),
#           NA otherwise. Each arm's sink is ITS OWN model under the conditions that define
#           it; draws are independent between arms. Same structure as attribution.rds.
# Output: doublechecks/attribution.rds, attribution_eqparams.rds with HIKET_ATTR_PARAMS=eq,
#         or attribution_recalibrated.rds with HIKET_ATTR_MODE=recalibrated
# Setup side effects (input bundle, sanity PNG) are removed at the end.
# =============================================================================

setwd("/Users/ilmenichetti/Library/CloudStorage/OneDrive-Valtion/HIKET/SOC_modeling")
a      <- commandArgs(trailingOnly = TRUE)
N_DRAW <- if (length(a) >= 1) as.integer(a[[1]]) else 50L
MODELS <- if (length(a) >= 2) a[-1] else c("SP1", "TP2", "TP3", "Yasso07", "Yasso15", "Yasso20")
CORES  <- as.integer(Sys.getenv("HIKET_ATTR_CORES", "4"))
PROJ_YEARS <- 60L; RECYCLE_YEARS <- 20L; PREINIT_YEAR <- 1917L; N_PRE <- 68L
Sys.setenv(HIKET_EQUILIBRIUM_INIT = "0")         # production arm
set.seed(2025)

snap <- function() c(list.files("Data/model_inputs", full.names = TRUE),
                     list.files("Calibration_real_data_transient/diagnostics", recursive = TRUE,
                                full.names = TRUE, include.dirs = TRUE))
before <- snap()
suppressMessages({ source("manuscript/figures/run_ids_eqinit.R"); source("manuscript/figures/obs_basis.R") })
EQP <- identical(Sys.getenv("HIKET_ATTR_PARAMS", ""), "eq")
if (EQP) { RID <- RID_EQ; cat("PARAMETER SOURCE: balanced-start calibration (eqinit posteriors, sigma_init = 1)\n") }
PTAG <- if (EQP) "_eqinit" else ""
cat("RUN_IDs:", paste(names(RID[MODELS]), RID[MODELS], sep = "=", collapse = " | "), "\n")

# --- the real calibration script up to the MCMC launch: its own engine functions ----
setup_model <- function(M) {
  src <- readLines(sprintf("Calibration_real_data_transient/run_%s_transient_calibration.R", M), warn = FALSE)
  cut <- grep("^t_run <- system.time", src)[1]
  e <- new.env(parent = globalenv())
  suppressMessages(capture.output(
    source(textConnection(paste(src[seq_len(cut - 1)], collapse = "\n")), local = e)))
  st <- grep("ll_fn <- make_likelihood", src)[1]; op <- 0L; en <- NA_integer_
  for (i in seq(st, length(src))) {
    ch <- strsplit(src[i], "")[[1]]; op <- op + sum(ch == "(") - sum(ch == ")")
    if (op == 0L) { en <- i; break }
  }
  ml <- as.list(str2lang(sub("^\\s*ll_fn\\s*<-\\s*", "", paste(src[seq(st, en)], collapse = "\n"))))[-1]
  for (k in c("assemble_params", "compute_xi", "compute_xi_mean", "steady_state", "run_model"))
    e[[paste0(".", k)]] <- eval(ml[[k]], envir = e)
  e$.assemble <- e$.assemble_params
  e$.n_ss <- eval(ml[["steady_state_n"]], envir = e)
  e
}

xi_rep <- function(x, n) if (is.list(x)) lapply(x, function(v) rep(unname(v)[1], n)) else rep(unname(x)[1], n)
xi_cat <- function(x, y) if (is.list(x)) Map(c, x, y) else c(x, y)
LIT <- c(paste0("nwl_", c("A", "W", "E", "N")), paste0("fwl_", c("A", "W", "E", "N")), paste0("cwl_", c("A", "W", "E", "N")))

# Litter rows at a fraction f of the way from the 1917 to the 1985 level (f = 1: the
# 1985-89 level). Unscaled by sigma_input -- the run function applies it, as in the
# calibration. Mirrors *_transient_init exactly.
litter_rows <- function(lm, template, f, sinit, years) {
  r <- template[rep(1L, length(years)), , drop = FALSE]; r$year <- years; rownames(r) <- NULL
  w <- sinit * (1 - f) + f
  if ("J_total" %in% names(r)) r$J_total <- lm$J_t0_mean * w
  else for (s in c("nwl", "fwl", "cwl")) {
    v <- lm[[paste0(s, "_t0_mean")]]
    for (k in c("A", "W", "E", "N")) r[[paste0(s, "_", k)]] <- unname(v[[paste0(s, "_", k)]]) * w
  }
  r
}

attribute_plot <- function(e, mp, pid) {
  clim <- e$climate_by_plot[[pid]]; inp <- e$inputs_by_plot[[pid]]; lm <- e$litter_means[[pid]]
  ny <- nrow(inp); rpy <- nrow(clim) %/% ny
  run   <- function(U, xi, C) e$.run_model(U, mp, C, xi)
  xi_h  <- e$.compute_xi(clim, mp)
  xi_ss <- e$.compute_xi_mean(hiket_ss_clim(pid, clim, min(e$.n_ss, nrow(clim))), mp)   # pre-run = reference climate
  sinit <- unname(mp[["sigma_init"]])
  C_init <- e$.steady_state(mp, lm, xi_ss)                            # production start (1985)
  mp_eq <- mp; mp_eq["sigma_init"] <- 1
  C_eq  <- e$.steady_state(mp_eq, lm, xi_ss)                          # equilibrium start
  lm0 <- lm; lm0$preinit_shape <- rep(0, N_PRE)
  C_1917 <- e$.steady_state(mp, lm0, xi_ss)                           # 1917 equilibrium
  # --- pre-run, replayed with the run function --------------------------------------
  shape <- if (!is.null(lm$preinit_shape) && length(lm$preinit_shape) == N_PRE) lm$preinit_shape
           else (seq_len(N_PRE) - 1L) / (N_PRE - 1L)
  py0 <- PREINIT_YEAR + seq_len(N_PRE) - 1L
  Upre <- litter_rows(lm, inp, shape, sinit, py0)
  if ("precip" %in% names(Upre)) Upre$precip <- lm$precip_mean
  R_pre <- run(Upre, xi_rep(xi_ss, N_PRE), C_1917)
  chk_pre <- abs(tail(R_pre$total_soc, 1) - sum(C_init)) / sum(C_init)
  F_pre <- diff(c(sum(C_1917), R_pre$total_soc))
  # --- 1985-2084 forcing (projection as the predictive stage) -----------------------
  yp <- max(inp$year) + seq_len(PROJ_YEARS)
  cr <- tail(clim, RECYCLE_YEARS * rpy)
  cp <- cr[((seq_len(PROJ_YEARS * rpy) - 1L) %% (RECYCLE_YEARS * rpy)) + 1L, , drop = FALSE]
  cp$year <- rep(yp, each = rpy); if (rpy == 12L && "month" %in% names(cp)) cp$month <- rep(1:12, PROJ_YEARS)
  rownames(cp) <- NULL
  ip <- inp[rep(ny, PROJ_YEARS), , drop = FALSE]; ip$year <- yp; rownames(ip) <- NULL
  U  <- rbind(inp, ip); N <- nrow(U)
  xi <- xi_cat(xi_h, e$.compute_xi(cp, mp)); xr <- xi_rep(xi_ss, N)
  Ur <- litter_rows(lm, U, rep(1, N), sinit, U$year); if ("precip" %in% names(U)) Ur$precip <- U$precip
  sinkof <- function(R, C0) diff(c(sum(C0), R$total_soc))
  F_tot <- sinkof(run(U,  xi, C_init), C_init)
  F_eq  <- sinkof(run(U,  xi, C_eq),   C_eq)
  F_ux0 <- sinkof(run(U,  xr, C_eq),   C_eq)                          # actual litter, reference climate
  F_u0x <- sinkof(run(Ur, xi, C_eq),   C_eq)                          # reference litter, actual climate
  F_00  <- sinkof(run(Ur, xr, C_eq),   C_eq)                          # both reference: must be ~0
  hist <- F_tot - F_eq
  inpt <- F_ux0 - F_00                                                # run 3 (- run 5 = 0)
  clmt <- F_eq - F_ux0                                                # run 2 - run 3
  inter <- F_eq - F_ux0 - F_u0x + F_00                                # litter x climate interaction
  # inherited x climate (diagnostic only): the inherited part measured under the REFERENCE
  # climate (run 1' = calibrated start, observed litter, reference climate, minus run 3)
  # versus under the observed climate (1 - 2). Linear in litter, so the litter does not matter.
  F_c0  <- sinkof(run(U, xr, C_init), C_init)                         # run 1'
  inter_hc <- hist - (F_c0 - F_ux0)
  out <- rbind(cbind(year = py0, F = F_pre, history = F_pre, inputs = 0, climate = 0, interaction = 0, inter_hist_clim = 0),
               cbind(year = U$year, F = F_tot, history = hist, inputs = inpt, climate = clmt, interaction = inter,
                     inter_hist_clim = inter_hc))
  attr(out, "chk") <- c(prerun = chk_pre, eq_drift = max(abs(F_00)),
                        sum = max(abs(F_tot - F_00 - (hist + inpt + clmt))))
  out
}

# --- RECALIBRATED MODE: one arm = one calibration, run under its own conditions ----------
RECAL <- identical(Sys.getenv("HIKET_ATTR_MODE", ""), "recalibrated")
ARMS <- list(S1  = list(sfx = "",                 eq = FALSE, clim = FALSE, lit = FALSE, need = TRUE),
             S2  = list(sfx = "_eqinit",          eq = TRUE,  clim = FALSE, lit = FALSE, need = TRUE),
             S3  = list(sfx = "_eqinit_climheld", eq = TRUE,  clim = TRUE,  lit = FALSE, need = TRUE),
             S4  = list(sfx = "_eqinit_litheld",  eq = TRUE,  clim = FALSE, lit = TRUE,  need = FALSE),
             S1p = list(sfx = "_climheld",        eq = FALSE, clim = TRUE,  lit = FALSE, need = FALSE))
arm_rid <- function(M, sfx) {
  fs <- list.files("Calibration_real_data_transient/runs",
                   pattern = sprintf("^%s%s_posterior_[0-9]{8}_[0-9]{6}\\.rds$", M, sfx))
  if (!length(fs)) return(NA_character_)
  sub(sprintf("^%s%s_posterior_(.+)\\.rds$", M, sfx), "\\1", sort(fs, decreasing = TRUE)[1])
}
# The arm's own sink 1985-2084 on one plot: start, climate and litter as the arm defines
# them (the same constructions as attribute_plot: xr = pre-run modifier repeated, Ur =
# litter at the 1985-89 level). sigma_init = 1 for a balanced start.
arm_plot <- function(e, mp, pid, A) {
  clim <- e$climate_by_plot[[pid]]; inp <- e$inputs_by_plot[[pid]]; lm <- e$litter_means[[pid]]
  ny <- nrow(inp); rpy <- nrow(clim) %/% ny
  if (A$eq) mp["sigma_init"] <- 1
  xi_h  <- e$.compute_xi(clim, mp)
  xi_ss <- e$.compute_xi_mean(hiket_ss_clim(pid, clim, min(e$.n_ss, nrow(clim))), mp)
  C0 <- e$.steady_state(mp, lm, xi_ss)
  yp <- max(inp$year) + seq_len(PROJ_YEARS)
  cr <- tail(clim, RECYCLE_YEARS * rpy)
  cp <- cr[((seq_len(PROJ_YEARS * rpy) - 1L) %% (RECYCLE_YEARS * rpy)) + 1L, , drop = FALSE]
  cp$year <- rep(yp, each = rpy); if (rpy == 12L && "month" %in% names(cp)) cp$month <- rep(1:12, PROJ_YEARS)
  rownames(cp) <- NULL
  ip <- inp[rep(ny, PROJ_YEARS), , drop = FALSE]; ip$year <- yp; rownames(ip) <- NULL
  U  <- rbind(inp, ip); N <- nrow(U)
  xi <- if (A$clim) xi_rep(xi_ss, N) else xi_cat(xi_h, e$.compute_xi(cp, mp))
  if (A$lit) { Ur <- litter_rows(lm, U, rep(1, N), 1, U$year); if ("precip" %in% names(U)) Ur$precip <- U$precip; U <- Ur }
  R <- e$.run_model(U, mp, C0, xi)
  cbind(year = U$year, F = diff(c(sum(C0), R$total_soc)))
}
arm_sinks <- function(e, M, A, rid, bal) {
  post <- readRDS(sprintf("Calibration_real_data_transient/runs/%s%s_posterior_%s.rds", M, A$sfx, rid))
  if (!"sigma_init" %in% colnames(post)) post <- cbind(post, sigma_init = 1)
  bun <- readRDS(sprintf("Data/model_inputs/%s%s_inputs_%s.rds", M, A$sfx, rid))
  hiket_set_prerun_clim(bun$climate_ss_by_plot, as.character(bun$plots_real))
  if (!is.null(bun$climate_ss_by_plot)) e$litter_means <- bun$litter_means
  ks <- sort(sample(nrow(post), N_DRAW))
  # CHECK (balanced arms): the same draw with litter AND climate held at their references
  # must have no sink on any plot -- the balanced start is in balance under this calibration.
  A0 <- modifyList(A, list(clim = TRUE, lit = TRUE)); drift <- 0
  S <- sapply(ks, function(k) {
    mp <- e$.assemble(setNames(as.numeric(post[k, ]), colnames(post)))
    R  <- parallel::mclapply(bal, function(p) tryCatch(arm_plot(e, mp, p, A), error = function(z) NULL),
                             mc.cores = CORES)
    if (A$eq) {
      D <- parallel::mclapply(bal, function(p) tryCatch(max(abs(arm_plot(e, mp, p, A0)[, "F"])),
                                                       error = function(z) NA_real_), mc.cores = CORES)
      drift <<- max(drift, unlist(D), na.rm = TRUE)
    }
    R  <- R[!vapply(R, is.null, logical(1))]
    Reduce(`+`, lapply(R, function(r) r[, "F"])) / length(R)
  })
  list(S = S, ks = ks, drift = if (A$eq) drift else NA_real_)
}

out <- list(rid = RID[MODELS], params = if (EQP) "eqinit" else "production", n_draw = N_DRAW, proj_years = PROJ_YEARS, version = "level_v3_sequential", models = list())
if (RECAL) {
  out <- list(version = "recalibrated", n_draw = N_DRAW, proj_years = PROJ_YEARS, rid = list(), models = list())
  for (M in MODELS) {
    t0 <- Sys.time(); e <- setup_model(M); bal <- as.character(balanced_plots(e$obs_meta))
    rids <- vapply(ARMS, function(A) arm_rid(M, A$sfx), character(1))
    miss <- names(rids)[is.na(rids) & vapply(ARMS, `[[`, TRUE, "need")]
    if (length(miss)) stop(M, ": recalibrated decomposition needs arms ", paste(miss, collapse = ", "),
                           " (", paste(vapply(ARMS[miss], `[[`, "", "sfx"), collapse = ", "), ")")
    cat(M, "arms:", paste(names(rids), ifelse(is.na(rids), "--", rids), sep = "=", collapse = " | "), "\n")
    SS <- lapply(names(ARMS), function(a) if (is.na(rids[[a]])) NULL else arm_sinks(e, M, ARMS[[a]], rids[[a]], bal))
    names(SS) <- names(ARMS)
    yr <- 1984L + seq_len(nrow(SS$S1$S))                 # 1985 .. 2084
    NAm <- matrix(NA_real_, nrow(SS$S1$S), N_DRAW)
    S1 <- SS$S1$S; S2 <- SS$S2$S; S3 <- SS$S3$S
    S4 <- if (is.null(SS$S4)) NAm else SS$S4$S; S1p <- if (is.null(SS$S1p)) NAm else SS$S1p$S
    arr <- simplify2array(lapply(seq_len(N_DRAW), function(d)
      cbind(year = yr, F = S1[, d], history = S1[, d] - S2[, d], inputs = S3[, d], climate = S2[, d] - S3[, d],
            interaction = S2[, d] - S3[, d] - S4[, d],
            inter_hist_clim = (S1[, d] - S2[, d]) - (S1p[, d] - S3[, d]))))
    chk <- c(eq_drift_S2 = SS$S2$drift, eq_drift_S3 = SS$S3$drift,
             sum = max(abs(arr[, "F", ] - (arr[, "history", ] + arr[, "inputs", ] + arr[, "climate", ]))))
    cat(sprintf("%-8s checks: balance drift S2 %.1e, S3 %.1e tC/ha/yr | sum %.1e\n", M,
                chk[["eq_drift_S2"]], chk[["eq_drift_S3"]], chk[["sum"]]))
    if (chk[["eq_drift_S2"]] > 1e-6 || chk[["eq_drift_S3"]] > 1e-6 || chk[["sum"]] > 1e-9)
      stop(M, ": recalibrated decomposition check failed")
    out$models[[M]] <- list(arr = arr, n_plots = length(bal), arms = rids, chk = chk,
                            draws = lapply(SS, function(z) if (is.null(z)) NULL else z$ks))
    out$rid[[M]] <- rids
    cat(sprintf("%-8s %d draws per arm x %d plots in %.1f min | 1985-2024: S1 %+.1f  S2 %+.1f  S3 %+.1f tC/ha\n",
                M, N_DRAW, length(bal), as.numeric(difftime(Sys.time(), t0, units = "mins")),
                mean(colSums(S1[yr <= 2024, ])), mean(colSums(S2[yr <= 2024, ])), mean(colSums(S3[yr <= 2024, ]))))
  }
  new <- setdiff(snap(), before)
  unlink(new[!dir.exists(new)]); unlink(new[dir.exists(new)], recursive = TRUE)
  saveRDS(out, "doublechecks/attribution_recalibrated.rds")
  cat("saved doublechecks/attribution_recalibrated.rds\n"); quit(save = "no")
}
for (M in MODELS) {
  t0 <- Sys.time()
  e  <- setup_model(M)
  bal  <- as.character(balanced_plots(e$obs_meta))
  post <- readRDS(sprintf("Calibration_real_data_transient/runs/%s%s_posterior_%s.rds", M, PTAG, RID[[M]]))
  if (EQP) post <- cbind(post, sigma_init = 1)
  # Pre-run climate = reference climate: the RUN's own (input bundle), not this shell's.
  bun <- readRDS(sprintf("Data/model_inputs/%s%s_inputs_%s.rds", M, PTAG, RID[[M]]))
  hiket_set_prerun_clim(bun$climate_ss_by_plot, as.character(bun$plots_real))
  if (!is.null(bun$climate_ss_by_plot)) e$litter_means <- bun$litter_means   # its precip_mean too
  cat(M, "pre-run / reference climate:", hiket_prerun_label(), "\n")
  ks   <- sort(sample(nrow(post), N_DRAW))
  A <- lapply(ks, function(k) {
    mp <- e$.assemble(setNames(as.numeric(post[k, ]), colnames(post)))
    R  <- parallel::mclapply(bal, function(p) tryCatch(attribute_plot(e, mp, p), error = function(z) NULL),
                             mc.cores = CORES)
    R  <- R[!vapply(R, is.null, logical(1))]
    ck <- do.call(rbind, lapply(R, attr, "chk"))
    list(m = Reduce(`+`, R) / length(R), chk = apply(ck, 2, max), n = length(R))
  })
  chk <- apply(do.call(rbind, lapply(A, `[[`, "chk")), 2, max)
  cat(sprintf("%-8s checks: pre-run vs initial state %.1e (rel) | equilibrium drift %.1e | sum %.1e | plots %d/%d\n",
              M, chk[["prerun"]], chk[["eq_drift"]], chk[["sum"]], min(vapply(A, `[[`, 0L, "n")), length(bal)))
  if (chk[["prerun"]] > 1e-8 || chk[["eq_drift"]] > 1e-6 || chk[["sum"]] > 1e-9) stop(M, ": decomposition check failed")
  out$models[[M]] <- list(arr = simplify2array(lapply(A, `[[`, "m")), n_plots = length(bal), draws = ks, chk = chk)
  cat(sprintf("%-8s %d draws x %d plots in %.1f min\n", M, N_DRAW, length(bal),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
new <- setdiff(snap(), before)                    # setup side effects
unlink(new[!dir.exists(new)]); unlink(new[dir.exists(new)], recursive = TRUE)
cat(sprintf("removed %d setup side-effect files\n", length(new)))
OUTF <- if (EQP) "doublechecks/attribution_eqparams.rds" else "doublechecks/attribution.rds"
saveRDS(out, OUTF)
cat("saved", OUTF, "\n")
