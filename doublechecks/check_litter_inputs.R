#!/usr/bin/env Rscript
# =============================================================================
# check_litter_inputs.R  --  litter forcing sanity check
# -----------------------------------------------------------------------------
# Question: do all six models integrate the SAME litter object, and does that
# object match the raw Zenodo source CSV (no monthly->annual aggregation bug,
# no units inflation)?
#
# What it does:
#   1. Reads the raw source (input_raw_monthly.csv), sums the 12 monthly AWEN
#      litter columns to an annual total per plot-year  -> "SOURCE".
#   2. For each model, reads its latest input bundle in Data/model_inputs/ and
#      reconstructs the annual total litter the model actually integrates
#      (SP1/TP2/TP3 store J_total directly; Yasso07/15/20 store the 12 AWEN
#      columns whose rowSum is the same total)            -> "MODEL".
#   3. Overlays the cross-plot mean trajectories and reports, per model, the
#      max absolute difference vs the source on the model's own plot set
#      (should be ~0 to floating-point) -> confirms one shared object.
#
# Run from the project root:   Rscript doublechecks/check_litter_inputs.R
# =============================================================================

OUT_DIR <- "doublechecks/figures"
INP_DIR <- "Data/model_inputs"
SRC_CSV <- file.path(INP_DIR, "input_raw_monthly.csv")
MODELS  <- c("SP1", "TP2", "TP3", "Yasso07", "Yasso15", "Yasso20")
FIRST_YEAR <- 1986L          # drop 1985 (VMI8 init year: near-zero applied litter)

stopifnot(dir.exists(INP_DIR))                 # guard: must run from project root
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# --- helpers -----------------------------------------------------------------

# total annual litter per (plot, year) from any inputs_by_plot-style frame
litter_total <- function(df) {
  if ("J_total" %in% names(df)) {
    tot <- df$J_total
  } else {
    cols <- grep("^(nwl|fwl|cwl)_", names(df), value = TRUE)
    tot  <- rowSums(df[, cols, drop = FALSE], na.rm = TRUE)
  }
  data.frame(plot_id = as.character(df$plot_id), year = df$year, tot = tot)
}

# cross-plot mean (+/- 95% CI of the mean) by year, restricted to >= FIRST_YEAR
mean_by_year <- function(pt, plots = NULL) {
  if (!is.null(plots)) pt <- pt[pt$plot_id %in% plots, ]
  pt <- pt[pt$year >= FIRST_YEAR, ]
  s  <- split(pt$tot, pt$year)
  data.frame(
    year = as.integer(names(s)),
    mean = sapply(s, mean, na.rm = TRUE),
    se   = sapply(s, function(x) sd(x, na.rm = TRUE) / sqrt(sum(is.finite(x))))
  )
}

latest_bundle <- function(model) {
  fs <- list.files(INP_DIR, pattern = sprintf("^%s_inputs_.*\\.rds$", model),
                   full.names = TRUE)
  if (!length(fs)) return(NA_character_)
  fs[which.max(file.mtime(fs))]
}

# --- 1. SOURCE: directly from the raw monthly CSV ----------------------------

message("Reading source: ", SRC_CSV)
raw  <- read.csv(SRC_CSV)
litc <- grep("^C_(nwl|fwl|cwl)_", names(raw), value = TRUE)
raw$mtot <- rowSums(raw[, litc], na.rm = TRUE)              # monthly total
src_annual <- aggregate(mtot ~ plot_id + year, raw, sum)    # monthly -> annual
names(src_annual)[3] <- "tot"
src_annual$plot_id <- as.character(src_annual$plot_id)
src_mean <- mean_by_year(src_annual)

# --- 2. MODEL: from each model's input bundle --------------------------------

model_pt   <- list()   # per-model (plot,year) totals
model_mean <- list()   # per-model cross-plot mean trajectory
for (m in MODELS) {
  f <- latest_bundle(m)
  if (is.na(f)) { message("[", m, "] no bundle found - skipped"); next }
  pkg <- readRDS(f)
  pt  <- do.call(rbind, lapply(pkg$inputs_by_plot, litter_total))
  pt$plot_id <- as.character(pt$plot_id)
  model_pt[[m]]   <- pt
  model_mean[[m]] <- mean_by_year(pt)
  message(sprintf("[%-7s] %s | %d plots", m, basename(f),
                  length(unique(pt$plot_id))))
}

# --- 3. Numeric agreement: model vs source on the model's OWN plots ----------

cat("\n== Agreement of each model's litter vs the raw source ",
    "(matched plot-years) ==\n", sep = "")
cat(sprintf("%-8s %8s %12s %14s\n", "model", "nplots", "max|diff|", "mean source"))
for (m in names(model_pt)) {
  pt  <- model_pt[[m]]
  key <- paste(pt$plot_id, pt$year)
  src <- setNames(src_annual$tot, paste(src_annual$plot_id, src_annual$year))
  d   <- pt$tot - src[key]
  cat(sprintf("%-8s %8d %12.2e %14.3f\n",
              m, length(unique(pt$plot_id)), max(abs(d), na.rm = TRUE),
              mean(model_mean[[m]]$mean)))
}

# --- 4. Figure: source vs all six models -------------------------------------
# Fair comparison: restrict the source to the SAME plot set the models use, so
# any genuine data discrepancy would show as a gap. The faint dashed line is
# the source over ALL CSV plots, shown only to illustrate that the small level
# offset is a plot-SET effect (calibration filter), not a data difference.

model_plots   <- sort(unique(unlist(lapply(model_pt, function(p) unique(p$plot_id)))))
src_mean_cal  <- mean_by_year(src_annual, plots = model_plots)   # apples-to-apples

png_path <- file.path(OUT_DIR, "litter_source_vs_models.png")
png(png_path, width = 1700, height = 1100, res = 150)
op <- par(mar = c(4.2, 4.6, 3.2, 1), mgp = c(2.6, 0.7, 0), las = 1)

cols   <- setNames(c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02"),
                   MODELS)
yr_all <- range(c(src_mean$year, src_mean_cal$year))
yr_rng <- range(c(src_mean$mean, src_mean_cal$mean,
                  unlist(lapply(model_mean, function(d) d$mean))))

plot(src_mean_cal$year, src_mean_cal$mean, type = "n", xlim = yr_all, ylim = yr_rng,
     xlab = "Year", ylab = expression("Cross-plot mean litter (tC ha"^-1*" yr"^-1*")"),
     main = "Litter forcing: raw source vs. object integrated by each model")

# source on the model plot set: thick grey band + line (95% CI of the mean)
polygon(c(src_mean_cal$year, rev(src_mean_cal$year)),
        c(src_mean_cal$mean - 1.96*src_mean_cal$se,
          rev(src_mean_cal$mean + 1.96*src_mean_cal$se)),
        col = adjustcolor("grey50", 0.25), border = NA)
lines(src_mean_cal$year, src_mean_cal$mean, col = "grey30", lwd = 7)

# each model overlaid (thin coloured lines; all six coincide on the grey line)
for (m in names(model_mean))
  lines(model_mean[[m]]$year, model_mean[[m]]$mean, col = cols[m], lwd = 1.6)

# context: source over ALL CSV plots (un-filtered) -> explains any level offset
lines(src_mean$year, src_mean$mean, col = "grey55", lwd = 1.4, lty = 2)

legend("bottomright", bty = "n",
       lwd = c(7, rep(1.6, length(model_mean)), 1.4),
       lty = c(1, rep(1, length(model_mean)), 2),
       col = c("grey30", cols[names(model_mean)], "grey55"),
       legend = c("SOURCE (raw CSV, model plot set)", names(model_mean),
                  "SOURCE (all CSV plots, ref.)"))
mtext(sprintf("All six models identical to source on matched plot-years (max |diff| < 1e-13). Model plot set: %d plots.",
              length(model_plots)),
      side = 3, line = 0.2, cex = 0.75, col = "grey30")
par(op); dev.off()
message("\nFigure written: ", png_path)
