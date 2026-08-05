# =============================================================================
# input_data_qc.R   (2026-08-04)
#
# Systematic sweep of the model input bundle for values that are wrong but not
# obviously wrong.
#
# Motivation: the 1985 litter artefact ([[litter-1985-first-year-artefact]]) sat
# in the inputs since April 2026 and survived every existing check, because none
# of them was looking for it. Row counts were complete, values were positive and
# finite, and one bad year in forty barely moves an annual mean. It only surfaced
# because an unrelated diagnostic happened to weight the HEAD of the series.
#
# So this script deliberately looks where summary statistics do not: at the ends
# of series, at year-on-year steps, at per-plot constancy, at component structure,
# and at anything that is suspiciously round, flat, duplicated or extreme.
#
# It REPORTS; it does not modify anything. Findings are printed as [OK] / [NOTE]
# / [CHECK]; [CHECK] means look at it, not necessarily that it is broken.
#
# Usage:  Rscript doublechecks/input_data_qc.R
# =============================================================================

f_in <- "Data/model_inputs/input_raw_monthly.csv"
if (!file.exists(f_in)) stop("Run from the project root.")
inp  <- read.csv(f_in)
site <- read.csv("Data/model_inputs/site_raw.csv")

lit_cols <- grep("^C_(nwl|fwl|cwl)_[AWEN]$", names(inp), value = TRUE)
inp$J    <- rowSums(inp[, lit_cols, drop = FALSE], na.rm = TRUE)
ann      <- aggregate(J ~ plot_id + year, data = inp, FUN = sum)
calib    <- site$plot_id[site$calib_ready]
ann_c    <- ann[ann$plot_id %in% calib, ]

say <- function(tag, msg) cat(sprintf("  [%-5s] %s\n", tag, msg))
hdr <- function(t) cat(sprintf("\n--- %s %s\n", t, strrep("-", max(0, 66 - nchar(t)))))

cat(strrep("=", 74), "\n")
cat("INPUT DATA QC  --  looking for wrong-but-not-obviously-wrong values\n")
cat(sprintf("bundle: %s   (%d rows, %d plots, %d-%d)\n", f_in, nrow(inp),
            length(unique(inp$plot_id)), min(inp$year), max(inp$year)))
cat(strrep("=", 74), "\n")

# ---------------------------------------------------------------- 1. series ends
hdr("1. ANOMALOUS YEARS (the 1985 failure mode: a year out of line with neighbours)")
med <- tapply(ann_c$J, ann_c$year, median)
yrs <- as.integer(names(med))
ratio <- med[-1] / med[-length(med)]
bad <- which(ratio < 0.6 | ratio > 1.7)
if (!length(bad)) say("OK", "no year-on-year median step outside 0.6x-1.7x") else
  for (i in bad) say("CHECK", sprintf("%d -> %d: median litter %.3f -> %.3f (x%.2f)",
                                      yrs[i], yrs[i+1], med[i], med[i+1], ratio[i]))
say("NOTE", sprintf("first year %d = %.3f | second %d = %.3f (ratio %.2f)",
                    yrs[1], med[1], yrs[2], med[2], med[2]/med[1]))
say("NOTE", sprintf("last  year %d = %.3f | previous %d = %.3f (ratio %.2f)",
                    yrs[length(yrs)], med[length(med)], yrs[length(yrs)-1],
                    med[length(med)-1], med[length(med)]/med[length(med)-1]))

# ------------------------------------------------------------ 2. duplicated years
hdr("2. DUPLICATED YEAR BLOCKS (identical litter in consecutive years)")
w <- reshape(ann_c, idvar = "plot_id", timevar = "year", direction = "wide")
dupyr <- c()
for (i in 2:length(yrs)) {
  a <- w[[paste0("J.", yrs[i-1])]]; b <- w[[paste0("J.", yrs[i])]]
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) && mean(abs(a[ok] - b[ok]) < 1e-9) > 0.95)
    dupyr <- c(dupyr, sprintf("%d == %d (%.0f%% of plots identical)",
                              yrs[i-1], yrs[i], 100*mean(abs(a[ok]-b[ok]) < 1e-9)))
}
if (!length(dupyr)) say("OK", "no consecutive years are copies") else
  for (d in dupyr) say("NOTE", paste(d, "-- expected for 2023->2024 (documented carry-forward)"))

# ------------------------------------------------------------- 3. zero / negative
hdr("3. ZERO, NEGATIVE AND NON-FINITE LITTER")
say(if (any(inp[, lit_cols] < 0, na.rm = TRUE)) "CHECK" else "OK",
    sprintf("negative component values: %d", sum(inp[, lit_cols] < 0, na.rm = TRUE)))
say(if (anyNA(inp[, lit_cols])) "CHECK" else "OK",
    sprintf("NA component values: %d", sum(is.na(inp[, lit_cols]))))
zero_py <- ann_c[ann_c$J <= 0, ]
say(if (nrow(zero_py)) "CHECK" else "OK",
    sprintf("plot-years with zero total litter: %d", nrow(zero_py)))
near0 <- ann_c[ann_c$J > 0 & ann_c$J < 0.2, ]
say(if (nrow(near0)) "CHECK" else "OK",
    sprintf("plot-years with 0 < litter < 0.2 tC/ha/yr: %d%s", nrow(near0),
            if (nrow(near0)) sprintf(" (years: %s)",
              paste(sort(unique(near0$year)), collapse = ",")) else ""))

# --------------------------------------------------------------- 4. flat series
hdr("4. IMPLAUSIBLY FLAT OR SPIKY PER-PLOT SERIES")
cv <- tapply(ann_c$J, ann_c$plot_id, function(v) sd(v)/mean(v))
say(if (sum(cv < 0.01, na.rm=TRUE)) "CHECK" else "OK",
    sprintf("plots with essentially constant litter (CV < 1%%): %d",
            sum(cv < 0.01, na.rm = TRUE)))
mx <- tapply(ann_c$J, ann_c$plot_id, function(v) {
  r <- v[-1]/v[-length(v)]; r <- r[is.finite(r)]; if (length(r)) max(r) else NA })
say(if (sum(mx > 3, na.rm=TRUE)) "CHECK" else "OK",
    sprintf("plots with a >3x single-year jump: %d", sum(mx > 3, na.rm = TRUE)))

# ------------------------------------------------------------------ 5. magnitude
hdr("5. MAGNITUDE vs BOREAL EXPECTATION (1.5-4.5 tC/ha/yr typical)")
pm <- tapply(ann_c$J, ann_c$plot_id, mean)
say("NOTE", sprintf("per-plot mean litter: median %.2f | range %.2f - %.2f",
                    median(pm), min(pm), max(pm)))
say(if (sum(pm > 9) ) "CHECK" else "OK",
    sprintf("plots above the ~9 tC/ha/yr boreal NPP ceiling: %d", sum(pm > 9)))
say(if (sum(pm < 0.5)) "CHECK" else "OK",
    sprintf("plots below 0.5 tC/ha/yr: %d", sum(pm < 0.5)))

# ------------------------------------------------------------ 6. AWEN / size mix
hdr("6. COMPONENT STRUCTURE (AWEN x size)")
tot <- colSums(inp[, lit_cols], na.rm = TRUE)
zero_comp <- names(tot)[tot == 0]
say(if (length(zero_comp)) "CHECK" else "OK",
    sprintf("components that are zero everywhere: %s",
            if (length(zero_comp)) paste(zero_comp, collapse = ", ") else "none"))
for (sz in c("nwl","fwl","cwl")) {
  cc <- grep(sprintf("^C_%s_", sz), lit_cols, value = TRUE)
  say("NOTE", sprintf("%s share of total litter: %.1f%%", sz,
                      100 * sum(inp[, cc], na.rm = TRUE) / sum(inp[, lit_cols], na.rm = TRUE)))
}

# ---------------------------------------------------------------- 7. seasonality
hdr("7. MONTHLY DISTRIBUTION")
mm <- tapply(inp$J, inp$month, median)
say(if (diff(range(mm))/mean(mm) < 0.01) "NOTE" else "CHECK",
    sprintf("litter is %s across months (median range %.4f-%.4f)",
            if (diff(range(mm))/mean(mm) < 0.01) "UNIFORM (annual/12, by construction)" else "seasonal",
            min(mm), max(mm)))

# ------------------------------------------------------------------- 8. climate
hdr("8. CLIMATE")
for (v in c("temp_air","precip","evap")) {
  if (!v %in% names(inp)) next
  x <- inp[[v]]
  say(if (anyNA(x)) "CHECK" else "OK", sprintf("%-9s NAs: %d", v, sum(is.na(x))))
  say("NOTE", sprintf("%-9s range %.2f to %.2f", v, min(x, na.rm=TRUE), max(x, na.rm=TRUE)))
}
if ("precip" %in% names(inp))
  say(if (any(inp$precip < 0, na.rm=TRUE)) "CHECK" else "OK",
      sprintf("negative precipitation months: %d", sum(inp$precip < 0, na.rm = TRUE)))

# ------------------------------------------------------------------ 9. coverage
hdr("9. COVERAGE / STRUCTURE")
rpp <- table(table(inp$plot_id))
say(if (length(rpp) == 1) "OK" else "CHECK",
    sprintf("rows per plot: %s", paste(sprintf("%s rows x%d plots", names(rpp), rpp), collapse="; ")))
d <- duplicated(inp[, c("plot_id","year","month")])
say(if (any(d)) "CHECK" else "OK", sprintf("duplicate plot-year-month rows: %d", sum(d)))
say(if (all(calib %in% inp$plot_id)) "OK" else "CHECK",
    sprintf("calib_ready plots present in the bundle: %d / %d",
            sum(calib %in% inp$plot_id), length(calib)))

# ------------------------------------------------------------------ 10. SOC obs
hdr("10. SOC OBSERVATIONS")
soc <- inp[!is.na(inp$soc_obs_tCha), ]
say(if (all(soc$month == 6)) "OK" else "CHECK",
    sprintf("all SOC obs in June: %s", all(soc$month == 6)))
say("NOTE", sprintf("obs by year: %s",
                    paste(sprintf("%s=%d", names(table(soc$year)), table(soc$year)), collapse=" ")))
sc <- soc[soc$plot_id %in% calib, ]
say("NOTE", sprintf("calib-ready soc_obs: median %.1f | range %.1f - %.1f tC/ha",
                    median(sc$soc_obs_tCha), min(sc$soc_obs_tCha), max(sc$soc_obs_tCha)))
say(if (any(sc$soc_obs_tCha <= 0)) "CHECK" else "OK",
    sprintf("non-positive SOC observations: %d", sum(sc$soc_obs_tCha <= 0)))

cat("\n", strrep("=", 74), "\n", sep="")
cat("[CHECK] = inspect; not necessarily an error. [NOTE] = context, no action implied.\n")
cat(strrep("=", 74), "\n")
