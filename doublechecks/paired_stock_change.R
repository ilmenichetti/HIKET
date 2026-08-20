# Stock-change rates on the CALIBRATION BASIS (decision 2026-08-17):
#   same plots the models were fitted to, balanced set (obs in all 3 campaigns),
#   unweighted, whole profile, TRUE observation years.
# The balanced set matters: on the pairwise 2006->2024 set (n=409) the observed
# rate is +0.209; on the balanced set (n=310) it is +0.117. Same data, 1.8x.

# RUN_IDs auto-selected (were hard-coded to 20260813/14 until 2026-08-20, which would
# have reported the PREVIOUS run's stock change for the current one). Pin an older run
# for comparison with HIKET_PSC_RID="SP1=<id>,...".
runs <- vapply(c("SP1","TP2","TP3","Yasso07","Yasso15","Yasso20"), function(m) {
  fs <- list.files("Calibration_real_data_transient/runs",
                   pattern = sprintf("^%s_posterior_predictive_[0-9]{8}_[0-9]{6}\\.rds$", m))
  if (!length(fs)) stop("paired_stock_change.R: no posterior_predictive for ", m, call. = FALSE)
  sub(sprintf("^%s_posterior_predictive_(.+)\\.rds$", m), "\\1", sort(fs, decreasing = TRUE)[1])
}, character(1))
local({
  ov <- Sys.getenv("HIKET_PSC_RID", "")
  if (nzchar(ov)) { kv <- strsplit(strsplit(ov, ",")[[1]], "=")
    for (e in kv) runs[[trimws(e[1])]] <<- trimws(e[2]) }
})
cat("RUN_IDs:", paste(names(runs), runs, sep="=", collapse=" | "), "\n")

om  <- readRDS(sprintf("Data/model_inputs/Yasso15_inputs_%s.rds", runs[["Yasso15"]]))$obs_meta
BAL <- as.integer(names(om)[vapply(om, function(z) length(z$soc_obs) >= 3L, logical(1))])
cat("balanced plot set: n =", length(BAL), "\n")

camp <- function(y) ifelse(y <= 2000, 1L, ifelse(y <= 2015, 2L, 3L))

out <- list()
for (m in names(runs)) {
  f <- sprintf("Calibration_real_data_transient/runs/%s_posterior_predictive_%s.rds", m, runs[[m]])
  d <- as.data.frame(readRDS(f)$residuals_df)[, c("plot_id", "year", "soc_obs_tCha", "soc_mean")]
  d <- d[d$plot_id %in% BAL, ]
  d$c <- camp(d$year)
  d <- d[!duplicated(d[, c("plot_id", "c")]), ]

  for (pr in list(c(1L, 3L), c(2L, 3L), c(1L, 2L))) {
    a <- d[d$c == pr[1], ]; b <- d[d$c == pr[2], ]
    k <- intersect(a$plot_id, b$plot_id)
    a <- a[match(k, a$plot_id), ]; b <- b[match(k, b$plot_id), ]
    dt <- b$year - a$year; ok <- dt > 0
    out[[length(out) + 1]] <- data.frame(
      model = m, pair = paste0("c", pr[1], "->c", pr[2]), n = sum(ok),
      dt = round(mean(dt[ok]), 1),
      obs = round(mean((b$soc_obs_tCha - a$soc_obs_tCha)[ok] / dt[ok]), 3),
      mod = round(mean((b$soc_mean     - a$soc_mean    )[ok] / dt[ok]), 3))
  }
}
res <- do.call(rbind, out)
for (p in unique(res$pair)) {
  r <- res[res$pair == p, ]
  cat("\n==== ", p, "  (n=", r$n[1], ", dt=", r$dt[1], " yr, obs ", sprintf("%+.3f", r$obs[1]),
      ") ====\n", sep = "")
  r$verdict <- ifelse(r$mod > 0, "sink", "SOURCE")
  r$ratio   <- round(r$mod / r$obs, 2)
  print(r[, c("model", "mod", "ratio", "verdict")], row.names = FALSE)
}
cat("\nrates tC/ha/yr; ratio = model/observed.\n")
