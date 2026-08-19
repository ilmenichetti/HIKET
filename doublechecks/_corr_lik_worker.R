# Worker for test_correlated_likelihood.R. Sources the real calibration script up
# to the MCMC launch under whatever HIKET_* environment it inherits, then prints
# the log-likelihood at a deterministic set of parameter vectors.
#
# Run ONLY as a subprocess: the engine's source() of correlated_likelihood.R is
# nested inside source(local = e), and a nested source() defaults to globalenv(),
# so two configurations in one session overwrite each other's switches.
M   <- commandArgs(trailingOnly = TRUE)[[1]]
src <- readLines(sprintf("Calibration_real_data_transient/run_%s_transient_calibration.R", M),
                 warn = FALSE)
cut <- grep("^t_run <- system.time", src)[1]
e   <- new.env(parent = globalenv())

msg <- character()
ok <- tryCatch({
  msg <<- capture.output(invisible(capture.output(
    source(textConnection(paste(src[seq_len(cut - 1)], collapse = "\n")), local = e))),
    type = "message")
  TRUE
}, error = function(cnd) {
  cat("ERROR ", conditionMessage(cnd), "\n", sep = "")
  FALSE
})
# the setup lines the parent displays; on failure everything, so the guard text shows
for (l in if (ok) grep("ERROR MODEL|tau_R", msg, value = TRUE) else msg)
  cat("SETUP ", l, "\n", sep = "")
if (!ok) quit(status = 1L)

ll <- get("ll_fn", e); bx <- get("best_x", e)
set.seed(2025)                                   # identical draws in every worker
xs <- c(list(bx), lapply(1:4, function(i) bx + rnorm(length(bx), 0, 0.05)))
for (k in seq_along(xs)) cat(sprintf("RESULT %d %.12f\n", k, ll(xs[[k]])))
cat(sprintf("TIMING %.4f\n", system.time(for (i in 1:5) ll(bx))[["elapsed"]] / 5))
