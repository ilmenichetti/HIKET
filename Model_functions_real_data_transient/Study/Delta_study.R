# Size modifier: hS(d) = MIN(1, (1 + delta1*d + delta2*d^2)^(-|r|))
# Published defaults: delta1=-1.708, r=-0.307
# Diameter range: 0-30 cm (Finnish stumps and coarse woody debris)
d <- seq(0, 30, by = 0.1)

delta1_def <- -1.7084
r_def      <- -0.3068
delta2_vals <- c(-1.17, -0.5, 0.50, 0.86, 1.27)
cols        <- c("firebrick", "tomato", "steelblue", "darkblue", "purple")

size_mod_raw <- function(d, delta1, delta2, r) {
  base <- 1 + delta1 * d + delta2 * d^2
  pmin(1, base^(-abs(r)))   # no NA masking — let it produce NaN naturally
}


plot(NA, xlim = c(0, 30), ylim = c(-0.5, 1.2),
     xlab = "Woody litter diameter (cm)",
     ylab = "Size modifier hS(d)",
     main = sprintf("Yasso size modifier  |  delta1 = %.3f,  |r| = %.3f",
                    delta1_def, abs(r_def)))
abline(h = c(0, 1), lty = c(1, 3), col = c("grey30", "grey70"), lwd = c(1, 1))
abline(v = c(2, 5, 10, 20), lty = 3, col = "grey85")
mtext(c("twigs", "branches", "stems", "logs"),
      at = c(2, 5, 10, 20), side = 3, cex = 0.7, col = "grey50", line = 0.2)

# Shade invalid zone (size_mod < 0)
polygon(c(0, 30, 30, 0), c(0, 0, -0.5, -0.5),
        col = adjustcolor("firebrick", 0.08), border = NA)
text(15, -0.25, "invalid zone (NaN in Fortran)", col = "firebrick", cex = 0.8)

for (i in seq_along(delta2_vals)) {
  d2 <- delta2_vals[i]
  base <- 1 + delta1_def * d + d2 * d^2
  sm   <- pmin(1, ifelse(base > 0, base^(-abs(r_def)), base))  # show raw value when negative
  lines(d, sm, col = cols[i], lwd = 2)
  
  # Mark zero crossing
  if (any(base <= 0)) {
    d_zero <- d[which(base <= 0)[1]]
    abline(v = d_zero, col = cols[i], lty = 2, lwd = 1)
  }
}

legend("topright",
       legend = sprintf("delta2 = %+.2f%s", delta2_vals,
                        ifelse(delta2_vals == 0.86, " (Yasso07 default)",
                               ifelse(delta2_vals == 1.27, " (Yasso15 default)",
                                      ifelse(delta2_vals == -1.17, " (Yasso15 repository)", "")))),
       col = cols, lwd = 2, bty = "n", cex = 0.85)
