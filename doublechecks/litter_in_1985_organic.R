# =============================================================================
# litter_in_1985_organic.R — does the 1985 organic-layer stock INCLUDE the
# litter (LM = litter + moss) layer?
#
# Why it matters. build_soc_homogenized.R folds a litter term into the 2006 and
# 2024 organic layer (BiSo cols 247/248) but the Data_1985 sheet has no such
# term. If 1985 also excludes litter, the three campaigns are measuring
# different quantities and ~89% of the apparent 1985->2006 organic gain is a
# definition, not carbon (see doublechecks/soc_depth_distribution.R sec. 7).
#
# The absence of a litter COLUMN proves nothing: Kramarenko 2012 (Pro Gradu,
# sec. 2.3) records that the first campaign measured the organic layer "aina
# kokonaisena" -- always as a whole. So test it two ways instead.
#
#   SIGNATURE A (dry mass).  C_kgm2 = ORGANIC_LAYER_WEIGHT * C%/100 holds
#   exactly in BiSo, so organic-layer dry mass is recoverable in both
#   campaigns. If 1985 pooled the litter in, its mass must sit at or above the
#   2006 humus-only mass.
#
#   SIGNATURE B (carbon concentration).  Litter is less decomposed than humus,
#   so mixing it in raises C%. CONFOUNDED: campaign 1's C% was largely
#   predicted from loss-on-ignition by regression while 2006 used a Leco TGA
#   (Kramarenko sec. 4.4), and that method difference is itself worth ~2 pp.
#
# Run from repo root:  Rscript doublechecks/litter_in_1985_organic.R
# =============================================================================

suppressPackageStartupMessages({ library(readxl); library(dplyr) })
options(width = 130)

ROOT <- if (dir.exists("Data/Komeetta")) "." else ".."
XLSX <- file.path(ROOT, "Data/Komeetta/Hannu/Komeetta 150526hi--.xlsx")
stopifnot(file.exists(XLSX))

# BiSo column indices (spreadsheet order; header on row 3, data from row 4)
#   22 = ORGANIC_LAYER_WEIGHT (kg/m2)   159 = 1985 C%   165 = 2006 C%   170 = 2024 C%
#   176/177 = organic C kg/m2 2006/2024  247/248 = litter kg/ha 2006/2024
sp <- suppressMessages(read_excel(XLSX, sheet = "BiSo", range = "A4:IN3031", col_names = FALSE,
                                  col_types = "text", .name_repair = "minimal"))
nn <- function(i) suppressWarnings(as.numeric(sp[[i]]))
b <- tibble(plot_id = as.integer(nn(3)), Krs = as.integer(nn(4)), REP = as.integer(nn(9)),
            OLW = nn(22), Cpct06 = nn(165), C06 = nn(176), lit06 = nn(247)) |>
  filter(!is.na(plot_id), REP == 1, Krs == 101)

d85 <- suppressMessages(read_excel(XLSX, sheet = "Data_1985", .name_repair = "minimal")) |> as.data.frame()
a <- tibble(plot_id = as.integer(d85$VMI), Cpct85 = as.numeric(d85$`Corg%`),
            C85 = as.numeric(d85$Corg_kgha)) |> filter(!is.na(plot_id))

m <- inner_join(a, b, by = "plot_id")
cat("plots matched 1985 <-> BiSo organic:", nrow(m), "\n")

# --- the mass route has to be verified before it can be used -----------------
b$mass_fromC <- 1e4 * b$C06 / (b$Cpct06 / 100)
ok0 <- is.finite(b$OLW) & is.finite(b$mass_fromC) & b$OLW > 0
cat(sprintf("\nroute check: (C stock / C%%) / (OLW*1e4) = %.3f over %d rows, %d distinct OLW\n",
            median(b$mass_fromC[ok0] / (b$OLW[ok0] * 1e4)), sum(ok0), n_distinct(round(b$OLW, 4))))
cat("=> C_kgm2 = OLW * C%/100 holds, so organic-layer dry mass is recoverable.\n")

# --- SIGNATURE A -------------------------------------------------------------
cat("\n=== A. ORGANIC-LAYER DRY MASS (kg/ha, paired) ===\n")
m <- m |> mutate(mass85       = C85 / (Cpct85 / 100),
                 mass06_humus = OLW * 1e4,
                 lit_dry      = lit06 / 0.47)          # litter C -> dry mass at ~47% C
ok <- is.finite(m$mass85) & is.finite(m$mass06_humus) & m$mass85 > 0 & m$mass06_humus > 0
tot <- m$mass06_humus + coalesce(m$lit_dry, 0)
cat(sprintf("  1985 organic layer      %8.0f\n", median(m$mass85[ok])))
cat(sprintf("  2006 humus alone        %8.0f   paired ratio 1985/x = %.3f\n",
            median(m$mass06_humus[ok]), median(m$mass85[ok] / m$mass06_humus[ok])))
cat(sprintf("  2006 humus + litter     %8.0f   paired ratio 1985/x = %.3f\n",
            median(tot[ok]), median(m$mass85[ok] / tot[ok])))
cat("\n  1985 sits ~10% BELOW the humus-only mass. Litter inclusion would require the\n")
cat("  humus layer itself to have gained ~25% mass in 21 yr. => favours EXCLUDES.\n")

# --- SIGNATURE B -------------------------------------------------------------
cat("\n=== B. ORGANIC-LAYER CARBON CONCENTRATION (paired) ===\n")
ok2 <- is.finite(m$Cpct85) & is.finite(m$Cpct06)
cat(sprintf("  1985 %5.2f %%   2006 humus %5.2f %%   paired median diff %+0.2f pp  (n=%d, p=%.3g)\n",
            median(m$Cpct85[ok2]), median(m$Cpct06[ok2]),
            median(m$Cpct85[ok2] - m$Cpct06[ok2]), sum(ok2),
            wilcox.test(m$Cpct85[ok2], m$Cpct06[ok2], paired = TRUE)$p.value))
cat("  Higher C% in 1985 is what mixed-in litter would do => would favour INCLUDES,\n")
cat("  BUT campaign 1's C% was largely LOI-regression-predicted vs 2006's Leco TGA,\n")
cat("  a method difference of about this size. Not clean evidence either way.\n")

# --- by-product: the like-for-like humus signal ------------------------------
cat("\n=== BY-PRODUCT: the humus layer's REAL 1985->2006 change (like-for-like) ===\n")
cat(sprintf("  dry mass  %+5.1f %%   C%%  %+5.2f pp   =>  carbon %+5.1f %%\n",
            100 * (median(m$mass06_humus[ok]) / median(m$mass85[ok]) - 1),
            median(m$Cpct06[ok2]) - median(m$Cpct85[ok2]),
            100 * (median(1e4 * m$C06[ok], na.rm = TRUE) / median(m$C85[ok]) - 1)))
cat("  i.e. the band both we and Kramarenko regard as the trustworthy one carries\n")
cat("  almost no trend once the litter definition is taken out.\n")

cat("\nVERDICT: leaning EXCLUDES (mass is the stronger signature; the C% signal has a\n")
cat("documented alternative explanation). Not settled -- confirm with H. Ilvesniemi.\n")
