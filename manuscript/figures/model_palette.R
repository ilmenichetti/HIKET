# === Shared manuscript palettes (source this from every figure script) ===
# ONE place so colour coding is uniform across all figures.

# --- per-model palette: ggthemes "Temperature Diverging" in COMPLEXITY order ---
# Diverging by construction, which does double duty: (1) encodes the complexity
# ramp SP1->Yasso20; (2) the two hue-sides CLASS the two model families --
# green-ish simple trio (SP1/TP2/TP3) vs warm gold->rust Yasso trio (07/15/20).
# Saturated yellow centre => no washed-out middle. Use wherever colour = MODEL.
MODEL_ORDER <- c("SP1","TP2","TP3","Yasso07","Yasso15","Yasso20")
MODEL_COL   <- setNames(c("#529985","#76A26A","#B4BC53","#E5C749","#E5A84E","#C26B51"),
                        MODEL_ORDER)

# --- campaign palette: single-hue sequential in TIME order ---
# Colour = the SOC CAMPAIGN (an observation), so it stays in the same red family
# the manuscript already uses for observed values (firebrick in F2, #AA3333 in F5),
# ramped light (oldest) -> dark (newest) so the time order is readable without a
# legend. Use wherever colour = CAMPAIGN rather than model.
CAMPAIGN_ORDER <- c("VMI8", "Biosoil", "Komeetta")
CAMPAIGN_COL   <- setNames(c("#EFB59C", "#C85A3C", "#6E2414"), CAMPAIGN_ORDER)
CAMPAIGN_LAB   <- setNames(c("VMI8 1985", "Biosoil 2006", "Komeetta 2024"), CAMPAIGN_ORDER)

# --- basal-area class palette: ggthemes "Green-Gold", 5 quintiles ---
# Use wherever SCATTERPLOT points are coloured by stand basal area (F6, S2):
# gold = low basal area -> dark green = high.
BASAL_COL <- setNames(c("#F4D166","#A7BE5A","#61A656","#39884C","#146C36"),
                      paste0("q", 1:5))

# Quintile breaks + legend labels for a basal-area vector (shared by F6 & S2 so the
# classes match). Returns list(brk, classify, labels).
basal_classes <- function(ba) {
  brk <- unique(stats::quantile(ba, seq(0, 1, 0.2), na.rm = TRUE))
  nb  <- length(brk)
  labs <- c(sprintf("<%.0f", brk[2]),
            vapply(2:(nb - 2), function(i) sprintf("%.0f-%.0f", brk[i], brk[i + 1]), character(1)),
            sprintf(">%.0f", brk[nb - 1]))
  list(brk = brk, nb = nb,
       classify = function(x) cut(x, breaks = brk, labels = names(BASAL_COL)[seq_len(nb - 1)],
                                  include.lowest = TRUE),
       labels = labs)
}
