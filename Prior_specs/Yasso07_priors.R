# Yasso07 prior specification — HIKET pipeline
# To regenerate: run Priors_model_matching.R interactively (Yasso07 is not
# auto-written; sigma_ppm values are from published posterior limits, not a
# posterior file).
# Prior centres (transfer fractions): YASSO07_DEFAULT_PARAMS (Tuomi et al. 2009)
# Prior centres (climate/size): GUI MAP, y07par_gui.csv (Tuomi et al. 2011 EMS)
# Prior widths (climate): Tuomi et al. 2009 (Ecol. Modelling) Table 3
# Prior widths (woody size): Tuomi et al. 2011 (Ecol. Modelling) Table 4
#   Published "±" limits read as 1σ (conservative — not divided by 1.96); for
#   log-transform params the relative SD = (1σ limit)/(paper MAP) via delta method.
# Transfer fraction widths: logit SD 0.4 (Tier-2 common weak prior — Finnish
#   data drives structure; see Prior_specs/PRIOR_HOMOGENIZATION_PLAN.md §4.2)

# Complete physical-space prior centres for all free parameters.
# r stored as positive (abs() applied; Tuomi reference stores it negative).
YASSO07_FREE_DEFAULTS <- c(
  # Transfer fractions (12)
  p_WA = 0.4888527989,
  p_EA = 0.01905768365,
  p_NA = 0.9696374536,
  p_AW = 0.9872559905,
  p_EW = 0.002843263559,
  p_NW = 0.003396461252,
  p_AE = 1.399370376e-05,
  p_WE = 1.796692413e-05,
  p_NE = 0.01218125969,
  p_AN = 0.002777846763,
  p_WN = 0.01269555371,
  p_EN = 0.9713827968,
  # Climate response
  beta1  = 0.09873183817,
  beta2  = -0.001571640489,
  gamma  = -1.271691799,
  # Woody size modifier
  delta1 = -1.708411336,
  delta2 = 0.8585553765,
  r      = 0.3068014085,
  # Auxiliary uncertainty parameters
  sigma_init  = 1.00,
  sigma_input = 1.00
)

# sigma_ppm in unconstrained (transformed) space. All free params listed
# explicitly (decision #5: explicit per-fraction listing for traceability).
# Climate/woody widths from Tuomi 2009 T3 / 2011 T4, "±" read as 1σ (§4.1).
YASSO07_SIGMA_PPM <- c(
  # Transfer fractions (12) — Tier-2 common logit SD 0.4
  p_WA = 0.4, p_EA = 0.4, p_NA = 0.4, p_AW = 0.4,
  p_EW = 0.4, p_NW = 0.4, p_AE = 0.4, p_WE = 0.4,
  p_NE = 0.4, p_AN = 0.4, p_WN = 0.4, p_EN = 0.4,
  # Climate response (Tuomi 2009 Table 3)
  beta1       = 0.26,     # log; rel SD 0.020/0.076 (paper MAP)
  beta2       = 0.00065,  # unconstrained; T3 −8.9 ±6.5 ×10⁻⁴
  gamma       = 0.20,     # unconstrained; T3 −1.27 ±0.20
  # Woody size modifier (Tuomi 2011 Table 4)
  delta1      = 0.16,     # unconstrained; T4 −1.71 ±0.16
  delta2      = 0.12,     # log; rel SD 0.10/0.86
  r           = 0.042,    # log; rel SD 0.013/0.306
  # Auxiliary uncertainty
  sigma_init  = 0.50,
  sigma_input = 0.50
)
