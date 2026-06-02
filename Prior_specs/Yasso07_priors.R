# Yasso07 prior specification — HIKET pipeline
# To regenerate: run Priors_model_matching.R interactively (Yasso07 is not
# auto-written; sigma_ppm values are from literature CIs, not a posterior file).
# Prior centres (transfer fractions): YASSO07_DEFAULT_PARAMS (Tuomi et al. 2009)
# Prior centres (climate/size): published MAP, y07par_gui.csv
# Prior widths (climate/size): literature CIs, Tuomi 2008-2010
# Transfer fraction widths: 1.0 (weakly informative — Finnish data drives structure)

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

# sigma_ppm overrides: only non-1.0 values listed.
# Transfer fractions default to 1.0 (weakly informative).
YASSO07_SIGMA_PPM <- c(
  beta1       = 0.20,
  beta2       = 0.05,
  gamma       = 0.30,
  delta1      = 0.15,
  delta2      = 0.10,
  r           = 0.015,
  sigma_init  = 0.50,
  sigma_input = 0.50
)
