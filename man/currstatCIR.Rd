% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/currstatCIR.R
\name{currstatCIR}
\alias{currstatCIR}
\title{Estimate a survival function under current status sampling}
\usage{
currstatCIR(
  time,
  event,
  W,
  direction = "increasing",
  SL_control = list(SL.library = c("SL.mean"), V = 5, method = "method.NNLS"),
  HAL_control = list(n_bins = c(5, 10), grid_type = c("equal_range", "equal_mass"), V =
    5),
  deriv_method = "m-spline",
  missing_method = "extended",
  eval_region,
  n_eval_pts = 101
)
}
\description{
Estimate a survival function under current status sampling
}
