#' Wrapper for various f_w stacked algorithms
#'
#' @param time Observed time
#' @param event Indicator of event (vs censoring)
#' @param X Covariate matrix
#' @param bin_size Size of time bin on which to discretize for estimation
#' @param censored Logical, indicates whether to run regression on censored observations (vs uncensored)
#' @param algorithm Which binary classification algorithm to use
#' @param V CV fold number, required for tuned algorithms (xgboost, ranger)
#' @param entry Variable indicating time of entry into the study (truncation variable) if applicable
#' @param time_basis How to treat time (continuous or dummy)
#' @param tuning_params Tuning parameters for binary classification
#' @param SL.library SuperLearner library
#'
#' @return An fitted pooled binary regression for the truncation distribution
#'
#' @export
#'
#' @examples
f_w_stack <- function(time,
                      event,
                      X,
                      entry,
                      censored,
                      bin_size = NULL,
                      algorithm,
                      V = NULL,
                      SL.library = NULL,
                      time_basis,
                      tuning_params = NULL,
                      direction = "forward"){

  if (algorithm == "xgboost"){ # do xgboost if speed not a concern
    fit <- f_w_stack_xgboost(time = time,
                             event = event,
                             X = X,
                             censored = censored,
                             bin_size = bin_size,
                             V = V,
                             time_basis = time_basis,
                             tuning_params = tuning_params,
                             entry = entry,
                             direction = direction)
  } else if (algorithm == "SuperLearner"){
    fit <- f_w_stack_SuperLearner(time = time,
                                  event = event,
                                  X = X,
                                  censored = censored,
                                  bin_size = bin_size,
                                  time_basis = time_basis,
                                  SL.library = SL.library,
                                  V = V,
                                  entry = entry,
                                  direction = direction)
  }
  return(fit)
}
