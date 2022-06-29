#' Wrapper for various f_y stacked algorithms
#'
#' @param time Observed time
#' @param event Indicator of event (vs censoring)
#' @param X Covariate matrix
#' @param censored Logical, indicates whether to run regression on censored observations (vs uncensored)
#' @param bin_size Size of time bin on which to discretize for estimation
#' @param algorithm Which binary classification algorithm to use
#' @param V CV fold number, required for tuned algorithms (xgboost, ranger)
#' @param time_basis How to treat time (continuous or dummy)
#' @param tuning_params Tuning parameters for binary classification
#'
#' @return An fitted pooled binary regression for the CDF
#'
#' @export
#'
#' @examples
f_y_stack <- function(time,
                      event,
                      X,
                      censored,
                      bin_size = NULL,
                      algorithm,
                      V = NULL,
                      SL.library = NULL,
                      time_basis,
                      tuning_params = NULL){

  if (algorithm == "xgboost"){ # do xgboost if speed not a concern
    fit <- f_y_stack_xgboost(time = time,
                             event = event,
                             X = X,
                             censored = censored,
                             bin_size = bin_size,
                             V = V,
                             time_basis = time_basis,
                             tuning_params = tuning_params)
  } else if (algorithm == "SuperLearner"){
    fit <- f_y_stack_SuperLearner(time = time,
                                  event = event,
                                  X = X,
                                  censored = censored,
                                  bin_size = bin_size,
                                  time_basis = time_basis,
                                  SL.library = SL.library,
                                  V = V)
  }
  return(fit)
}
