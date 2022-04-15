#' Wrapper for various p_delta algorithms
#'
#' @param event Indicator of event (vs censoring)
#' @param X Covariate matrix
#' @param algorithm Which binary classification algorithm to use
#' @param V CV fold number, required for tuned algorithms (xgboost, ranger)
#' @param tuning_params Tuning parameters for binary classification
#' @param SL.library SuperLearner library
#'
#' @return An fitted binary regression for (complement of) probability of censoring
#'
#' @export
#'
#' @examples
p_delta <- function(event,
                    X,
                    algorithm,
                    V = NULL,
                    SL.library = NULL,
                    tuning_params = NULL){

  if (algorithm == "xgboost"){ # do xgboost if speed not a concern
    fit <- p_delta_xgboost(event = event,
                           X = X,
                           V = V,
                           tuning_params = tuning_params)
  } else if (algorithm == "SuperLearner"){
    fit <- p_delta_SuperLearner(event = event,
                                X = X,
                                SL.library = SL.library,
                                V = V)
  }
  return(fit)
}
