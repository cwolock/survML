#' Wrapper for various p_delta algorithms
#'
#' @param event \code{n x 1} numeric vector of status indicators of
#' whether an event was observed. Defaults to a vector of 1s, i.e. no censoring.
#' @param X \code{n x p} data.frame of observed covariate values
#' on which to train the estimator.
#' @param SL.library Library of algorithms to include in the binary classification
#' Super Learner. Should have the same structure as the \code{SL.library}
#' argument to the \code{SuperLearner} function in the \code{SuperLearner} package.
#' @param V Number of cross validation folds on which to train the Super Learner
#' classifier. Defaults to 10.
#'
#' @return An fitted binary regression for (complement of)
#' probability of censoring
p_delta <- function(event,
                    X,
                    V = 10,
                    SL.library){

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
