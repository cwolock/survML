#' Wrapper for various f_w stacked algorithms
#'
#' @param time \code{n x 1} numeric vector of observed
#' follow-up times If there is censoring, these are the minimum of the
#' event and censoring times.
#' @param event \code{n x 1} numeric vector of status indicators of
#' whether an event was observed. Defaults to a vector of 1s, i.e. no censoring.
#' @param entry Study entry variable, if applicable. Defaults to \code{NULL},
#' indicating that there is no truncation.
#' @param X \code{n x p} data.frame of observed covariate values
#' on which to train the estimator.
#' @param censored Logical, indicates whether to run regression on censored
#' observations (\code{event == 0}) vs. uncensored (\code{event == 1}).
#' @param bin_size Size of time bin on which to discretize for estimation
#' of cumulative probability functions. Can be a number between 0 and 1,
#' indicating the size of quantile grid (e.g. \code{0.1} estimates
#' the cumulative probability functions on a grid based on deciles of
#' observed \code{time}s). If \code{NULL}, creates a grid of
#' all observed \code{time}s.
#' @param time_basis How to treat time for training the binary
#' classifier. Options are \code{"continuous"} and \code{"dummy"}, meaning
#' an indicator variable is included for each time in the time grid.
#' @param SL.library Library of algorithms to include in the binary classification
#' Super Learner. Should have the same structure as the \code{SL.library}
#' argument to the \code{SuperLearner} function in the \code{SuperLearner} package.
#' @param V Number of cross validation folds on which to train the Super Learner
#' classifier. Defaults to 10.
#' @param obsWeights Optional observation weights. These weights are passed
#' directly to \code{SuperLearner}, which in turn passes them directly to the
#' prediction algorithms.
#'
#' @return An fitted pooled binary regression for the truncation distribution
#'
#' @noRd
f_w_stack <- function(time,
                      event,
                      X,
                      entry,
                      censored,
                      bin_size,
                      bin_variable,
                      learner = "SuperLearner",
                      SL_control,
                      xgb_control,
                      time_basis){

  if (learner == "SuperLearner"){
    fit <- f_w_stack_SuperLearner(time = time,
                                  event = event,
                                  X = X,
                                  censored = censored,
                                  bin_size = bin_size,
                                  bin_variable = bin_variable,
                                  SL_control = SL_control,
                                  time_basis = time_basis,
                                  entry = entry)

  }

  return(fit)
}
