#' Calculate the mean squared error of a binary probability estimator
#'
#' @param fit Fitted estimator of the cdf
#' @param test_event Event indicators to test against
#' @param test_X Covariates corresponding to \code{test_times}
#' @param weights TBD, maybe want to weight different time points?
#'
#' @return A scalar estimate of MSE
#'
#' @noRd
calculate_MSE <- function(fit, test_event, test_X, weights = NULL){

  preds <- predict(fit, newX = test_X)
  MSE <- mean((preds - test_event)^2, na.rm = T)

  return(MSE)
}
