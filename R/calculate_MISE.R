#' Calculate the mean integrated squared error of a cdf estimator
#'
#' @param fit Fitted estimator of the cdf
#' @param time_grid Grid of time points on which to evaluate
#' @param test_time Observed times to evaluate against
#' @param test_event Event indicators to test against
#' @param test_X Covariates corresponding to \code{test_times}
#' @param censored Logical, indicates whether to run regression on censored observations (vs uncensored)
#' @param weights TBD, maybe want to weight different time points?
#'
#' @return A scalar estimate of MISE
#'
#' @noRd
calculate_MISE <- function(fit, time_grid, test_time, test_event, test_X, censored, weights = NULL){

  if (censored){
    test_time <- test_time[!as.logical(test_event)]
    test_X <- test_X[!as.logical(test_event),]
  } else{
    test_time <- test_time[as.logical(test_event)]
    test_X <- test_X[as.logical(test_event),]
  }

  preds <- predict(fit, newX = test_X, newtimes = time_grid)
  MSEs <- rep(NA, length(time_grid))
  for (i in 1:length(time_grid)){
    truth <- ifelse(test_time <= time_grid[i], 1, 0)
    MSE <- mean((preds[,i] - truth)^2, na.rm = TRUE)
    MSEs[i] <- MSE
  }
  MISE <- mean(MSEs)

  return(MISE)
}
