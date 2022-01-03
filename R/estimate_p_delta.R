#' Estimate the conditional probability of Delta = 1 given X
#'
#' @param event Numeric vector of status indicators of whether an event was observed
#' @param X Data frame of observed covariate values on which to train the prediction algorithm
#' @param rescale A logical value indicating whether or not to rescale covariates to have equal variance.
#'        Scaling has an impact on some algorithms, such as LOESS.
#'
#' @return A return object
#' @noRd
#'
#' @examples
estimate_p_delta <- function(event, X, rescale = TRUE){
  p_delta_fit <- p_delta_logit_reg(event = event, X = X)
  print(p_delta_fit)
  p_delta_pred <- predict(p_delta_fit, newX = data.frame(X[1,]))
  p_delta_fit2 <- p_delta_nw(event = event, X = X, span = 0.5)
  print(p_delta_fit2)
  p_delta_pred2 <- predict(p_delta_fit2, newX = data.frame(X[1,]))
  p_delta_fit3 <- p_delta_mean(event = event, X = X)
  print(p_delta_fit3)
  p_delta_pred3 <- predict(p_delta_fit3, newX = data.frame(X[1,]))
  return(c(p_delta_pred, p_delta_pred2, p_delta_pred3))
}
