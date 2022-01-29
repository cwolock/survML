#' Estimate the conditional probability of Delta = 1 given X
#'
#' @param event Numeric vector of status indicators of whether an event was observed
#' @param X Data frame of observed covariate values on which to train the prediction algorithm
#' @param test_event Event indicators to test against
#' @param test_X Covariates corresponding to \code{test_event}
#' @param rescale A logical value indicating whether or not to rescale covariates to have equal variance.
#'        Scaling has an impact on some algorithms, such as LOESS.
#'
#' @return Optimal fit for estimating p_delta
#' @noRd
estimate_p_delta <- function(event, X, test_event, test_X, rescale = TRUE){

  # bws <- seq(0.1, 1, by = 0.1)
  # MSEs <- rep(NA, length(bws))
  #
  # for (i in 1:length(bws)){
  #   fit <- p_delta_nw(event = event,
  #                     X = X,
  #                     bw = bws[i],
  #                     kernel_type = "gaussian",
  #                     kernel_order = 2)
  #   MSEs[i] <- calculate_MSE(fit = fit,
  #                            test_event = test_event,
  #                            test_X = test_X)
  # }
  #
  # opt_bw <- bws[which.min(MSEs)]
  #
  # opt_fit <- p_delta_nw(event = event,
  #                       X = X,
  #                       bw = opt_bw,
  #                       kernel_type = "gaussian",
  #                       kernel_order = 2)
  SL.library <- c("SL.mean", "SL.glm", "SL.gam", "SL.randomForest", "SL.gbm")
  opt_fit <- SuperLearner::SuperLearner(Y = event,
                           X = X,
                           family = "binomial",
                           SL.library = SL.library,
                           method = "method.NNloglik",
                           verbose = FALSE)

  # opt_fit <- p_delta_rf(event = event,
  #                       X = X)

  return(opt_fit)
}
