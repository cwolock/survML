#' Estimate the conditional cdf of Y given Delta and X
#'
#' @param time Observed time
#' @param event Numeric vector of status indicators of whether an event was observed
#' @param X Data frame of observed covariate values on which to train the prediction algorithm
#' @param test_time Observed times to evaluate against
#' @param test_event Event indicators to test against
#' @param test_X Covariates corresponding to \code{test_event}
#' @param time_grid_eval Grid of time points on which to evaluate (e.g. approximate a MISE integral)
#' @param censored Logical, indicates whether to run regression on censored observations (vs uncensored)
#' @param rescale A logical value indicating whether or not to rescale covariates to have equal variance.
#'        Scaling has an impact on some algorithms, such as LOESS.
#'
#' @return Optimal fit for estimating f_y
#' @noRd
estimate_f_y <- function(time, event, X, test_time, test_event, test_X, time_grid_eval, censored, rescale = TRUE){
#
#   bws <- expand.grid(x = seq(0.1, 1, by = 0.1), y = c(seq(0.5, 1, by = 0.1), seq(2, 5, by = 1)))
#
#   MISEs <- rep(NA, nrow(bws))
#
#   for (i in 1:nrow(bws)){
#     bw_x <- bws[i,1]
#     bw_y <- bws[i,2]
#     # fit <- f_y_smoothllkern(time = time,
#     #                     event = event,
#     #                     X = X,
#     #                     censored = censored,
#     #                     bw = bw_x,
#     #                     bwy = bw_y,
#     #                     kernel_type = "gaussian",
#     #                     kernel_order = 2)
#     fit <- f_y_qrf(time = time,
#                             event = event,
#                             X = X,
#                             censored = censored)
#     MISEs[i] <- calculate_MISE(fit = fit,
#                                time_grid = time_grid_eval,
#                                test_time = test_time,
#                                test_event = test_event,
#                                test_X = test_X,
#                                censored = censored)
#   }


  # pick optimal tuning parameters
  #opt_bw<- bws[which.min(MISEs),]

  # opt_fit <- f_y_smoothllkern(time = time,
  #                         event = event,
  #                         X = X,
  #                         censored = censored,
  #                         bw = opt_bw[1,1],
  #                         bwy = opt_bw[1,2],
  #                         kernel_type = "gaussian",
  #                         kernel_order = 2)
  bin_sizes <- c(0.025, 0.05, 0.1, 0.2)
  MISEs <- rep(NA, length(bin_sizes))
  #
  for (i in 1:length(bin_sizes)){
    fit <- f_y_discSL(time = time,
                       event = event,
                       X = X,
                       censored = censored,
                       bin_size = bin_sizes[i])
    MISEs[i] <- calculate_MISE(fit = fit,
                               time_grid = time_grid_eval,
                               test_time = test_time,
                               test_event = test_event,
                               test_X = test_X,
                               censored = censored)
  }


  # pick optimal tuning parameters
  opt_bin_size <- bin_sizes[which.min(MISEs)]
  opt_fit <- f_y_discSL(time = time,
                 event = event,
                 X = X,
                 censored = censored,
                 bin_size = opt_bin_size)

  return(opt_fit)
}
