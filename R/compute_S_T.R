#' Compute the plug-in product-limit survival function estimator
#'
#' @param cdf_uncens Matrix of predictions of the cdf of the uncensored times (F_Y_1) on chosen time grid
#' @param cdf_cens Predictions of the cdf of the censored times (F_Y_0) on chosen time grid
#' @param p_uncens Prediction of the probability of being uncensored
#' @param newtimes Times at which to make the prediction
#' @param time_grid Grid of time points over which to discretize the product integral
#'
#' @return A vector of estimates of the survival function over \code{time_grid}
#'
#' @noRd
compute_S_T <- function(cdf_uncens, cdf_cens, p_uncens, newtimes, time_grid){

  estimate_S_T <- function(t){
    curr_length <- sum(time_grid <= t)

    # get S_Y estimates up to t
    S_Y_1_curr <- cdf_uncens[1:curr_length]
    S_Y_0_curr <- cdf_cens[1:curr_length]

    dF_Y_1_pred <- c(S_Y_1_curr[1], diff(S_Y_1_curr))

    S_Y_1_pred_left <- c(1, 1-S_Y_1_curr[-length(S_Y_1_curr)])# probability of being "at risk" at time t
    ### CHECK TO MAKE SURE THIS IS CORRECT WITH THE DISCRETIZATION OF TIME
    S_Y_0_pred_left <- c(1, 1-S_Y_0_curr[-length(S_Y_0_curr)])# probability of being "at risk" at time t

    low_left <- S_Y_1_pred_left * p_uncens
    low_right <- S_Y_0_pred_left * (1 - p_uncens)
    # product form
    S_T_est <- prod(1 - p_uncens * dF_Y_1_pred/(low_left + low_right))

    return(S_T_est)
  }

  S_T_ests <- apply(X = as.matrix(newtimes),
                    MARGIN = 1,
                    FUN = estimate_S_T)

  return(S_T_ests)
}
