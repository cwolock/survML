#' Compute the plug-in product-limit survival function estimator
#'
#' @param cdf_uncens Matrix of predictions of the cdf of the uncensored times (F_Y_1) on chosen time grid
#' @param cdf_cens Predictions of the cdf of the censored times (F_Y_0) on chosen time grid
#' @param p_uncens Prediction of the probability of being uncensored
#' @param newtimes Times at which to make the prediction
#' @param time_grid Grid of time points over which to discretize the product integral
#' @param denom_method Method of computing the denominator
#'
#' @return A vector of estimates of the survival function over \code{time_grid}
#'
#' @noRd
compute_prodint_reverse <- function(F_T_preds,
                            F_W_preds,
                            newtimes,
                            time_grid){

  estimate_S_T <- function(t){
    curr_length <- sum(time_grid < t)
    dF_T <- c(F_T_preds[1], diff(F_T_preds))
    S_W_preds_left <- c(1, 1 - F_W_preds[-length(F_W_preds)])
    F_T_curr <- F_T_preds[curr_length:length(F_T_preds)] # prob of being at risk at time t
    dF_T_curr <- dF_T[curr_length:length(F_T_preds)] # change in cdf at time t
    S_W_curr <- S_W_preds_left[curr_length:length(S_W_preds_left)]
    S_T_complement <- prod(1 - dF_T_curr/(F_T_curr*S_W_curr))
    S_T_est <- 1 - S_T_complement
    return(S_T_est)
  }

  S_T_ests <- apply(X = as.matrix(newtimes),
                    MARGIN = 1,
                    FUN = estimate_S_T)

  return(S_T_ests)
}
