#' Construct one-step estimator of MSE for predicting restricted survival time
#'
#' @param time n \times 1 vector of observed follow-up times
#' @param event n \times 1 vector of event indicators (1 = event, 0 = right censoring)
#' @param S_hat n \times J matrix of conditional event survival function estimates.
#'  Each row corresponds to an observation, and each column to one of the J times in
#'  the approximation grid.
#' @param G_hat n \times J matrix of conditional censoring survival function estimates.
#'  Each row corresponds to an observation, and each column to one of the J times in
#'  the approximation grid.
#' @param approx_times Time grid of length J to approximate integrals taken with respect to the
#'  conditional cumulative hazard.
#' @param tau restriction time
#' @param preds n \times 1 vector of predictions
#'
#' @return An estimate of the MSE
#'
#' @noRd
estimate_rmst_mse <- function(time,
                              event,
                              approx_times,
                              preds,
                              S_hat,
                              G_hat,
                              tau){

  n <- length(time)

  KM_IFs <- calc_KM_IF(time = time,
                       event = event,
                       S_hat = S_hat,
                       G_hat = G_hat,
                       approx_times = approx_times)

  S_hat_k <- S_hat[,approx_times <= tau]
  KM_IFs_k <- KM_IFs[,approx_times <= tau]
  pi_k <- S_hat_k + KM_IFs_k

  approx_times_k <- approx_times[approx_times <= tau]
  approx_times_k <- approx_times_k[-1]

  calc_phi_01_combined <- function(i){
    -sum((preds[i] - approx_times_k)^2 * diff(pi_k[i,])) +
      (preds[i] - tau)^2*pi_k[i,which(approx_times == tau)]
  }

  phi_01_combined <- unlist(lapply(1:n, FUN = calc_phi_01_combined))

  calc_phi_01 <- function(i){
    -sum((preds[i] - approx_times_k)^2 * diff(KM_IFs_k[i,])) + (preds[i] - tau)^2*KM_IFs_k[i,which(approx_times == tau)]
  }

  phi_01 <- unlist(lapply(1:n, FUN = calc_phi_01))

  calc_phi_tilde_01 <- function(i){
    -sum((preds[i] - approx_times_k)^2 * diff(S_hat_k[i,])) + (preds[i] - tau)^2*S_hat_k[i,which(approx_times == tau)]
  }
  phi_tilde_01 <- unlist(lapply(1:n, FUN = calc_phi_tilde_01))

  plug_in <- mean(phi_tilde_01)

  phi_tilde_01 <- phi_tilde_01 - plug_in

  if_func <- phi_01 + phi_tilde_01

  one_step <- mean(phi_01_combined)
  return(list(one_step = one_step,
              plug_in = plug_in,
              EIF = if_func))
}
