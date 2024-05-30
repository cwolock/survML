#' Construct one-step estimator of landmark Brier score
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
#' @param tau landmark time
#' @param preds n \times 1 vector of predictions
#'
#' @return An estimate of the Brier score
#'
#' @noRd
estimate_rsquared <- function(time,
                              event,
                              approx_times,
                              tau,
                              preds,
                              S_hat,
                              G_hat){

  KM_IFs <- calc_KM_IF(time = time,
                       event = event,
                       S_hat = S_hat,
                       G_hat = G_hat,
                       approx_times = approx_times)

  k <- min(which(approx_times >= tau))
  S_hat_k <- S_hat[,k]
  G_hat_k <- G_hat[,k]
  KM_IFs <- KM_IFs[,k]

  # calculate two components of EIF
  phi0 <- 2*preds*KM_IFs - KM_IFs
  phi_tilde_0 <- 2*preds*S_hat_k - preds^2 - S_hat_k - mean(2*preds*S_hat_k - preds^2 - S_hat_k)

  EIF <- phi0 + phi_tilde_0

  one_step <- mean(2*preds*S_hat_k - preds^2 - S_hat_k) + mean(EIF)
  plug_in <- mean(2*preds*S_hat_k - preds^2 - S_hat_k)

  V1_os <- one_step
  V1_plugin <- plug_in
  V2_os <- mean(S_hat_k + KM_IFs) - (mean(S_hat_k + KM_IFs))^2
  V2_plugin <- mean(S_hat_k) - (mean(S_hat_k))^2
  one_step <- one_step/V2_os
  plug_in <- plug_in/V2_plugin

  return(list(one_step = one_step,
              plug_in = plug_in,
              EIF = EIF,
              numerator = V1_os,
              denominator = V2_os))
}
