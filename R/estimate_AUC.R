#' Construct one-step estimator of the landmark time AUC
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
#' @param robust Whether or not to perform the DR debiasing approach
#'
#' @return An estimate of the AUC
#'
#' @noRd
estimate_AUC<- function(time,
                        event,
                        approx_times,
                        tau,
                        preds,
                        S_hat,
                        G_hat,
                        robust){

  n <- length(time)

  KM_IFs <- calc_KM_IF(time = time,
                       event = event,
                       S_hat = S_hat,
                       G_hat = G_hat,
                       approx_times = approx_times)

  k <- min(which(approx_times >= tau))
  S_hat_k <- S_hat[,k]
  G_hat_k <- G_hat[,k]
  KM_IFs <- KM_IFs[,k]

  calc_phi_01 <- function(j){
    fx <- preds[j]
    varphi_x <- KM_IFs[j]
    int <- mean(ifelse(preds > fx, 1, 0) * (1 - S_hat_k) - ifelse(preds <= fx, 1, 0) * S_hat_k)
    return(varphi_x*int)
  }
  calc_phi_02 <- function(j){
    varphi_x <- KM_IFs[j]
    int <- mean((1 - S_hat_k) - S_hat_k)
    return(varphi_x * int)
  }

  calc_phi_tilde_01 <- function(j){
    fx <- preds[j]
    Sx <- S_hat_k[j]
    int <- mean(ifelse(preds > fx, 1, 0) * (1 - S_hat_k) * (Sx) +
                  ifelse(preds <= fx, 1, 0) * (S_hat_k) * (1 - Sx))
    return(int)
  }
  calc_phi_tilde_02 <- function(j){
    Sx <- S_hat_k[j]
    int <- mean((1 - S_hat_k) * Sx + S_hat_k * (1 - Sx))
    return(int)
  }

  calc_phi_01_combined <- function(j){
    fx <- preds[j]
    Sx <- S_hat_k[j]
    varphi_x <- KM_IFs[j]
    int <- mean(ifelse(preds > fx, 1, 0) * (1 - S_hat_k - KM_IFs) * (Sx + varphi_x) +
                  ifelse(preds <= fx, 1, 0) * (S_hat_k + KM_IFs) * (1 - Sx - varphi_x))
  }
  calc_phi_02_combined <- function(j){
    Sx <- S_hat_k[j]
    varphi_x <- KM_IFs[j]
    int <- mean((1 - S_hat_k - KM_IFs) * (Sx + varphi_x) +
                 (S_hat_k + KM_IFs) * (1 - Sx - varphi_x))
  }

  phi_01 <- unlist(lapply(1:n, FUN = calc_phi_01))
  phi_tilde_01_uncentered <- unlist(lapply(1:n, FUN = calc_phi_tilde_01))
  phi_01_combined <- unlist(lapply(1:n, FUN = calc_phi_01_combined))
  phi_02_combined <- unlist(lapply(1:n, FUN = calc_phi_02_combined))
  phi_02 <- unlist(lapply(1:n, FUN = calc_phi_02))
  phi_tilde_02_uncentered <- unlist(lapply(1:n, FUN = calc_phi_tilde_02))

  phi_tilde_01 <- phi_tilde_01_uncentered - mean(phi_tilde_01_uncentered)
  phi_tilde_02 <- phi_tilde_02_uncentered - mean(phi_tilde_02_uncentered)

  V_1 <- mean(phi_tilde_01_uncentered)/2
  V_2 <- mean(phi_tilde_02_uncentered)/2

  if_func_1 <- phi_01 + phi_tilde_01
  if_func_2 <- phi_02 + phi_tilde_02

  V_1_os <- V_1 + mean(if_func_1)
  V_1_alternative <- mean(phi_01_combined)/2
  V_2_os <- V_2 + mean(if_func_2)
  V_2_alternative <- mean(phi_02_combined)/2

  one_step <- V_1_alternative/V_2_alternative

  EIF <- (phi_01 + phi_tilde_01)/V_2 - V_1/(V_2^2)*(phi_02 + phi_tilde_02)
  plug_in <- V_1/V_2

  V_1_final <- ifelse(robust, V_1_alternative, V_1_os)
  V_2_final <- ifelse(robust, V_2_alternative, V_2_os)

  return(list(one_step = one_step,
              plug_in = plug_in,
              EIF = EIF,
              numerator = V_1_final,
              denominator = V_2_final))
}
