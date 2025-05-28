#' Construct one-step estimator of the concordance index
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
#' @return An estimate of the C-index
#'
#' @noRd
estimate_cindex <- function(time,
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

  k <- sum(approx_times <= tau)

  calc_phi_01 <- function(j){
    fx <- preds[j]
    varphi_x <- KM_IFs_k[j,]
    exceed_probs1 <- -rowSums(sweep(S_hat_k[,-k], MARGIN=2, diff(varphi_x), `*`))
    exceed_probs2 <- -rowSums(sweep(t(diff(t(S_hat_k))), MARGIN=2, varphi_x[-k], `*`))
    int <- mean(ifelse(fx >= preds, 1, 0)* exceed_probs1 + ifelse(preds > fx, 1, 0)* exceed_probs2)
    return(int)
  }

  calc_phi_tilde_01 <- function(j){
    fx <- preds[j]
    Sx <- S_hat_k[j,]
    exceed_probs1 <- -rowSums(sweep(S_hat_k[,-k], MARGIN=2, diff(Sx), `*`))
    exceed_probs2 <- -rowSums(sweep(t(diff(t(S_hat_k))), MARGIN=2, Sx[-k], `*`))
    int <- mean(ifelse(fx >= preds, 1, 0)* exceed_probs1 + ifelse(preds > fx, 1, 0)* exceed_probs2)
    return(int)
  }

  calc_phi_01_combined <- function(j){
    fx <- preds[j]
    Sx <- S_hat_k[j,]
    varphi_x <- KM_IFs_k[j,]
    pi_x <- Sx + varphi_x
    pi_k <- S_hat_k + KM_IFs_k
    exceed_probs1 <- -rowSums(sweep(pi_k[,-k], MARGIN=2, diff(pi_x), `*`))
    exceed_probs2 <- -rowSums(sweep(t(diff(t(pi_k))), MARGIN=2, pi_x[-k], `*`))
    int <- mean(ifelse(fx >= preds, 1, 0)* exceed_probs1 + ifelse(preds > fx, 1, 0)* exceed_probs2)
    return(int)
  }

  calc_phi_02_combined <- function(j){
    Sx <- S_hat_k[j,]
    varphi_x <- KM_IFs_k[j,]
    pi_x <- Sx + varphi_x
    pi_k <- S_hat_k + KM_IFs_k
    exceed_probs1 <- -rowSums(sweep(pi_k[,-k], MARGIN=2, diff(pi_x), `*`))
    exceed_probs2 <- -rowSums(sweep(t(diff(t(pi_k))), MARGIN=2, pi_x[-k], `*`))
    int <- mean(exceed_probs1 + exceed_probs2)
    return(int)
  }

  calc_phi_02 <- function(j){
    varphi_x <- KM_IFs_k[j,]
    exceed_probs1 <- -rowSums(sweep(S_hat_k[,-k], MARGIN=2, diff(varphi_x), `*`))
    exceed_probs2 <- -rowSums(sweep(t(diff(t(S_hat_k))), MARGIN=2, varphi_x[-k], `*`))
    int <- mean(exceed_probs1 + exceed_probs2)
    return(int)
  }

  calc_phi_tilde_02 <- function(j){
    Sx <- S_hat_k[j,]
    exceed_probs1 <- -rowSums(sweep(S_hat_k[,-k], MARGIN=2, diff(Sx), `*`))
    exceed_probs2 <- -rowSums(sweep(t(diff(t(S_hat_k))), MARGIN=2, Sx[-k], `*`))
    int <- mean(exceed_probs1 + exceed_probs2)
    return(int)
  }

  phi_01 <- unlist(lapply(1:n, FUN = calc_phi_01))
  phi_tilde_01_uncentered <- unlist(lapply(1:n, FUN = calc_phi_tilde_01))
  phi_02 <- unlist(lapply(1:n, FUN = calc_phi_02))
  phi_tilde_02_uncentered <- unlist(lapply(1:n, FUN = calc_phi_tilde_02))

  phi_tilde_01 <- phi_tilde_01_uncentered - mean(phi_tilde_01_uncentered)
  phi_tilde_02 <- phi_tilde_02_uncentered - mean(phi_tilde_02_uncentered)

  phi_01_combined <- unlist(lapply(1:n, FUN = calc_phi_01_combined))
  phi_02_combined <- unlist(lapply(1:n, FUN = calc_phi_02_combined))

  V_1 <- mean(phi_tilde_01_uncentered)/2
  V_2 <- mean(phi_tilde_02_uncentered)/2

  V_1_alternative <- mean(phi_01_combined)/2
  V_2_alternative <- mean(phi_02_combined)/2

  one_step <- V_1_alternative/V_2_alternative

  EIF <- (phi_01 + phi_tilde_01)/V_2 - V_1/(V_2^2)*(phi_02 + phi_tilde_02)
  plug_in <- V_1/V_2

  return(list(one_step = one_step,
              plug_in = plug_in,
              EIF = EIF,
              numerator = V_1_alternative,
              denominator = V_2_alternative))
}
