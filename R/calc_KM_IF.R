#' Compute the Kaplan-Meier influence function for use in one-step debiasing
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
#'
#' @return An n \times J matrix of KM influence function estimates.
#'
#' @noRd
calc_KM_IF <- function(time,
                       event,
                       S_hat,
                       G_hat,
                       approx_times){

  n <- length(time)
  # calculate integral for KM influence function
  int.vals <- t(sapply(1:n, function(i) {
    vals <- diff(1/S_hat[i,])* 1/ G_hat[i,-ncol(G_hat)]
    if(any(approx_times[-1] > time[i])){
      vals[approx_times[-1] > time[i]] <- 0
    }
    c(0,cumsum(vals))
  }))
  S_hat_Y <- sapply(1:n, function(i) stepfun(approx_times, c(1,S_hat[i,]), right = FALSE)(time[i]))
  G_hat_Y <- sapply(1:n, function(i) stepfun(approx_times, c(1,G_hat[i,]), right = TRUE)(time[i]))

  calc_one <- function(i){
    t0 <- approx_times[i]
    S_hat_i <- S_hat[,i]
    inner.func.1 <- ifelse(time <= t0 & event == 1, 1/(S_hat_Y * G_hat_Y), 0 )
    inner.func.2 <- int.vals[,i]
    KM_IF <- -S_hat_i * ( inner.func.1 - inner.func.2)
    return(KM_IF)
  }

  KM_IFs <- matrix(unlist(lapply(1:length(approx_times), FUN = calc_one)),
                   nrow = n)
  return(KM_IFs)
}
