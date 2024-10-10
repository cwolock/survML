#' Generate K-fold cross-fit survival predictions for downstream use
#'
#' @return data frame giving results
#'
#' @export
DR_pseudo_outcome_regression <- function(time,
                                         event,
                                         X,
                                         newX,
                                         S_hat,
                                         G_hat,
                                         newtimes,
                                         outcome,
                                         approx_times,
                                         SL.library,
                                         V){
  DR_predictions <- matrix(NA, nrow = nrow(newX), ncol = length(newtimes))
  for (i in 1:length(newtimes)){
    tau <- newtimes[i]
    Delta_t <- event * (time <= tau) + (time > tau)
    Y_t <- pmin(time, tau)

    # calculate m_0 at a specific time t
    calc_one <- function(t, outcome){
      if (outcome == "survival_probability"){
        if (approx_times[t] >= tau){
          m <-  rep(1, nrow(S_hat))
        } else{
          m <- S_hat[,which(approx_times == tau)]/S_hat[,t]
        }
      } else if (outcome == "restricted_survival_time"){
        if (approx_times[t] >= tau){
          m <- rep(tau, nrow(S_hat))
        } else{
          int.vals <- t(sapply(1:nrow(S_hat), function(i) {
            indices <- which(approx_times <= tau & approx_times >= t)
            S_hat_ind <- S_hat[,indices]
            approx_times_ind <- approx_times[indices]
            vals <- diff(S_hat_ind[i,])*approx_times_ind[-length(approx_times_ind)]
            return(sum(vals))
          }))
          m1 <- -int.vals
          m2 <- tau * S_hat[,which(approx_times == tau)]
          m <- (m1 + m2)/S_hat[,t]
        }
      }
      return(m)
    }
    # calculate m_0 over the whole grid
    ms <- matrix(unlist(lapply(1:length(approx_times), FUN = calc_one, outcome = outcome)),
                 nrow = length(time))

    m_Y <- apply(X = matrix(1:length(Y_t)), MARGIN = 1,
                 FUN = function(x) ms[x,which.min(abs(approx_times - Y_t[x]))])

    G_hat_Y <- apply(X = matrix(1:length(Y_t)), MARGIN = 1,
                     FUN = function(x) G_hat[x,which.min(abs(approx_times - Y_t[x]))])

    term1 <- Delta_t*(Y_t >= tau) / G_hat_Y

    term2 <- (1 - Delta_t) * m_Y / G_hat_Y

    int.vals <- t(sapply(1:length(time), function(j) {
      vals <- diff(1/G_hat[j,])* ms[j,-ncol(ms)]
      if(any(approx_times[-1] > Y_t[j])){
        vals[approx_times[-1] > Y_t[j]] <- 0
      }
      sum(vals)
    }))

    term3 <- t(int.vals)

    DR_pseudo_outcome <- term1 + term2 - term3

    SL_fit <- SuperLearner::SuperLearner(Y = DR_pseudo_outcome,
                                         X = X,
                                         family = stats::gaussian(),
                                         SL.library = SL.library,
                                         method = "method.NNLS",
                                         cvControl = list(V = V),
                                         verbose = FALSE)
    SL_preds <- predict(SL_fit, newdata = newX)$pred
    DR_predictions[,i] <- SL_preds
  }
  return(DR_predictions)
}
