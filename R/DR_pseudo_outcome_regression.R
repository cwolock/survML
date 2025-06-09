#' Doubly-robust pseudo-outcome regression
#'
#' Generate estimates of conditional survival probability or conditional restrcited mean survival time
#' using doubly-robust pseudo-outcome regression with SuperLearner
#'
#' @param time \code{n x 1} numeric vector of observed
#' follow-up times. If there is censoring, these are the minimum of the
#' event and censoring times.
#' @param event \code{n x 1} numeric vector of status indicators of
#' whether an event was observed.
#' @param X \code{n x p} data.frame of observed covariate values
#' @param newX \code{m x p} data.frame of new observed covariate
#' values at which to obtain \code{m} predictions for the estimated algorithm.
#' Must have the same names and structure as \code{X}.
#' @param approx_times Numeric vector of length J2 giving times at which to
#' approximate integral appearing in the pseudo-outcomes
#' @param S_hat \code{n x J2} matrix of conditional event time survival function estimates
#' @param G_hat \code{n x J2} matrix of conditional censoring time survival function estimates
#' @param newtimes Numeric vector of times at which to generate oracle prediction function estimates. For outcome \code{"survival_probability"},
#' this is the times at which the survival function is to be estimated. For outcome \code{"restricted_survival_time"}, this is simply the restriction
#' time.
#' @param outcome Outcome type, either \code{"survival_probability"} or \code{"restricted_survival_time"}
#' @param SL.library Super Learner library
#' @param V Number of cross-validation folds, to be passed to \code{SuperLearner}
#'
#' @return Matrix of predictions corresponding to \code{newX} and \code{newtimes}.
#'
#' @examples
#' # This is a small simulation example
#' set.seed(123)
#' n <- 250
#' X <- data.frame(X1 = rnorm(n), X2 = rbinom(n, size = 1, prob = 0.5))
#'
#' T <- rexp(n, rate = exp(-2 + X[,1] - X[,2] + .5 *  X[,1] * X[,2]))
#' C <- rexp(n, exp(-2 -.5 * X[,1] - .25 * X[,2] + .5 * X[,1] * X[,2]))
#' C[C > 15] <- 15
#'
#' time <- pmin(T, C)
#' event <- as.numeric(T <= C)
#'
#' # Note that this a very small Super Learner library, for computational purposes.
#' SL.library <- c("SL.mean", "SL.glm")
#'
#' approx_times <- c(0, sort(unique(time)))
#'
#' # estimate conditional survival functions at approx_times
#' fit <- stackG(time = time,
#'               event = event,
#'               X = X,
#'               newX = X,
#'               newtimes = approx_times,
#'               direction = "prospective",
#'               bin_size = 0.1,
#'               time_basis = "continuous",
#'               surv_form = "PI",
#'               learner = "SuperLearner",
#'               time_grid_approx = approx_times,
#'               SL_control = list(SL.library = SL.library,
#'                                 V = 3))
#'
#' # use DR pseudo-outcome regression to (robustly) estimate survival at t = 5
#' DR_preds <- DR_pseudo_outcome_regression(time = time,
#'                                         event = event,
#'                                         X = X,
#'                                         newX = X,
#'                                         newtimes = 5,
#'                                         approx_times = approx_times,
#'                                         S_hat = fit$S_T_preds,
#'                                         G_hat = fit$S_C_preds,
#'                                         outcome = "survival_probability",
#'                                         SL.library = SL.library,
#'                                         V = 3)
#' DR_preds
#'
#' @export
DR_pseudo_outcome_regression <- function(time,
                                         event,
                                         X,
                                         newX,
                                         approx_times,
                                         S_hat,
                                         G_hat,
                                         newtimes,
                                         outcome,
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
          m <- S_hat[,which.min(abs(approx_times - tau))]/S_hat[,t]
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

    if (outcome == "survival_probability"){
      term1 <- Delta_t*(Y_t >= tau) / G_hat_Y
    } else if (outcome == "restricted_survival_time"){
      term1 <- Delta_t*pmin(Y_t, tau) / G_hat_Y
    }

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
    SL_preds <- stats::predict(SL_fit, newdata = newX)$pred
    DR_predictions[,i] <- SL_preds
  }
  return(DR_predictions)
}
