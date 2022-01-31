#' Estimate a conditional survival function
#'
#' @param time Observed time
#' @param event Indicator of event (vs censoring)
#' @param X Covariate matrix
#' @param newX Values of covariates at which to make a prediction
#' @param newtimes Times at which to make the survival function prediction
#' @param time_grid_approx Grid of time points on which to approximate PL integral
#' @param time_grid_eval Grid of time points on which to evaluate (e.g. approximate a MISE integral)
#' @param weights TBD, possible weighting of time when calculating MISE
#' @param test_time Observed times to evaluate against
#' @param test_event Event indicators to test against
#' @param test_X Covariates corresponding to \code{test_times} and \code{test_event}
#' @param SL.library SuperLearner library
#'
#' @return An object of class \code{conSurv}
#'
#' @export
#'
#' @examples
conSurv <- function(time,
                    event,
                    X,
                    newX,
                    newtimes,
                    time_grid_approx,
                    time_grid_eval = NULL,
                    weights = NULL,
                    test_time = NULL,
                    test_event = NULL,
                    test_X = NULL,
                    SL.library){

  # determine optimal models (currently using oracle tuning)
  P_Delta_opt <- estimate_p_delta(event = event,
                                            X = X,
                                            test_event = test_event,
                                            test_X = test_X,
                                  SL.library = SL.library)
  S_Y_1_opt <-estimate_f_y(time = time,
                                      event = event,
                                      X = X,
                                      test_time = test_time,
                                      test_event = test_event,
                                      test_X = test_X,
                                      time_grid_eval = time_grid_eval,
                                      censored = FALSE,
                           SL.library = SL.library)
  S_Y_0_opt <- estimate_f_y(time = time,
                                      event = event,
                                      X = X,
                                      test_time = test_time,
                                      test_event = test_event,
                                      test_X = test_X,
                                      time_grid_eval = time_grid_eval,
                                      censored = TRUE,
                            SL.library = SL.library)

  # fit optimal models
  #P_Delta_opt_preds <- predict(P_Delta_opt, newX = newX) # this is for my wrapped algorithms
  P_Delta_opt_preds <- predict(P_Delta_opt, newdata = newX)$pred # this is for SuperLearner

  S_Y_1_opt_preds <- predict(S_Y_1_opt,
                             newX = newX,
                             newtimes = time_grid_approx) # this was previously newtimes, which was possibly causing a HUGE issue
  S_Y_0_opt_preds <- predict(S_Y_0_opt,
                             newX = newX,
                             newtimes = time_grid_approx)

  estimate_S_T <- function(i){
    # get S_Y estimates up to t
    S_Y_1_curr <- S_Y_1_opt_preds[i,]
    S_Y_0_curr <- S_Y_0_opt_preds[i,]
    pi_curr <- P_Delta_opt_preds[i]

    S_T_ests <- compute_prodint(cdf_uncens = S_Y_1_curr,
                            cdf_cens = S_Y_0_curr,
                            p_uncens = pi_curr,
                            newtimes = newtimes,
                            time_grid = time_grid_approx)

    return(S_T_ests)
  }

  S_T_preds <- t(apply(X = as.matrix(seq(1, nrow(newX))),
                    MARGIN = 1,
                    FUN = estimate_S_T))

  estimate_S_C <- function(i){
    # get S_Y estimates up to t
    S_Y_1_curr <- S_Y_1_opt_preds[i,]
    S_Y_0_curr <- S_Y_0_opt_preds[i,]
    pi_curr <- P_Delta_opt_preds[i]

    # switch roles of T and C for estimating S_C
    S_T_ests <- compute_prodint(cdf_uncens = S_Y_0_curr,
                            cdf_cens = S_Y_1_curr,
                            p_uncens = 1-pi_curr,
                            newtimes = newtimes,
                            time_grid = time_grid_approx)

    return(S_T_ests)
  }

  S_C_preds <- t(apply(X = as.matrix(seq(1, nrow(newX))),
                       MARGIN = 1,
                       FUN = estimate_S_C))

  res <- list(S_T_preds = S_T_preds,
              S_C_preds = S_C_preds,
              P_Delta_preds = P_Delta_opt_preds,
              F_Y_1_preds = S_Y_1_opt_preds,
              F_Y_0_preds = S_Y_0_opt_preds,
              P_Delta_algo = P_Delta_opt,
              F_Y_1_algo = S_Y_1_opt,
              F_Y_0_algo = S_Y_0_opt)
  class(res) <- "conSurv"
  return(res)
}
