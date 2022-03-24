#' Estimate a conditional survival function by binary regression on a grid
#'
#' @param time Observed time
#' @param event Indicator of event (vs censoring)
#' @param X Covariate matrix
#' @param newX Values of covariates at which to make a prediction
#' @param newtimes Times at which to make the survival function prediction
#' @param time_grid_approx Grid of time points on which to approximate PL integral
#' @param test_time Observed times to evaluate against
#' @param test_event Event indicators to test against
#' @param test_X Covariates corresponding to \code{test_times} and \code{test_event}
#' @param V CV fold number, required for tuned algorithms (xgboost, ranger)
#' @param entry Variable indicating time of entry into the study (truncation variable)
#' @param algorithm Which binary classification algorithm to use (xgboost, ranger, gam with smoothing splines)
#'
#' @return An object of class \code{gridSurv}
#'
#' @export
#'
#' @examples
gridSurv_reverse <- function(time,
                             X,
                             entry,
                             newX,
                             newtimes,
                             time_grid_approx,
                             bin_size = NULL,
                             algorithm = "xgboost",
                             V = 10,
                             time_basis,
                             surv_form = "PI",
                             tuning_params = NULL){


  if (algorithm == "xgboost"){ # do xgboost if speed not a concern
    # determine optimal models
    S_Y_opt <- f_y_stack_xgboost(time = time,
                                 event = NULL,
                                 X = X,
                                 censored = NULL,
                                 bin_size = bin_size,
                                 V = V,
                                 time_basis = time_basis,
                                 tuning_params = tuning_params)
    S_Y_opt_preds <- predict(S_Y_opt,
                             newX = newX,
                             newtimes = time_grid_approx)

    F_W_opt <-f_w_stack_xgboost(time = time,
                                event = NULL,
                                X = X,
                                censored = NULL,
                                bin_size = bin_size,
                                V = V,
                                entry = entry,
                                time_basis = time_basis,
                                direction = "reverse",
                                tuning_params = tuning_params)
    F_W_opt_preds <- predict(F_W_opt,
                             newX = newX,
                             newtimes = time_grid_approx)
  } else if (algorithm == "earth"){
    S_Y_opt <- f_y_stack_earth(time = time,
                               event = NULL,
                               X = X,
                               censored = NULL,
                               bin_size = bin_size,
                               time_basis = time_basis)
    S_Y_opt_preds <- predict(S_Y_opt,
                             newX = newX,
                             newtimes = time_grid_approx)

    F_W_opt <-f_w_stack_earth(time = time,
                              event = NULL,
                              X = X,
                              censored = NULL,
                              bin_size = bin_size,
                              entry = entry,
                              time_basis = time_basis,
                              direction = "reverse")
    F_W_opt_preds <- predict(F_W_opt,
                             newX = newX,
                             newtimes = time_grid_approx)
  }

  estimate_S_T <- function(i){
    # get S_Y estimates up to t
    S_Y_curr <- S_Y_opt_preds[i,]
    F_W_curr <- F_W_opt_preds[i,]
    if (surv_form == "PI"){
      S_T_ests <-compute_prodint_reverse(F_T_preds = S_Y_curr,
                                 F_W_preds = F_W_curr,
                                 newtimes = newtimes,
                                 time_grid = time_grid_approx)
    } else if (surv_form == "exp"){
      S_T_ests <-compute_exponential_reverse(F_T_preds = S_Y_curr,
                                         F_W_preds = F_W_curr,
                                         newtimes = newtimes,
                                         time_grid = time_grid_approx)
    }

    return(S_T_ests)
  }

  S_T_preds <- t(apply(X = as.matrix(seq(1, nrow(newX))),
                       MARGIN = 1,
                       FUN = estimate_S_T))

  res <- list(S_T_preds = S_T_preds,
              F_T_preds = S_Y_opt_preds,
              F_W_preds = F_W_opt_preds,
              F_T_opt = S_Y_opt,
              F_W_opt = F_W_opt)
  class(res) <- "gridSurv_reverse"
  return(res)
}
