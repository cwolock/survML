#' Estimate a conditional survival function by binary regression on a grid
#'
#' @param time Observed time
#' @param event Indicator of event (vs censoring)
#' @param X Covariate matrix
#' @param newX Values of covariates at which to make a prediction
#' @param newtimes Times at which to make the survival function prediction
#' @param time_grid_approx Grid of time points on which to approximate PL integral or cumulative hazard. Defaults to observed times
#' @param bin_size Size of time bin on which to discretize for estimation
#' @param denom_method Stratified or marginal denominator
#' @param algorithm Which binary classification algorithm to use
#' @param V CV fold number, required for tuned algorithms (xgboost, ranger)
#' @param entry Variable indicating time of entry into the study (truncation variable) if applicable
#' @param time_basis How to treat time (continuous or dummy)
#' @param surv_form Product integral or exponential mapping
#' @param tuning_params Tuning parameters for binary classification
#' @param direction Prospective or retrospective setting
#'
#' @return An object of class \code{survMLc}
#'
#' @export
#'
#' @examples
survMLc <- function(time,
                    event,
                    X,
                    newX,
                    newtimes,
                    time_grid_approx = sort(unique(time)),
                    bin_size = NULL,
                    denom_method = "stratified",
                    algorithm = "xgboost",
                    V = 10,
                    entry = NULL,
                    time_basis,
                    surv_form = "PI",
                    tuning_params = NULL,
                    direction = "prospective",
                    SL.library = NULL,
                    tau = NULL){

  if (direction == "retrospective"){
    if (is.null(tau)){
      tau <- max(time)
    }
    time <- tau - time
    entry <- tau - time
    event <- rep(1, length(time))
    newtimes <- tau - newtimes
    approx_times <- tau - approx_times
    denom_method <- "marginal"
    P_Delta_opt_preds <- rep(1, nrow(newX))
  }

  # if there is a censoring probability to estimate, i.e. if there is censoring
  if (sum(event == 0) != 0){
    P_Delta_opt <- p_delta(event = event,
                           X = X,
                           V = V,
                           tuning_params = tuning_params,
                           algorithm = algorithm,
                           SL.library = SL.library)
    P_Delta_opt_preds <- predict(P_Delta_opt, newX = newX) # this is for my wrapped algorithms

    if (denom_method == "stratified"){
      S_Y_0_opt <- f_y_stack(time = time,
                             event = event,
                             X = X,
                             censored = TRUE,
                             bin_size = bin_size,
                             V = V,
                             time_basis = time_basis,
                             tuning_params = tuning_params,
                             algorithm = algorithm,
                             SL.library = SL.library)
      S_Y_0_opt_preds <- predict(S_Y_0_opt,
                                 newX = newX,
                                 newtimes = time_grid_approx)
    } else{
      S_Y_opt <- f_y_stack(time = time,
                           event = event,
                           X = X,
                           censored = NULL,
                           bin_size = bin_size,
                           V = V,
                           time_basis = time_basis,
                           tuning_params = tuning_params,
                           algorithm = algorithm,
                           SL.library = SL.library)
      S_Y_opt_preds <- predict(S_Y_opt,
                               newX = newX,
                               newtimes = time_grid_approx)
    }
  }

  S_Y_1_opt <- f_y_stack(time = time,
                         event = event,
                         X = X,
                         censored = FALSE,
                         bin_size = bin_size,
                         V = V,
                         time_basis = time_basis,
                         tuning_params = tuning_params,
                         algorithm = algorithm,
                         SL.library = SL.library)

  if (!is.null(entry)){ # if a truncation variable is given
    if (denom_method == "stratified"){
      F_W_1_opt <- f_w_stack(time = time,
                             event = event,
                             X = X,
                             censored = FALSE,
                             bin_size = bin_size,
                             V = V,
                             entry = entry,
                             time_basis = time_basis,
                             tuning_params = tuning_params,
                             algorithm = algorithm,
                             SL.library = SL.library)
      F_W_1_opt_preds <- predict(F_W_1_opt,
                                 newX = newX,
                                 newtimes = time_grid_approx)
      F_W_0_opt <- f_w_stack(time = time,
                             event = event,
                             X = X,
                             censored = TRUE,
                             bin_size = bin_size,
                             V = V,
                             entry = entry,
                             time_basis = time_basis,
                             tuning_params = tuning_params,
                             algorithm = algorithm,
                             SL.library = SL.library)
      F_W_0_opt_preds <- predict(F_W_0_opt,
                                 newX = newX,
                                 newtimes = time_grid_approx)
    } else{ # retrospective setting is automatically marginal denominator
      F_W_opt <-f_w_stack(time = time,
                          event = event,
                          X = X,
                          censored = NULL,
                          bin_size = bin_size,
                          V = V,
                          entry = entry,
                          time_basis = time_basis,
                          tuning_params = tuning_params,
                          algorithm = algorithm,
                          SL.library = SL.library)
      F_W_opt_preds <- predict(F_W_opt,
                               newX = newX,
                               newtimes = time_grid_approx)
    }
  }




  S_Y_1_opt_preds <- predict(S_Y_1_opt,
                             newX = newX,
                             newtimes = time_grid_approx)

  if (direction == "retrospective"){
    S_Y_opt_preds <- S_Y_1_opt_preds
  }

  estimate_S_T <- function(i){
    # get S_Y estimates up to t
    S_Y_1_curr <- S_Y_1_opt_preds[i,]

    pi_curr <- P_Delta_opt_preds[i]

    if (denom_method == "stratified"){
      S_Y_0_curr <- S_Y_0_opt_preds[i,]
      if (!is.null(entry)){ # truncation
        F_W_0_curr <- F_W_0_opt_preds[i,]
        F_W_1_curr <- F_W_1_opt_preds[i,]
        if (surv_form == "PI"){
          S_T_ests <-compute_prodint(cdf_uncens = S_Y_1_curr,
                                     cdf_cens = S_Y_0_curr,
                                     entry_uncens = F_W_1_curr,
                                     entry_cens = F_W_0_curr,
                                     p_uncens = pi_curr,
                                     newtimes = newtimes,
                                     time_grid = time_grid_approx,
                                     denom_method = denom_method,
                                     truncation = TRUE)
        } else if (surv_form == "exp"){
          S_T_ests <-compute_exponential(cdf_uncens = S_Y_1_curr,
                                         cdf_cens = S_Y_0_curr,
                                         entry_uncens = F_W_1_curr,
                                         entry_cens = F_W_0_curr,
                                         p_uncens = pi_curr,
                                         newtimes = newtimes,
                                         time_grid = time_grid_approx,
                                         denom_method = denom_method,
                                         truncation = TRUE)
        }
      } else{ # no truncation
        if (surv_form == "PI"){
          S_T_ests <- compute_prodint(cdf_uncens = S_Y_1_curr,
                                      cdf_cens = S_Y_0_curr,
                                      p_uncens = pi_curr,
                                      newtimes = newtimes,
                                      time_grid = time_grid_approx,
                                      denom_method = denom_method,
                                      truncation = FALSE)
        } else if (surv_form == "exp"){
          S_T_ests <- compute_exponential(cdf_uncens = S_Y_1_curr,
                                          cdf_cens = S_Y_0_curr,
                                          p_uncens = pi_curr,
                                          newtimes = newtimes,
                                          time_grid = time_grid_approx,
                                          denom_method = denom_method,
                                          truncation = FALSE)
        }
      }
    } else{ # marginal denominator, or retrospective
      S_Y_curr <- S_Y_opt_preds[i,]
      if (!is.null(entry)){ # truncation
        F_W_curr <- F_W_opt_preds[i,]
        if (surv_form == "PI"){
          S_T_ests <-compute_prodint(cdf_uncens = S_Y_1_curr,
                                     cdf_marg = S_Y_curr,
                                     p_uncens = pi_curr,
                                     entry_marg = F_W_curr,
                                     newtimes = newtimes,
                                     time_grid = time_grid_approx,
                                     denom_method = denom_method,
                                     truncation = TRUE)
        } else if (surv_form == "exp"){
          S_T_ests <-compute_exponential(cdf_uncens = S_Y_1_curr,
                                         cdf_marg = S_Y_curr,
                                         p_uncens = pi_curr,
                                         entry_marg = F_W_curr,
                                         newtimes = newtimes,
                                         time_grid = time_grid_approx,
                                         denom_method = denom_method,
                                         truncation = TRUE)
        }
      } else{
        if (surv_form == "PI"){
          S_T_ests <- compute_prodint(cdf_uncens = S_Y_1_curr,
                                      cdf_marg = S_Y_curr,
                                      p_uncens = pi_curr,
                                      newtimes = newtimes,
                                      time_grid = time_grid_approx,
                                      denom_method = denom_method,
                                      truncation = FALSE)
        } else if (surv_form == "exp"){
          S_T_ests <- compute_exponential(cdf_uncens = S_Y_1_curr,
                                          cdf_marg = S_Y_curr,
                                          p_uncens = pi_curr,
                                          newtimes = newtimes,
                                          time_grid = time_grid_approx,
                                          denom_method = denom_method,
                                          truncation = FALSE)
        }
      }
    }
    return(S_T_ests)
  }

  S_T_preds <- t(apply(X = as.matrix(seq(1, nrow(newX))),
                       MARGIN = 1,
                       FUN = estimate_S_T))

  if (direction == "retrospective"){
    S_T_preds <- 1 - S_T_preds
  }

  res <- list(S_T_preds = S_T_preds)
  class(res) <- "survMLc"
  return(res)
}
