#' Estimate a conditional survival function by binary regression on a grid
#'
#' @param time Observed time
#' @param event Indicator of event (vs censoring)
#' @param X Covariate matrix
#' @param newX Values of covariates at which to make a prediction
#' @param newtimes Times at which to make the survival function prediction
#' @param time_grid_approx Grid of time points on which to approximate PL integral
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
                     time_grid_approx,
                     bin_size = NULL,
                     denom_method = "stratified",
                     algorithm = "xgboost",
                     V = 10,
                     entry = NULL,
                     time_basis,
                     surv_form = "PI",
                     tuning_params = NULL,
                     direction = "prospective"){

  if (direction == "prospective"){
    if (algorithm == "xgboost"){ # do xgboost if speed not a concern
      # determine optimal models
      P_Delta_opt <- p_delta_xgboost(event = event,
                                     X = X,
                                     V = V,
                                     tuning_params = tuning_params)
      S_Y_1_opt <- f_y_stack_xgboost(time = time,
                                     event = event,
                                     X = X,
                                     censored = FALSE,
                                     bin_size = bin_size,
                                     V = V,
                                     time_basis = time_basis,
                                     tuning_params = tuning_params)

      if (denom_method == "stratified"){
        S_Y_0_opt <- f_y_stack_xgboost(time = time,
                                       event = event,
                                       X = X,
                                       censored = TRUE,
                                       bin_size = bin_size,
                                       V = V,
                                       time_basis = time_basis,
                                       tuning_params = tuning_params)
        S_Y_0_opt_preds <- predict(S_Y_0_opt,
                                   newX = newX,
                                   newtimes = time_grid_approx)
      } else{
        S_Y_opt <- f_y_stack_xgboost(time = time,
                                     event = event,
                                     X = X,
                                     censored = NULL,
                                     bin_size = bin_size,
                                     V = V,
                                     time_basis = time_basis,
                                     tuning_params = tuning_params)
        S_Y_opt_preds <- predict(S_Y_opt,
                                 newX = newX,
                                 newtimes = time_grid_approx)
      }
      if (!is.null(entry)){ # if a truncation variable is given
        if (denom_method == "stratified"){
          F_W_1_opt <- f_w_stack_xgboost(time = time,
                                         event = event,
                                         X = X,
                                         censored = FALSE,
                                         bin_size = bin_size,
                                         V = V,
                                         entry = entry,
                                         time_basis = time_basis,
                                         tuning_params = tuning_params)
          F_W_1_opt_preds <- predict(F_W_1_opt,
                                     newX = newX,
                                     newtimes = time_grid_approx)
          F_W_0_opt <- f_w_stack_xgboost(time = time,
                                         event = event,
                                         X = X,
                                         censored = TRUE,
                                         bin_size = bin_size,
                                         V = V,
                                         entry = entry,
                                         time_basis = time_basis,
                                         tuning_params = tuning_params)
          F_W_0_opt_preds <- predict(F_W_0_opt,
                                     newX = newX,
                                     newtimes = time_grid_approx)


        } else{
          F_W_opt <-f_w_stack_xgboost(time = time,
                                      event = event,
                                      X = X,
                                      censored = NULL,
                                      bin_size = bin_size,
                                      V = V,
                                      entry = entry,
                                      time_basis = time_basis,
                                      tuning_params = tuning_params)
          F_W_opt_preds <- predict(F_W_opt,
                                   newX = newX,
                                   newtimes = time_grid_approx)
        }
      }
    } else if (algorithm == "ranger"){ # if speed is a concern, use ranger
      # determine optimal models
      P_Delta_opt <- p_delta_ranger(event = event,
                                    X = X,
                                    V = V)
      S_Y_1_opt <- f_y_stack_ranger(time = time,
                                    event = event,
                                    X = X,
                                    censored = FALSE,
                                    bin_size = bin_size,
                                    V = V,
                                    time_basis = time_basis)

      if (denom_method == "stratified"){
        S_Y_0_opt <- f_y_stack_ranger(time = time,
                                      event = event,
                                      X = X,
                                      censored = TRUE,
                                      bin_size = bin_size,
                                      V = V,
                                      time_basis = time_basis)
        S_Y_0_opt_preds <- predict(S_Y_0_opt,
                                   newX = newX,
                                   newtimes = time_grid_approx)
      } else{
        S_Y_opt <- f_y_stack_ranger(time = time,
                                    event = event,
                                    X = X,
                                    censored = NULL,
                                    bin_size = bin_size,
                                    V = V,
                                    time_basis = time_basis)
        S_Y_opt_preds <- predict(S_Y_opt,
                                 newX = newX,
                                 newtimes = time_grid_approx)
      }
      if (!is.null(entry)){ # if a truncation variable is given
        if (denom_method == "stratified"){
          F_W_1_opt <- f_w_stackCVranger(time = time,
                                         event = event,
                                         X = X,
                                         censored = FALSE,
                                         bin_size = bin_size,
                                         V = V,
                                         entry = entry,
                                         time_basis = time_basis)
          F_W_1_opt_preds <- predict(F_W_1_opt,
                                     newX = newX,
                                     newtimes = time_grid_approx)
          F_W_0_opt <- f_w_stackCVranger(time = time,
                                         event = event,
                                         X = X,
                                         censored = TRUE,
                                         bin_size = bin_size,
                                         V = V,
                                         entry = entry,
                                         time_basis = time_basis)
          F_W_0_opt_preds <- predict(F_W_0_opt,
                                     newX = newX,
                                     newtimes = time_grid_approx)
        } else{
          F_W_opt <-f_w_stackCVranger(time = time,
                                      event = event,
                                      X = X,
                                      censored = NULL,
                                      bin_size = bin_size,
                                      V = V,
                                      entry = entry,
                                      time_basis = time_basis)
          F_W_opt_preds <- predict(F_W_opt,
                                   newX = newX,
                                   newtimes = time_grid_approx)
        }
      }
    } else if (algorithm == "gam"){
      # determine optimal models
      SL.library <- c("SL.mean", "SL.glm", "SL.gam", "SL.earth", "SL.ranger")
      P_Delta_opt <- p_delta_SuperLearner(event = event,
                                          X = X,
                                          SL.library = SL.library)
      S_Y_1_opt <- f_y_stack_gam(time = time,
                                 event = event,
                                 X = X,
                                 censored = FALSE,
                                 bin_size = bin_size,
                                 time_basis = time_basis)

      if (denom_method == "stratified"){
        S_Y_0_opt <- f_y_stack_gam(time = time,
                                   event = event,
                                   X = X,
                                   censored = TRUE,
                                   bin_size = bin_size,
                                   time_basis = time_basis)
        S_Y_0_opt_preds <- predict(S_Y_0_opt,
                                   newX = newX,
                                   newtimes = time_grid_approx)
      } else{
        S_Y_opt <- f_y_stack_gam(time = time,
                                 event = event,
                                 X = X,
                                 censored = NULL,
                                 bin_size = bin_size,
                                 time_basis = time_basis)
        S_Y_opt_preds <- predict(S_Y_opt,
                                 newX = newX,
                                 newtimes = time_grid_approx)
      }
      if (!is.null(entry)){ # if a truncation variable is given
        if (denom_method == "stratified"){
          F_W_1_opt <- f_w_stack_gam(time = time,
                                     event = event,
                                     X = X,
                                     censored = FALSE,
                                     bin_size = bin_size,
                                     entry = entry,
                                     time_basis = time_basis)
          F_W_1_opt_preds <- predict(F_W_1_opt,
                                     newX = newX,
                                     newtimes = time_grid_approx)
          F_W_0_opt <- f_w_stack_gam(time = time,
                                     event = event,
                                     X = X,
                                     censored = TRUE,
                                     bin_size = bin_size,
                                     entry = entry,
                                     time_basis = time_basis)
          F_W_0_opt_preds <- predict(F_W_0_opt,
                                     newX = newX,
                                     newtimes = time_grid_approx)


        } else{
          F_W_opt <-f_w_stack_gam(time = time,
                                  event = event,
                                  X = X,
                                  censored = NULL,
                                  bin_size = bin_size,
                                  entry = entry,
                                  time_basis = time_basis)
          F_W_opt_preds <- predict(F_W_opt,
                                   newX = newX,
                                   newtimes = time_grid_approx)
        }
      }
    } else if (algorithm == "earth"){
      # determine optimal models
      SL.library <- c("SL.mean", "SL.glm", "SL.gam", "SL.earth", "SL.ranger")
      P_Delta_opt <- p_delta_SuperLearner(event = event,
                                          X = X,
                                          SL.library = SL.library)
      S_Y_1_opt <- f_y_stack_earth(time = time,
                                   event = event,
                                   X = X,
                                   censored = FALSE,
                                   bin_size = bin_size,
                                   time_basis = time_basis)

      if (denom_method == "stratified"){
        S_Y_0_opt <- f_y_stack_earth(time = time,
                                     event = event,
                                     X = X,
                                     censored = TRUE,
                                     bin_size = bin_size,
                                     time_basis = time_basis)
        S_Y_0_opt_preds <- predict(S_Y_0_opt,
                                   newX = newX,
                                   newtimes = time_grid_approx)
      } else{
        S_Y_opt <- f_y_stack_earth(time = time,
                                   event = event,
                                   X = X,
                                   censored = NULL,
                                   bin_size = bin_size,
                                   time_basis = time_basis)
        S_Y_opt_preds <- predict(S_Y_opt,
                                 newX = newX,
                                 newtimes = time_grid_approx)
      }
      if (!is.null(entry)){ # if a truncation variable is given
        if (denom_method == "stratified"){
          F_W_1_opt <- f_w_stack_earth(time = time,
                                       event = event,
                                       X = X,
                                       censored = FALSE,
                                       bin_size = bin_size,
                                       entry = entry,
                                       time_basis = time_basis)
          F_W_1_opt_preds <- predict(F_W_1_opt,
                                     newX = newX,
                                     newtimes = time_grid_approx)
          F_W_0_opt <- f_w_stack_earth(time = time,
                                       event = event,
                                       X = X,
                                       censored = TRUE,
                                       bin_size = bin_size,
                                       entry = entry,
                                       time_basis = time_basis)
          F_W_0_opt_preds <- predict(F_W_0_opt,
                                     newX = newX,
                                     newtimes = time_grid_approx)


        } else{
          F_W_opt <-f_w_stack_earth(time = time,
                                    event = event,
                                    X = X,
                                    censored = NULL,
                                    bin_size = bin_size,
                                    entry = entry,
                                    time_basis = time_basis)
          F_W_opt_preds <- predict(F_W_opt,
                                   newX = newX,
                                   newtimes = time_grid_approx)
        }
      }
    }

    # fit optimal models
    P_Delta_opt_preds <- predict(P_Delta_opt, newX = newX) # this is for my wrapped algorithms

    S_Y_1_opt_preds <- predict(S_Y_1_opt,
                               newX = newX,
                               newtimes = time_grid_approx)

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
      } else{ # marginal denominator
        S_Y_curr <- S_Y_opt_preds[i,]
        if (!is.null(entry)){ # trunation
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
  } else if (direction == "retrospective"){
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

  }

  S_T_preds <- t(apply(X = as.matrix(seq(1, nrow(newX))),
                       MARGIN = 1,
                       FUN = estimate_S_T))

  res <- list(S_T_preds = S_T_preds)
  class(res) <- "survMLc"
  return(res)
}
