#' Estimate a conditional survival function by binary regression on a grid
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
#' @param stack Logical, indicating whether to fit a single bin
#' @param V CV fold number
#' @param entry Variable indicating time of entry into the study (truncation variable)
#'
#' @return An object of class \code{gridSurv}
#'
#' @export
#'
#' @examples
gridSurv <- function(time,
                     event,
                     X,
                     newX,
                     newtimes,
                     time_grid_approx,
                     bin_size = NULL,
                     stack = TRUE,
                     denom_method = "stratified",
                     fast = FALSE,
                     V = 10,
                     entry = NULL){

  if (!fast){ # do xgboost if speed not a concern
    # determine optimal models
    P_Delta_opt <- p_delta_xgboost(event = event,
                                   X = X,
                                   V = V)
    S_Y_1_opt <- f_y_stackCVcdf(time = time,
                                event = event,
                                X = X,
                                censored = FALSE,
                                bin_size = bin_size,
                                V = V)

    if (denom_method == "stratified"){
      S_Y_0_opt <- f_y_stackCVcdf(time = time,
                                  event = event,
                                  X = X,
                                  censored = TRUE,
                                  bin_size = bin_size,
                                  V = V)
      S_Y_0_opt_preds <- predict(S_Y_0_opt,
                                 newX = newX,
                                 newtimes = time_grid_approx)
    } else{
      S_Y_opt <- f_y_stackCVcdf(time = time,
                                event = event,
                                X = X,
                                censored = NULL,
                                bin_size = bin_size,
                                V = V)
      S_Y_opt_preds <- predict(S_Y_opt,
                               newX = newX,
                               newtimes = time_grid_approx)
    }
    if (!is.null(entry)){ # if a truncation variable is given
      if (denom_method == "stratified"){
        F_W_1_opt <- f_w_stackCVcdf(time = time,
                                    event = event,
                                    X = X,
                                    censored = FALSE,
                                    bin_size = bin_size,
                                    V = V,
                                    entry = entry)
        F_W_1_opt_preds <- predict(F_W_1_opt,
                                   newX = newX,
                                   newtimes = newtimes)
        F_W_0_opt <- f_w_stackCVcdf(time = time,
                                    event = event,
                                    X = X,
                                    censored = TRUE,
                                    bin_size = bin_size,
                                    V = V,
                                    entry = entry)
        F_W_0_opt_preds <- predict(F_W_0_opt,
                                   newX = newX,
                                   newtimes = newtimes)


      } else{
        F_W_opt <-f_w_stackCVcdf(time = time,
                                 event = event,
                                 X = X,
                                 censored = NULL,
                                 bin_size = bin_size,
                                 V = V,
                                 entry = entry)
        F_W_opt_preds <- predict(F_W_opt,
                                   newX = newX,
                                   newtimes = newtimes)
      }
    }
  } else{ # if speed is a concern, use ranger
    # determine optimal models
    P_Delta_opt <- p_delta_ranger(event = event,
                                   X = X,
                                   V = V)
    S_Y_1_opt <- f_y_stackCVranger(time = time,
                                event = event,
                                X = X,
                                censored = FALSE,
                                bin_size = bin_size,
                                V = V)

    if (denom_method == "stratified"){
      S_Y_0_opt <- f_y_stackCVranger(time = time,
                                  event = event,
                                  X = X,
                                  censored = TRUE,
                                  bin_size = bin_size,
                                  V = V)
      S_Y_0_opt_preds <- predict(S_Y_0_opt,
                                 newX = newX,
                                 newtimes = time_grid_approx)
    } else{
      S_Y_opt <- f_y_stackCVranger(time = time,
                                event = event,
                                X = X,
                                censored = NULL,
                                bin_size = bin_size,
                                V = V)
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
                                    entry = entry)
        F_W_1_opt_preds <- predict(F_W_1_opt,
                                   newX = newX,
                                   newtimes = newtimes)
        F_W_0_opt <- f_w_stackCVranger(time = time,
                                    event = event,
                                    X = X,
                                    censored = TRUE,
                                    bin_size = bin_size,
                                    V = V,
                                    entry = entry)
        F_W_0_opt_preds <- predict(F_W_0_opt,
                                   newX = newX,
                                   newtimes = newtimes)
      } else{
        F_W_opt <-f_w_stackCVranger(time = time,
                                 event = event,
                                 X = X,
                                 censored = NULL,
                                 bin_size = bin_size,
                                 V = V,
                                 entry = entry)
        F_W_opt_preds <- predict(F_W_opt,
                                   newX = newX,
                                   newtimes = newtimes)
      }
    }
  }

  # fit optimal models
  P_Delta_opt_preds <- predict(P_Delta_opt, newX = newX) # this is for my wrapped algorithms
  #P_Delta_opt_preds <- predict(P_Delta_opt, newdata = newX)$pred # this is for SuperLearner

  S_Y_1_opt_preds <- predict(S_Y_1_opt,
                             newX = newX,
                             newtimes = time_grid_approx) # this was previously newtimes, which was possibly causing a HUGE issue

  estimate_S_T <- function(i){
    # get S_Y estimates up to t
    S_Y_1_curr <- S_Y_1_opt_preds[i,]

    pi_curr <- P_Delta_opt_preds[i]

    if (denom_method == "stratified"){
      S_Y_0_curr <- S_Y_0_opt_preds[i,]
      if (!is.null(entry)){ # truncation
        F_W_0_curr <- F_W_0_opt_preds[i,]
        F_W_1_curr <- F_W_1_opt_preds[i,]
        # S_T_ests <-compute_prodint(cdf_uncens = S_Y_1_curr,
        #                            cdf_cens = S_Y_0_curr,
        #                            #cdf_marg = S_Y_curr,
        #                            entry_uncens = F_W_1_curr,
        #                            entry_cens = F_W_0_curr,
        #                            p_uncens = pi_curr,
        #                            newtimes = newtimes,
        #                            time_grid = time_grid_approx,
        #                            denom_method = denom_method,
        #                            truncation = TRUE)
        S_T_ests <-compute_exponential(cdf_uncens = S_Y_1_curr,
                                   cdf_cens = S_Y_0_curr,
                                   #cdf_marg = S_Y_curr,
                                   entry_uncens = F_W_1_curr,
                                   entry_cens = F_W_0_curr,
                                   p_uncens = pi_curr,
                                   newtimes = newtimes,
                                   time_grid = time_grid_approx,
                                   denom_method = denom_method,
                                   truncation = TRUE)
      } else{ # no truncation
        S_T_ests <- compute_prodint(cdf_uncens = S_Y_1_curr,
                                    cdf_cens = S_Y_0_curr,
                                    #cdf_marg = S_Y_curr,
                                    p_uncens = pi_curr,
                                    newtimes = newtimes,
                                    time_grid = time_grid_approx,
                                    denom_method = denom_method,
                                    truncation = FALSE)
      }

    } else{ # marginal denominator
      S_Y_curr <- S_Y_opt_preds[i,]
      if (!is.null(entry)){
        F_W_curr <- F_W_opt_preds[i,]
        # S_T_ests <-compute_prodint(cdf_uncens = S_Y_1_curr,
        #                            #cdf_cens = S_Y_0_curr,
        #                            cdf_marg = S_Y_curr,
        #                            p_uncens = pi_curr,
        #                            entry_marg = F_W_curr,
        #                            newtimes = newtimes,
        #                            time_grid = time_grid_approx,
        #                            denom_method = denom_method,
        #                            truncation = TRUE)
        S_T_ests <-compute_exponential(cdf_uncens = S_Y_1_curr,
                                   #cdf_cens = S_Y_0_curr,
                                   cdf_marg = S_Y_curr,
                                   p_uncens = pi_curr,
                                   entry_marg = F_W_curr,
                                   newtimes = newtimes,
                                   time_grid = time_grid_approx,
                                   denom_method = denom_method,
                                   truncation = TRUE)
      } else{
        S_T_ests <- compute_prodint(cdf_uncens = S_Y_1_curr,
                                    #cdf_cens = S_Y_0_curr,
                                    cdf_marg = S_Y_curr,
                                    p_uncens = pi_curr,
                                    newtimes = newtimes,
                                    time_grid = time_grid_approx,
                                    denom_method = denom_method,
                                    truncation = FALSE)
      }

    }

    return(S_T_ests)
  }

  S_T_preds <- t(apply(X = as.matrix(seq(1, nrow(newX))),
                       MARGIN = 1,
                       FUN = estimate_S_T))
#
#   estimate_S_C <- function(i){
#     # get S_Y estimates up to t
#     S_Y_1_curr <- S_Y_1_opt_preds[i,]
#     S_Y_0_curr <- S_Y_0_opt_preds[i,]
#     #S_Y_curr <- S_Y_opt_preds[i,]
#     pi_curr <- P_Delta_opt_preds[i]
#
#
#
#
#     # switch roles of T and C for estimating S_C
#     S_T_ests <- compute_prodint(cdf_uncens = S_Y_0_curr,
#                                 cdf_cens = S_Y_1_curr,
#                                 #cdf_marg = S_Y_curr,
#                                 p_uncens = 1-pi_curr,
#                                 newtimes = newtimes,
#                                 time_grid = time_grid_approx,
#                                 denom_method = denom_method)
#
#     return(S_T_ests)
#   }
#
#   S_C_preds <- t(apply(X = as.matrix(seq(1, nrow(newX))),
#                        MARGIN = 1,
#                        FUN = estimate_S_C))

  res <- list(S_T_preds = S_T_preds,
              S_C_preds = NA)#,
              #P_Delta_preds = P_Delta_opt_preds,
              #F_Y_1_preds = S_Y_1_opt_preds,
              #F_Y_0_preds = S_Y_0_opt_preds,
              #F_W_1_preds = F_W_1_opt_preds,
              #F_W_0_preds = F_W_0_opt_preds,
              #P_Delta_algo = P_Delta_opt,
              #F_Y_1_algo = S_Y_1_opt,
              #F_Y_0_algo = S_Y_0_opt,
              #F_W_1_algo = F_W_1_opt,
              #F_W_0_algo = F_W_0_opt)
  class(res) <- "gridSurv"
  return(res)
}
