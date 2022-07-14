#' Estimate a conditional survival function using cumulative probability estimator
#'
#' @param time \code{n x 1} numeric vector of observed
#' follow-up times If there is censoring, these are the minimum of the
#' event and censoring times.
#' @param event \code{n x 1} numeric vector of status indicators of
#' whether an event was observed. Defaults to a vector of 1s, i.e. no censoring.
#' @param entry Study entry variable, if applicable. Defaults to \code{NULL},
#' indicating that there is no truncation.
#' @param X \code{n x p} data.frame of observed covariate values
#' on which to train the estimator.
#' @param newX \code{m x p} data.frame of new observed covariate
#' values at which to obtain \code{m} predictions for the estimated algorithm.
#' Must have the same names and structure as \code{X}.
#' @param newtimes \code{k x 1} numeric vector of times at which to obtain \code{k}
#' predicted conditional survivals.
#' @param direction Whether the data come from a prospective or retrospective study.
#' This determines whether the data are treated as subject to left truncation and
#' right censoring (\code{"prospective"}) or right truncation alone
#' (\code{"retrospective"}).
#' @param bin_size Size of time bin on which to discretize for estimation
#' of cumulative probability functions. Can be a number between 0 and 1,
#' indicating the size of quantile grid (e.g. \code{0.1} estimates
#' the cumulative probability functions on a grid based on deciles of
#' observed \code{time}s). If \code{NULL}, creates a grid of
#' all observed \code{time}s.
#' @param time_basis How to treat time for training the binary
#' classifier. Options are \code{"continuous"} and \code{"dummy"}, meaning
#' an indicator variable is included for each time in the time grid.
#' @param time_grid_approx Numeric vector of times at which to
#' approximate product integral or cumulative hazard interval.
#' Defaults to \code{times} argument.
#' @param denom_method Denominator form for the hazard identification. Can
#' be either \code{"stratified"} (estimate the at-risk probability
#' within strata of the \code{event} variable) or \code{"marginal"} (estimate
#' the at-risk probability using all observations).
#' @param surv_form Mapping from hazard estimate to survival estimate.
#' Can be either \code{"PI"} (product integral mapping) or \code{"exp"}
#' (exponentiated cumulative hazard estimate).
#' @param SL.library Library of algorithms to include in the binary classification
#' Super Learner. Should have the same structure as the \code{SL.library}
#' argument to the \code{SuperLearner} function in the \code{SuperLearner} package.
#' @param V Number of cross validation folds on which to train the Super Learner
#' classifier. Defaults to 10.
#' @param tau The maximum time of interest in a study, used for
#' retrospective conditional survival estimation. Rather than dealing
#' with right truncation separately than left truncation, it is simpler to
#' estimate the survival function of \code{tau - time}. Defaults to code{NULL},
#' in which case the maximum study entry time is chosen as the
#' reference point.
#' @param obsWeights Optional observation weights. These weights are passed
#' directly to \code{SuperLearner}, which in turn passes them directly to the
#' prediction algorithms.
#'
#' @return A named list of class \code{survMLc}
#'
#' @export
#'
#' @examples
#'
#' # This is a small simulation example
#' set.seed(92)
#' n <- 100
#' X <- data.frame(X1 = rnorm(n), X2 = rbinom(n, size = 1, prob = 0.5))
#'
#' S0 <- function(t, x){
#'   pexp(t, rate = exp(-2 + x[,1] - x[,2] + .5 * x[,1] * x[,2]), lower.tail = FALSE)
#' }
#' T <- rexp(n, rate = exp(-2 + X[,1] - X[,2] + .5 *  X[,1] * X[,2]))
#'
#' G0 <- function(t, x) {
#'   as.numeric(t < 15) *.9*pexp(t,
#'                               rate = exp(-2 -.5*x[,1]-.25*x[,2]+.5*x[,1]*x[,2]),
#'                               lower.tail=FALSE)
#' }
#' C <- rexp(n, exp(-2 -.5 * X[,1] - .25 * X[,2] + .5 * X[,1] * X[,2]))
#' C[C > 15] <- 15
#'
#' entry <- runif(n, 0, 15)
#'
#' time <- pmin(T, C)
#' event <- as.numeric(T <= C)
#'
#' sampled <- which(time >= entry)
#' X <- X[sampled,]
#' time <- time[sampled]
#' event <- event[sampled]
#' entry <- entry[sampled]
#'
#' SL.library <- c("SL.mean", "SL.glm", "SL.gam", "SL.earth")
#'
#' fit <- survMLc(time = time,
#'                event = event,
#'                entry = entry,
#'                X = X,
#'                newX = X,
#'                newtimes = seq(0, 15, .1),
#'                direction = "prospective",
#'                bin_size = 0.02,
#'                time_basis = "continuous",
#'                time_grid_approx = sort(unique(time)),
#'                denom_method = "stratified",
#'                surv_form = "exp",
#'                SL.library = SL.library,
#'                V = 5)
#'
#' plot(fit$S_T_preds[1,], S0(t =  seq(0, 15, .1), X[1,]))
#' abline(0,1,col='red')
survMLc <- function(time,
                    event = rep(1, length(time)),
                    entry = NULL,
                    X,
                    newX,
                    newtimes,
                    direction = "prospective",
                    bin_size = NULL,
                    time_basis,
                    time_grid_approx = sort(unique(time)),
                    denom_method = "stratified",
                    surv_form = "PI",
                    SL.library,
                    V = 10,
                    obsWeights = NULL,
                    tau = NULL){
  P_Delta_opt <- NULL
  S_Y_opt <- NULL
  S_Y_1_opt <- NULL
  S_Y_0_opt <- NULL
  F_W_1_opt <- NULL
  F_W_0_opt <- NULL
  F_W_opt <- NULL

  if (direction == "retrospective"){
    if (is.null(tau)){
      tau <- max(entry)
    }
    time <- tau - time
    entry <- tau - entry
    event <- rep(1, length(time))
    newtimes <- tau - newtimes
    time_grid_approx <- sort(tau - time_grid_approx)
    denom_method <- "marginal"
    P_Delta_opt_preds <- rep(1, nrow(newX))
  }

  # if there is a censoring probability to estimate, i.e. if there is censoring
  if (sum(event == 0) != 0){
    P_Delta_opt <- p_delta(event = event,
                           X = X,
                           V = V,
                           SL.library = SL.library,
                           obsWeights = obsWeights)
    P_Delta_opt_preds <- stats::predict(P_Delta_opt, newX = newX) # this is for my wrapped algorithms

    if (denom_method == "stratified"){
      S_Y_0_opt <- f_y_stack(time = time,
                             event = event,
                             X = X,
                             censored = TRUE,
                             bin_size = bin_size,
                             V = V,
                             time_basis = time_basis,
                             SL.library = SL.library,
                             obsWeights = obsWeights)
      S_Y_0_opt_preds <- stats::predict(S_Y_0_opt,
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
                           SL.library = SL.library,
                           obsWeights = obsWeights)
      S_Y_opt_preds <- stats::predict(S_Y_opt,
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
                         SL.library = SL.library,
                         obsWeights = obsWeights)

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
                             SL.library = SL.library,
                             obsWeights = obsWeights)
      F_W_1_opt_preds <- stats::predict(F_W_1_opt,
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
                             SL.library = SL.library,
                             obsWeights = obsWeights)
      F_W_0_opt_preds <- stats::predict(F_W_0_opt,
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
                          SL.library = SL.library,
                          obsWeights = obsWeights)
      F_W_opt_preds <- stats::predict(F_W_opt,
                               newX = newX,
                               newtimes = time_grid_approx)
    }
  }




  S_Y_1_opt_preds <- stats::predict(S_Y_1_opt,
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

  res <- list(S_T_preds = S_T_preds,
              fits = list(P_Delta = P_Delta_opt,
                          S_Y_1 = S_Y_1_opt,
                          S_Y_0 = S_Y_0_opt,
                          S_Y = S_Y_opt,
                          F_W_1 = F_W_1_opt,
                          F_W_0 = F_W_0_opt,
                          F_W = F_W_opt))
  class(res) <- "survMLc"
  return(res)
}
