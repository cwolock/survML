#' Estimate a conditional survival function using global survival stacking
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
#' @param surv_form Mapping from hazard estimate to survival estimate.
#' Can be either \code{"PI"} (product integral mapping) or \code{"exp"}
#' (exponentiated cumulative hazard estimate).
#' @param learner Which binary regression algorithm to use. Currently, only
#' \code{SuperLearner} is supported, but more learners will be added.
#' See below for algorithm-specific arguments.
#' @param SL_control Named list of parameters controlling the Super Learner fitting
#' process. These parameters are passed directly to the \code{SuperLearner} function.
#' Parameters include \code{SL.library} (library of algorithms to include in the
#' binary classification Super Learner), \code{V} (Number of cross validation folds on
#' which to train the Super Learner classifier, defaults to 10), \code{method} (Method for
#' estimating coefficients for the Super Learner, defaults to \code{"method.NNLS"}), code{stratifyCV}
#' (logical indicating whether to stratify by outcome in \code{SuperLearner}'s cross-validation
#' scheme), and \code{obsWeights}
#' (observation weights, passed directly to prediction algorithms by \code{SuperLearner}).
#' @param tau The maximum time of interest in a study, used for
#' retrospective conditional survival estimation. Rather than dealing
#' with right truncation separately than left truncation, it is simpler to
#' estimate the survival function of \code{tau - time}. Defaults to code{NULL},
#' in which case the maximum study entry time is chosen as the
#' reference point.
#'
#' @return A named list of class \code{stackG}, with the following components:
#' \item{S_T_preds}{An \code{m x k} matrix of estimated event time survival probabilities at the
#' \code{m} covariate vector values and \code{k} times provided by the user in
#' \code{newX} and \code{newtimes}, respectively.}
#' \item{S_C_preds}{An \code{m x k} matrix of estimated censoring time survival probabilities at the
#' \code{m} covariate vector values and \code{k} times provided by the user in
#' \code{newX} and \code{newtimes}, respectively.}
#' \item{time_grid_approx}{The approximation grid for the product integral or cumulative hazard integral,
#' (user-specified).}
#' \item{direction}{Whether the data come from a prospective or retrospective study (user-specified).}
#' \item{tau}{The maximum time of interest in a study, used for
#' retrospective conditional survival estimation (user-specified).}
#' \item{surv_form}{Exponential or product-integral form (user-specified).}
#' \item{time_basis}{Whether time is included in the regression as \code{continuous} or
#' \code{dummy} (user-specified).}
#' \item{SL_control}{Named list of parameters controlling the Super Learner fitting
#' process (user-specified).}
#' \item{fits}{A named list of fitted regression objects corresponding to the constituent regressions needed for
#' global survival stacking. Includes \code{P_Delta} (probability of event given covariates),
#' \code{F_Y_1} (conditional cdf of follow-up times given covariates among uncensored),
#' \code{F_Y_0} (conditional cdf of follow-up times given covariates among censored),
#' \code{G_W_1} (conditional distribution of entry times given covariates and follow-up time among uncensored),
#' \code{G_W_0} (conditional distribution of entry times given covariates and follow-up time among uncensored).
#' Each of these objects includes estimated coefficients from the \code{SuperLearner} fit, as well as the
#' time grid used to create the stacked dataset (where applicable).}
#'
#' @export
#'
#' @examples
#'
#' # This is a small simulation example
#' set.seed(123)
#' n <- 500
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
#' # Note that this a very small Super Learner library, for computational purposes.
#' SL.library <- c("SL.mean", "SL.glm")
#'
#' fit <- stackG(time = time,
#'               event = event,
#'               entry = entry,
#'               X = X,
#'               newX = X,
#'               newtimes = seq(0, 15, .1),
#'               direction = "prospective",
#'               bin_size = 0.1,
#'               time_basis = "continuous",
#'               time_grid_approx = sort(unique(time)),
#'               surv_form = "exp",
#'               learner = "SuperLearner",
#'               SL_control = list(SL.library = SL.library,
#'                                 V = 5))
#'
#' plot(fit$S_T_preds[1,], S0(t =  seq(0, 15, .1), X[1,]))
#' abline(0,1,col='red')
#'
#' @references Wolock C.J., Gilbert P.B., Simon N., and Carone, M. (2022).
#'   "A framework for leveraging machine learning tools to estimate personalized survival curves."
stackG <- function(time,
                   event = rep(1, length(time)),
                   entry = NULL,
                   X,
                   newX = NULL,
                   newtimes = NULL,
                   direction = "prospective",
                   time_grid_fit = NULL,
                   bin_size = NULL,
                   time_basis,
                   time_grid_approx = sort(unique(time)),
                   surv_form = "PI",
                   learner = "SuperLearner",
                   SL_control = list(SL.library = c("SL.mean"),
                                     V = 10,
                                     method = "method.NNLS",
                                     stratifyCV = FALSE),
                   tau = NULL){
  P_Delta_opt <- NULL
  S_Y_opt <- NULL
  F_Y_1_opt <- NULL
  F_Y_0_opt <- NULL
  G_W_1_opt <- NULL
  G_W_0_opt <- NULL
  F_W_opt <- NULL

  tau <- NULL

  if (is.null(newX)){
    newX <- X
  }

  if (is.null(newtimes)){
    newtimes <- time_grid_approx
  }

  if (direction == "retrospective"){
    if (is.null(tau)){
      tau <- max(entry)
    }
    time <- tau - time
    entry <- tau - entry
    event <- rep(1, length(time))
    newtimes <- tau - newtimes
    time_grid_approx <- sort(tau - time_grid_approx)
    P_Delta_opt_preds <- rep(1, nrow(newX))
    F_Y_0_opt_preds <- matrix(0, nrow = nrow(newX), ncol = length(time_grid_approx))
    G_W_0_opt_preds <- matrix(0, nrow = nrow(newX), ncol = length(time_grid_approx))
  }

  # if there is a censoring probability to estimate, i.e. if there is censoring
  if (sum(event == 0) != 0){
    P_Delta_opt <- p_delta(event = event,
                           X = X,
                           learner = learner,
                           SL_control = SL_control)
    P_Delta_opt_preds <- stats::predict(P_Delta_opt, newX = newX) # this is for my wrapped algorithms

    F_Y_0_opt <- f_y_stack(time = time,
                           event = event,
                           X = X,
                           censored = TRUE,
                           time_grid = time_grid_fit$F_Y_0_grid,
                           bin_size = bin_size,
                           learner = learner,
                           SL_control = SL_control,
                           time_basis = time_basis)
    F_Y_0_opt_preds <- stats::predict(F_Y_0_opt,
                                      newX = newX,
                                      newtimes = time_grid_approx)
  }

  F_Y_1_opt <- f_y_stack(time = time,
                         event = event,
                         X = X,
                         censored = FALSE,
                         time_grid = time_grid_fit$F_Y_1_grid,
                         bin_size = bin_size,
                         learner = learner,
                         SL_control = SL_control,
                         time_basis = time_basis)

  if (!is.null(entry)){ # if a truncation variable is given
    G_W_1_opt <- f_w_stack(time = time,
                           event = event,
                           X = X,
                           censored = FALSE,
                           time_grid = time_grid_fit$G_W_1_grid,
                           bin_size = bin_size,
                           learner = learner,
                           SL_control = SL_control,
                           entry = entry,
                           time_basis = time_basis)
    G_W_1_opt_preds <- stats::predict(G_W_1_opt,
                                      newX = newX,
                                      newtimes = time_grid_approx)
    if (sum(event == 0) != 0){ # if there's censoring
      G_W_0_opt <- f_w_stack(time = time,
                             event = event,
                             X = X,
                             censored = TRUE,
                             time_grid = time_grid_fit$G_W_0_grid,
                             bin_size = bin_size,
                             learner = learner,
                             SL_control = SL_control,
                             entry = entry,
                             time_basis = time_basis)
      G_W_0_opt_preds <- stats::predict(G_W_0_opt,
                                        newX = newX,
                                        newtimes = time_grid_approx)
    } else{
      G_W_0_opt_preds <- matrix(1, nrow = nrow(newX), ncol = length(time_grid_approx))
    }
  } else{
    G_W_0_opt_preds <- matrix(1, nrow = nrow(newX), ncol = length(time_grid_approx))
    G_W_1_opt_preds <- matrix(1, nrow = nrow(newX), ncol = length(time_grid_approx))
  }

  F_Y_1_opt_preds <- stats::predict(F_Y_1_opt,
                                    newX = newX,
                                    newtimes = time_grid_approx)

  estimate_S_T <- function(i){
    # get S_Y estimates up to t
    F_Y_1_curr <- F_Y_1_opt_preds[i,]
    pi_curr <- P_Delta_opt_preds[i]
    F_Y_0_curr <- F_Y_0_opt_preds[i,]
    G_W_0_curr <- G_W_0_opt_preds[i,]
    G_W_1_curr <- G_W_1_opt_preds[i,]
    if (surv_form == "PI"){
      S_T_ests <-compute_prodint(cdf_uncens = F_Y_1_curr,
                                 cdf_cens = F_Y_0_curr,
                                 entry_uncens = G_W_1_curr,
                                 entry_cens = G_W_0_curr,
                                 p_uncens = pi_curr,
                                 newtimes = newtimes,
                                 time_grid = time_grid_approx)
      S_C_ests <-compute_prodint(cdf_uncens = F_Y_0_curr,
                                 cdf_cens = F_Y_1_curr,
                                 entry_uncens = G_W_0_curr,
                                 entry_cens = G_W_1_curr,
                                 p_uncens = 1 - pi_curr,
                                 newtimes = newtimes,
                                 time_grid = time_grid_approx)
    } else if (surv_form == "exp"){
      S_T_ests <-compute_exponential(cdf_uncens = F_Y_1_curr,
                                     cdf_cens = F_Y_0_curr,
                                     entry_uncens = G_W_1_curr,
                                     entry_cens = G_W_0_curr,
                                     p_uncens = pi_curr,
                                     newtimes = newtimes,
                                     time_grid = time_grid_approx)
      S_C_ests <-compute_exponential(cdf_uncens = F_Y_0_curr,
                                     cdf_cens = F_Y_1_curr,
                                     entry_uncens = G_W_0_curr,
                                     entry_cens = G_W_1_curr,
                                     p_uncens = 1 - pi_curr,
                                     newtimes = newtimes,
                                     time_grid = time_grid_approx)
    }
    return(list(S_T_ests = S_T_ests, S_C_ests = S_C_ests))
  }



  preds <- t(matrix(unlist(apply(X = as.matrix(seq(1, nrow(newX))),
                                 MARGIN = 1,
                                 FUN = estimate_S_T)), nrow = 2*length(newtimes)))

  S_T_preds <- preds[,1:length(newtimes)]
  S_C_preds <- preds[,(length(newtimes) + 1):(2*length(newtimes))]

  if (direction == "retrospective"){
    S_T_preds <- 1 - S_T_preds
    S_C_preds <- NULL
  }

  res <- list(S_T_preds = S_T_preds,
              S_C_preds = S_C_preds,
              time_grid_approx = time_grid_approx,
              direction = direction,
              tau = tau,
              surv_form = surv_form,
              time_basis = time_basis,
              learner = learner,
              SL_control = SL_control,
              fits = list(P_Delta = P_Delta_opt,
                          F_Y_1 = F_Y_1_opt,
                          F_Y_0 = F_Y_0_opt,
                          G_W_1 = G_W_1_opt,
                          G_W_0 = G_W_0_opt))
  class(res) <- "stackG"
  return(res)
}


#' @noRd
#' @export
predict.stackG <- function(object,
                           newX,
                           newtimes,
                           surv_form = object$surv_form,
                           time_grid_approx = object$time_grid_approx,
                           ...){

  if (object$direction == "retrospective"){
    newtimes <- object$tau - newtimes
  }

  if (!is.null(object$fits$P_Delta)){
    P_Delta_opt_preds <- stats::predict(object$fits$P_Delta, newX = newX)
  } else{
    P_Delta_opt_preds <- rep(1, nrow(newX))
  }
  if (!is.null(object$fits$G_W_1)){
    G_W_1_opt_preds <- stats::predict(object$fits$G_W_1,
                                      newX = newX,
                                      newtimes = time_grid_approx)
  } else{
    G_W_1_opt_preds <- matrix(1, nrow = nrow(newX), ncol = length(time_grid_approx))
  }
  if (!is.null(object$fits$G_W_0)){
    G_W_0_opt_preds <- stats::predict(object$fits$G_W_0,
                                      newX = newX,
                                      newtimes = time_grid_approx)
  } else{
    G_W_0_opt_preds <- matrix(1, nrow = nrow(newX), ncol = length(time_grid_approx))
  }
  if (!is.null(object$fits$F_Y_0)){
    F_Y_0_opt_preds <- stats::predict(object$fits$F_Y_0,
                                      newX = newX,
                                      newtimes = time_grid_approx)
  } else{
    F_Y_0_opt_preds <- matrix(1, nrow = nrow(newX), ncol = length(time_grid_approx))
  }

  F_Y_1_opt_preds <- stats::predict(object$fits$F_Y_1,
                                    newX = newX,
                                    newtimes = time_grid_approx)

  estimate_S_T <- function(i){
    # get S_Y estimates up to t
    F_Y_1_curr <- F_Y_1_opt_preds[i,]
    pi_curr <- P_Delta_opt_preds[i]
    F_Y_0_curr <- F_Y_0_opt_preds[i,]
    G_W_0_curr <- G_W_0_opt_preds[i,]
    G_W_1_curr <- G_W_1_opt_preds[i,]
    if (surv_form == "PI"){
      S_T_ests <-compute_prodint(cdf_uncens = F_Y_1_curr,
                                 cdf_cens = F_Y_0_curr,
                                 entry_uncens = G_W_1_curr,
                                 entry_cens = G_W_0_curr,
                                 p_uncens = pi_curr,
                                 newtimes = newtimes,
                                 time_grid = time_grid_approx)
      S_C_ests <-compute_prodint(cdf_uncens = F_Y_0_curr,
                                 cdf_cens = F_Y_1_curr,
                                 entry_uncens = G_W_0_curr,
                                 entry_cens = G_W_1_curr,
                                 p_uncens = 1 - pi_curr,
                                 newtimes = newtimes,
                                 time_grid = time_grid_approx)
    } else if (surv_form == "exp"){
      S_T_ests <-compute_exponential(cdf_uncens = F_Y_1_curr,
                                     cdf_cens = F_Y_0_curr,
                                     entry_uncens = G_W_1_curr,
                                     entry_cens = G_W_0_curr,
                                     p_uncens = pi_curr,
                                     newtimes = newtimes,
                                     time_grid = time_grid_approx)
      S_C_ests <-compute_exponential(cdf_uncens = F_Y_0_curr,
                                     cdf_cens = F_Y_1_curr,
                                     entry_uncens = G_W_0_curr,
                                     entry_cens = G_W_1_curr,
                                     p_uncens = 1 - pi_curr,
                                     newtimes = newtimes,
                                     time_grid = time_grid_approx)
    }

    return(list(S_T_ests = S_T_ests, S_C_ests = S_C_ests))
  }

  preds <- t(matrix(unlist(apply(X = as.matrix(seq(1, nrow(newX))),
                                 MARGIN = 1,
                                 FUN = estimate_S_T)), nrow = 2*length(newtimes)))

  S_T_preds <- preds[,1:length(newtimes)]
  S_C_preds <- preds[,(length(newtimes) + 1):(2*length(newtimes))]

  if (object$direction == "retrospective"){
    S_T_preds <- 1 - S_T_preds
    S_C_preds <- NULL
  }

  res <- list(S_T_preds = S_T_preds,
              S_C_preds = S_C_preds,
              surv_form = surv_form,
              time_grid_approx = time_grid_approx)
  return(res)

}
