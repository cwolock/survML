#' Estimate a conditional survival function via local survival stacking
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
#' @param bin_size Size of bins for the discretization of time.
#' A value between 0 and 1 indicating the size of observed event time quantiles
#' on which to grid times (e.g. 0.02 creates a grid of 50 times evenly spaced on the
#' quantile scaled). If NULL, defaults to every observed event time.
#' @param time_basis How to treat time for training the binary
#' classifier. Options are \code{"continuous"} and \code{"dummy"}, meaning
#' an indicator variable is included for each time in the time grid.
#' @param SL_control Named list of parameters controlling the Super Learner fitting
#' process. These parameters are passed directly to the \code{SuperLearner} function.
#' Parameters include \code{SL.library} (library of algorithms to include in the
#' binary classification Super Learner), \code{V} (Number of cross validation folds on
#' which to train the Super Learner classifier, defaults to 10), \code{method} (Method for
#' estimating coefficients for the Super Learner, defaults to \code{"method.NNLS"}), and \code{obsWeights}
#' (observation weights, passed directly to prediction algorithms by \code{SuperLearner}).
#' @param tau The maximum time of interest in a study, used for
#' retrospective conditional survival estimation. Rather than dealing
#' with right truncation separately than left truncation, it is simpler to
#' estimate the survival function of \code{tau - time}. Defaults to code{NULL},
#' in which case the maximum study entry time is chosen as the
#' reference point.
#'
#' @return A named list of class \code{stackL}.
#' \item{S_T_preds}{An \code{m x k} matrix of estimated event time survival probabilities at the
#' \code{m} covariate vector values and \code{k} times provided by the user in
#' \code{newX} and \code{newtimes}, respectively.}
#' \item{fit}{The Super Learner fit for binary classification on the stacked
#' dataset.}
#'
#' @export
#'
#' @examples
#'
#' # This is a small simulation example
#' set.seed(123)
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
#' fit <- stackL(time = time,
#'                event = event,
#'                entry = entry,
#'                X = X,
#'                newX = X,
#'                newtimes = seq(0, 15, .1),
#'                direction = "prospective",
#'                bin_size = 0.02,
#'                time_basis = "continuous",
#'                SL_control = list(SL.library = SL.library,
#'                                  V = 5))
#'
#' plot(fit$S_T_preds[1,], S0(t =  seq(0, 15, .1), X[1,]))
#' abline(0,1,col='red')
stackL <- function(time,
                   event = rep(1, length(time)),
                   entry = NULL,
                   X,
                   newX,
                   newtimes,
                   direction = "prospective",
                   bin_size = NULL,
                   time_basis = "continuous",
                   SL_control = list(SL.library = c("SL.mean"),
                                     V = 10,
                                     method = "method.NNLS"),
                   tau = NULL){

  if (is.null(newX)){
    newX <- X
  }

  if (is.null(newtimes)){
    newtimes <- time
  }

  if (direction == "retrospective"){
    if (is.null(tau)){
      tau <- max(entry)
    }
    time <- tau - time
    newtimes <- tau - newtimes
    entry <- tau - entry
    event <- rep(1, length(time))
  }

  X <- as.matrix(X)
  time <- as.matrix(time)
  event <- as.matrix(event)
  dat <- data.frame(X, time, event)

  # if user gives bin size, set time grid based on quantiles. otherwise, every observed time
  if (!is.null(bin_size)){
    time_grid <- sort(unique(stats::quantile(dat$time[dat$event == 1], probs = seq(0, 1, by = bin_size))))
    time_grid[1] <- 0
  } else{
    time_grid <- sort(unique(dat$time[dat$event == 1]))
    time_grid <- c(0, time_grid)
  }

  # this truncated time grid does not include the first time, since our discretization
  # convention pushes events with t= < time < t + 1 to time t
  trunc_time_grid <- time_grid#[-length(trunc_time_grid)]

  if (!is.null(SL_control$obsWeights)){
    stackX <- as.matrix(data.frame(X, obsWeights = SL_control$obsWeights))
  } else{
    stackX <- X
  }

  # create stacked dataset
  stacked <- survML:::stack_haz(time = time,
                                event = event,
                                X = stackX,
                                time_grid = time_grid,
                                entry = entry,
                                time_basis = "continuous")

  # change t to dummy variable
  if (time_basis == "dummy"){
    stacked$t <- factor(stacked$t)
    dummy_mat <- model.matrix(~-1 + t, data=stacked)
    risk_set_names <- paste0("risk_set_", seq(1, (length(trunc_time_grid))))
    colnames(dummy_mat) <- risk_set_names
    stacked$t <- NULL
    stacked <- cbind(dummy_mat, stacked)
  }

  long_obsWeights <- stacked$obsWeights
  stacked$obsWeights <- NULL
  .Y <- stacked[,ncol(stacked)]
  .X <- data.frame(stacked[,-ncol(stacked)])
  # fit Super Learner
  if (is.null(SL_control$method)){
    SL_control$method <- "method.NNLS"
  }
  if (is.null(SL_control$V)){
    SL_control$V <- 10
  }
  if (is.null(SL_control$SL.library)){
    SL_control$SL.library <- c("SL.mean")
  }
  fit <- SuperLearner::SuperLearner(Y = .Y,
                                    X = .X,
                                    SL.library = SL_control$SL.library,
                                    family = stats::binomial(),
                                    method = SL_control$method,
                                    verbose = FALSE,
                                    obsWeights = long_obsWeights,
                                    cvControl = list(V = SL_control$V))

  # create function to get discrete hazard predictions
  if (time_basis == "continuous"){
    get_hazard_preds <- function(index){
      new_stacked <- data.frame(t = trunc_time_grid[index], newX)
      preds <- stats::predict(fit, newdata=new_stacked)$pred
      return(preds)
    }
  } else if (time_basis == "dummy"){
    get_hazard_preds <- function(index){
      dummies <- matrix(0, ncol = length(trunc_time_grid), nrow = nrow(newX))
      dummies[,index] <- 1
      new_stacked <- cbind(dummies, newX)
      risk_set_names <- paste0("risk_set_", seq(1, (length(trunc_time_grid))))
      colnames(new_stacked)[1:length(trunc_time_grid)] <- risk_set_names
      new_stacked <- data.frame(new_stacked)
      preds <- stats::predict(fit, newdata=new_stacked)$pred
      return(preds)
    }
  }

  # don't estimate hazard at t =0
  #hazard_preds <- apply(X = matrix(time_grid), FUN = get_hazard_preds, MARGIN = 1)
  hazard_preds <- apply(X = matrix(1:length(trunc_time_grid)),
                        FUN = get_hazard_preds,
                        MARGIN = 1)

  get_surv_preds <- function(t){
    if (sum(trunc_time_grid <= t) != 0){ # if you don't fall before the first time in the grid
      final_index <- max(which(trunc_time_grid <= t))
      haz <- as.matrix(hazard_preds[,1:final_index])
      anti_haz <- 1 - haz
      surv <- apply(anti_haz, MARGIN = 1, prod)
    } else{
      surv <- rep(1, nrow(hazard_preds))
    }
    return(surv)
  }

  surv_preds <- apply(X = matrix(newtimes), FUN = get_surv_preds, MARGIN = 1)

  if (direction == "retrospective"){
    surv_preds <- 1 - surv_preds
  }

  res <- list(S_T_preds = surv_preds,
              direction = direction,
              time_basis = time_basis,
              time_grid = time_grid,
              tau = tau,
              fit = fit)
  class(res) <- "stackL"
  return(res)
}

#' @noRd
#' @export

predict.stackL <- function(object,
                           newX,
                           newtimes){

  trunc_time_grid <- object$time_grid

  # create function to get discrete hazard predictions
  if (object$time_basis == "continuous"){
    get_hazard_preds <- function(index){
      new_stacked <- data.frame(t = trunc_time_grid[index], newX)
      preds <- stats::predict(object$fit, newdata=new_stacked)$pred
      return(preds)
    }
  } else if (object$time_basis == "dummy"){
    get_hazard_preds <- function(index){
      dummies <- matrix(0, ncol = length(trunc_time_grid), nrow = nrow(newX))
      dummies[,index] <- 1
      new_stacked <- cbind(dummies, newX)
      risk_set_names <- paste0("risk_set_", seq(1, (length(trunc_time_grid))))
      colnames(new_stacked)[1:length(trunc_time_grid)] <- risk_set_names
      new_stacked <- data.frame(new_stacked)
      preds <- stats::predict(object$fit, newdata=new_stacked)$pred
      return(preds)
    }
  }

  # don't estimate hazard at t =0
  #hazard_preds <- apply(X = matrix(time_grid), FUN = get_hazard_preds, MARGIN = 1)
  hazard_preds <- apply(X = matrix(1:length(trunc_time_grid)),
                        FUN = get_hazard_preds,
                        MARGIN = 1)

  get_surv_preds <- function(t){
    if (sum(trunc_time_grid <= t) != 0){ # if you don't fall before the first time in the grid
      final_index <- max(which(trunc_time_grid <= t))
      haz <- as.matrix(hazard_preds[,1:final_index])
      anti_haz <- 1 - haz
      surv <- apply(anti_haz, MARGIN = 1, prod)
    } else{
      surv <- rep(1, nrow(hazard_preds))
    }
    return(surv)
  }

  surv_preds <- apply(X = matrix(newtimes), FUN = get_surv_preds, MARGIN = 1)


  if (object$direction == "retrospective"){
    surv_preds <- 1 - surv_preds
  }

  return(list(S_T_preds = surv_preds))
}
