#' Stacked binary regression with SuperLearner, using the cdf
#'
#' @param time Observed time
#' @param event Indicator of event (vs censoring)
#' @param X Covariate matrix
#' @param censored Logical, indicates whether to run regression on censored observations (vs uncensored)
#' @param bin_size Size of quantiles over which to make the stacking bins
#' @param isotonize Logical, indicates whether or not to isotonize cdf estimates using PAVA
#' @param time_basis How to treat time
#' @param SL_control SuperLearner control parameters
#'
#' @return An object of class \code{f_y_stack_SuperLearner}
#' @noRd
f_y_stack_SuperLearner <- function(time,
                                   event,
                                   X,
                                   censored,
                                   time_grid,
                                   bin_size,
                                   isotonize = TRUE,
                                   SL_control,
                                   time_basis){

  if (!is.null(censored)){
    if (censored == TRUE){
      time <- time[!as.logical(event)]
      X <- X[!as.logical(event),]
      obsWeights <- SL_control$obsWeights[!as.logical(event)]
    } else if (censored == FALSE){
      time <- time[as.logical(event)]
      X <- X[as.logical(event),]
      obsWeights <- SL_control$obsWeights[as.logical(event)]
    }
  } else{
    time <- time
    X <- X
    obsWeights <- SL_control$obsWeights
  }

  if (is.null(time_grid)){
    bin_variable <- time

    # if user gives bin size, set time grid based on quantiles. otherwise, every observed time
    if (!is.null(bin_size)){
      time_grid <- sort(unique(stats::quantile(bin_variable, type = 1, probs = seq(0, 1, by = bin_size))))
      time_grid <- c(0, time_grid) # 013123 changed this to try to get better predictions at time 0
      #time_grid[1] <- 0 # manually set first point to 0, instead of first observed time
    } else{
      time_grid <- sort(unique(bin_variable))
      time_grid <- c(0, time_grid)
    }

  }

  cv_folds <- split(sample(1:length(time)), rep(1:SL_control$V, length = length(time)))

  X <- as.matrix(X)
  time <- as.matrix(time)
  dat <- data.frame(X, time)

  ids <- seq(1:length(time))

  if (!is.null(obsWeights)){
    stackX <- as.matrix(data.frame(X,
                                   obsWeights = obsWeights,
                                   ids = ids))
  } else{
    stackX <- as.matrix(data.frame(X,
                                   ids = ids))
  }

  stacked <- stack_cdf(time = time,
                       X = stackX,
                       time_grid = time_grid)

  # change t to dummy variable
  if (time_basis == "dummy"){
    stacked$t <- factor(stacked$t)
    dummy_mat <- stats::model.matrix(~-1 + t, data=stacked)
    risk_set_names <- paste0("risk_set_", seq(1, (length(time_grid))))
    colnames(dummy_mat) <- risk_set_names
    stacked$t <- NULL
    stacked <- cbind(dummy_mat, stacked)
  }

  long_obsWeights <- stacked$obsWeights
  stacked_ids <- stacked$ids
  stacked$obsWeights <- NULL
  stacked$ids <- NULL
  .Y <- stacked[,ncol(stacked)]
  .X <- stacked[,-ncol(stacked)]

  get_validRows <- function(fold_sample_ids){
    validRows <- which(stacked_ids %in% fold_sample_ids)
    return(validRows)
  }

  validRows <- lapply(cv_folds, get_validRows)

  if (is.null(SL_control$method)){
    SL_control$method <- "method.NNLS"
  }
  if (is.null(SL_control$V)){
    SL_control$V <- 10
  }
  if (is.null(SL_control$SL.library)){
    SL_control$SL.library <- c("SL.mean")
  }
  if (is.null(SL_control$stratifyCV)){
    SL_control$stratifyCV <- FALSE
  }

  fit <- SuperLearner::SuperLearner(Y = .Y,
                                    X = .X,
                                    SL.library = SL_control$SL.library,
                                    family = stats::binomial(),
                                    method = SL_control$method,
                                    verbose = FALSE,
                                    obsWeights = long_obsWeights,
                                    cvControl = list(V = SL_control$V,
                                                     validRows = validRows,
                                                     stratifyCV = SL_control$stratifyCV))

  fit <- list(reg.object = fit,
              time_grid = time_grid,
              isotonize = isotonize,
              time_basis = time_basis)
  class(fit) <- c("f_y_stack_SuperLearner")
  return(fit)
}

#' Prediction function for stacked SuperLearner CDF
#'
#' @param fit Fitted regression object
#' @param newX Values of covariates at which to make a prediction
#' @param newtimes
#'
#' @return Matrix of predictions
#' @noRd
predict.f_y_stack_SuperLearner <- function(fit, newX, newtimes){

  time_grid <- fit$time_grid

  if (fit$time_basis == "continuous"){
    get_stacked_pred <- function(t){
      new_stacked <- data.frame(t = t, newX)
      preds <- stats::predict(fit$reg.object, newdata=new_stacked)$pred
      return(preds)
    }

    predictions <- matrix(apply(X = matrix(newtimes), FUN = get_stacked_pred, MARGIN = 1),
                          nrow = nrow(newX))
  } else{
    get_preds <- function(t){
      dummies <- matrix(0, ncol = length(time_grid), nrow = nrow(newX))
      index <- max(which(time_grid <= t))
      dummies[,index] <- 1
      new_stacked <- cbind(dummies, newX)
      risk_set_names <- paste0("risk_set_", seq(1, (length(time_grid))))
      colnames(new_stacked)[1:length(time_grid)] <- risk_set_names
      new_stacked <- data.frame(new_stacked)
      preds <- stats::predict(fit$reg.object, newdata=new_stacked)$pred
      return(preds)
    }

    predictions <- matrix(apply(X = matrix(newtimes), FUN = get_preds, MARGIN = 1),
                          nrow = nrow(newX))
  }

  if (fit$isotonize){
    iso.cdf.ests <- t(apply(predictions, MARGIN = 1, FUN = Iso::pava))
  } else{
    iso.cdf.ests <- predictions
  }

  return(iso.cdf.ests)
}
