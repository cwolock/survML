#' Stacked binary regression with SuperLearner, using the cdf
#'
#' @param time Observed time
#' @param event Indicator of event (vs censoring)
#' @param X Covariate matrix
#' @param censored Logical, indicates whether to run regression on censored observations (vs uncensored)
#' @param bin_size Size of quantiles over which to make the stacking bins
#' @param isotonize Logical, indicates whether or not to isotonize cdf estimates using PAVA
#' @param V Number of CV folds
#' @param time_basis How to treat time
#' @param SL.library SuperLearner library
#'
#' @return An object of class \code{f_y_stack_xgboost}
#' @noRd
f_y_stack_SuperLearner <- function(time,
                                   event,
                                   X,
                                   censored,
                                   bin_size,
                                   isotonize = TRUE,
                                   V,
                                   time_basis = "continuous",
                                   SL.library){

  if (!is.null(censored)){
    if (censored == TRUE){
      time <- time[!as.logical(event)]
      X <- X[!as.logical(event),]
    } else if (censored == FALSE){
      time <- time[as.logical(event)]
      X <- X[as.logical(event),]
    }
  } else{
    time <- time
    X <- X
  }

  cv_folds <- split(sample(1:length(time)), rep(1:V, length = length(time)))

  X <- as.matrix(X)
  time <- as.matrix(time)
  dat <- data.frame(X, time)

  # if user gives bin size, set time grid based on quantiles. otherwise, every observed time
  if (!is.null(bin_size)){
    time_grid <- stats::quantile(dat$time, probs = seq(0, 1, by = bin_size))
    time_grid[1] <- 0 # manually set first point to 0, instead of first observed time
  } else{
    time_grid <- sort(unique(time))
    time_grid <- c(0, time_grid)
  }

  stacked_list <- stack_cdf(time = time,
                            X = X,
                            time_grid = time_grid,
                            time_basis = time_basis,
                            ids = TRUE)
  stacked <- stacked_list$stacked
  stacked_ids <- stacked_list$ids

  get_validRows <- function(fold_sample_ids){
    validRows <- which(stacked_ids %in% fold_sample_ids)
    return(validRows)
  }

  validRows <- lapply(cv_folds, get_validRows)

  .Y <- stacked[,ncol(stacked)]
  .X <- stacked[,-ncol(stacked)]

  fit <- SuperLearner::SuperLearner(Y = .Y,
                                    X = .X,
                                    SL.library = SL.library,
                                    family = stats::binomial(),
                                    method = 'method.NNLS',
                                    verbose = FALSE,
                                    cvControl = list(V = V,
                                                     validRows = validRows))

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
  trunc_time_grid <- time_grid[-length(time_grid)]

  if (fit$time_basis == "continuous"){
    get_stacked_pred <- function(t){
      new_stacked <- data.frame(t = t, newX)
      preds <- stats::predict(fit$reg.object, newdata=new_stacked)$pred
      return(preds)
    }

    predictions <- apply(X = matrix(newtimes), FUN = get_stacked_pred, MARGIN = 1)
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

    predictions <- apply(X = matrix(newtimes), FUN = get_preds, MARGIN = 1)
  }

  if (fit$isotonize){
    iso.cdf.ests <- t(apply(predictions, MARGIN = 1, FUN = Iso::pava))
  } else{
    iso.cdf.ests <- predictions
  }

  return(iso.cdf.ests)
}
