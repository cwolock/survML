#' Stacked binary regression with homemade cross validation, xgboost
#'
#' @param time Observed time
#' @param event Indicator of event (vs censoring)
#' @param entry Truncation variable, time of entry into the study
#' @param X Covariate matrix
#' @param censored Logical, indicates whether to run regression on censored observations (vs uncensored)
#' @param bin_size Size of quantiles over which to make the stacking bins
#' @param V Number of CV folds
#' @param time_basis How to treat time
#'
#' @return An object of class \code{f_w_stack_xgboost}
#' @noRd
f_w_stack_xgboost <- function(time, event, entry, X, censored, bin_size, V,
                              time_basis = "continuous", direction = "forward",
                              tuning_params = NULL){

  if (!is.null(censored)){
    if (censored == TRUE){
      time <- time[!as.logical(event)]
      entry <- entry[!as.logical(event)]
      X <- X[!as.logical(event),]
    } else if (censored == FALSE){
      time <- time[as.logical(event)]
      X <- X[as.logical(event),]
      entry <- entry[as.logical(event)]
    }
  } else{
    time <- time
    entry <- entry
    X <- X
  }

  cv_folds <- split(sample(1:length(time)), rep(1:V, length = length(time)))

  X <- as.matrix(X)
  time <- as.matrix(time)
  entry <- as.matrix(entry)
  dat <- data.frame(X, time, entry)

  if (!is.null(bin_size)){
    #time_grid <- quantile(dat$time, probs = seq(0, 1, by = bin_size))

    time_grid <- quantile(time, probs = seq(0, 1, by = bin_size))
    if (direction == "reverse"){
      time_grid <- c(time_grid, max(entry)) # manually set first point to 0, instead of first observed time
    } else{
      time_grid[1] <- 0 # manually set first point to 0, instead of first observed time
    }
  } else{
    time_grid <- sort(unique(time))
    if (direction == "reverse"){
      time_grid <- c(time_grid, max(entry)) # manually set first point to 0, instead of first observed time
    } else{
      time_grid <- c(0, time_grid)
    }
  }
  # time_grid <- quantile(dat$time, probs = seq(0, 1, by = bin_size))
  # time_grid[1] <- 0 # manually set first point to 0, instead of first observed time

  if (is.null(tuning_params)){
    tune <- list(ntrees = c(50, 100, 250, 500), max_depth = c(1,2,3),
                 eta = c(0.1), subsample = c(0.5))
  } else{
    tune <- tuning_params
  }

  param_grid <- expand.grid(ntrees = tune$ntrees,
                            max_depth = tune$max_depth,
                            eta = tune$eta,
                            subsample = tune$subsample)

  get_CV_risk <- function(i){
    ntrees <- param_grid$ntrees[i]
    max_depth <- param_grid$max_depth[i]
    eta <- param_grid$eta[i]
    subsample <- param_grid$subsample[i]
    risks <- rep(NA, V)
    for (j in 1:V){
      train_X <- X[-cv_folds[[j]],,drop=FALSE]
      train_time <- time[-cv_folds[[j]]]
      train_entry <- entry[-cv_folds[[j]]]
      train_stack <- survML:::stack_entry(time = train_time, entry = train_entry, X = train_X, time_grid = time_grid,
                                           direction = direction)$stacked
      xgmat <- xgboost::xgb.DMatrix(data = as.matrix(train_stack[,-ncol(train_stack)]),
                                    label = as.matrix(train_stack[,ncol(train_stack)]))
      fit <- xgboost::xgboost(data = xgmat, objective="binary:logistic", nrounds = ntrees,
                       max_depth = max_depth, eta = eta,
                       verbose = FALSE, nthread = 1,
                       save_period = NULL, eval_metric = "logloss",
                       subsample = subsample)
      test_X <- X[cv_folds[[j]],,drop=FALSE]
      test_time <- time[cv_folds[[j]]]
      test_entry <- entry[cv_folds[[j]]]
      test_stack <- survML:::stack_entry(time = test_time, entry = test_entry, X = test_X, time_grid = time_grid,
                                          direction = direction)$stacked
      preds <- predict(fit, newdata = as.matrix(test_stack[,-ncol(test_stack)]))
      preds[preds == 1] <- 0.99 # this is a hack, but come back to it later
      truth <- test_stack[,ncol(test_stack)]
      log_loss <- lapply(1:length(preds), function(x) { # using log loss right now
        -truth[x] * log(preds[x]) - (1-truth[x])*log(1 - preds[x])
      })
      log_loss <- unlist(log_loss)
      sum_log_loss <- sum(log_loss)
      risks[j] <- sum_log_loss
    }
    return(sum(risks))
  }

  CV_risks <- unlist(lapply(1:nrow(param_grid), get_CV_risk))

  opt_param_index <- which.min(CV_risks)
  opt_ntrees <- param_grid$ntrees[opt_param_index]
  opt_max_depth <- param_grid$max_depth[opt_param_index]
  opt_eta <- param_grid$eta[opt_param_index]
  opt_subsample <- param_grid$subsample[opt_param_index]
  opt_params <- list(ntrees = opt_ntrees, max_depth = opt_max_depth, eta = opt_eta, subsample = opt_subsample)
  stacked <- survML:::stack_entry(time = time, entry = entry, X = X, time_grid = time_grid, direction = direction)$stacked
  .Y <- stacked[,ncol(stacked)]
  .X <- as.matrix(stacked[,-ncol(stacked)])
  xgmat <- xgboost::xgb.DMatrix(data = .X, label = .Y)
  fit <- xgboost::xgboost(data = xgmat, objective="binary:logistic", nrounds = opt_ntrees,
                          max_depth = opt_max_depth, eta = opt_eta,
                          verbose = FALSE, nthread = 1,
                          save_period = NULL, eval_metric = "logloss", subsample = opt_subsample)

  print(censored)
  print(CV_risks)
  print(fit$params)
  print(fit$niter)
  fit <- list(reg.object = fit, time_grid = time_grid, time_basis = time_basis)
  class(fit) <- c("f_w_stack_xgboost")
  return(fit)
}

#' Prediction function for stacked xgboost
#'
#' @param fit Fitted regression object
#' @param newX Values of covariates at which to make a prediction
#' @param newtimes
#'
#' @return Matrix of predictions
#' @noRd
predict.f_w_stack_xgboost <- function(fit, newX, newtimes){

  time_grid <- fit$time_grid
  trunc_time_grid <- time_grid[-length(time_grid)]

  get_stacked_pred <- function(t){
    new_stacked <- as.matrix(data.frame(t = t, newX))
    preds <- predict(fit$reg.object, newdata=new_stacked)
    return(preds)
  }

  predictions <- apply(X = matrix(newtimes), FUN = get_stacked_pred, MARGIN = 1)

  return(predictions)
}

#' Stacked binary regression using SuperLearner
#'
#' @param time Observed time
#' @param event Indicator of event (vs censoring)
#' @param entry Truncation variable, time of entry into the study
#' @param X Covariate matrix
#' @param censored Logical, indicates whether to run regression on censored observations (vs uncensored)
#' @param bin_size Size of quantiles over which to make the stacking bins
#' @param V Number of CV folds
#' @param time_basis How to treat time
#' @param SL.library SuperLearner library
#'
#' @return An object of class \code{f_w_stack_SuperLearner}
#' @noRd
f_w_stack_SuperLearner <- function(time, event, entry, X, censored, bin_size, V,
                              time_basis = "continuous", direction = "forward",
                              SL.library){

  if (!is.null(censored)){
    if (censored == TRUE){
      time <- time[!as.logical(event)]
      entry <- entry[!as.logical(event)]
      X <- X[!as.logical(event),]
    } else if (censored == FALSE){
      time <- time[as.logical(event)]
      X <- X[as.logical(event),]
      entry <- entry[as.logical(event)]
    }
  } else{
    time <- time
    entry <- entry
    X <- X
  }

  cv_folds <- split(sample(1:length(time)), rep(1:V, length = length(time)))

  X <- as.matrix(X)
  time <- as.matrix(time)
  entry <- as.matrix(entry)
  dat <- data.frame(X, time, entry)

  if (!is.null(bin_size)){
    #time_grid <- quantile(dat$time, probs = seq(0, 1, by = bin_size))

    time_grid <- quantile(time, probs = seq(0, 1, by = bin_size))
    if (direction == "reverse"){
      time_grid <- c(time_grid, max(entry)) # manually set first point to 0, instead of first observed time
    } else{
      time_grid[1] <- 0 # manually set first point to 0, instead of first observed time
    }
  } else{
    time_grid <- sort(unique(time))
    if (direction == "reverse"){
      time_grid <- c(time_grid, max(entry)) # manually set first point to 0, instead of first observed time
    } else{
      time_grid <- c(0, time_grid)
    }
  }
  # time_grid <- quantile(dat$time, probs = seq(0, 1, by = bin_size))
  # time_grid[1] <- 0 # manually set first point to 0, instead of first observed time

  stacked_list <- survML:::stack_entry(time = time,
                                       entry = entry,
                                       X = X,
                                       time_grid = time_grid,
                                       direction = direction,
                                       ids = TRUE)
  stacked <- stacked_list$stacked
  stacked_ids <- stacked_list$ids
  .Y <- stacked[,ncol(stacked)]
  .X <- stacked[,-ncol(stacked)]

  get_validRows <- function(fold_sample_ids){
    validRows <- which(stacked_ids %in% fold_sample_ids)
    return(validRows)
  }

  validRows <- lapply(cv_folds, get_validRows)

  fit <- SuperLearner::SuperLearner(Y = .Y,
                                    X = .X,
                                    SL.library = SL.library,
                                    family = binomial(),
                                    method = 'method.NNLS',
                                    verbose = TRUE,
                                    cvControl = list(V = V,
                                                     validRows = validRows))

  fit <- list(reg.object = fit, time_grid = time_grid, time_basis = time_basis)
  class(fit) <- c("f_w_stack_SuperLearner")
  return(fit)
}

#' Prediction function for stacked SuperLearner
#'
#' @param fit Fitted regression object
#' @param newX Values of covariates at which to make a prediction
#' @param newtimes
#'
#' @return Matrix of predictions
#' @noRd
predict.f_w_stack_SuperLearner <- function(fit, newX, newtimes){

  time_grid <- fit$time_grid
  trunc_time_grid <- time_grid[-length(time_grid)]

  get_stacked_pred <- function(t){
    new_stacked <- data.frame(t = t, newX)
    preds <- predict(fit$reg.object, newdata=new_stacked)$pred
    return(preds)
  }

  predictions <- apply(X = matrix(newtimes), FUN = get_stacked_pred, MARGIN = 1)

  return(predictions)
}



