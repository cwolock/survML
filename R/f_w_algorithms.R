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
f_w_stack_xgboost <- function(time, event, entry, X, censored, bin_size, V, time_basis = "continuous"){

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

  time_grid <- quantile(dat$time, probs = seq(0, 1, by = bin_size))
  time_grid[1] <- 0 # manually set first point to 0, instead of first observed time

  tune <- list(ntrees = c(100, 200, 300, 500, 1000), max_depth = c(1,2,3),
              eta = c(0.05))

  param_grid <- expand.grid(ntrees = tune$ntrees,
                            max_depth = tune$max_depth,
                            eta = tune$eta)

  get_CV_risk <- function(i){
    ntrees <- param_grid$ntrees[i]
    max_depth <- param_grid$max_depth[i]
    eta <- param_grid$eta[i]
    risks <- rep(NA, V)
    for (j in 1:V){
      train_X <- X[-cv_folds[[j]],,drop=FALSE]
      train_time <- time[-cv_folds[[j]]]
      train_entry <- entry[-cv_folds[[j]]]
      train_stack <- conSurv:::stack_entry(time = train_time, entry = train_entry, X = train_X, time_grid = time_grid)
      xgmat <- xgboost::xgb.DMatrix(data = train_stack[,-ncol(train_stack)], label = train_stack[,ncol(train_stack)])
      fit <- xgboost::xgboost(data = xgmat, objective="binary:logistic", nrounds = ntrees,
                       max_depth = max_depth, eta = eta,
                       verbose = FALSE, nthread = 1,
                       save_period = NULL, eval_metric = "logloss")
      test_X <- X[cv_folds[[j]],,drop=FALSE]
      test_time <- time[cv_folds[[j]]]
      test_entry <- entry[cv_folds[[j]]]
      test_stack <- conSurv:::stack_entry(time = test_time, entry = test_entry, X = test_X, time_grid = time_grid)
      preds <- predict(fit, newdata = test_stack[,-ncol(test_stack)])
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
  opt_params <- list(ntrees = opt_ntrees, max_depth = opt_max_depth, eta = opt_eta)
  stacked <- conSurv:::stack_entry(time = time, entry = entry, X = X, time_grid = time_grid)
  Y <- stacked[,ncol(stacked)]
  X <- as.matrix(stacked[,-ncol(stacked)])
  xgmat <- xgboost::xgb.DMatrix(data = X, label = Y)
  fit <- xgboost::xgboost(data = xgmat, objective="binary:logistic", nrounds = opt_ntrees,
                          max_depth = opt_max_depth, eta = opt_eta,
                          verbose = FALSE, nthread = 1,
                          save_period = NULL, eval_metric = "logloss")

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


#' Stacked binary regression using ranger with homemade cross validation
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
#' @return An object of class \code{f_w_stack_ranger}
#' @noRd
f_w_stack_ranger <- function(time, event, entry, X, censored, bin_size, V, time_basis = "continuous"){

  if (!is.null(censored)){
    if (censored == TRUE){
      time <- time[!as.logical(event)]
      entry <- entry[!as.logical(event)]
      X <- X[!as.logical(event),]
    } else if (censored == FALSE){
      time <- time[as.logical(event)]
      entry <- entry[as.logical(event)]
      X <- X[as.logical(event),]
    }
  } else{
    time <- time
    X <- X
    entry <- entry
  }

  cv_folds <- split(sample(1:length(time)), rep(1:V, length = length(time)))

  X <- as.matrix(X)
  time <- as.matrix(time)
  entry <- as.matrix(entry)
  dat <- data.frame(X, time, entry)

  time_grid <- quantile(dat$entry, probs = seq(0, 1, by = bin_size))
  time_grid[1] <- 0 # manually set first point to 0, instead of first observed time

  tune <- list(num.trees = c(250, 500, 1000), max.depth = c(1,2,3,4), mtry = c(1,2,3))

  param_grid <- expand.grid(num.trees = tune$num.trees,
                            max.depth = tune$max.depth,
                            mtry = tune$mtry)

  get_CV_risk <- function(i){
    num.trees <- param_grid$num.trees[i]
    max.depth <- param_grid$max.depth[i]
    mtry <- param_grid$mtry[i]
    risks <- rep(NA, V)
    for (j in 1:V){
      train_X <- X[-cv_folds[[j]],]
      train_time <- time[-cv_folds[[j]]]
      train_entry <- entry[-cv_folds[[j]]]
      train_stack <- conSurv:::stack_entry(time = train_time, entry = train_entry, X = train_X, time_grid = time_grid)
      fit <- ranger::ranger(formula = event_indicators ~ .,
                            data = train_stack,
                            num.trees = num.trees,
                            max.depth = max.depth,
                            mtry = mtry,
                            probability = TRUE)
      test_X <- X[cv_folds[[j]],]
      test_time <- time[cv_folds[[j]]]
      test_entry <- entry[cv_folds[[j]]]
      test_stack <- conSurv:::stack_entry(time = test_time, entry = test_entry, X = test_X, time_grid = time_grid)
      preds <- predict(fit, data = test_stack[,-ncol(test_stack)])$predictions
      preds <- preds[,2] # I think this is choosing the correct column but the output is not labeled...
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
  opt_num.trees <- param_grid$num.trees[opt_param_index]
  opt_max.depth <- param_grid$max.depth[opt_param_index]
  opt_mtry <- param_grid$mtry[opt_param_index]
  opt_params <- list(ntrees = opt_num.trees, max_depth = opt_max.depth, mtry = opt_mtry)
  stacked <- conSurv:::stack_entry(time = time, entry = entry, X = X, time_grid = time_grid)
  Y <- stacked[,ncol(stacked)]
  X <- as.matrix(stacked[,-ncol(stacked)])
  fit <- ranger::ranger(formula = event_indicators ~ .,
                        data = stacked,
                        num.trees = opt_num.trees,
                        max.depth = opt_max.depth,
                        mtry = opt_mtry,
                        probability = TRUE)

  print(censored)
  print(CV_risks)
  print(opt_max.depth)
  print(fit$num.trees)
  print(opt_mtry)
  fit <- list(reg.object = fit, time_grid = time_grid, time_basis = time_basis)
  class(fit) <- c("f_w_stack_ranger")
  return(fit)
}

#' Prediction function for stacked CV cdf, ranger
#'
#' @param fit Fitted regression object
#' @param newX Values of covariates at which to make a prediction
#' @param newtimes
#'
#' @return Matrix of predictions
#' @noRd
predict.f_w_stack_ranger <- function(fit, newX, newtimes){

  time_grid <- fit$time_grid
  trunc_time_grid <- time_grid[-length(time_grid)]

  get_stacked_pred <- function(t){
    new_stacked <- as.matrix(data.frame(t = t, newX))
    preds <- predict(fit$reg.object, data=new_stacked)$predictions
    preds <- preds[,2]
    return(preds)
  }

  predictions <- apply(X = matrix(newtimes), FUN = get_stacked_pred, MARGIN = 1)

  return(predictions)
}

#' Stacked binary regression with gam
#'
#' @param time Observed time
#' @param event Indicator of event (vs censoring)
#' @param entry Truncation variable, time of entry into the study
#' @param X Covariate matrix
#' @param censored Logical, indicates whether to run regression on censored observations (vs uncensored)
#' @param bin_size Size of quantiles over which to make the stacking bins
#' @param time_basis How to treat time
#'
#' @return An object of class \code{f_w_stack_gam}
#' @noRd
f_w_stack_gam <- function(time, event, entry, X, censored, bin_size, deg.gam = 2, cts.num = 4,time_basis = "continuous"){

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

  X <- as.matrix(X)
  time <- as.matrix(time)
  entry <- as.matrix(entry)
  dat <- data.frame(X, time, entry)

  # should my grid be quantiles of entry, or of Y??
  time_grid <- quantile(dat$time, probs = seq(0, 1, by = bin_size))
  time_grid[1] <- 0 # manually set first point to 0, instead of first observed time

  stacked <- conSurv:::stack_entry(time = time, entry = entry, X = X, time_grid = time_grid)
  Y <- stacked[,ncol(stacked)]
  X <- as.matrix(stacked[,-ncol(stacked)])
  cts.x <- apply(X, 2, function(x) (length(unique(x)) > cts.num))
  X <- as.data.frame(X)
  if (sum(!cts.x) > 0) {
    gam.model <- as.formula(paste("Y~",
                                  paste(paste("s(", colnames(X[, cts.x, drop = FALSE]), ",", deg.gam,")", sep=""),
                                        collapse = "+"), "+", paste(colnames(X[, !cts.x, drop=FALSE]), collapse = "+")))
  } else {
    gam.model <- as.formula(paste("Y~",
                                  paste(paste("s(", colnames(X[, cts.x, drop = FALSE]), ",", deg.gam, ")", sep=""),
                                        collapse = "+")))
  }
  # fix for when all variables are binomial
  if (sum(!cts.x) == length(cts.x)) {
    gam.model <- as.formula(paste("Y~", paste(colnames(X), collapse = "+"), sep = ""))
  }
  fit.gam <- gam::gam(gam.model, data = X, family = binomial(), control = gam::gam.control(maxit = 50, bf.maxit = 50))

  fit <- list(reg.object = fit.gam, time_grid = time_grid, time_basis = time_basis)
  class(fit) <- c("f_w_stack_gam")
  return(fit)
}

#' Prediction function for stacked gam
#'
#' @param fit Fitted regression object
#' @param newX Values of covariates at which to make a prediction
#' @param newtimes
#'
#' @return Matrix of predictions
#' @noRd
predict.f_w_stack_gam <- function(fit, newX, newtimes){

  time_grid <- fit$time_grid
  trunc_time_grid <- time_grid[-length(time_grid)]

  get_stacked_pred <- function(t){
    new_stacked <- data.frame(t = t, newX)
    preds <- gam::predict.Gam(fit$reg.object, newdata=new_stacked, type = "response")
    return(preds)
  }

  predictions <- apply(X = matrix(newtimes), FUN = get_stacked_pred, MARGIN = 1)

  return(predictions)
}

#' Stacked binary regression with earth
#'
#' @param time Observed time
#' @param event Indicator of event (vs censoring)
#' @param entry Truncation variable, time of entry into the study
#' @param X Covariate matrix
#' @param censored Logical, indicates whether to run regression on censored observations (vs uncensored)
#' @param bin_size Size of quantiles over which to make the stacking bins
#' @param time_basis How to treat time
#'
#' @return An object of class \code{f_w_stack_gam}
#' @noRd
f_w_stack_earth <- function(time, event, entry, X, censored, bin_size,time_basis = "continuous", degree = 2,
                            penalty = 3, nk = max(21, 2*ncol(X) + 1), pmethod = "backward",
                            nfold = 0, ncross = 1, minspan = 0, endspan = 0){

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

  X <- as.matrix(X)
  time <- as.matrix(time)
  entry <- as.matrix(entry)
  dat <- data.frame(X, time, entry)

  # should my grid be quantiles of entry, or of Y??
  time_grid <- quantile(dat$time, probs = seq(0, 1, by = bin_size))
  time_grid[1] <- 0 # manually set first point to 0, instead of first observed time

  stacked <- conSurv:::stack_entry(time = time, entry = entry, X = X, time_grid = time_grid)
  Y <- stacked[,ncol(stacked)]
  X <- as.matrix(stacked[,-ncol(stacked)])
  X <- as.data.frame(X)

  fit.earth <- earth::earth(x = X, y = Y, degree = degree, nk = nk, penalty = penalty, pmethod = pmethod,
                            nfold = nfold, ncross = ncross, minspan = minspan, endspan = endspan,
                            glm = list(family = binomial))

  fit <- list(reg.object = fit.earth, time_grid = time_grid, time_basis = time_basis)
  class(fit) <- c("f_w_stack_earth")
  return(fit)
}

#' Prediction function for stacked earth
#'
#' @param fit Fitted regression object
#' @param newX Values of covariates at which to make a prediction
#' @param newtimes
#'
#' @return Matrix of predictions
#' @noRd
predict.f_w_stack_earth <- function(fit, newX, newtimes){

  time_grid <- fit$time_grid
  trunc_time_grid <- time_grid[-length(time_grid)]

  get_stacked_pred <- function(t){
    new_stacked <- data.frame(t = t, newX)
    preds <- predict(fit$reg.object, newdata=new_stacked, type = "response")
    return(preds)
  }

  predictions <- apply(X = matrix(newtimes), FUN = get_stacked_pred, MARGIN = 1)

  return(predictions)
}

