#' Stacked binary regression with homemade cross validation, using the cdf
#'
#' @param time Observed time
#' @param event Indicator of event (vs censoring)
#' @param X Covariate matrix
#' @param censored Logical, indicates whether to run regression on censored observations (vs uncensored)
#' @param bin_size Size of quantiles over which to make the stacking bins
#' @param isotonize Logical, indicates whether or not to isotonize cdf estimates using PAVA
#' @param V Number of CV folds
#' @param time_basis How to treat time
#'
#' @return An object of class \code{f_y_stackCVcdf}
#' @noRd
f_y_stack_xgboost <- function(time, event, X, censored, bin_size, isotonize = TRUE, V, time_basis = "continuous"){

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
    time_grid <- quantile(dat$time, probs = seq(0, 1, by = bin_size))
    time_grid[1] <- 0 # manually set first point to 0, instead of first observed time
  } else{
    time_grid <- sort(unique(time))
    time_grid <- c(0, time_grid)
  }
  # time_grid <- quantile(dat$time, probs = seq(0, 1, by = bin_size))
  # time_grid[1] <- 0 # manually set first point to 0, instead of first observed time

  tune <- list(ntrees = c(50, 100, 250, 500), max_depth = c(1,2,3),
              eta = c(0.1))

  param_grid <- expand.grid(ntrees = tune$ntrees,
                            max_depth = tune$max_depth,
                            eta = tune$eta)

  if (time_basis == "continuous"){
    get_CV_risk <- function(i){
      ntrees <- param_grid$ntrees[i]
      max_depth <- param_grid$max_depth[i]
      eta <- param_grid$eta[i]
      risks <- rep(NA, V)
      for (j in 1:V){
        train_X <- X[-cv_folds[[j]],]
        train_time <- time[-cv_folds[[j]]]
        train_stack <- conSurv:::stack(time = train_time, X = train_X, time_grid = time_grid)
        xgmat <- xgboost::xgb.DMatrix(data = train_stack[,-ncol(train_stack)], label = train_stack[,ncol(train_stack)])
        fit <- xgboost::xgboost(data = xgmat, objective="binary:logistic", nrounds = ntrees,
                                max_depth = max_depth, eta = eta,
                                verbose = FALSE, nthread = 1,
                                save_period = NULL, eval_metric = "logloss",
                                subsample = 0.5)
        test_X <- X[cv_folds[[j]],]
        test_time <- time[cv_folds[[j]]]
        test_stack <- conSurv:::stack(time = test_time, X = test_X, time_grid = time_grid)
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
    stacked <- conSurv:::stack(time = time, X = X, time_grid = time_grid)
    Y <- stacked[,ncol(stacked)]
    X <- as.matrix(stacked[,-ncol(stacked)])
    xgmat <- xgboost::xgb.DMatrix(data = X, label = Y)
    fit <- xgboost::xgboost(data = xgmat, objective="binary:logistic", nrounds = opt_ntrees,
                            max_depth = opt_max_depth, eta = opt_eta,
                            verbose = FALSE, nthread = 1,
                            save_period = NULL, eval_metric = "logloss")
  } else{
    get_CV_risk <- function(i){
      ntrees <- param_grid$ntrees[i]
      max_depth <- param_grid$max_depth[i]
      eta <- param_grid$eta[i]
      risks <- rep(NA, V)
      for (j in 1:V){
        train_X <- X[-cv_folds[[j]],]
        train_time <- time[-cv_folds[[j]]]
        train_stack <- conSurv:::stack_dummy(time = train_time, X = train_X, time_grid = time_grid)
        xgmat <- xgboost::xgb.DMatrix(data = as.matrix(train_stack[,-ncol(train_stack)]), label = as.matrix(train_stack[,ncol(train_stack)]))
        fit <- xgboost::xgboost(data = xgmat, objective="binary:logistic", nrounds = ntrees,
                                max_depth = max_depth, eta = eta,
                                verbose = FALSE, nthread = 1,
                                save_period = NULL, eval_metric = "logloss")
        test_X <- X[cv_folds[[j]],]
        test_time <- time[cv_folds[[j]]]
        test_stack <- conSurv:::stack_dummy(time = test_time, X = test_X, time_grid = time_grid)
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
    opt_params <- list(ntrees = opt_ntrees, max_depth = opt_max_depth, eta = opt_eta)
    stacked <- conSurv:::stack_dummy(time = time, X = X, time_grid = time_grid)
    Y <- stacked[,ncol(stacked)]
    X <- as.matrix(stacked[,-ncol(stacked)])
    xgmat <- xgboost::xgb.DMatrix(data = X, label = Y)
    fit <- xgboost::xgboost(data = xgmat, objective="binary:logistic", nrounds = opt_ntrees,
                            max_depth = opt_max_depth, eta = opt_eta,
                            verbose = FALSE, nthread = 1,
                            save_period = NULL, eval_metric = "logloss")
  }

  print(censored)
  print(CV_risks)
  print(fit$params)
  print(fit$niter)
  fit <- list(reg.object = fit, time_grid = time_grid, isotonize = isotonize, time_basis = time_basis)
  class(fit) <- c("f_y_stack_xgboost")
  return(fit)
}

#' Prediction function for stacked xgboost CDF
#'
#' @param fit Fitted regression object
#' @param newX Values of covariates at which to make a prediction
#' @param newtimes
#'
#' @return Matrix of predictions
#' @noRd
predict.f_y_stack_xgboost <- function(fit, newX, newtimes){

  time_grid <- fit$time_grid
  trunc_time_grid <- time_grid[-length(time_grid)]

  if (fit$time_basis == "continuous"){
    get_stacked_pred <- function(t){
      new_stacked <- as.matrix(data.frame(t = t, newX))
      preds <- predict(fit$reg.object, newdata=new_stacked)
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
      new_stacked <- as.matrix(new_stacked)
      preds <- predict(fit$reg.object, newdata=new_stacked)
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


#' Stacked binary regression using ranger with homemade cross validation
#'
#' @param time Observed time
#' @param event Indicator of event (vs censoring)
#' @param X Covariate matrix
#' @param censored Logical, indicates whether to run regression on censored observations (vs uncensored)
#' @param bin_size Size of quantiles over which to make the stacking bins
#' @param isotonize Logical, indicates whether or not to isotonize cdf estimates using PAVA
#' @param V Number of CV folds
#' @param time_basis How to treat time
#'
#' @return An object of class \code{f_y_stack_ranger}
#' @noRd
f_y_stack_ranger <- function(time, event, X, censored, bin_size, isotonize = TRUE, V, time_basis = "continuous"){

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
    time_grid <- quantile(dat$time, probs = seq(0, 1, by = bin_size))
    time_grid[1] <- 0 # manually set first point to 0, instead of first observed time
  } else{
    time_grid <- sort(unique(time))
    time_grid <- c(0, time_grid)
  }
  # time_grid <- quantile(dat$time, probs = seq(0, 1, by = bin_size))
  # time_grid[1] <- 0 # manually set first point to 0, instead of first observed time

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
      train_stack <- conSurv:::stack(time = train_time, X = train_X, time_grid = time_grid)
      fit <- ranger::ranger(formula = event_indicators ~ .,
                            data = train_stack,
                            num.trees = num.trees,
                            max.depth = max.depth,
                            mtry = mtry,
                            probability = TRUE)
      test_X <- X[cv_folds[[j]],]
      test_time <- time[cv_folds[[j]]]
      test_stack <- conSurv:::stack(time = test_time, X = test_X, time_grid = time_grid)
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
  stacked <- conSurv:::stack(time = time, X = X, time_grid = time_grid)
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
  fit <- list(reg.object = fit, time_grid = time_grid, isotonize = isotonize, time_basis = time_basis)
  class(fit) <- c("f_y_stackCVranger")
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
predict.f_y_stack_ranger <- function(fit, newX, newtimes){

  time_grid <- fit$time_grid
  trunc_time_grid <- time_grid[-length(time_grid)]

  get_stacked_pred <- function(t){
    new_stacked <- as.matrix(data.frame(t = t, newX))
    preds <- predict(fit$reg.object, data=new_stacked)$predictions
    preds <- preds[,2]
    return(preds)
  }

  predictions <- apply(X = matrix(newtimes), FUN = get_stacked_pred, MARGIN = 1)

  if (fit$isotonize){
    iso.cdf.ests <- t(apply(predictions, MARGIN = 1, FUN = Iso::pava))
  } else{
    iso.cdf.ests <- predictions
  }

  return(iso.cdf.ests)
}


#' Stacked binary cdf regression with gam
#'
#' @param time Observed time
#' @param event Indicator of event (vs censoring)
#' @param X Covariate matrix
#' @param censored Logical, indicates whether to run regression on censored observations (vs uncensored)
#' @param bin_size Size of quantiles over which to make the stacking bins
#' @param isotonize Logical, indicates whether or not to isotonize cdf estimates using PAVA
#' @param time_basis How to treat time
#' @param cts.num If a variable has more than this many unique values, consider it continuous
#'
#' @return An object of class \code{f_y_stack_gam}
#' @noRd
f_y_stack_gam <- function(time, event, X, censored, bin_size, deg.gam = 2,isotonize = TRUE, time_basis = "continuous", cts.num = 4){

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

  # if user gives bin size, set time grid based on quantiles. otherwise, every observed time
  if (!is.null(bin_size)){
    time_grid <- quantile(time, probs = seq(0, 1, by = bin_size))
    time_grid[1] <- 0 # manually set first point to 0, instead of first observed time
  } else{
    time_grid <- sort(unique(time))
    time_grid <- c(0, time_grid)
  }

  # time_grid <- quantile(time, probs = seq(0, 1, by = bin_size))
  # time_grid[1] <- 0 # manually set first point to 0, instead of first observed time

  if (time_basis == "continuous"){
    stacked <- conSurv:::stack(time = time, X = X, time_grid = time_grid)
  } else{
    stacked <- conSurv:::stack_dummy(time = time, X = X, time_grid = time_grid)
  }

  Y <- stacked[,ncol(stacked)]
  X <- stacked[,-ncol(stacked)]
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

  fit <- list(reg.object = fit.gam, time_grid = time_grid, isotonize = isotonize, time_basis = time_basis)
  class(fit) <- c("f_y_stack_gam")
  return(fit)
}

#' Prediction function for stacked GAM cdf estimation
#'
#' @param fit Fitted regression object
#' @param newX Values of covariates at which to make a prediction
#' @param newtimes
#'
#' @return Matrix of predictions
#' @noRd
predict.f_y_stack_gam <- function(fit, newX, newtimes){

  time_grid <- fit$time_grid
  trunc_time_grid <- time_grid[-length(time_grid)]

  if (fit$time_basis == "continuous"){
    get_stacked_pred <- function(t){
      new_stacked <- data.frame(t = t, newX)
      preds <- gam::predict.Gam(fit$reg.object, newdata=new_stacked, type = "response")
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
      new_stacked <- new_stacked
      preds <- gam::predict.Gam(fit$reg.object, newdata=new_stacked, type = "response")
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

#' Stacked binary cdf regression with earth
#'
#' @param time Observed time
#' @param event Indicator of event (vs censoring)
#' @param X Covariate matrix
#' @param censored Logical, indicates whether to run regression on censored observations (vs uncensored)
#' @param bin_size Size of quantiles over which to make the stacking bins
#' @param isotonize Logical, indicates whether or not to isotonize cdf estimates using PAVA
#' @param time_basis How to treat time
#' @param cts.num If a variable has more than this many unique values, consider it continuous
#'
#' @return An object of class \code{f_y_stack_earth}
#' @noRd
f_y_stack_earth <- function(time, event, X, censored, bin_size, isotonize = TRUE, time_basis = "continuous",
                            V = 5){

  #degree = 2
  penalty = 3
  nk = max(21, 2*ncol(X) + 1)
  pmethod = "backward"
  nfold = 0
  ncross = 1
  minspan = 0
  endspan = 0

  tune <- list(degree = c(1,2,3), nprune = seq.int(1, nk, length.out = 5))

  param_grid <- expand.grid(degree = tune$degree,
                            nprune = tune$nprune)

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

  # if user gives bin size, set time grid based on quantiles. otherwise, every observed time
  if (!is.null(bin_size)){
    time_grid <- quantile(time, probs = seq(0, 1, by = bin_size))
    time_grid[1] <- 0 # manually set first point to 0, instead of first observed time
  } else{
    time_grid <- sort(unique(time))
    time_grid <- c(0, time_grid)
  }

  cv_folds <- split(sample(1:length(time)), rep(1:V, length = length(time)))

  if (time_basis == "continuous"){
    get_CV_risk <- function(i){
      degree = param_grid$degree[i]
      nprune = param_grid$nprune[i]
      risks <- rep(NA, V)
      for (j in 1:V){
        train_X <- X[-cv_folds[[j]],]
        train_time <- time[-cv_folds[[j]]]
        train_stack <- conSurv:::stack(time = train_time, X = train_X, time_grid = time_grid)
        .Y <- train_stack[,ncol(train_stack)]
        .X <- train_stack[,-ncol(train_stack)]
        .X <- as.data.frame(.X)
        fit <- earth::earth(x = .X, y = .Y, degree = degree, nk = nk, penalty = penalty, pmethod = pmethod,
                            nfold = nfold, ncross = ncross, minspan = minspan, endspan = endspan, nprune = nprune,
                            glm = list(family = binomial))
        test_X <- X[cv_folds[[j]],]
        test_time <- time[cv_folds[[j]]]
        test_stack <- conSurv:::stack(time = test_time, X = test_X, time_grid = time_grid)
        preds <- predict(fit, newdata = test_stack[,-ncol(test_stack)], type = "response")
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
  } else{
    get_CV_risk <- function(i){
      degree = param_grid$degree[i]
      nprune = param_grid$nprune[i]
      risks <- rep(NA, V)
      for (j in 1:V){
        train_X <- X[-cv_folds[[j]],]
        train_time <- time[-cv_folds[[j]]]
        train_stack <- conSurv:::stack_dummy(time = train_time, X = train_X, time_grid = time_grid)
        .Y <- train_stack[,ncol(train_stack)]
        .X <- train_stack[,-ncol(train_stack)]
        .X <- as.data.frame(X)
        fit <- earth::earth(x = .X, y = .Y, degree = degree, nk = nk, penalty = penalty, pmethod = pmethod,
                            nfold = nfold, ncross = ncross, minspan = minspan, endspan = endspan, nprune = nprune,
                            glm = list(family = binomial))
        test_X <- X[cv_folds[[j]],]
        test_time <- time[cv_folds[[j]]]
        test_stack <- conSurv:::stack_dummy(time = test_time, X = test_X, time_grid = time_grid)
        preds <- predict(fit, newdata = as.matrix(test_stack[,-ncol(test_stack)]), type = "response")
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
  }

  CV_risks <- unlist(lapply(1:nrow(param_grid), get_CV_risk))
  opt_param_index <- which.min(CV_risks)
  opt_degree <- param_grid$degree[opt_param_index]
  opt_nprune <- param_grid$nprune[opt_param_index]

  if (time_basis == "continuous"){
    stacked <- conSurv:::stack(time = time, X = X, time_grid = time_grid)
  } else{
    stacked <- conSurv:::stack_dummy(time = time, X = X, time_grid = time_grid)
  }

  Y <- stacked[,ncol(stacked)]
  X <- stacked[,-ncol(stacked)]
  X <- as.data.frame(X)

  fit.earth <- earth::earth(x = X, y = Y, degree = opt_degree, nk = nk, penalty = penalty, pmethod = pmethod,
                            nfold = nfold, ncross = ncross, minspan = minspan, endspan = endspan, nprune = opt_nprune,
                            glm = list(family = binomial))

  fit <- list(reg.object = fit.earth, time_grid = time_grid, isotonize = isotonize, time_basis = time_basis)
  class(fit) <- c("f_y_stack_earth")
  return(fit)
}

#' Prediction function for stacked earth cdf estimation
#'
#' @param fit Fitted regression object
#' @param newX Values of covariates at which to make a prediction
#' @param newtimes
#'
#' @return Matrix of predictions
#' @noRd
predict.f_y_stack_earth <- function(fit, newX, newtimes){

  time_grid <- fit$time_grid
  trunc_time_grid <- time_grid[-length(time_grid)]

  if (fit$time_basis == "continuous"){
    get_stacked_pred <- function(t){
      new_stacked <- data.frame(t = t, newX)
      preds <- predict(fit$reg.object, newdata = new_stacked, type = "response")
      return(preds)
    }

    predictions <- apply(X = matrix(newtimes), FUN = get_stacked_pred, MARGIN = 1)
    #predictions <- apply(X = matrix(fit$time_grid), FUN = get_stacked_pred, MARGIN = 1)
  } else{
    get_preds <- function(t){
      dummies <- matrix(0, ncol = length(time_grid), nrow = nrow(newX))
      index <- max(which(time_grid <= t))
      dummies[,index] <- 1
      new_stacked <- cbind(dummies, newX)
      risk_set_names <- paste0("risk_set_", seq(1, (length(time_grid))))
      colnames(new_stacked)[1:length(time_grid)] <- risk_set_names
      new_stacked <- new_stacked
      preds <- predict(fit$reg.object, newdata = new_stacked, type = "response")
      return(preds)
    }

    predictions <- apply(X = matrix(newtimes), FUN = get_preds, MARGIN = 1)
  }

  #predictions[,1] <- 0

  if (fit$isotonize){
    iso.cdf.ests <- t(apply(predictions, MARGIN = 1, FUN = Iso::pava))
  } else{
    iso.cdf.ests <- predictions
  }

  # cdf.ests <- matrix(unlist(lapply(1:nrow(iso.cdf.ests),
  #                           function(v){
  #                             approx(x = fit$time_grid, y = iso.cdf.ests[v,], xout = newtimes, method = "linear",
  #                                    rule = 2)$y
  #                           })),
  #                    nrow = nrow(newX),
  #                    ncol = length(newtimes),
  #                    byrow = TRUE)

  return(iso.cdf.ests)
}
