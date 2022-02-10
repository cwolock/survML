#' Smoothed Nadaraya-Watson estimator
#'
#' @param time Observed time
#' @param event Indicator of event (vs censoring)
#' @param X Covariate matrix
#' @param censored Logical, indicates whether to run regression on censored observations (vs uncensored)
#' @param bw bw
#' @param bwy bwy
#' @param kernel_type Kernel_type
#' @param kernel_order Kernel_order
#'
#' @return An object of class \code{f_y_smoothnw}
#' @noRd
f_y_smoothnw <- function(time, event, X, censored, bw, bwy, kernel_type = "gaussian", kernel_order = 2){

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

  fmla <- as.formula(paste("ind ~ ", paste(colnames(X), collapse = "+")))
  ind <- time
  dat <- data.frame(cbind(ind = ind, X))
  bws <- eval(bquote(np::npcdistbw(formula = .(fmla),
                                   data = dat,
                                   regtype = "lc",
                                   ckertype = kernel_type,
                                   ckerorder = kernel_order)))
  nw_fit <- eval(bquote(np::npcdist(formula = .(fmla),
                                    data = dat,
                                    bws = bws, #c(bwy, rep(bw, ncol(X))),
                                    regtype = "lc",
                                    ckertype = kernel_type,
                                    ckerorder = kernel_order)))

  fit <- list(reg.object = nw_fit)
  class(fit) <- c("f_y_smoothnw")
  return(fit)
}

#' Prediction function for smoothed Nadaraya-Watson estimator
#'
#' @param fit Fitted regression object
#' @param newX Values of covariates at which to make a prediction
#' @param newtimes
#'
#' @return Matrix of predictions
#' @noRd
predict.f_y_smoothnw <- function(fit, newX, newtimes){
  predict_nw <- function(t){
    pred <- predict(fit$reg.object, newdata = cbind(ind = rep(t, nrow(newX)), newX))
    return(pred)
  }
  predictions <- apply(X = as.matrix(newtimes),
                       MARGIN = 1,
                       FUN = predict_nw)
  return(predictions)
}

#' Smoothed locally linear kernel regression estimator
#'
#' @param time Observed time
#' @param event Indicator of event (vs censoring)
#' @param X Covariate matrix
#' @param censored Logical, indicates whether to run regression on censored observations (vs uncensored)
#' @param bw bw
#' @param bwy bwy
#' @param kernel_type Kernel_type
#' @param kernel_order Kernel_order
#'
#' @return An object of class \code{f_y_smoothlllkern}
#' @noRd
f_y_smoothllkern <- function(time, event, X, censored, bw, bwy, kernel_type = "gaussian", kernel_order = 2){

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

  fmla <- as.formula(paste("ind ~ ", paste(colnames(X), collapse = "+")))
  ind <- time
  dat <- data.frame(cbind(ind = ind, X))
  bws <- eval(bquote(np::npcdistbw(formula = .(fmla),
                                   data = dat,
                                   regtype = "ll",
                                   ckertype = kernel_type,
                                   ckerorder = kernel_order)))
  ll_fit <- eval(bquote(np::npcdist(formula = .(fmla),
                                    data = dat,
                                    bws = bws,#c(bwy, rep(bw, ncol(X))),
                                    regtype = "ll",
                                    ckertype = kernel_type,
                                    ckerorder = kernel_order)))

  fit <- list(reg.object = ll_fit)
  class(fit) <- c("f_y_smoothllkern")
  return(fit)
}

#' Prediction function for smoothed locally linear kernel estimator
#'
#' @param fit Fitted regression object
#' @param newX Values of covariates at which to make a prediction
#' @param newtimes
#'
#' @return Matrix of predictions
#' @noRd
predict.f_y_smoothllkern <- function(fit, newX, newtimes){
  predict_ll <- function(t){
    pred <- predict(fit$reg.object, newdata = cbind(ind = rep(t, nrow(newX)), newX))
    return(pred)
  }
  predictions <- apply(X = as.matrix(newtimes),
                       MARGIN = 1,
                       FUN = predict_ll)
  return(predictions)
}

#' Stacked binary regression with Super Learner, using the cdf
#'
#' @param time Observed time
#' @param event Indicator of event (vs censoring)
#' @param X Covariate matrix
#' @param censored Logical, indicates whether to run regression on censored observations (vs uncensored)
#' @param bin_size Size of quantiles over which to make the stacking bins
#' @param isotonize Logical, indicates whether or not to isotonize cdf estimates using PAVA
#' @param SL.library Super Learner library
#' @param time_basis How to treat time
#'
#' @return An object of class \code{f_y_stackSLcdf}
#' @noRd
f_y_stackSLcdf <- function(time, event, X, censored, bin_size, isotonize = TRUE, SL.library, time_basis = "continuous"){

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

  X <- as.matrix(X)
  time <- as.matrix(time)
  dat <- data.frame(X, time)

  time_grid <- quantile(dat$time, probs = seq(0, 1, by = bin_size))
  time_grid[1] <- 0 # manually set first point to 0, instead of first observed time
  trunc_time_grid <- time_grid[-length(time_grid)] # do I need to truncate if treating time as continuous? look at this later
  # alternatively, time grid could just be observed times, or a fixed (non-data-dependent) grid


  # we will treat time as continuous
  if (time_basis == "continuous"){
    ncol_stacked <- ncol(X) + 2 # covariates, time, binary outcome
  } else{
    ncol_stacked <- ncol(X) + length(trunc_time_grid) + 1 # covariates, time, binary outcome
  }
  stacked <- matrix(NA, ncol = ncol_stacked, nrow = 1)
  for (i in 1:(length(trunc_time_grid)-1)){ # changed this to not do anything in last time bin
    event_indicators <- matrix(ifelse(time <= time_grid[i + 1], 1, 0))
    if (time_basis == "continuous"){
      t <- time_grid[i + 1]
      newdata <- as.matrix(cbind(t, X, event_indicators))
    }
    else{
      dummies <- matrix(0, ncol = length(trunc_time_grid), nrow = nrow(X))
      dummies[,i] <- 1
      newdata <- as.matrix(cbind(dummies, X, event_indicators))
    }
    stacked <- rbind(stacked, newdata)
  }

  stacked <- stacked[-1,]
  if (time_basis == "continuous"){
    colnames(stacked)[ncol(stacked)] <- "event_indicators"
  } else{
    risk_set_names <- paste0("risk_set_", seq(1, (length(trunc_time_grid))))
    colnames(stacked)[1:(length(trunc_time_grid))] <- risk_set_names
    colnames(stacked)[ncol(stacked)] <- "event_indicators"
  }

  stacked <- data.frame(stacked)
  Y <- stacked$event_indicators
  X <- stacked[,-ncol(stacked)]


  tune = list(ntrees = c(1000, 5000), max_depth = c(2,3, 4), minobspernode = 10,
              shrinkage = c(0.05, 0.1, 0.5))
  xgb_grid = SuperLearner::create.SL.xgboost(tune = tune, detailed_names = TRUE)
  fit <- SuperLearner::SuperLearner(Y = Y,
                                    X = X,
                                    family = binomial(),
                                    SL.library = xgb_grid$names,#SL.library,
                                    method = "method.NNloglik",
                                    verbose = FALSE)
  print(fit)

  fit <- list(reg.object = fit, time_grid = time_grid, isotonize = isotonize, time_basis = time_basis)
  class(fit) <- c("f_y_stackSLcdf")
  return(fit)
}

#' Prediction function for stacked SL cdf
#'
#' @param fit Fitted regression object
#' @param newX Values of covariates at which to make a prediction
#' @param newtimes
#'
#' @return Matrix of predictions
#' @noRd
predict.f_y_stackSLcdf <- function(fit, newX, newtimes){

  time_grid <- fit$time_grid
  trunc_time_grid <- time_grid[-length(time_grid)]

  get_stacked_pred <- function(t){
    new_stacked <- data.frame(t = t, newX)
    preds <- predict(fit$reg.object, newdata=new_stacked)$pred
    return(preds)
  }


  # get_stacked_pred <- function(t){
  #   new_stacked <- newX
  #   dummy_stacked <- matrix(NA, ncol = length(trunc_time_grid), nrow = 1)
  #   for (i in 1:n_bins){
  #     dummies <- matrix(0, ncol = length(trunc_time_grid), nrow = nrow(newX))
  #     if (t > trunc_time_grid[i]){
  #       dummies[,i] <- 1
  #     }
  #     dummy_stacked <- rbind(dummy_stacked, dummies)
  #   }
  #   dummy_stacked <- dummy_stacked[-1,]
  #   new_stacked <- cbind(dummy_stacked, new_stacked)
  #   risk_set_names <- paste0("risk_set_", seq(1, (length(trunc_time_grid))))
  #   colnames(new_stacked) <- c(risk_set_names, colnames(newX))
  #   new_stacked <- data.frame(new_stacked)
  #   haz_preds <- predict(fit, newdata=new_stacked)$pred
  #   haz_preds <- matrix(haz_preds, nrow = nrow(newX))
  #   surv_preds <- 1 - haz_preds
  #   total_surv_preds <- apply(surv_preds, prod, MARGIN = 1)
  #   return(total_surv_preds)
  # }

  predictions <- apply(X = matrix(newtimes), FUN = get_stacked_pred, MARGIN = 1)

  if (fit$isotonize){
    # isotonize??
    iso.cdf.ests <- t(apply(predictions, MARGIN = 1, FUN = Iso::pava))
  } else{
    iso.cdf.ests <- predictions
  }

  return(iso.cdf.ests)
}

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
f_y_stackCVcdf <- function(time, event, X, censored, bin_size, isotonize = TRUE, V, time_basis = "continuous"){

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

  time_grid <- quantile(dat$time, probs = seq(0, 1, by = bin_size))
  time_grid[1] <- 0 # manually set first point to 0, instead of first observed time

  tune <- list(ntrees = c(500, 1000, 2000), max_depth = c(1,2,3),
              eta = c(0.01, 0.1))

  param_grid <- expand.grid(ntrees = tune$ntrees,
                            max_depth = tune$max_depth,
                            eta = tune$eta)

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
                       save_period = NULL, eval_metric = "logloss")
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

  print(censored)
  print(CV_risks)
  print(fit$params)
  print(fit$niter)
  fit <- list(reg.object = fit, time_grid = time_grid, isotonize = isotonize, time_basis = time_basis)
  class(fit) <- c("f_y_stackCVcdf")
  return(fit)
}

#' Prediction function for stacked CV cdf
#'
#' @param fit Fitted regression object
#' @param newX Values of covariates at which to make a prediction
#' @param newtimes
#'
#' @return Matrix of predictions
#' @noRd
predict.f_y_stackCVcdf <- function(fit, newX, newtimes){

  time_grid <- fit$time_grid
  trunc_time_grid <- time_grid[-length(time_grid)]

  get_stacked_pred <- function(t){
    new_stacked <- as.matrix(data.frame(t = t, newX))
    preds <- predict(fit$reg.object, newdata=new_stacked)
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
