#' Stacked binary regression using SuperLearner
#'
#' @param time Observed time
#' @param event Indicator of event (vs censoring)
#' @param entry Truncation variable, time of entry into the study
#' @param X Covariate matrix
#' @param censored Logical, indicates whether to run regression on censored observations (vs uncensored)
#' @param bin_size Size of quantiles over which to make the stacking bins
#' @param time_basis How to treat time
#'
#' @return An object of class \code{f_w_stack_SuperLearner}
#' @noRd
f_w_stack_SuperLearner <- function(time,
                                   event,
                                   entry,
                                   X,
                                   censored,
                                   bin_size,
                                   SL_control,
                                   time_basis,
                                   direction){

  if (!is.null(censored)){
    if (censored == TRUE){
      time <- time[!as.logical(event)]
      entry <- entry[!as.logical(event)]
      X <- X[!as.logical(event),]
      obsWeights <- SL_control$obsWeights[!as.logical(event)]
    } else if (censored == FALSE){
      time <- time[as.logical(event)]
      X <- X[as.logical(event),]
      entry <- entry[as.logical(event)]
      obsWeights <- SL_control$obsWeights[as.logical(event)]
    }
  } else{
    time <- time
    entry <- entry
    X <- X
    obsWeights <- SL_control$obsWeights
  }

  cv_folds <- split(sample(1:length(time)), rep(1:SL_control$V, length = length(time)))

  X <- as.matrix(X)
  time <- as.matrix(time)
  entry <- as.matrix(entry)
  dat <- data.frame(X, time, entry)


  if (!is.null(bin_size)){
    #time_grid <- quantile(dat$time, probs = seq(0, 1, by = bin_size))
    time_grid <- sort(unique(stats::quantile(time, probs = seq(0, 1, by = bin_size))))
    time_grid <- c(0, time_grid) # 013123 changed this to try to get better predictions at time 0
    #time_grid[1] <- 0 # manually set first point to 0, instead of first observed time
  } else{
    time_grid <- sort(unique(time))
    time_grid <- c(0, time_grid)

  }

  ids <- seq(1:length(time))

  if (!is.null(obsWeights)){
    stackX <- as.matrix(data.frame(X,
                                   obsWeights = obsWeights,
                                   ids = ids))
  } else{
    stackX <- as.matrix(data.frame(X,
                                   ids = ids))
  }

  stacked <- stack_entry(time = time,
                         entry = entry,
                         X = stackX,
                         time_grid = time_grid,
                         time_basis = "continuous")

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

  if (fit$time_basis == "continuous"){
    get_stacked_pred <- function(t){
      new_stacked <- data.frame(t = t, newX)
      preds <- stats::predict(fit$reg.object, newdata=new_stacked)$pred
      return(preds)
    }

    predictions <- apply(X = matrix(newtimes), FUN = get_stacked_pred, MARGIN = 1)
  } else if (fit$time_basis == "dummy"){
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

  return(predictions)
}


#' Stacked binary regression with xgboost, using the cdf
#'
#' @param time Observed time
#' @param event Indicator of event (vs censoring)
#' @param X Covariate matrix
#' @param censored Logical, indicates whether to run regression on censored observations (vs uncensored)
#' @param bin_size Size of quantiles over which to make the stacking bins
#' @param isotonize Logical, indicates whether or not to isotonize cdf estimates using PAVA
#' @param time_basis How to treat time
#' @param xgb_control Control parameters for xgboost
#'
#' @return An object of class \code{f_y_stack_xgboost}
#' @noRd
f_w_stack_xgboost <- function(time,
                              event,
                              entry,
                              X,
                              censored,
                              bin_size,
                              isotonize = TRUE,
                              xgb_control,
                              time_basis = "continuous"){

  if (!is.null(censored)){
    if (censored == TRUE){
      time <- time[!as.logical(event)]
      entry <- entry[!as.logical(event)]
      X <- X[!as.logical(event),]
      obsWeights <- xgb_control$obsWeights[!as.logical(event)]
    } else if (censored == FALSE){
      time <- time[as.logical(event)]
      X <- X[as.logical(event),]
      entry <- entry[as.logical(event)]
      obsWeights <- xgb_control$obsWeights[as.logical(event)]
    }
  } else{
    time <- time
    entry <- entry
    X <- X
    obsWeights <- xgb_control$obsWeights
  }

  cv_folds <- split(sample(1:length(time)), rep(1:xgb_control$V, length = length(time)))

  X <- as.matrix(X)
  time <- as.matrix(time)
  entry <- as.matrix(entry)
  dat <- data.frame(X, time, entry)

  if (!is.null(bin_size)){
    #time_grid <- quantile(dat$time, probs = seq(0, 1, by = bin_size))
    time_grid <- sort(unique(stats::quantile(time, probs = seq(0, 1, by = bin_size))))
    time_grid <- c(0, time_grid) # 013123 changed this to try to get better predictions at time 0
    #time_grid[1] <- 0 # manually set first point to 0, instead of first observed time
  } else{
    time_grid <- sort(unique(time))
    time_grid <- c(0, time_grid)

  }

  ids <- seq(1:length(time))

  if (!is.null(obsWeights)){
    stackX <- as.matrix(data.frame(X,
                                   obsWeights = obsWeights,
                                   ids = ids))
  } else{
    stackX <- as.matrix(data.frame(X,
                                   ids = ids))
  }


  if (is.null(xgb_control$tuning_params)){
    xgb_control$tuning_params <- list(ntrees = c(50, 100, 250, 500), max_depth = c(1,2,3),
                                      eta = c(0.1), subsample = c(0.5))
  }
  if (is.null(xgb_control$objective)){
    xgb_control$objective <- "binary:logistic"
  }
  if (is.null(xgb_control$eval_metric)){
    xgb_control$eval_metric <- "logloss"
  }
  if (is.null(xgb_control$V)){
    xgb_control$V <- 10
  }

  if (is.null(xgb_control$tune)){
    xgb_control$tune <- TRUE
  }

  if (xgb_control$tune){
    param_grid <- expand.grid(ntrees = xgb_control$tuning_params$ntrees,
                              max_depth = xgb_control$tuning_params$max_depth,
                              eta = xgb_control$tuning_params$eta,
                              subsample = xgb_control$tuning_params$subsample)

    get_CV_risk <- function(i){
      ntrees <- param_grid$ntrees[i]
      max_depth <- param_grid$max_depth[i]
      eta <- param_grid$eta[i]
      subsample <- param_grid$subsample[i]
      risks <- rep(NA, xgb_control$V)
      for (j in 1:xgb_control$V){
        train_stackX <- stackX[-cv_folds[[j]],]
        train_time <- time[-cv_folds[[j]]]
        train_entry <- entry[-cv_folds[[j]]]
        train_stack <- stack_entry(time = train_time,
                                   entry = train_entry,
                                   X = train_stackX,
                                   time_grid = time_grid,
                                   time_basis = "continuous")

        # change t to dummy variable
        if (time_basis == "dummy"){
          train_stack$t <- factor(train_stack$t)
          dummy_mat <- stats::model.matrix(~-1 + t, data=train_stack)
          risk_set_names <- paste0("risk_set_", seq(1, (length(time_grid))))
          colnames(dummy_mat) <- risk_set_names
          train_stack$t <- NULL
          train_stack <- cbind(dummy_mat, train_stack)
        }

        # long_obsWeights <- train_stack$obsWeights
        # stacked_ids <- train_stack$ids
        train_stack$obsWeights <- NULL
        train_stack$ids <- NULL

        xgmat <- xgboost::xgb.DMatrix(data = as.matrix(train_stack[,-ncol(train_stack)]),
                                      label = as.matrix(train_stack[,ncol(train_stack)]))
        fit <- xgboost::xgboost(data = xgmat,
                                objective=xgb_control$objective,
                                nrounds = ntrees,
                                max_depth = max_depth,
                                eta = eta,
                                verbose = FALSE,
                                nthread = 1,
                                save_period = NULL,
                                eval_metric = xgb_control$eval_metric,
                                subsample = subsample)
        test_stackX <- stackX[cv_folds[[j]],]
        test_time <- time[cv_folds[[j]]]
        test_entry <- entry[cv_folds[[j]]]
        test_stack <- stack_entry(time = test_time,
                                  entry = test_entry,
                                  X = test_stackX,
                                  time_grid = time_grid,
                                  time_basis = "continuous")

        # change t to dummy variable
        if (time_basis == "dummy"){
          test_stack$t <- factor(test_stack$t)
          dummy_mat <- stats::model.matrix(~-1 + t, data=test_stack)
          risk_set_names <- paste0("risk_set_", seq(1, (length(time_grid))))
          colnames(dummy_mat) <- risk_set_names
          test_stack$t <- NULL
          test_stack <- cbind(dummy_mat, test_stack)
        }

        # long_obsWeights <- train_stack$obsWeights
        # stacked_ids <- train_stack$ids
        test_stack$obsWeights <- NULL
        test_stack$ids <- NULL

        preds <- stats::predict(fit, newdata = as.matrix(test_stack[,-ncol(test_stack)]))
        preds[preds == 1] <- 0.99 # this is a hack, but come back to it later
        truth <- test_stack[,ncol(test_stack)]
        log_loss <- lapply(1:length(preds), function(x) { # using log loss right now
          -truth[x] * log(preds[x]) - (1-truth[x])*log(1 - preds[x])
        })
        log_loss <- unlist(log_loss)
        risks[j] <- sum(log_loss)
      }
      return(sum(risks))
    }

    CV_risks <- unlist(lapply(1:nrow(param_grid), get_CV_risk))

    opt_param_index <- which.min(CV_risks)
    opt_ntrees <- param_grid$ntrees[opt_param_index]
    opt_max_depth <- param_grid$max_depth[opt_param_index]
    opt_eta <- param_grid$eta[opt_param_index]
    opt_subsample <- param_grid$subsample[opt_param_index]
    opt_params <- list(ntrees = opt_ntrees,
                       max_depth = opt_max_depth,
                       eta = opt_eta,
                       subsample = opt_subsample)

    CV_mat <- cbind(param_grid, CV_risks)
  } else{
    opt_params <- xgb_control$tuning_params
    CV_mat <- NA
  }



  stacked <- stack_entry(time = time,
                         entry = entry,
                         X = stackX,
                         time_grid = time_grid,
                         time_basis = "continuous")

  # change t to dummy variable
  if (time_basis == "dummy"){
    stacked$t <- factor(stacked$t)
    dummy_mat <- stats::model.matrix(~-1 + t, data=stacked)
    risk_set_names <- paste0("risk_set_", seq(1, (length(time_grid))))
    colnames(dummy_mat) <- risk_set_names
    stacked$t <- NULL
    stacked <- cbind(dummy_mat, stacked)
  }
  #
  # long_obsWeights <- stacked$obsWeights
  # stacked_ids <- stacked$ids
  stacked$obsWeights <- NULL
  stacked$ids <- NULL

  .Y <- stacked[,ncol(stacked)]
  .X <- as.matrix(stacked[,-ncol(stacked)])
  xgmat <- xgboost::xgb.DMatrix(data = .X, label = .Y)
  fit <- xgboost::xgboost(data = xgmat,
                          objective=xgb_control$objective,
                          nrounds = opt_ntrees,
                          max_depth = opt_max_depth,
                          eta = opt_eta,
                          verbose = FALSE,
                          nthread = 1,
                          save_period = NULL,
                          eval_metric = xgb_control$eval_metric,
                          subsample = opt_subsample)

  fit <- list(reg.object = fit,
              time_grid = time_grid,
              time_basis = time_basis,
              CV_mat = CV_mat)
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

  if (fit$time_basis == "continuous"){
    get_stacked_pred <- function(t){
      new_stacked <- as.matrix(data.frame(t = t, newX))
      preds <- stats::predict(fit$reg.object, newdata=new_stacked)
      return(preds)
    }

    predictions <- apply(X = matrix(newtimes), FUN = get_stacked_pred, MARGIN = 1)
  } else if (fit$time_basis == "dummy"){
    get_preds <- function(t){
      dummies <- matrix(0, ncol = length(time_grid), nrow = nrow(newX))
      index <- max(which(time_grid <= t))
      dummies[,index] <- 1
      new_stacked <- cbind(dummies, newX)
      risk_set_names <- paste0("risk_set_", seq(1, (length(time_grid))))
      colnames(new_stacked)[1:length(time_grid)] <- risk_set_names
      new_stacked <- as.matrix(data.frame(new_stacked))
      preds <- stats::predict(fit$reg.object, newdata=new_stacked)
      return(preds)
    }

    predictions <- apply(X = matrix(newtimes), FUN = get_preds, MARGIN = 1)
  }

  return(predictions)
}
