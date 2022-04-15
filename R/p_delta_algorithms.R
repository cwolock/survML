#' Binary SuperLearner
#'
#' @param event Indicator of event (vs censoring)
#' @param X Covariate matrix
#' @param SL.library SuperLearner library
#' @param V number of CV folds
#'
#' @return An object of class \code{p_delta_SuperLearner}
#' @noRd
p_delta_SuperLearner <- function(event, X, SL.library, V = 10){

  X <- as.data.frame(X)
  opt_fit <- SuperLearner::SuperLearner(Y = event,
                                        X = X,
                                        family = "binomial",
                                        SL.library = SL.library,
                                        method = "method.NNLS",
                                        verbose = FALSE)

  fit <- list(reg.object = opt_fit)
  class(fit) <- c("p_delta_SuperLearner")
  return(fit)
}

#' Prediction function for p delta SuperLearner
#'
#' @param fit Fitted regression object
#' @param newX Values of covariates at which to make a prediction
#'
#' @return Matrix of predictions
#' @noRd
predict.p_delta_SuperLearner <- function(fit, newX){
  X <- as.data.frame(newX)
  preds <- predict(fit$reg.object, newdata = newX)$pred
  return(preds)
}

#' Binary xgboost regression with homemade cross validation
#'
#' @param event Indicator of event (vs censoring)
#' @param X Covariate matrix
#' @param V Number of CV folds
#'
#' @return An object of class \code{f_y_stackCVcdf}
#' @noRd
p_delta_xgboost <- function(event, X, V = 10, tuning_params = NULL){

  cv_folds <- split(sample(1:length(event)), rep(1:V, length = length(event)))

  X <- as.matrix(X)
  event <- as.matrix(event)
  dat <- data.frame(X, event)

  if (is.null(tuning_params)){
    tune <- list(ntrees = c(50, 100, 250, 500), max_depth = c(1,2,3),
                 eta = c(0.1))
  } else{
    tune <- tuning_params
  }

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
      train_event <- event[-cv_folds[[j]]]
      xgmat <- xgboost::xgb.DMatrix(data = train_X, label = train_event)
      fit <- xgboost::xgboost(data = xgmat, objective="binary:logistic", nrounds = ntrees,
                              max_depth = max_depth, eta = eta,
                              verbose = FALSE, nthread = 1,
                              save_period = NULL, eval_metric = "logloss")
      test_X <- X[cv_folds[[j]],]
      test_event <- event[cv_folds[[j]]]
      preds <- predict(fit, newdata = test_X)
      preds[preds == 1] <- 0.99 # this is a hack, but come back to it later
      truth <- test_event
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
  xgmat <- xgboost::xgb.DMatrix(data = X, label = event)
  fit <- xgboost::xgboost(data = xgmat, objective="binary:logistic", nrounds = opt_ntrees,
                          max_depth = opt_max_depth, eta = opt_eta,
                          verbose = FALSE,
                          save_period = NULL, eval_metric = "logloss")

  fit <- list(reg.object = fit)
  class(fit) <- c("p_delta_xgboost")
  return(fit)
}

#' Prediction function for p delta xgboost
#'
#' @param fit Fitted regression object
#' @param newX Values of covariates at which to make a prediction
#'
#' @return Matrix of predictions
#' @noRd
predict.p_delta_xgboost <- function(fit, newX){
  X <- as.matrix(newX)
  preds <- predict(fit$reg.object, newdata=X)
  return(preds)
}
