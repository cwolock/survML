#' Binary SuperLearner
#'
#' @param event Indicator of event (vs censoring)
#' @param X Covariate matrix
#' @param SL_control Super Learner control parameters
#'
#' @return An object of class \code{p_delta_SuperLearner}
#' @noRd
p_delta_SuperLearner <- function(event,
                                 X,
                                 SL_control){

  X <- as.data.frame(X)

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

  opt_fit <- SuperLearner::SuperLearner(Y = event,
                                        X = X,
                                        family = stats::binomial(),
                                        SL.library = SL_control$SL.library,
                                        method = SL_control$method,
                                        verbose = FALSE,
                                        cvControl = list(V = SL_control$V,
                                                         stratifyCV = SL_control$stratifyCV),
                                        obsWeights = SL_control$obsWeights)

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
predict.p_delta_SuperLearner <- function(fit,
                                         newX){
  X <- as.data.frame(newX)
  preds <- stats::predict(fit$reg.object, newdata = newX)$pred
  return(preds)
}


#' Binary xgboost regression with homemade cross validation
#'
#' @param event Indicator of event (vs censoring)
#' @param X Covariate matrix
#' @param xgb_control
#'
#' @return An object of class \code{f_y_stackCVcdf}
#' @noRd
p_delta_xgboost <- function(event, X, xgb_control){

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
    xgb_control$V <- 5
  }

  if (is.null(xgb_control$tune)){
    xgb_control$tune <- TRUE
  }


  cv_folds <- split(sample(1:length(event)), rep(1:xgb_control$V, length = length(event)))

  X <- as.matrix(X)
  event <- as.matrix(event)
  dat <- data.frame(X, event)

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
        train_X <- X[-cv_folds[[j]],]
        train_event <- event[-cv_folds[[j]]]
        xgmat <- xgboost::xgb.DMatrix(data = train_X, label = train_event)
        fit <- xgboost::xgboost(data = xgmat, objective=xgb_control$objective, nrounds = ntrees,
                                max_depth = max_depth, eta = eta,
                                verbose = FALSE, nthread = 1,
                                save_period = NULL, eval_metric = xgb_control$eval_metric,
                                subsample = subsample)
        test_X <- X[cv_folds[[j]],]
        test_event <- event[cv_folds[[j]]]
        preds <- stats::predict(fit, newdata = test_X)
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
    opt_subsample <- param_grid$subsample[opt_param_index]
    opt_params <- list(ntrees = opt_ntrees,
                       max_depth = opt_max_depth,
                       eta = opt_eta,
                       subsample = opt_subsample)

    CV_mat <- cbind(param_grid, CV_risks)
  } else{
    opt_params <- xgb_control$param_grid
    CV_mat <- NA
  }

  xgmat <- xgboost::xgb.DMatrix(data = X, label = event)
  fit <- xgboost::xgboost(data = xgmat, objective=xgb_control$objective, nrounds = opt_ntrees,
                          max_depth = opt_max_depth, eta = opt_eta,
                          verbose = FALSE,
                          save_period = NULL, eval_metric = xgb_control$eval_metric,
                          subsample = opt_subsample)

  fit <- list(reg.object = fit,
              CV_mat = CV_mat)
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
  preds <- stats::predict(fit$reg.object, newdata=X)
  return(preds)
}

