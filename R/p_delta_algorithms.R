#' Binary xgboost regression with homemade cross validation
#'
#' @param event Indicator of event (vs censoring)
#' @param X Covariate matrix
#' @param V Number of CV folds
#'
#' @return An object of class \code{f_y_stackCVcdf}
#' @noRd
p_delta_xgboost <- function(event, X, V){

  cv_folds <- split(sample(1:length(event)), rep(1:V, length = length(event)))

  X <- as.matrix(X)
  event <- as.matrix(event)
  dat <- data.frame(X, event)

  tune <- list(ntrees = c(250, 500, 1000, 2500), max_depth = c(1,2,3),
               eta = c(0.01))

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
                          verbose = FALSE, nthread = 1,
                          save_period = NULL, eval_metric = "logloss")

  print(CV_risks)
  print(fit$params)
  print(fit$niter)
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

#' Binary xgboost regression with homemade cross validation
#'
#' @param event Indicator of event (vs censoring)
#' @param X Covariate matrix
#' @param V Number of CV folds
#'
#' @return An object of class \code{f_y_stackCVcdf}
#' @noRd
p_delta_ranger <- function(event, X, V){

  cv_folds <- split(sample(1:length(event)), rep(1:V, length = length(event)))

  X <- as.matrix(X)
  event <- as.matrix(event)
  dat <- data.frame(X, event)

  tune <- list(num.trees = c(250, 500, 1000, 2000), max.depth = c(1,2,3,4), mtry = c(1,2))

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
      train_event <- event[-cv_folds[[j]]]
      train <- data.frame(X = train_X, event = train_event)
      fit <- ranger::ranger(formula = event ~ .,
                            data = train,
                            num.trees = num.trees,
                            max.depth = max.depth,
                            mtry = mtry,
                            probability = TRUE)
      test_X <- X[cv_folds[[j]],]
      test_event <- event[cv_folds[[j]]]
      test <- data.frame(X = test_X, event = test_event)
      preds <- predict(fit, data = test[,-ncol(test)])$predictions
      preds <- preds[,2] # I think this is choosing the correct column but the output is not labeled...
      preds[preds == 1] <- 0.99 # this is a hack, but come back to it later
      truth <- test[,ncol(test)]
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
  fit <- ranger::ranger(formula = event ~ .,
                        data = dat,
                        num.trees = opt_num.trees,
                        max.depth = opt_max.depth,
                        mtry = opt_mtry,
                        probability = TRUE)

  print(CV_risks)
  print(opt_max.depth)
  print(fit$num.trees)
  print(opt_mtry)

  fit <- list(reg.object = fit)
  class(fit) <- c("p_delta_ranger")
  return(fit)
}

#' Prediction function for p delta xgboost
#'
#' @param fit Fitted regression object
#' @param newX Values of covariates at which to make a prediction
#'
#' @return Matrix of predictions
#' @noRd
predict.p_delta_ranger <- function(fit, newX){
  X <- as.matrix(newX)
  preds <- predict(fit$reg.object, data=X)$predictions[,2]
  return(preds)
}

#' Nadaraya-Watson estimator
#'
#' @param event Event indicator
#' @param X Covariates
#' @param bw bw
#' @param kernel_type Kernel_type
#' @param kernel_order Kernel_order
#'
#' @return An object of class \code{p_delta_nw}
#' @noRd
p_delta_nw <- function(event, X, bw, kernel_type = "gaussian", kernel_order = 2){
  fmla <- as.formula(paste("event ~ ", paste(colnames(X), collapse = "+")))
  dat <- data.frame(cbind(event = event, X))
  fit.nw <- np::npreg(fmla,
                   data = dat,
                   bws = rep(bw, ncol(X)),
                   regtype = "lc",
                   ckertype = kernel_type,
                   ckerorder = kernel_order)
  fit <- list(reg.object = fit.nw)
  class(fit) <- c("p_delta_nw")
  return(fit)
}

#' Prediction function for Nadaraya-Watson estimator
#'
#' @param fit Fitted regression object
#' @param newX Values of covariates at which to make a prediction
#'
#' @return Matrix of predictions
#' @noRd
predict.p_delta_nw <- function(fit, newX){
  predictions <- predict(fit$reg.object, newdata = newX)
  return(predictions)
}

#' Logistic regression
#'
#' @param event Event indicator
#' @param X Covariates
#'
#' @return An object of class \code{p_delta_logit_reg}
#' @noRd
p_delta_logit_reg <- function(event, X){
  fit.logit_reg <- stats::glm(event ~ .,
                    family = binomial(link = "logit"),
                    data = cbind(event = event, X))
  fit <- list(reg.object = fit.logit_reg)
  class(fit) <- c("p_delta_logit_reg")
  return(fit)
}

#' Prediction function for logistic regression
#'
#' @param fit Fitted regression object
#' @param newX Values of covariates at which to make a prediction
#'
#' @return Matrix of predictions
#' @noRd
predict.p_delta_logit_reg <- function(fit, newX){
  predictions <- predict(fit$reg.object, newdata = newX)
  predictions <- apply(X = as.matrix(predictions),
                       MARGIN = 1,
                       FUN = function(x) exp(x)/(1 + exp(x)))
  return(predictions)
}

#' Mean
#'
#' @param event Event indicator
#' @param X Covariates
#'
#' @return An object of class \code{p_delta_mean}
#' @noRd
p_delta_mean <- function(event, X){
  fit.mean <- stats::lm(event ~ 1,
                        data = cbind(event = event, X))
  fit <- list(reg.object = fit.mean)
  class(fit) <- c("p_delta_mean")
  return(fit)
}

#' Prediction function for mean
#'
#' @param fit Fitted regression object
#' @param newX Values of covariates at which to make a prediction
#'
#' @return Matrix of predictions
#' @noRd
predict.p_delta_mean <- function(fit, newX){
  predictions <- predict(fit$reg.object, newdata = newX)
  return(predictions)
}

#' Ranger random forest
#'
#' @param event Event indicator
#' @param X Covariates
#' @param mtry Number of varaibles sampled as candidates as each split
#' @param ntree Number of trees to grow
#'
#' @return An object of class \code{p_delta_ranger}
#' @noRd
# p_delta_ranger <- function(event, X, mtry = floor(sqrt(ncol(X))), num.trees = 500){
#
#   event <- as.factor(event)
#
#   # Ranger does not seem to work with X as a matrix, so we explicitly convert to
#   # data.frame rather than cbind. newX can remain as-is though.
#   if (is.matrix(X)) {
#     X <- data.frame(X)
#   }
#
#   fit <- ranger::ranger(event ~ .,
#                         data = cbind(event = event, X),
#                         num.trees = num.trees,
#                         mtry = mtry,
#                         min.node.size = 1,
#                         replace = TRUE,
#                         sample.fraction = 1,
#                         write.forest = TRUE,
#                         probability = TRUE,
#                         #num.threads = num.threads,
#                         verbose = TRUE)
#
#   fit <- list(reg.object = fit)
#   class(fit) <- c("p_delta_ranger")
#   return(fit)
# }

#' Prediction function for ranger
#'
#' @param fit Fitted regression object
#' @param newX Values of covariates at which to make a prediction
#'
#' @return Matrix of predictions
#' @noRd
# predict.p_delta_ranger <- function(fit, newX){
#   pred <- predict(fit$reg.object, data = newX)$predictions
#   pred <- pred[, "1"]
#   return(pred)
# }

#' Random forest
#'
#' @param event Event indicator
#' @param X Covariates
#' @param mtry Number of variables sampled as candidates as each split
#' @param ntree Number of trees to grow
#'
#' @return An object of class \code{p_delta_rf}
#' @noRd
p_delta_rf <- function(event, X, mtry = floor(sqrt(ncol(X))), ntree = 500){

  event <- as.factor(event)

  fit.rf <- randomForest::randomForest(y = event,
                                       x = X,
                                       ntree = ntree,
                                       keep.forest = TRUE,
                                       mtry = mtry,
                                       nodesize = 1,
                                       importance = FALSE)
  fit <- list(reg.object = fit.rf)
  class(fit) <- c("p_delta_rf")
  return(fit)
}

#' Prediction function for random forest
#'
#' @param fit Fitted regression object
#' @param newX Values of covariates at which to make a prediction
#'
#' @return Matrix of predictions
#' @noRd
predict.p_delta_rf <- function(fit, newX){
  pred <- predict(fit$reg.object, newdata = newX, type = 'vote')[,2]
  return(pred)
}
