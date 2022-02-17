#' Estimate a conditional survival function via stacking
#'
#' @param time Observed time
#' @param event Indicator of event (vs censoring)
#' @param X Covariate matrix
#' @param newX Values of covariates at which to make a prediction
#' @param newtimes Times at which to make the survival function prediction
#' @param bin_size Quantiles on which to grid times. If NULL, defaults to every observed time
#' @param SL.library SuperLearner library
#'
#' @return An object of class \code{stack}
#'
#' @export
#'
#' @examples
stackSurv <- function(time,
                      event,
                      X,
                      newX,
                      newtimes,
                      bin_size = NULL,
                      test_time = NULL,
                      test_event = NULL,
                      test_X = NULL,
                      SL.library = NULL,
                      V = 10){

  X <- as.matrix(X)
  time <- as.matrix(time)
  event <- as.matrix(event)
  dat <- data.frame(X, time, event)

  # if user gives bin size, set time grid based on quantiles. otherwise, every observed time
  if (!is.null(bin_size)){
    time_grid <- quantile(dat$time, probs = seq(0, 1, by = bin_size))
    time_grid[1] <- 0 # manually set first point to 0, instead of first observed time
  } else{
    time_grid <- sort(unique(time))
    time_grid <- c(0, time_grid)
  }

  trunc_time_grid <- time_grid[-length(time_grid)]

  tune <- list(ntrees = c(250, 500, 1000, 2500), max_depth = c(1,2,3),
               eta = c(0.01))

  param_grid <- expand.grid(ntrees = tune$ntrees,
                            max_depth = tune$max_depth,
                            eta = tune$eta)
  cv_folds <- split(sample(1:length(time)), rep(1:V, length = length(time)))
  stacked <- conSurv:::stack_haz(time = time, event = event, X = X, time_grid = time_grid)
  # I guess for stacking, I can do cross validation on stacked dataset, rather than on individuals? shouldn't matter too
  # much I'd think
  get_CV_risk <- function(i){
    ntrees <- param_grid$ntrees[i]
    max_depth <- param_grid$max_depth[i]
    eta <- param_grid$eta[i]
    risks <- rep(NA, V)
    for (j in 1:V){
      train <- stacked[-cv_folds[[j]],]
      xgmat <- xgboost::xgb.DMatrix(data = train[,-ncol(train)], label = train[,ncol(train)])
      fit <- xgboost::xgboost(data = xgmat, objective="binary:logistic", nrounds = ntrees,
                              max_depth = max_depth, eta = eta,
                              verbose = FALSE, nthread = 1,
                              save_period = NULL, eval_metric = "logloss")
      test <- stacked[cv_folds[[j]],]
      preds <- predict(fit, newdata = test[,-ncol(test)])
      preds[preds == 1] <- 0.99 # this is a hack, but come back to it later
      truth <- test[,ncol(test)]
      log_loss <- lapply(1:length(preds), function(x) { # using log loss right now
        -truth[x] * log(preds[x]) - (1-truth[x])*log(1 - preds[x])
      })
      log_loss <- unlist(log_loss)
      sum_log_loss <- sum(log_loss)
      risks[j] <- sum_log_loss


      # train_X <- X[-cv_folds[[j]],]
      # train_time <- time[-cv_folds[[j]]]
      # train_event <- event[-cv_folds[[j]]]
      # train_stack <- conSurv:::stack_haz(time = train_time, event = train_event, X = train_X, time_grid = time_grid)
      # xgmat <- xgboost::xgb.DMatrix(data = train_stack[,-ncol(train_stack)], label = train_stack[,ncol(train_stack)])
      # fit <- xgboost::xgboost(data = xgmat, objective="binary:logistic", nrounds = ntrees,
      #                         max_depth = max_depth, eta = eta,
      #                         verbose = FALSE, nthread = 1,
      #                         save_period = NULL, eval_metric = "logloss")
      # test_X <- X[cv_folds[[j]],]
      # test_time <- time[cv_folds[[j]]]
      # test_event <- event[cv_folds[[j]]]
      # test_stack <- conSurv:::stack_haz(time = test_time, event = test_event, X = test_X, time_grid = time_grid)
      # preds <- predict(fit, newdata = test_stack[,-ncol(test_stack)])
      # preds[preds == 1] <- 0.99 # this is a hack, but come back to it later
      # truth <- test_stack[,ncol(test_stack)]
      # log_loss <- lapply(1:length(preds), function(x) { # using log loss right now
      #   -truth[x] * log(preds[x]) - (1-truth[x])*log(1 - preds[x])
      # })
      # log_loss <- unlist(log_loss)
      # sum_log_loss <- sum(log_loss)
      # risks[j] <- sum_log_loss
    }
    return(sum(risks))
  }

  CV_risks <- unlist(lapply(1:nrow(param_grid), get_CV_risk))

  opt_param_index <- which.min(CV_risks)
  opt_ntrees <- param_grid$ntrees[opt_param_index]
  opt_max_depth <- param_grid$max_depth[opt_param_index]
  opt_eta <- param_grid$eta[opt_param_index]
  opt_params <- list(ntrees = opt_ntrees, max_depth = opt_max_depth, eta = opt_eta)
  #stacked <- conSurv:::stack_haz(time = time, event = event, X = X, time_grid = time_grid)
  Y <- stacked[,ncol(stacked)]
  X <- as.matrix(stacked[,-ncol(stacked)])
  xgmat <- xgboost::xgb.DMatrix(data = X, label = Y)
  fit <- xgboost::xgboost(data = xgmat, objective="binary:logistic", nrounds = opt_ntrees,
                          max_depth = opt_max_depth, eta = opt_eta,
                          verbose = FALSE, nthread = 1,
                          save_period = NULL, eval_metric = "logloss")

  get_hazard_preds <- function(t){
    new_stacked <- as.matrix(data.frame(t = t, newX))
    preds <- predict(fit, newdata=new_stacked)
    return(preds)
  }

  hazard_preds <- apply(X = matrix(newtimes), FUN = get_hazard_preds, MARGIN = 1)

  get_surv_preds <- function(t){
    final_index <- max(which(newtimes <= t))
    haz <- as.matrix(hazard_preds[,1:final_index])
    anti_haz <- 1 - haz
    surv <- apply(anti_haz, MARGIN = 1, prod)
    return(surv)
  }

  surv_preds <- apply(X = matrix(newtimes), FUN = get_surv_preds, MARGIN = 1)
  print(CV_risks)
  print(fit$params)
  print(fit$niter)

  # this version does dummy variables for time - include this as an option
  # ncol_stacked <- ncol(X) + length(trunc_time_grid) + 1 # covariates, risk set dummies, binary outcome
  # stacked <- matrix(NA, ncol = ncol_stacked, nrow = 1)
  # for (i in 1:(length(trunc_time_grid))){
  #   risk_set <- dat[dat$time > time_grid[i],] # should this be <= rather than <??
  #   risk_set_covariates <- risk_set[,1:ncol(X)]
  #   event_indicators <- matrix(ifelse(risk_set$time <= time_grid[i + 1 ] & risk_set$event == 1, 1, 0))
  #   dummies <- matrix(0, ncol = length(trunc_time_grid), nrow = nrow(risk_set))
  #   dummies[,i] <- 1
  #   newdata <- as.matrix(cbind(dummies, risk_set_covariates, event_indicators))
  #   stacked <- rbind(stacked, newdata)
  # }
  #
  # stacked <- stacked[-1,]
  # risk_set_names <- paste0("risk_set_", seq(1, (length(trunc_time_grid))))
  # colnames(stacked)[1:(length(trunc_time_grid))] <- risk_set_names
  # stacked <- data.frame(stacked)
  # Y <- stacked$event_indicators
  # X <- stacked[,-ncol(stacked)]
  #
  # fit <- SuperLearner::SuperLearner(Y = Y,
  #                                   X = X,
  #                                   family = binomial(),
  #                                   SL.library = SL.library,
  #                                   method = "method.NNLS",
  #                                   verbose = FALSE)
  #
  #
  # get_stacked_pred <- function(t){
  #   n_bins <- sum(t > trunc_time_grid)
  #   new_stacked <- matrix(rep(t(newX), n_bins), ncol = ncol(newX) , byrow = TRUE)
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

  #predictions <- apply(X = matrix(newtimes), FUN = get_stacked_pred, MARGIN = 1)

  res <- list(S_T_preds = surv_preds,
              fit = fit)
  class(res) <- "stackSurv"
  return(res)
}
