#' Estimate a conditional survival function via stacking
#'
#' @param time Observed time
#' @param event Indicator of event (vs censoring)
#' @param X Covariate matrix
#' @param newX Values of covariates at which to make a prediction
#' @param newtimes Times at which to make the survival function prediction
#' @param bin_size Quantiles on which to grid times. If NULL, defaults to every observed time
#' @param V Number of cross validation folds
#' @param time_basis How to treat time (continuous or dummy)
#' @param algorithm Binary classification algorithm to use
#' @param entry Study entry variable, if applicable
#' @param direction Prospective or retrospective study
#' @param tuning_params Tuning parameters to be passed to algorithm for tuning
#'
#' @return An object of class \code{survMLs}
#'
#' @export
#'
#' @examples
survMLs <- function(time,
                    event,
                    X,
                    newX,
                    newtimes,
                    bin_size = NULL,
                    V = 10,
                    time_basis = "continuous",
                    algorithm = "xgboost",
                    entry = NULL,
                    direction = "forward",
                    tuning_params = NULL){

  cts.num <- 4
  deg.gam <- 2
  if (direction == "prospective"){
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

    if (algorithm == "xgboost"){
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
      cv_folds <- split(sample(1:length(time)), rep(1:V, length = length(time)))

      get_CV_risk <- function(i){
        ntrees <- param_grid$ntrees[i]
        max_depth <- param_grid$max_depth[i]
        eta <- param_grid$eta[i]
        subsample <- param_grid$subsample[i]
        risks <- rep(NA, V)
        for (j in 1:V){
          train <- stacked[-cv_folds[[j]],]
          xgmat <- xgboost::xgb.DMatrix(data = as.matrix(train[,-ncol(train)]), label = as.matrix(train[,ncol(train)]))
          #ratio <- min(c(subsamp_size/nrow(train), 1))
          fit <- xgboost::xgboost(data = xgmat, objective="binary:logistic", nrounds = ntrees,
                                  max_depth = max_depth, eta = eta,
                                  verbose = FALSE, nthread = 1,
                                  save_period = NULL, eval_metric = "logloss",
                                  subsample = subsample)
          test <- as.matrix(stacked[cv_folds[[j]],])
          preds <- predict(fit, newdata = test[,-ncol(test)])
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

      if (time_basis == "continuous"){
        stacked <- survML:::stack_haz(time = time,
                                      event = event,
                                      X = X,
                                      time_grid = time_grid,
                                      entry = entry)
        # I guess for stacking, I can do cross validation on stacked dataset, rather than on individuals? shouldn't matter too
        # much I'd think

        CV_risks <- unlist(lapply(1:nrow(param_grid), get_CV_risk))

        opt_param_index <- which.min(CV_risks)
        opt_ntrees <- param_grid$ntrees[opt_param_index]
        opt_max_depth <- param_grid$max_depth[opt_param_index]
        opt_eta <- param_grid$eta[opt_param_index]
        opt_subsample <- param_grid$subsample[opt_param_index]
        opt_params <- list(ntrees = opt_ntrees, max_depth = opt_max_depth, eta = opt_eta, subsampe = opt_subsample)
        #stacked <- survML:::stack_haz(time = time, event = event, X = X, time_grid = time_grid)
        Y2 <- stacked[,ncol(stacked)]
        X2 <- as.matrix(stacked[,-ncol(stacked)])
        xgmat <- xgboost::xgb.DMatrix(data = X2, label = Y2)
        #ratio <- min(c(opt_subsamp_size/nrow(stacked), 1))
        fit <- xgboost::xgboost(data = xgmat, objective="binary:logistic", nrounds = opt_ntrees,
                                max_depth = opt_max_depth, eta = opt_eta,
                                verbose = FALSE, nthread = 1,
                                save_period = NULL, eval_metric = "logloss", subsample = opt_subsample)

        get_hazard_preds <- function(t){
          new_stacked <- as.matrix(data.frame(t = t, newX))
          preds <- predict(fit, newdata=new_stacked)
          return(preds)
        }

        hazard_preds <- apply(X = matrix(time_grid[-1]), FUN = get_hazard_preds, MARGIN = 1) # don't estimate hazard at t =0

        get_surv_preds <- function(t){
          if (sum(time_grid[-1] <= t) != 0){ # if you don't fall before the first time in the grid
            final_index <- max(which(time_grid[-1] <= t))
            haz <- as.matrix(hazard_preds[,1:final_index])
            anti_haz <- 1 - haz
            surv <- apply(anti_haz, MARGIN = 1, prod)
          } else{
            surv <- rep(1, nrow(hazard_preds))
          }
          return(surv)
        }

        surv_preds <- apply(X = matrix(newtimes), FUN = get_surv_preds, MARGIN = 1)
        print(CV_risks)
        print(fit$params)
        print(fit$niter)
      } else if (time_basis == "dummy"){
        # time dummies
        ncol_stacked <- ncol(X) + length(trunc_time_grid) + 1 # covariates, risk set dummies, binary outcome
        stacked <- matrix(NA, ncol = ncol_stacked, nrow = 1)
        for (i in 1:(length(trunc_time_grid))){
          if (is.null(entry)){
            risk_set <- dat[dat$time > time_grid[i],] # should this be <= rather than <??
          } else{
            risk_set <- dat[dat$time > time_grid[i] & entry <= time_grid[i+1],] # should this be <= rather than <??
          }
          risk_set_covariates <- risk_set[,1:ncol(X)]
          event_indicators <- matrix(ifelse(risk_set$time <= time_grid[i + 1 ] & risk_set$event == 1, 1, 0))
          dummies <- matrix(0, ncol = length(trunc_time_grid), nrow = nrow(risk_set))
          dummies[,i] <- 1
          newdata <- as.matrix(cbind(dummies, risk_set_covariates, event_indicators))
          stacked <- rbind(stacked, newdata)
        }

        stacked <- stacked[-1,]
        risk_set_names <- paste0("risk_set_", seq(1, (length(trunc_time_grid))))
        colnames(stacked)[1:(length(trunc_time_grid))] <- risk_set_names
        stacked <- data.frame(stacked)
        CV_risks <- unlist(lapply(1:nrow(param_grid), get_CV_risk))
        opt_param_index <- which.min(CV_risks)
        opt_ntrees <- param_grid$ntrees[opt_param_index]
        opt_max_depth <- param_grid$max_depth[opt_param_index]
        opt_eta <- param_grid$eta[opt_param_index]
        opt_params <- list(ntrees = opt_ntrees, max_depth = opt_max_depth, eta = opt_eta)
        #stacked <- survML:::stack_haz(time = time, event = event, X = X, time_grid = time_grid)
        Y1 <- stacked[,ncol(stacked)]
        X1 <- as.matrix(stacked[,-ncol(stacked)])
        xgmat <- xgboost::xgb.DMatrix(data = X1, label = Y1)
        fit <- xgboost::xgboost(data = xgmat, objective="binary:logistic", nrounds = opt_ntrees,
                                max_depth = opt_max_depth, eta = opt_eta,
                                verbose = FALSE, nthread = 1,
                                save_period = NULL, eval_metric = "logloss")
        get_hazard_preds <- function(t){
          dummies <- matrix(0, ncol = length(trunc_time_grid), nrow = nrow(newX))
          index <- max(which(trunc_time_grid <= t))
          dummies[,index] <- 1
          new_stacked <- cbind(dummies, newX)
          risk_set_names <- paste0("risk_set_", seq(1, (length(trunc_time_grid))))
          colnames(new_stacked)[1:length(trunc_time_grid)] <- risk_set_names
          new_stacked <- as.matrix(new_stacked)
          preds <- predict(fit, newdata=new_stacked)
          return(preds)
        }

        hazard_preds <- apply(X = matrix(trunc_time_grid), FUN = get_hazard_preds, MARGIN = 1)

        get_surv_preds <- function(t){
          final_index <- max(which(trunc_time_grid <= t))
          haz <- as.matrix(hazard_preds[,1:final_index])
          anti_haz <- 1 - haz
          surv <- apply(anti_haz, MARGIN = 1, prod)
          return(surv)
        }

        surv_preds <- apply(X = matrix(newtimes), FUN = get_surv_preds, MARGIN = 1)
        print(CV_risks)
        print(fit$params)
        print(fit$niter)
      }
    } else if (algorithm == "gam"){
      if (time_basis == "continuous"){
        # continuous time
        # make sure the two ways of stacking agree
        stacked <- survML:::stack_haz(time = time, event = event, X = X, time_grid = time_grid, entry = entry)
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

        get_hazard_preds <- function(t){
          new_stacked <- data.frame(t = t, newX)
          preds <- gam::predict.Gam(fit.gam, newdata=new_stacked, type = "response")
          return(preds)
        }

        hazard_preds <- apply(X = matrix(time_grid[-1]), FUN = get_hazard_preds, MARGIN = 1) # don't estimate hazard at t =0

        get_surv_preds <- function(t){
          if (sum(time_grid[-1] <= t) != 0){ # if you don't fall before the first time in the grid
            final_index <- max(which(time_grid[-1] <= t))
            haz <- as.matrix(hazard_preds[,1:final_index])
            anti_haz <- 1 - haz
            surv <- apply(anti_haz, MARGIN = 1, prod)
          } else{
            surv <- rep(1, nrow(hazard_preds))
          }
          return(surv)
        }

        surv_preds <- apply(X = matrix(newtimes), FUN = get_surv_preds, MARGIN = 1)
      } else if (time_basis == "dummy"){
        # # Super learner version of stacking
        # # dummies
        # # this version does dummy variables for time - include this as an option
        ncol_stacked <- ncol(X) + length(trunc_time_grid) + 1 # covariates, risk set dummies, binary outcome
        stacked <- matrix(NA, ncol = ncol_stacked, nrow = 1)
        for (i in 1:(length(trunc_time_grid))){
          if (is.null(entry)){
            risk_set <- dat[dat$time > time_grid[i],] # should this be <= rather than <??
          } else{
            risk_set <- dat[dat$time > time_grid[i] & entry <= time_grid[i+1],] # should this be <= rather than <??
          }
          risk_set_covariates <- risk_set[,1:ncol(X)]
          event_indicators <- matrix(ifelse(risk_set$time <= time_grid[i + 1 ] & risk_set$event == 1, 1, 0))
          dummies <- matrix(0, ncol = length(trunc_time_grid), nrow = nrow(risk_set))
          dummies[,i] <- 1
          newdata <- as.matrix(cbind(dummies, risk_set_covariates, event_indicators))
          stacked <- rbind(stacked, newdata)
        }

        stacked <- stacked[-1,]
        risk_set_names <- paste0("risk_set_", seq(1, (length(trunc_time_grid))))
        colnames(stacked)[1:(length(trunc_time_grid))] <- risk_set_names
        stacked <- data.frame(stacked)
        Y <- stacked$event_indicators
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


        get_hazard_preds <- function(t){
          dummies <- matrix(0, ncol = length(trunc_time_grid), nrow = nrow(newX))
          index <- max(which(trunc_time_grid <= t))
          dummies[,index] <- 1
          new_stacked <- cbind(dummies, newX)
          risk_set_names <- paste0("risk_set_", seq(1, (length(trunc_time_grid))))
          colnames(new_stacked)[1:length(trunc_time_grid)] <- risk_set_names
          preds <- gam::predict.Gam(fit.gam, newdata=new_stacked, type = "response")
          return(preds)
        }

        hazard_preds <- apply(X = matrix(trunc_time_grid), FUN = get_hazard_preds, MARGIN = 1)

        get_surv_preds <- function(t){
          final_index <- max(which(trunc_time_grid <= t))
          haz <- as.matrix(hazard_preds[,1:final_index])
          anti_haz <- 1 - haz
          surv <- apply(anti_haz, MARGIN = 1, prod)
          return(surv)
        }

        surv_preds <- apply(X = matrix(newtimes), FUN = get_surv_preds, MARGIN = 1)
      }
    } else if (algorithm == "earth"){

      degree = 2
      penalty = 3
      nk = max(21, 2*ncol(X) + 1)
      pmethod = "backward"
      nfold = 0
      ncross = 1
      minspan = 0
      endspan = 0

      if (time_basis == "continuous"){
        # continuous time

        # make sure the two ways of stacking agree
        stacked <- survML:::stack_haz(time = time, event = event, X = X, time_grid = time_grid, entry = entry)
        Y <- stacked[,ncol(stacked)]
        X <- stacked[,-ncol(stacked)]
        X <- as.data.frame(X)

        fit.earth <- earth::earth(x = X, y = Y, degree = degree, nk = nk, penalty = penalty, pmethod = pmethod,
                                  nfold = nfold, ncross = ncross, minspan = minspan, endspan = endspan,
                                  glm = list(family = binomial))

        get_hazard_preds <- function(t){
          new_stacked <- data.frame(t = t, newX)
          preds <- predict(fit.earth, newdata=new_stacked, type = "response")
          return(preds)
        }

        hazard_preds <- apply(X = matrix(time_grid[-1]), FUN = get_hazard_preds, MARGIN = 1) # don't estimate hazard at t =0

        get_surv_preds <- function(t){
          if (sum(time_grid[-1] <= t) != 0){ # if you don't fall before the first time in the grid
            final_index <- max(which(time_grid[-1] <= t))
            haz <- as.matrix(hazard_preds[,1:final_index])
            anti_haz <- 1 - haz
            surv <- apply(anti_haz, MARGIN = 1, prod)
          } else{
            surv <- rep(1, nrow(hazard_preds))
          }
          return(surv)
        }

        surv_preds <- apply(X = matrix(newtimes), FUN = get_surv_preds, MARGIN = 1)
      } else if (time_basis == "dummy"){
        # # Super learner version of stacking
        # # dummies
        # # this version does dummy variables for time - include this as an option
        ncol_stacked <- ncol(X) + length(trunc_time_grid) + 1 # covariates, risk set dummies, binary outcome
        stacked <- matrix(NA, ncol = ncol_stacked, nrow = 1)
        for (i in 1:(length(trunc_time_grid))){
          if (is.null(entry)){
            risk_set <- dat[dat$time > time_grid[i],] # should this be <= rather than <??
          } else{
            risk_set <- dat[dat$time > time_grid[i] & entry <= time_grid[i+1],] # should this be <= rather than <??
          }
          risk_set_covariates <- risk_set[,1:ncol(X)]
          event_indicators <- matrix(ifelse(risk_set$time <= time_grid[i + 1 ] & risk_set$event == 1, 1, 0))
          dummies <- matrix(0, ncol = length(trunc_time_grid), nrow = nrow(risk_set))
          dummies[,i] <- 1
          newdata <- as.matrix(cbind(dummies, risk_set_covariates, event_indicators))
          stacked <- rbind(stacked, newdata)
        }

        stacked <- stacked[-1,]
        risk_set_names <- paste0("risk_set_", seq(1, (length(trunc_time_grid))))
        colnames(stacked)[1:(length(trunc_time_grid))] <- risk_set_names
        stacked <- data.frame(stacked)
        Y <- stacked$event_indicators
        X <- stacked[,-ncol(stacked)]

        X <- as.data.frame(X)

        fit.earth <- earth::earth(x = X, y = Y, degree = degree, nk = nk, penalty = penalty, pmethod = pmethod,
                                  nfold = nfold, ncross = ncross, minspan = minspan, endspan = endspan,
                                  glm = list(family = binomial))

        get_hazard_preds <- function(t){
          dummies <- matrix(0, ncol = length(trunc_time_grid), nrow = nrow(newX))
          index <- max(which(trunc_time_grid <= t))
          dummies[,index] <- 1
          new_stacked <- cbind(dummies, newX)
          risk_set_names <- paste0("risk_set_", seq(1, (length(trunc_time_grid))))
          colnames(new_stacked)[1:length(trunc_time_grid)] <- risk_set_names
          preds <- predict(fit.earth, newdata=new_stacked, type = "response")
          return(preds)
        }

        hazard_preds <- apply(X = matrix(trunc_time_grid), FUN = get_hazard_preds, MARGIN = 1)

        get_surv_preds <- function(t){
          final_index <- max(which(trunc_time_grid <= t))
          haz <- as.matrix(hazard_preds[,1:final_index])
          anti_haz <- 1 - haz
          surv <- apply(anti_haz, MARGIN = 1, prod)
          return(surv)
        }

        surv_preds <- apply(X = matrix(newtimes), FUN = get_surv_preds, MARGIN = 1)
      }
    }
  } else if (direction == "retrospective"){
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

    trunc_time_grid <- time_grid[-length(time_grid)]

    if (algorithm == "xgboost"){
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
      cv_folds <- split(sample(1:length(time)), rep(1:V, length = length(time)))

      get_CV_risk <- function(i){
        ntrees <- param_grid$ntrees[i]
        max_depth <- param_grid$max_depth[i]
        eta <- param_grid$eta[i]
        subsample <- param_grid$subsample[i]
        risks <- rep(NA, V)
        for (j in 1:V){
          train <- stacked[-cv_folds[[j]],]
          xgmat <- xgboost::xgb.DMatrix(data = as.matrix(train[,-ncol(train)]), label = as.matrix(train[,ncol(train)]))
          fit <- xgboost::xgboost(data = xgmat, objective="binary:logistic", nrounds = ntrees,
                                  max_depth = max_depth, eta = eta,
                                  verbose = FALSE, nthread = 1,
                                  save_period = NULL, eval_metric = "logloss", subsample = subsample)
          test <- as.matrix(stacked[cv_folds[[j]],])
          preds <- predict(fit, newdata = test[,-ncol(test)])
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

      if (time_basis == "continuous"){
        stacked <- survML:::stack_haz(time = time, event = rep(1,length(time)),
                                      X = X, time_grid = time_grid, entry = entry, direction = "reverse")
        # I guess for stacking, I can do cross validation on stacked dataset, rather than on individuals? shouldn't matter too
        # much I'd think

        CV_risks <- unlist(lapply(1:nrow(param_grid), get_CV_risk))

        opt_param_index <- which.min(CV_risks)
        opt_ntrees <- param_grid$ntrees[opt_param_index]
        opt_max_depth <- param_grid$max_depth[opt_param_index]
        opt_eta <- param_grid$eta[opt_param_index]
        opt_subsample <- param_grid$subsample[opt_param_index]
        opt_params <- list(ntrees = opt_ntrees, max_depth = opt_max_depth, eta = opt_eta,
                           subsample = opt_subsample)
        #stacked <- survML:::stack_haz(time = time, event = event, X = X, time_grid = time_grid)
        Y2 <- stacked[,ncol(stacked)]
        X2 <- as.matrix(stacked[,-ncol(stacked)])
        xgmat <- xgboost::xgb.DMatrix(data = X2, label = Y2)
        fit <- xgboost::xgboost(data = xgmat, objective="binary:logistic", nrounds = opt_ntrees,
                                max_depth = opt_max_depth, eta = opt_eta,
                                verbose = FALSE, nthread = 1,
                                save_period = NULL, eval_metric = "logloss", subsample = opt_subsample)

        get_hazard_preds <- function(t){
          new_stacked <- as.matrix(data.frame(t = t, newX))
          preds <- predict(fit, newdata=new_stacked)
          return(preds)
        }


        hazard_preds <- apply(X = matrix(time_grid[-length(time_grid)]), FUN = get_hazard_preds, MARGIN = 1)
        # don't estimate hazard at final time point (it's 0)

        get_surv_preds <- function(t){
          if (sum(time_grid[-length(time_grid)] >= t) != 0){ # if you don't fall before the first time in the grid
            final_index <- min(which(time_grid[-length(time_grid)] >= t))
            haz <- as.matrix(hazard_preds[,final_index:(length(time_grid)-1)])
            anti_haz <- 1 - haz
            surv <- apply(anti_haz, MARGIN = 1, prod)
            surv <- 1 - surv
          } else{
            surv <- rep(0, nrow(hazard_preds))
          }
          return(surv)
        }

        surv_preds <- apply(X = matrix(newtimes), FUN = get_surv_preds, MARGIN = 1)
        print(CV_risks)
        print(fit$params)
        print(fit$niter)
      } else if (time_basis == "dummy"){
        # time dummies
        ncol_stacked <- ncol(X) + length(trunc_time_grid) + 1 # covariates, risk set dummies, binary outcome
        stacked <- matrix(NA, ncol = ncol_stacked, nrow = 1)
        for (i in 1:(length(trunc_time_grid))){
          if (is.null(entry)){
            risk_set <- dat[dat$time > time_grid[i],] # should this be <= rather than <??
          } else{
            risk_set <- dat[dat$time > time_grid[i] & entry <= time_grid[i+1],] # should this be <= rather than <??
          }
          risk_set_covariates <- risk_set[,1:ncol(X)]
          event_indicators <- matrix(ifelse(risk_set$time <= time_grid[i + 1 ] & risk_set$event == 1, 1, 0))
          dummies <- matrix(0, ncol = length(trunc_time_grid), nrow = nrow(risk_set))
          dummies[,i] <- 1
          newdata <- as.matrix(cbind(dummies, risk_set_covariates, event_indicators))
          stacked <- rbind(stacked, newdata)
        }

        stacked <- stacked[-1,]
        risk_set_names <- paste0("risk_set_", seq(1, (length(trunc_time_grid))))
        colnames(stacked)[1:(length(trunc_time_grid))] <- risk_set_names
        stacked <- data.frame(stacked)
        CV_risks <- unlist(lapply(1:nrow(param_grid), get_CV_risk))
        opt_param_index <- which.min(CV_risks)
        opt_ntrees <- param_grid$ntrees[opt_param_index]
        opt_max_depth <- param_grid$max_depth[opt_param_index]
        opt_eta <- param_grid$eta[opt_param_index]
        opt_params <- list(ntrees = opt_ntrees, max_depth = opt_max_depth, eta = opt_eta)
        #stacked <- survML:::stack_haz(time = time, event = event, X = X, time_grid = time_grid)
        Y1 <- stacked[,ncol(stacked)]
        X1 <- as.matrix(stacked[,-ncol(stacked)])
        xgmat <- xgboost::xgb.DMatrix(data = X1, label = Y1)
        fit <- xgboost::xgboost(data = xgmat, objective="binary:logistic", nrounds = opt_ntrees,
                                max_depth = opt_max_depth, eta = opt_eta,
                                verbose = FALSE, nthread = 1,
                                save_period = NULL, eval_metric = "logloss")
        get_hazard_preds <- function(t){
          dummies <- matrix(0, ncol = length(trunc_time_grid), nrow = nrow(newX))
          index <- max(which(trunc_time_grid <= t))
          dummies[,index] <- 1
          new_stacked <- cbind(dummies, newX)
          risk_set_names <- paste0("risk_set_", seq(1, (length(trunc_time_grid))))
          colnames(new_stacked)[1:length(trunc_time_grid)] <- risk_set_names
          new_stacked <- as.matrix(new_stacked)
          preds <- predict(fit, newdata=new_stacked)
          return(preds)
        }

        hazard_preds <- apply(X = matrix(trunc_time_grid), FUN = get_hazard_preds, MARGIN = 1)

        get_surv_preds <- function(t){
          final_index <- max(which(trunc_time_grid <= t))
          haz <- as.matrix(hazard_preds[,1:final_index])
          anti_haz <- 1 - haz
          surv <- apply(anti_haz, MARGIN = 1, prod)
          return(surv)
        }

        surv_preds <- apply(X = matrix(newtimes), FUN = get_surv_preds, MARGIN = 1)
        print(CV_risks)
        print(fit$params)
        print(fit$niter)
      }
    } else if (algorithm == "gam"){
      if (time_basis == "continuous"){
        # continuous time

        # make sure the two ways of stacking agree
        stacked <- survML:::stack_haz(time = time, event = event, X = X, time_grid = time_grid, entry = entry)
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

        get_hazard_preds <- function(t){
          new_stacked <- data.frame(t = t, newX)
          preds <- gam::predict.Gam(fit.gam, newdata=new_stacked, type = "response")
          return(preds)
        }

        hazard_preds <- apply(X = matrix(time_grid[-1]), FUN = get_hazard_preds, MARGIN = 1) # don't estimate hazard at t =0

        get_surv_preds <- function(t){
          if (sum(time_grid[-1] <= t) != 0){ # if you don't fall before the first time in the grid
            final_index <- max(which(time_grid[-1] <= t))
            haz <- as.matrix(hazard_preds[,1:final_index])
            anti_haz <- 1 - haz
            surv <- apply(anti_haz, MARGIN = 1, prod)
          } else{
            surv <- rep(1, nrow(hazard_preds))
          }
          return(surv)
        }

        surv_preds <- apply(X = matrix(newtimes), FUN = get_surv_preds, MARGIN = 1)
      } else if (time_basis == "dummy"){
        # # Super learner version of stacking
        # # dummies
        # # this version does dummy variables for time - include this as an option
        ncol_stacked <- ncol(X) + length(trunc_time_grid) + 1 # covariates, risk set dummies, binary outcome
        stacked <- matrix(NA, ncol = ncol_stacked, nrow = 1)
        for (i in 1:(length(trunc_time_grid))){
          if (is.null(entry)){
            risk_set <- dat[dat$time > time_grid[i],] # should this be <= rather than <??
          } else{
            risk_set <- dat[dat$time > time_grid[i] & entry <= time_grid[i+1],] # should this be <= rather than <??
          }
          risk_set_covariates <- risk_set[,1:ncol(X)]
          event_indicators <- matrix(ifelse(risk_set$time <= time_grid[i + 1 ] & risk_set$event == 1, 1, 0))
          dummies <- matrix(0, ncol = length(trunc_time_grid), nrow = nrow(risk_set))
          dummies[,i] <- 1
          newdata <- as.matrix(cbind(dummies, risk_set_covariates, event_indicators))
          stacked <- rbind(stacked, newdata)
        }

        stacked <- stacked[-1,]
        risk_set_names <- paste0("risk_set_", seq(1, (length(trunc_time_grid))))
        colnames(stacked)[1:(length(trunc_time_grid))] <- risk_set_names
        stacked <- data.frame(stacked)
        Y <- stacked$event_indicators
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


        get_hazard_preds <- function(t){
          dummies <- matrix(0, ncol = length(trunc_time_grid), nrow = nrow(newX))
          index <- max(which(trunc_time_grid <= t))
          dummies[,index] <- 1
          new_stacked <- cbind(dummies, newX)
          risk_set_names <- paste0("risk_set_", seq(1, (length(trunc_time_grid))))
          colnames(new_stacked)[1:length(trunc_time_grid)] <- risk_set_names
          preds <- gam::predict.Gam(fit.gam, newdata=new_stacked, type = "response")
          return(preds)
        }

        hazard_preds <- apply(X = matrix(trunc_time_grid), FUN = get_hazard_preds, MARGIN = 1)

        get_surv_preds <- function(t){
          final_index <- max(which(trunc_time_grid <= t))
          haz <- as.matrix(hazard_preds[,1:final_index])
          anti_haz <- 1 - haz
          surv <- apply(anti_haz, MARGIN = 1, prod)
          return(surv)
        }

        surv_preds <- apply(X = matrix(newtimes), FUN = get_surv_preds, MARGIN = 1)
      }
    } else if (algorithm == "earth"){

      degree = 2
      penalty = 3
      nk = max(21, 2*ncol(X) + 1)
      pmethod = "backward"
      nfold = 0
      ncross = 1
      minspan = 0
      endspan = 0

      if (time_basis == "continuous"){
        # continuous time

        # make sure the two ways of stacking agree
        stacked <- survML:::stack_haz(time = time, event = rep(1,length(time)),
                                      X = X, time_grid = time_grid, entry = entry, direction = "reverse")
        Y <- stacked[,ncol(stacked)]
        X <- stacked[,-ncol(stacked)]
        X <- as.data.frame(X)

        fit.earth <- earth::earth(x = X, y = Y, degree = degree, nk = nk, penalty = penalty, pmethod = pmethod,
                                  nfold = nfold, ncross = ncross, minspan = minspan, endspan = endspan,
                                  glm = list(family = binomial))

        get_hazard_preds <- function(t){
          new_stacked <- data.frame(t = t, newX)
          preds <- predict(fit.earth, newdata=new_stacked, type = "response")
          return(preds)
        }

        hazard_preds <- apply(X = matrix(time_grid[-length(time_grid)]), FUN = get_hazard_preds, MARGIN = 1)
        # don't estimate hazard at final time point (it's 0)

        get_surv_preds <- function(t){
          if (sum(time_grid[-length(time_grid)] >= t) != 0){ # if you don't fall before the first time in the grid
            final_index <- min(which(time_grid[-length(time_grid)] >= t))
            haz <- as.matrix(hazard_preds[,final_index:(length(time_grid)-1)])
            anti_haz <- 1 - haz
            surv <- apply(anti_haz, MARGIN = 1, prod)
            surv <- 1 - surv
          } else{
            surv <- rep(0, nrow(hazard_preds))
          }
          return(surv)
        }

        surv_preds <- apply(X = matrix(newtimes), FUN = get_surv_preds, MARGIN = 1)
      } else if (time_basis == "dummy"){
        # # Super learner version of stacking
        # # dummies
        # # this version does dummy variables for time - include this as an option
        ncol_stacked <- ncol(X) + length(trunc_time_grid) + 1 # covariates, risk set dummies, binary outcome
        stacked <- matrix(NA, ncol = ncol_stacked, nrow = 1)
        for (i in 1:(length(trunc_time_grid))){
          if (is.null(entry)){
            risk_set <- dat[dat$time > time_grid[i],] # should this be <= rather than <??
          } else{
            risk_set <- dat[dat$time > time_grid[i] & entry <= time_grid[i+1],] # should this be <= rather than <??
          }
          risk_set_covariates <- risk_set[,1:ncol(X)]
          event_indicators <- matrix(ifelse(risk_set$time <= time_grid[i + 1 ] & risk_set$event == 1, 1, 0))
          dummies <- matrix(0, ncol = length(trunc_time_grid), nrow = nrow(risk_set))
          dummies[,i] <- 1
          newdata <- as.matrix(cbind(dummies, risk_set_covariates, event_indicators))
          stacked <- rbind(stacked, newdata)
        }

        stacked <- stacked[-1,]
        risk_set_names <- paste0("risk_set_", seq(1, (length(trunc_time_grid))))
        colnames(stacked)[1:(length(trunc_time_grid))] <- risk_set_names
        stacked <- data.frame(stacked)
        Y <- stacked$event_indicators
        X <- stacked[,-ncol(stacked)]

        X <- as.data.frame(X)

        fit.earth <- earth::earth(x = X, y = Y, degree = degree, nk = nk, penalty = penalty, pmethod = pmethod,
                                  nfold = nfold, ncross = ncross, minspan = minspan, endspan = endspan,
                                  glm = list(family = binomial))

        get_hazard_preds <- function(t){
          dummies <- matrix(0, ncol = length(trunc_time_grid), nrow = nrow(newX))
          index <- max(which(trunc_time_grid <= t))
          dummies[,index] <- 1
          new_stacked <- cbind(dummies, newX)
          risk_set_names <- paste0("risk_set_", seq(1, (length(trunc_time_grid))))
          colnames(new_stacked)[1:length(trunc_time_grid)] <- risk_set_names
          preds <- predict(fit.earth, newdata=new_stacked, type = "response")
          return(preds)
        }

        hazard_preds <- apply(X = matrix(trunc_time_grid), FUN = get_hazard_preds, MARGIN = 1)

        get_surv_preds <- function(t){
          final_index <- max(which(trunc_time_grid <= t))
          haz <- as.matrix(hazard_preds[,1:final_index])
          anti_haz <- 1 - haz
          surv <- apply(anti_haz, MARGIN = 1, prod)
          return(surv)
        }

        surv_preds <- apply(X = matrix(newtimes), FUN = get_surv_preds, MARGIN = 1)
      }
    }
  }

  res <- list(S_T_preds = surv_preds)
  class(res) <- "survMLs"
  return(res)
}
