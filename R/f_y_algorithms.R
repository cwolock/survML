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


  tune = list(ntrees = c(1000, 5000, 10000), max_depth = c(2,3, 4), minobspernode = 10,
              shrinkage = c(0.01, 0.1, 0.5))
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
