#' Empirical cdf
#'
#' @param time Observed time
#' @param event Indicator of event (vs censoring)
#' @param X Covariate matrix
#' @param censored Logical, indicates whether to run regression on censored observations (vs uncensored)
#'
#' @return An object of class \code{f_y_ecdf}
#' @noRd
f_y_ecdf <- function(time, event, X, censored){
  if (censored){
    time <- time[!as.logical(event)]
  } else{
    time <- time[as.logical(event)]
  }
  ecdf <- function(t){
    return(mean(time <= t))
  }
  time <- sort(time)
  estimated_cdf <- apply(X = as.matrix(time),
                         MARGIN = 1,
                         FUN = ecdf)
  fit.ecdf <- list(time = time, ecdf = estimated_cdf)
  fit <- list(reg.object = fit.ecdf)
  class(fit) <- c("f_y_ecdf")
  return(fit)
}

#' Prediction function for empirical cdf estimator
#'
#' @param fit Fitted regression object
#' @param newX Values of covariates at which to make a prediction
#' @param newtimes Times at which to make the cdf prediction
#'
#' @return Matrix of predictions
#' @noRd
predict.f_y_ecdf <- function(fit, newX, newtimes){
  predict_ecdf <- function(t){
    return(fit$reg.object$ecdf[max(which(fit$reg.object$time <= t))])
  }
  predictions <- apply(X = as.matrix(newtimes),
                       MARGIN = 1,
                       FUN = predict_ecdf)
  print(predictions)
  predictions <- t(replicate(nrow(newX), predictions))
  return(predictions)
}

#' Exponential regression
#'
#' @param time Observed time
#' @param event Indicator of event (vs censoring)
#' @param X Covariate matrix
#' @param censored Logical, indicates whether to run regression on censored observations (vs uncensored)
#'
#' @return An object of class \code{f_y_expreg}
#' @noRd
f_y_expreg <- function(time, event, X, censored) {
  if (censored){
    time <- time[!as.logical(event)]
    X <- X[!as.logical(event),]
  } else{
    time <- time[as.logical(event)]
    X <- X[as.logical(event),]
  }
  event <- rep(1, length(time))
  fit.expreg <- survival::survreg(survival::Surv(time, event) ~ .,
                                  data = cbind(time = time, event = event, X),
                                  dist = 'exponential')
  fit <- list(reg.object = fit.expreg)
  class(fit) <- c("f_y_expreg")
  return(fit)
}

#' Prediction function for exponential regression estimator
#'
#' @param fit Fitted regression object
#' @param newX Values of covariates at which to make a prediction
#' @param newtimes Times at which to make the cdf prediction
#'
#' @return Matrix of predictions
#' @noRd

predict.f_y_expreg <- function(object, newX, newtimes) {

  predictions <- predict(object$reg.object, newdata = newX, type = "quantile", p = seq(0, .999, by=.001))
  if (nrow(newX) == 1){ # problem here is, if it's just 1 row, when you do as.matrix you have to then transpose
    predictions <- t(as.matrix(predictions))
  }
  predictions <- t(sapply(1:nrow(predictions), function(j) {
    (1-stats::approx(predictions[j,], seq(0, .999, by=.001), xout = newtimes, method = "linear", rule = 2)$y)
  }))

  return(predictions)
}

#' Weibull regression
#'
#' @param time Observed time
#' @param event Indicator of event (vs censoring)
#' @param X Covariate matrix
#' @param censored Logical, indicates whether to run regression on censored observations (vs uncensored)
#'
#' @return An object of class \code{f_y_weibreg}
#' @noRd
f_y_weibreg <- function(time, event, X, censored) {
  if (censored){
    time <- time[!as.logical(event)]
    X <- X[!as.logical(event),]
  } else{
    time <- time[as.logical(event)]
    X <- X[as.logical(event),]
  }
  event <- rep(1, length(time))
  fit.expreg <- survival::survreg(survival::Surv(time, event) ~ .,
                                  data = cbind(time = time, event = event, X),
                                  dist = 'weibull')
  fit <- list(reg.object = fit.expreg)
  class(fit) <- c("f_y_weibreg")
  return(fit)
}

#' Prediction function for Weibull regression estimator
#'
#' @param fit Fitted regression object
#' @param newX Values of covariates at which to make a prediction
#' @param newtimes Times at which to make the cdf prediction
#'
#' @return Matrix of predictions
#' @noRd
predict.f_y_weibreg <- function(object, newX, newtimes) {

  predictions <- predict(object$reg.object, newdata = newX, type = "quantile", p = seq(0, .999, by=.001))
  if (nrow(newX) == 1){ # problem here is, if it's just 1 row, when you do as.matrix you have to then transpose
    predictions <- t(as.matrix(predictions))
  }
  predictions <- t(sapply(1:nrow(predictions), function(j) {
    (1-stats::approx(predictions[j,], seq(0, .999, by=.001), xout = newtimes, method = "linear", rule = 2)$y)
  }))

  return(predictions)
}

#' Log logistic regression
#'
#' @param time Observed time
#' @param event Indicator of event (vs censoring)
#' @param X Covariate matrix
#' @param censored Logical, indicates whether to run regression on censored observations (vs uncensored)
#'
#' @return An object of class \code{f_y_loglogreg}
#' @noRd
f_y_loglogreg <- function(time, event, X, censored) {
  if (censored){
    time <- time[!as.logical(event)]
    X <- X[!as.logical(event),]
  } else{
    time <- time[as.logical(event)]
    X <- X[as.logical(event),]
  }
  event <- rep(1, length(time))
  fit.expreg <- survival::survreg(survival::Surv(time, event) ~ .,
                                  data = cbind(time = time, event = event, X),
                                  dist = 'loglogistic')
  fit <- list(reg.object = fit.expreg)
  class(fit) <- c("f_y_weibreg")
  return(fit)
}

#' Prediction function for log logistic regression estimator
#'
#' @param fit Fitted regression object
#' @param newX Values of covariates at which to make a prediction
#' @param newtimes Times at which to make the cdf prediction
#'
#' @return Matrix of predictions
#' @noRd
predict.f_y_loglogreg <- function(object, newX, newtimes) {

  predictions <- predict(object$reg.object, newdata = newX, type = "quantile", p = seq(0, .999, by=.001))
  if (nrow(newX) == 1){ # problem here is, if it's just 1 row, when you do as.matrix you have to then transpose
    predictions <- t(as.matrix(predictions))
  }
  predictions <- t(sapply(1:nrow(predictions), function(j) {
    (1-stats::approx(predictions[j,], seq(0, .999, by=.001), xout = newtimes, method = "linear", rule = 2)$y)
  }))

  return(predictions)
}

#' Nadaraya-Watson estimator
#'
#' @param time Observed time
#' @param event Indicator of event (vs censoring)
#' @param X Covariate matrix
#' @param censored Logical, indicates whether to run regression on censored observations (vs uncensored)
#' @param time_grid Grid on which to fit the NW estimator
#' @param bw bw
#' @param kernel_type Kernel_type
#' @param kernel_order Kernel_order
#'
#' @return An object of class \code{f_y_nw}
#' @noRd
f_y_nw <- function(time, event, X, censored, time_grid, bw, kernel_type = "gaussian", kernel_order = 2){

  if (censored){
    time <- time[!as.logical(event)]
    X <- X[!as.logical(event),]
  } else{
    time <- time[as.logical(event)]
    X <- X[as.logical(event),]
  }

  fit_one <- function(t){
    fmla <- as.formula(paste("ind ~ ", paste(colnames(X), collapse = "+")))
    ind <- ifelse(time <= t, 1, 0)
    dat <- data.frame(cbind(ind = ind, X))
    nw_fit <- np::npreg(fmla,
                        data = dat,
                        bws = rep(bw, ncol(X)),
                        regtype = "lc",
                        ckertype = kernel_type,
                        ckerorder = kernel_order)
    return(nw_fit)
  }

  fits <- apply(X = as.matrix(time_grid), MARGIN = 1, FUN = fit_one)
  fit <- list(reg.object = fits, time_grid = time_grid)
  class(fit) <- c("f_y_nw")
  return(fit)
}

#' Prediction function for Nadaraya-Watson estimator
#'
#' @param fit Fitted regression object
#' @param newX Values of covariates at which to make a prediction
#' @param newtimes
#'
#' @return Matrix of predictions
#' @noRd
predict.f_y_nw <- function(fit, newX, newtimes){
  predict_nw <- function(t){
    if (t < min(fit$time_grid)){
      pred <- rep(0, nrow(newX))
    } else{
      t_floor <-  max(which(fit$time_grid <= t))
      print(t_floor)
      reg_floor <- fit$reg.object[[t_floor]]
      print(reg_floor)
      pred <- predict(reg_floor, newdata = newX)
    }
    return(pred)
  }
  predictions <- apply(X = as.matrix(newtimes),
                       MARGIN = 1,
                       FUN = predict_nw)
  return(predictions)
}

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

  if (censored){
    time <- time[!as.logical(event)]
    X <- X[!as.logical(event),]
  } else{
    time <- time[as.logical(event)]
    X <- X[as.logical(event),]
  }

  fmla <- as.formula(paste("ind ~ ", paste(colnames(X), collapse = "+")))
  ind <- time
  dat <- data.frame(cbind(ind = ind, X))
  nw_fit <- eval(bquote(np::npcdist(formula = .(fmla),
                          data = dat,
                          bws = c(bwy, rep(bw, ncol(X))),
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

  if (censored){
    time <- time[!as.logical(event)]
    X <- X[!as.logical(event),]
  } else{
    time <- time[as.logical(event)]
    X <- X[as.logical(event),]
  }

  fmla <- as.formula(paste("ind ~ ", paste(colnames(X), collapse = "+")))
  ind <- time
  dat <- data.frame(cbind(ind = ind, X))
  ll_fit <- eval(bquote(np::npcdist(formula = .(fmla),
                                    data = dat,
                                    bws = c(bwy, rep(bw, ncol(X))),
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

