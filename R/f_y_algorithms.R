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
    (stats::approx(predictions[j,], seq(0, .999, by=.001), xout = newtimes, method = "linear", rule = 2)$y)
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
    (stats::approx(predictions[j,], seq(0, .999, by=.001), xout = newtimes, method = "linear", rule = 2)$y)
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
    (stats::approx(predictions[j,], seq(0, .999, by=.001), xout = newtimes, method = "linear", rule = 2)$y)
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

#' Quantile regression forest
#'
#' @param time Observed time
#' @param event Indicator of event (vs censoring)
#' @param X Covariate matrix
#' @param censored Logical, indicates whether to run regression on censored observations (vs uncensored)
#' @param mtry Number of variables sampled as candidates as each split
#' @param ntree Number of trees to grow
#'
#' @return An object of class \code{f_y_qrf}
#' @noRd
f_y_qrf <- function(time, event, X, censored, mtry = floor(sqrt(ncol(X))), ntree = 500){

  if (censored){
    time <- time[!as.logical(event)]
    X <- X[!as.logical(event),]
  } else{
    time <- time[as.logical(event)]
    X <- X[as.logical(event),]
  }

  qrf_fit <- quantregForest::quantregForest(x = X, y = time)

  fit <- list(reg.object = qrf_fit)
  class(fit) <- c("f_y_qrf")
  return(fit)
}

#' Prediction function for quantile regression forest
#'
#' @param fit Fitted regression object
#' @param newX Values of covariates at which to make a prediction
#' @param newtimes
#'
#' @return Matrix of predictions
#' @noRd
predict.f_y_qrf <- function(fit, newX, newtimes){

  predict_qrf <- function(x){
    pred <- predict(fit$reg.object, newdata = newX, what = function(z) ecdf(z)(x))
    return(pred[1:nrow(newX)])
  }
  predictions <- apply(X = as.matrix(newtimes),
                       MARGIN = 1,
                       FUN = predict_qrf)
  return(predictions)
}

#' Linear quantile regression
#'
#' @param time Observed time
#' @param event Indicator of event (vs censoring)
#' @param X Covariate matrix
#' @param censored Logical, indicates whether to run regression on censored observations (vs uncensored)
#'
#' @return An object of class \code{f_y_quantreg}
#' @noRd
f_y_quantreg <- function(time, event, X, censored){

  if (censored){
    time <- time[!as.logical(event)]
    X <- X[!as.logical(event),]
  } else{
    time <- time[as.logical(event)]
    X <- X[as.logical(event),]
  }

  dat <- data.frame(cbind(time = time, X))
  fmla <- as.formula(paste("time ~ ", paste(colnames(X), collapse = "+")))

  qr_fit <- quantreg::rq(formula = fmla,
                         tau = seq(0.01, 0.99, by = 0.01),
                         data = dat)

  fit <- list(reg.object = qr_fit)
  class(fit) <- c("f_y_quantreg")
  return(fit)
}

#' Prediction function for linear quantile regression
#'
#' @param fit Fitted regression object
#' @param newX Values of covariates at which to make a prediction
#' @param newtimes
#'
#' @return Matrix of predictions
#' @noRd
predict.f_y_quantreg <- function(fit, newX, newtimes){

  predictions <- predict(fit$reg.object, newdata = newX)
  predictions <- t(sapply(1:nrow(predictions), function(j) {
    (stats::approx(predictions[j,], seq(0.01, .99, by=.01), xout = newtimes, method = "linear", rule = 2)$y)
  }))

  return(predictions)
}

#' Additive neural net quantile regression
#'
#' @param time Observed time
#' @param event Indicator of event (vs censoring)
#' @param X Covariate matrix
#' @param censored Logical, indicates whether to run regression on censored observations (vs uncensored)
#' @param n.hidden Number of hidden nodes, passed to \code{qrnn.fit}
#'
#' @return An object of class \code{f_y_qrnn}
#' @noRd
f_y_qrnn <- function(time, event, X, censored, n.hidden = 3){

  if (censored){
    time <- time[!as.logical(event)]
    X <- X[!as.logical(event),]
  } else{
    time <- time[as.logical(event)]
    X <- X[as.logical(event),]
  }

  X <- as.matrix(X)
  time <- as.matrix(time)
  tau <- seq(0.01, 0.99, by = 0.01)
  X.time.tau <- qrnn::composite.stack(X, time, tau)
  fit <- qrnn::qrnn.fit(cbind(X.time.tau$tau, X.time.tau$x),
                        X.time.tau$y,
                        tau = X.time.tau$tau,
                        n.hidden = n.hidden,
                        monotone = NULL,
                        additive = TRUE,
                        trace = FALSE)
  #fit <- qrnn::mcqrnn.fit(x = X, y = time, tau = tau, n.hidden= 2)

  fit <- list(reg.object = fit)
  class(fit) <- c("f_y_qrnn")
  return(fit)
}

#' Prediction function for quantile regression neural net
#'
#' @param fit Fitted regression object
#' @param newX Values of covariates at which to make a prediction
#' @param newtimes
#'
#' @return Matrix of predictions
#' @noRd
predict.f_y_qrnn <- function(fit, newX, newtimes){

  newX <- as.matrix(newX)

  tau <- seq(0.01, 0.99, by = 0.01)
  X.tau <- qrnn::composite.stack(x = newX, y = as.matrix(0), tau = tau)
  predictions <- matrix(qrnn::qrnn.predict(cbind(X.tau$tau, X.tau$x), fit$reg.object), ncol = length(tau))
  #predictions <- qrnn::mcqrnn.predict(newX, fit$reg.object)

  predictions <- t(sapply(1:nrow(predictions), function(j) {
    (stats::approx(predictions[j,], seq(0.01, .99, by=.01), xout = newtimes, method = "linear", rule = 2)$y)
  }))

  return(predictions)
}

#' Stacked binary regression with Super Learner
#'
#' @param time Observed time
#' @param event Indicator of event (vs censoring)
#' @param X Covariate matrix
#' @param censored Logical, indicates whether to run regression on censored observations (vs uncensored)
#' @param bin_size Size of quantiles over which to make the stacking bins
#'
#' @return An object of class \code{f_y_stackSL}
#' @noRd
f_y_stackSL <- function(time, event, X, censored, bin_size){

  if (censored){
    time <- time[!as.logical(event)]
    X <- X[!as.logical(event),]
  } else{
    time <- time[as.logical(event)]
    X <- X[as.logical(event),]
  }

  X <- as.matrix(X)
  time <- as.matrix(time)
  dat <- data.frame(X, time)

  time_grid <- quantile(dat$time, probs = seq(0, 1, by = bin_size))
  time_grid[1] <- 0 # manually set first point to 0, instead of first observed time
  trunc_time_grid <- time_grid[-length(time_grid)]
  # alternatively, time grid could just be observed times, or a fixed (non-data-dependent) grid

  ncol_stacked <- ncol(X) + length(trunc_time_grid) + 1 # covariates, risk set dummies, binary outcome
  stacked <- matrix(NA, ncol = ncol_stacked, nrow = 1)
  for (i in 1:(length(trunc_time_grid))){
    risk_set <- dat[dat$time > time_grid[i],]
    risk_set_covariates <- risk_set[,-ncol(risk_set)]
    event_indicators <- matrix(ifelse(risk_set$time < time_grid[i + 1], 1, 0))
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

  SL.library <- c("SL.mean", "SL.glm", "SL.gam", "SL.randomForest")
  fit <- SuperLearner::SuperLearner(Y = Y,
                                    X = X,
                                    family = "binomial",
                                    SL.library = SL.library,
                                    method = "method.NNloglik",
                                    verbose = FALSE)

  fit <- list(reg.object = fit, time_grid = time_grid)
  class(fit) <- c("f_y_stackSL")
  return(fit)
}

#' Prediction function for stacked SL
#'
#' @param fit Fitted regression object
#' @param newX Values of covariates at which to make a prediction
#' @param newtimes
#'
#' @return Matrix of predictions
#' @noRd
predict.f_y_stackSL <- function(fit, newX, newtimes){

  time_grid <- fit$time_grid
  trunc_time_grid <- time_grid[-length(time_grid)]

  get_stacked_pred <- function(t){
    n_bins <- sum(t > trunc_time_grid)
    new_stacked <- matrix(rep(t(newX), n_bins), ncol = ncol(newX) , byrow = TRUE)
    dummy_stacked <- matrix(NA, ncol = length(trunc_time_grid), nrow = 1)
    for (i in 1:n_bins){
      dummies <- matrix(0, ncol = length(trunc_time_grid), nrow = nrow(newX))
      if (t > trunc_time_grid[i]){
        dummies[,i] <- 1
      }
      dummy_stacked <- rbind(dummy_stacked, dummies)
    }
    dummy_stacked <- dummy_stacked[-1,]
    new_stacked <- cbind(dummy_stacked, new_stacked)
    risk_set_names <- paste0("risk_set_", seq(1, (length(trunc_time_grid))))
    colnames(new_stacked) <- c(risk_set_names, colnames(newX))
    new_stacked <- data.frame(new_stacked)
    haz_preds <- predict(fit$reg.object, newdata=new_stacked)$pred
    haz_preds <- matrix(haz_preds, nrow = nrow(newX))
    surv_preds <- 1 - haz_preds
    total_surv_preds <- apply(surv_preds, prod, MARGIN = 1)
    return(total_surv_preds)
  }

  predictions <- apply(X = matrix(newtimes), FUN = get_stacked_pred, MARGIN = 1)
  predictions <- 1 - predictions
  return(predictions)
}

#' Discrete hazard SuperLearner
#'
#' @param time Observed time
#' @param event Indicator of event (vs censoring)
#' @param X Covariate matrix
#' @param censored Logical, indicates whether to run regression on censored observations (vs uncensored)
#' @param bin_size Size of quantiles over which to make the stacking bins
#'
#' @return An object of class \code{f_y_discSL}
#' @noRd
f_y_discSL <- function(time, event, X, censored, bin_size){

  if (censored){
    time <- time[!as.logical(event)]
    X <- X[!as.logical(event),]
  } else{
    time <- time[as.logical(event)]
    X <- X[as.logical(event),]
  }

  # X <- as.matrix(X)
  # time <- as.matrix(time)
  dat <- data.frame(X, time)

  time_grid <- quantile(dat$time, probs = seq(0, 1, by = bin_size))
  time_grid[1] <- 0 # manually set first point to 0, instead of first observed time
  intervals <- cbind(time_grid[-length(time_grid)], time_grid[-1])
  n.intervals <- length(time_grid) - 1

  SL.library <- c("SL.mean", "SL.glm", "SL.gam")#, "SL.randomForest")
  sl.fits <- lapply(1:(n.intervals-1), function(j) {
    if(j == 1){
      samp <- time >= 0
    } else {
      samp <- time > intervals[j,1]
    }
    outcome <- as.numeric(time <= intervals[j,2])
    # } else {
    #   outcome <- rep(1, length(time))
    # }
    fit <- SuperLearner::SuperLearner(Y = outcome[samp],
                                      X = X[samp,],
                                      SL.library = SL.library,
                                      family = 'binomial',
                                      method = 'method.NNloglik',
                                      verbose = FALSE)
    fit
  })

  fit <- list(reg.object = sl.fits, time_grid = time_grid)
  class(fit) <- c("f_y_discSL")
  return(fit)
}

#' Prediction function for discrete hazard SL
#'
#' @param fit Fitted regression object
#' @param newX Values of covariates at which to make a prediction
#' @param newtimes
#'
#' @return Matrix of predictions
#' @noRd
predict.f_y_discSL <- function(fit, newX, newtimes){

  time_grid <- fit$time_grid
  n.intervals <- length(time_grid) - 1
  intervals <- cbind(time_grid[-length(time_grid)], time_grid[-1])

  new.time.bins <- findInterval(newtimes, time_grid, all.inside = TRUE)

  hazard.ests <- sapply(1:(n.intervals-1), function(j) {
    predict(fit$reg.object[[j]], newdata=newX)$pred
  })

  # set hazard in last time bin to 1
  hazard.ests <- cbind(hazard.ests, rep(1, nrow(newX)))

  hazard.ests <- 1-hazard.ests

  get_pred <- function(j){
    bin <- new.time.bins[j]
    preds <- apply(hazard.ests[,1:bin,drop=FALSE], MARGIN = 1, FUN = prod)
  }

  predictions <- sapply(1:length(newtimes), get_pred)

  predictions <- 1 - predictions
  return(predictions)
}

