#' List wrappers for p_delta prediction
#'
#' This function lists all prediction algorithms built into the \code{conSurv} package.
#'
#' @usage \code{p_delta_listWrappers()}
#'
#' @return Invisible character vector with the requested functions.
#'
#' @export

p_delta_listWrappers <- function() {
  everything <- sort(unclass(lsf.str(envir = asNamespace("conSurv"), all.names = T)))
  message("All prediction algorithm wrappers in conSurv:\n")
  funs <- everything[grepl(pattern = "^p_delta", everything)]
  funs <- funs[!funs %in% c("p_delta_listWrappers")]
  print(funs)
  invisible(everything)
}

#' Nadaraya-Watson estimator
#'
#' @param event Event indicator
#' @param X Covariates
#' @param span Span
#' @param kernel_type Kernel_type
#' @param kernel_order Kernel_order
#'
#' @return An object of class \code{p_delta_nw}
#' @noRd
p_delta_nw <- function(event, X, span, kernel_type = "gaussian", kernel_order = 2){
  fmla <- as.formula(paste("event ~ ", paste(colnames(X), collapse = "+")))
  dat <- data.frame(cbind(event = event, X))
  fit.nw <- np::npreg(fmla,
                   data = dat,
                   bws = rep(span, ncol(X)),
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
