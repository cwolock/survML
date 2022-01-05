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
p_delta_ranger <- function(event, X, mtry = floor(sqrt(ncol(X))), num.trees = 500){

  event <- as.factor(event)

  # Ranger does not seem to work with X as a matrix, so we explicitly convert to
  # data.frame rather than cbind. newX can remain as-is though.
  if (is.matrix(X)) {
    X <- data.frame(X)
  }

  fit <- ranger::ranger(event ~ .,
                        data = cbind(event = event, X),
                        num.trees = num.trees,
                        mtry = mtry,
                        min.node.size = 1,
                        replace = TRUE,
                        sample.fraction = 1,
                        write.forest = TRUE,
                        probability = TRUE,
                        #num.threads = num.threads,
                        verbose = TRUE)

  fit <- list(reg.object = fit)
  class(fit) <- c("p_delta_ranger")
  return(fit)
}

#' Prediction function for ranger
#'
#' @param fit Fitted regression object
#' @param newX Values of covariates at which to make a prediction
#'
#' @return Matrix of predictions
#' @noRd
predict.p_delta_ranger <- function(fit, newX){
  pred <- predict(fit$reg.object, data = newX)$predictions
  print(pred)
  pred <- pred[, "1"]
  return(pred)
}

#' Random forest
#'
#' @param event Event indicator
#' @param X Covariates
#' @param mtry Number of varaibles sampled as candidates as each split
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
