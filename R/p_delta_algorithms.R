#' Binary SuperLearner
#'
#' @param event Indicator of event (vs censoring)
#' @param X Covariate matrix
#' @param SL.library SuperLearner library
#' @param V number of CV folds
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
