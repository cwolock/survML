#' Generate cross-fitted conditional survival predictions
#'
#' @param time \code{n x 1} numeric vector of observed
#' follow-up times. If there is censoring, these are the minimum of the
#' event and censoring times.
#' @param event \code{n x 1} numeric vector of status indicators of
#' whether an event was observed.
#' @param X \code{n x p} data.frame of observed covariate values
#' @param newtimes Numeric vector of times on which to estimate the conditional survival functions
#' @param folds \code{n x 1} numeric vector of folds identifiers for cross-fitting
#' @param pred_generator Function to be used to estimate conditional survival function.
#' @param ... Additional arguments to be passed to \code{pred_generator}.
#'
#' @return Named list of cross-fitted conditional survival predictions
#'
#' @export
crossfit_surv_preds <- function(time,
                                event,
                                X,
                                newtimes,
                                folds,
                                pred_generator,
                                ...){
  .V <- length(unique(folds))
  CV_S_preds <- list()
  CV_S_preds_train <- list()
  CV_G_preds <- list()
  CV_G_preds_train <- list()

  for (j in 1:.V){

    if (.V == 1){ # if not actually cross fitting
      time_train <- time[folds == j]
      event_train <- event[folds == j]
      X_train <- X[folds == j,]
    } else{ # if actually cross fitting
      time_train <- time[folds != j]
      event_train <- event[folds != j]
      X_train <- X[folds != j,]
    }
    X_holdout <- X[folds == j,]
    preds <- pred_generator(time = time_train,
                            event = event_train,
                            X = X_train,
                            X_holdout = X_holdout,
                            newtimes = newtimes,
                            ...)
    CV_S_preds[[j]] <- preds$S_hat
    CV_S_preds_train[[j]] <- preds$S_hat_train
    CV_G_preds[[j]] <- preds$G_hat
    CV_G_preds_train[[j]] <- preds$G_hat_train
  }

  return(list(S_hat = CV_S_preds,
              S_hat_train = CV_S_preds_train,
              G_hat = CV_G_preds,
              G_hat_train = CV_G_preds_train))
}
