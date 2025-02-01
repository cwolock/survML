#' Generate cross-fitted oracle prediction function estimates
#'
#' @param time \code{n x 1} numeric vector of observed
#' follow-up times. If there is censoring, these are the minimum of the
#' event and censoring times.
#' @param event \code{n x 1} numeric vector of status indicators of
#' whether an event was observed.
#' @param X \code{n x p} data.frame of observed covariate values
#' @param folds \code{n x 1} numeric vector of folds identifiers for cross-fitting
#' @param nuisance_preds Named list of conditional event and censoring survival functions
#' that will be used to estimate the oracle prediction function.
#' @param pred_generator Function to be used to estimate oracle prediction function.
#' @param ... Additional arguments to be passed to \code{pred_generator}.
#'
#' @return Named list of cross-fitted oracle prediction estimates
#'
#' @export
crossfit_oracle_preds <- function(time,
                                  event,
                                  X,
                                  folds,
                                  nuisance_preds,
                                  pred_generator,
                                  ...){
  .V <- length(unique(folds))
  CV_oracle_preds <- list()
  CV_oracle_preds_train <- list()

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
    nuisance_preds_j <- sapply(nuisance_preds, "[", j)
    preds <- pred_generator(time = time_train,
                            event = event_train,
                            X = X_train,
                            X_holdout = X_holdout,
                            nuisance_preds = nuisance_preds_j,
                            ...)
    CV_oracle_preds[[j]] <- preds$f0_hat
    CV_oracle_preds_train[[j]] <- preds$f0_hat_train
  }

  return(list(f_hat = CV_oracle_preds,
              f_hat_train = CV_oracle_preds_train))
}
