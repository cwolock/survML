#' Generate K-fold cross-fit survival predictions for downstream use
#'
#' @return data frame giving results
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
    nuisance_preds_j <- sapply(nuisance_preds, "[[", j)
    preds <- pred_generator(time = time_train,
                            event = event_train,
                            X = X_train,
                            X_holdout = X_holdout,
                            nuisance_preds = nuisance_preds_j,
                            ...)
    CV_oracle_preds[[j]] <- preds$f0_hat
    CV_oracle_preds_train[[j]] <- preds$f0_hat_train
  }

  return(list(oracle_preds = CV_oracle_preds,
              oracle_preds_train = CV_oracle_preds_train))
}
