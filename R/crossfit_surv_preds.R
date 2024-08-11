#' Generate K-fold cross-fit survival predictions for downstream use
#'
#' @return data frame giving results
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

  return(list(S_preds = CV_S_preds,
              S_preds_train = CV_S_preds_train,
              G_preds = CV_G_preds,
              G_preds_train = CV_G_preds_train))
}
