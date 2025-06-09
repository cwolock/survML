#' Estimate conditional survival function nuisance parameters using survival stacking
#'
#' @param time \code{n x 1} numeric vector of observed
#' follow-up times. If there is censoring, these are the minimum of the
#' event and censoring times.
#' @param event \code{n x 1} numeric vector of status indicators of
#' whether an event was observed.
#' @param X \code{n x p} data.frame of observed covariate values
#' @param X_holdout \code{m x p} data.frame of new observed covariate
#' values at which to obtain \code{m} predictions for the estimated algorithm.
#' Must have the same names and structure as \code{X}.
#' @param approx_times Numeric vector of times at which to
#' approximate product integral or cumulative hazard interval. See [stackG] documentation.
#' @param SL.library Super Learner library
#' @param V Number of cross-validation folds, to be passed to \code{SuperLearner}
#' @param newtimes \code{k x 1} numeric vector of times at which to obtain \code{k}
#' predicted conditional survivals.
#' @param bin_size Size of time bin on which to discretize for estimation
#' of cumulative probability functions. Can be a number between 0 and 1,
#' indicating the size of quantile grid (e.g. \code{0.1} estimates
#' the cumulative probability functions on a grid based on deciles of
#' observed \code{time}s). If \code{NULL}, creates a grid of
#' all observed \code{time}s. See [stackG] documentation.
#'
#' @return A list containing elements \code{S_hat} (conditional event survival function, corresponding to \code{X_holdout} and \code{newtimes}),
#' \code{S_hat_train} (conditional event survival function, corresponding to \code{X} and \code{newtimes}),
#' \code{G_hat} (conditional censoring survival function, corresponding to \code{X_holdout} and \code{newtimes}),
#' and \code{G_hat_train} (conditional censoring survival function, corresponding to \code{X} and \code{newtimes})
#'
#' @seealso [stackG]
#'
#' @export
generate_nuisance_predictions_stackG <- function(time,
                                                 event,
                                                 X,
                                                 X_holdout,
                                                 newtimes,
                                                 SL.library = c("SL.mean", "SL.glm", "SL.earth", "SL.gam", "SL.ranger"),
                                                 V = 5,
                                                 bin_size = 0.05,
                                                 approx_times){

  surv_out <- survML::stackG(time = time,
                             event = event,
                             X = X,
                             newX = rbind(X_holdout, X),
                             newtimes = newtimes,
                             time_grid_approx = approx_times,
                             bin_size = bin_size,
                             time_basis = "continuous",
                             surv_form = "PI",
                             SL_control = list(SL.library = SL.library,
                                               V = V))
  S_hat <- surv_out$S_T_preds[1:nrow(X_holdout),]
  G_hat <- surv_out$S_C_preds[1:nrow(X_holdout),]
  S_hat_train <- surv_out$S_T_preds[(nrow(X_holdout)+1):(nrow(X_holdout)+nrow(X)),]
  G_hat_train <- surv_out$S_C_preds[(nrow(X_holdout)+1):(nrow(X_holdout)+nrow(X)),]
  return(list(S_hat = S_hat,
              G_hat = G_hat,
              S_hat_train = S_hat_train,
              G_hat_train = G_hat_train))
}

#' Estimate small oracle prediction function using Super Learner regression
#'
#' When the oracle prediction function is a conditional mean, the small oracle prediction function can be estimated by
#' regressing the large oracle prediction function on the small feature vector. This function performs this regression using
#' Super Learner.
#'
#' @param time \code{n x 1} numeric vector of observed
#' follow-up times. If there is censoring, these are the minimum of the
#' event and censoring times.
#' @param event \code{n x 1} numeric vector of status indicators of
#' whether an event was observed.
#' @param X \code{n x p} data.frame of observed covariate values
#' @param X_holdout \code{m x p} data.frame of new observed covariate
#' values at which to obtain \code{m} predictions for the estimated algorithm.
#' Must have the same names and structure as \code{X}.
#' @param approx_times Numeric vector of length J2 giving times at which to
#' approximate integral appearing in the pseudo-outcomes
#' @param landmark_times Numeric vector of length J1 giving
#' landmark times at which to estimate VIM (\code{"accuracy"}, \code{"AUC"}, \code{"Brier"}, \code{"R-squared"}).
#' @param restriction_time Maximum follow-up time for calculation of \code{"survival_time_MSE"}. Essentially, this time
#' should be chosen such that the conditional survival function is identified at this time for all covariate values \code{X} present in the data.
#' Choosing the restriction time such that roughly 10% of individuals remain at-risk at that time has been shown to work reasonably well in simulations.
#' @param nuisance_preds Named list of conditional survival function predictions with elements \code{"S_hat"}, \code{"S_hat_train"},
#' \code{"G_hat"}, and \code{"G_hat_train"}. This should match the output of \code{conditional_surv_generator}.
#' @param outcome Outcome type, either \code{"survival_probability"} or \code{"restricted_survival_time"}
#' @param SL.library Super Learner library
#' @param V Number of cross-validation folds, to be passed to \code{SuperLearner}
#' @param indx Numeric index of column(s) of \code{X} to be removed, i.e., not used in the oracle prediction function.
#'
#' @return A list containing elements \code{f0_hat} and \code{f0_hat_train}, the estimated small oracle prediction functions for
#' \code{X_holdout} and \code{X}, respectively.
#'
#' @export
generate_oracle_predictions_SL <- function(time,
                                           event,
                                           X,
                                           X_holdout,
                                           nuisance_preds,
                                           outcome,
                                           landmark_times,
                                           restriction_time,
                                           approx_times,
                                           SL.library = c("SL.mean", "SL.glm", "SL.earth", "SL.gam", "SL.ranger"),
                                           V = 5,
                                           indx){
  X_reduced_train <- X[,-indx,drop=FALSE]
  X_reduced_holdout <- X_holdout[,-indx,drop=FALSE]

  if (outcome == "survival_probability"){
    preds_j <- matrix(NA, nrow = nrow(X_reduced_holdout), ncol = length(landmark_times))
    preds_j_train <- matrix(NA, nrow = nrow(X_reduced_train), ncol = length(landmark_times))
  }
  # else if (outcome == "restricted_survival_time"){
  # preds_j <- matrix(NA, nrow = nrow(X_reduced_holdout), ncol = 1)
  # preds_j_train <- matrix(NA, nrow = nrow(X_reduced_train), ncol = 1)
  # }

  if (outcome == "survival_probability"){
    for (t in landmark_times){
      # if (is.null(nuisance_preds$f_hat_train)){
      # outcomes <- nuisance_preds$S_hat_train[,which(approx_times == t),drop=FALSE]
      # } else{
      outcomes <- nuisance_preds$f_hat_train[,which(landmark_times == t)]
      # }

      long_dat <- data.frame(f_hat = outcomes, X_reduced_train)
      long_new_dat <- data.frame(X_reduced_holdout)
      long_new_dat_train <- data.frame(X_reduced_train)
      reduced_fit <- SuperLearner::SuperLearner(Y = long_dat$f_hat,
                                                X = long_dat[,2:ncol(long_dat),drop=FALSE],
                                                family = stats::gaussian(),
                                                SL.library = SL.library,
                                                cvControl = list(V = V),
                                                method = "method.NNLS",
                                                verbose = FALSE)
      reduced_preds <- matrix(stats::predict(reduced_fit, newdata = long_new_dat)$pred,
                              nrow = nrow(X_reduced_holdout),
                              ncol = 1)
      reduced_preds_train <- matrix(stats::predict(reduced_fit, newdata = long_new_dat_train)$pred,
                                    nrow = nrow(X_reduced_train),
                                    ncol = 1)
      preds_j[,which(landmark_times == t)] <- reduced_preds
      preds_j_train[,which(landmark_times == t)] <- reduced_preds_train
    }
  } else if (outcome == "restricted_survival_time"){
    # if (is.null(nuisance_preds$f_hat_train)){
    # outcomes <- nuisance_preds$S_hat_train[,which(approx_times == t),drop=FALSE]
    # } else{
    outcomes <- nuisance_preds$f_hat_train
    # }

    long_dat <- data.frame(f_hat = outcomes, X_reduced_train)
    long_new_dat <- data.frame(X_reduced_holdout)
    long_new_dat_train <- data.frame(X_reduced_train)
    reduced_fit <- SuperLearner::SuperLearner(Y = long_dat$f_hat,
                                              X = long_dat[,2:ncol(long_dat),drop=FALSE],
                                              family = stats::gaussian(),
                                              SL.library = SL.library,
                                              cvControl = list(V = V),
                                              method = "method.NNLS",
                                              verbose = FALSE)
    preds_j <- stats::predict(reduced_fit, newdata = long_new_dat)$pred
    preds_j_train <- stats::predict(reduced_fit, newdata = long_new_dat_train)$pred
  }


  return(list(f0_hat = preds_j,
              f0_hat_train = preds_j_train))
}

#' Estimate full oracle prediction function using DR pseudo-outcome regression
#'
#' @param time \code{n x 1} numeric vector of observed
#' follow-up times. If there is censoring, these are the minimum of the
#' event and censoring times.
#' @param event \code{n x 1} numeric vector of status indicators of
#' whether an event was observed.
#' @param X \code{n x p} data.frame of observed covariate values
#' @param X_holdout \code{m x p} data.frame of new observed covariate
#' values at which to obtain \code{m} predictions for the estimated algorithm.
#' Must have the same names and structure as \code{X}.
#' @param approx_times Numeric vector of length J2 giving times at which to
#' approximate integral appearing in the pseudo-outcomes
#' @param landmark_times Numeric vector of length J1 giving
#' landmark times at which to estimate VIM (\code{"accuracy"}, \code{"AUC"}, \code{"Brier"}, \code{"R-squared"}).
#' @param restriction_time Maximum follow-up time for calculation of \code{"survival_time_MSE"}. Essentially, this time
#' should be chosen such that the conditional survival function is identified at this time for all covariate values \code{X} present in the data.
#' Choosing the restriction time such that roughly 10% of individuals remain at-risk at that time has been shown to work reasonably well in simulations.
#' @param nuisance_preds Named list of conditional survival function predictions with elements \code{"S_hat"}, \code{"S_hat_train"},
#' \code{"G_hat"}, and \code{"G_hat_train"}. This should match the output of \code{conditional_surv_generator}.
#' @param outcome Outcome type, either \code{"survival_probability"} or \code{"restricted_survival_time"}
#' @param SL.library Super Learner library
#' @param V Number of cross-validation folds, to be passed to \code{SuperLearner}
#' @param indx Numeric index of column(s) of \code{X} to be removed, i.e., not used in the oracle prediction function.
#'
#' @return A list containing elements \code{f0_hat} and \code{f0_hat_train}, the estimated oracle prediction functions for
#' \code{X_holdout} and \code{X}, respectively.
#'
#' @seealso [DR_pseudo_outcome_regression]
#'
#' @export
generate_oracle_predictions_DR <- function(time,
                                           event,
                                           X,
                                           X_holdout,
                                           nuisance_preds,
                                           outcome,
                                           landmark_times,
                                           restriction_time,
                                           approx_times,
                                           SL.library = c("SL.mean", "SL.glm", "SL.earth", "SL.gam", "SL.ranger"),
                                           V = 5,
                                           indx){

  S_hat <- nuisance_preds$S_hat_train
  G_hat <- nuisance_preds$G_hat_train

  if (sum(indx) != 0){ # only remove column if there is a column to remove
    X <- X[,-indx,drop=FALSE]
    X_holdout <- X_holdout[,-indx,drop=FALSE]
  }

  if (outcome == "survival_probability"){
    newtimes <- landmark_times
  } else if (outcome == "restricted_survival_time"){
    newtimes <- restriction_time
  }
  DR_predictions_combined <- DR_pseudo_outcome_regression(time = time,
                                                          event = event,
                                                          X = X,
                                                          newX = rbind(X_holdout, X),
                                                          S_hat = S_hat,
                                                          G_hat = G_hat,
                                                          newtimes = newtimes,
                                                          outcome = outcome,
                                                          approx_times = approx_times,
                                                          SL.library = SL.library,
                                                          V = V)

  DR_predictions <- DR_predictions_combined[1:nrow(X_holdout),]
  DR_predictions_train <- DR_predictions_combined[(nrow(X_holdout) + 1):nrow(DR_predictions_combined),]

  if (outcome == "survival_probability"){
    f0_hat <- 1 - DR_predictions
    f0_hat_train <- 1 - DR_predictions_train
    if (length(landmark_times) == 1){
      f0_hat <- matrix(f0_hat, ncol = 1)
      f0_hat_train <- matrix(f0_hat_train, ncol = 1)
    }
  } else if (outcome == "restricted_survival_time"){
    f0_hat <- DR_predictions
    f0_hat_train <- DR_predictions_train
  }

  return(list(f0_hat = f0_hat,
              f0_hat_train = f0_hat_train))

}

#' Estimate oracle prediction function using DR gradient boosting
#'
#' @param time \code{n x 1} numeric vector of observed
#' follow-up times. If there is censoring, these are the minimum of the
#' event and censoring times.
#' @param event \code{n x 1} numeric vector of status indicators of
#' whether an event was observed.
#' @param X \code{n x p} data.frame of observed covariate values
#' @param X_holdout \code{m x p} data.frame of new observed covariate
#' values at which to obtain \code{m} predictions for the estimated algorithm.
#' Must have the same names and structure as \code{X}.
#' @param approx_times Numeric vector of length J2 giving times at which to
#' approximate C-index integral.
#' @param restriction_time Maximum follow-up time for calculation of the C-index.
#' Essentially, this time should be chosen such that the conditional survival function is identified at
#' this time for all covariate values \code{X} present in the data. Choosing the restriction time such that roughly 10% of individuals remain at-risk
#' at that time has been shown to work reasonably well in simulations.
#' @param nuisance_preds Named list of conditional survival function predictions with elements \code{"S_hat"}, \code{"S_hat_train"},
#' \code{"G_hat"}, and \code{"G_hat_train"}. This should match the output of \code{conditional_surv_generator}.
#' @param V Number of cross-validation folds for selection of tuning parameters
#' @param tuning Logical, whether or not to use cross-validation to select tuning parameters
#' @param subsample_n Number of samples to use for boosting procedure. Using a subsample of the full sample can greatly reduce runtime
#' @param boosting_params Named list of parameter values for the boosting procedure. Elements of this list include \code{mstop} (number of
#' boosting iterations), \code{nu} (learning rate), \code{sigma} (smoothness parameter for sigmoid approximation, with smaller meaning
#' less smoothing), and \code{learner} (base learner, can take values \code{"glm"}, \code{"gam"}, or \code{"tree"})
#' @param indx Numeric index of column(s) of \code{X} to be removed, i.e., not used in the oracle prediction function.
#'
#' @return A list containing elements \code{f0_hat} and \code{f0_hat_train}, the estimated oracle prediction functions for
#' \code{X_holdout} and \code{X}, respectively.
#'
#' @seealso [boost_c_index]
#'
#' @export
generate_oracle_predictions_boost <- function(time,
                                              event,
                                              X,
                                              X_holdout,
                                              nuisance_preds,
                                              restriction_time,
                                              approx_times,
                                              V = 5,
                                              indx,
                                              tuning = FALSE,
                                              subsample_n = length(time),
                                              boosting_params = list(mstop = c(100),
                                                                     nu = c(0.1),
                                                                     sigma = c(0.01),
                                                                     learner = c("glm"))){

  if (sum(indx) != 0){ # only remove column if there is a column to remove
    X <- X[,-indx,drop=FALSE]
    X_holdout <- X_holdout[,-indx,drop=FALSE]
  }

  boost_results <- boost_c_index(time = time,
                                 event = event,
                                 X = X,
                                 newX = rbind(X_holdout, X),
                                 S_hat = nuisance_preds$S_hat_train,
                                 G_hat = nuisance_preds$G_hat_train,
                                 approx_times = approx_times[approx_times <= restriction_time],
                                 tuning = tuning,
                                 produce_fit = TRUE,
                                 subsample_n = subsample_n,
                                 V = V,
                                 boosting_params = boosting_params)

  f0_hat <- boost_results[1:nrow(X_holdout)]
  f0_hat_train <- boost_results[(nrow(X_holdout) + 1):length(boost_results)]

  return(list(f0_hat = f0_hat,
              f0_hat_train = f0_hat_train))

}

