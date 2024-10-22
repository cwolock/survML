#' Estimate nuisances using stackG
#'
#' @noRd
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
                             newtimes = approx_times,
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

#' Estimate residual oracle prediction function using SL regression
#'
#' @noRd
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
#' @noRd
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
  } else if (outcome == "restricted_survival_time"){
    f0_hat <- DR_predictions
    f0_hat_train <- DR_predictions_train
  }

  return(list(f0_hat = f0_hat,
              f0_hat_train = f0_hat_train))

}

#' Estimate full oracle prediction function using DR pseudo-outcome regression
#'
#' @noRd
generate_oracle_predictions_boost <- function(time,
                                              event,
                                              X,
                                              X_holdout,
                                              nuisance_preds,
                                              restriction_time,
                                              approx_times,
                                              V = 5,
                                              indx,
                                              tuning = "none",
                                              subsample_n = length(time),
                                              params = list(mstop = c(100),
                                                            nu = c(0.1),
                                                            sigma = c(0.01),
                                                            learner = c("glm"))){

  if (sum(indx) != 0){ # only remove column if there is a column to remove
    X <- X[,-indx,drop=FALSE]
    X_holdout <- X_holdout[,-indx,drop=FALSE]
  }

  boost_results <- boost_c_index_DR(time = time,
                                    event = event,
                                    X = X,
                                    S_hat = nuisance_preds$S_hat_train,
                                    G_hat = nuisance_preds$G_hat_train,
                                    approx_times = approx_times[approx_times <= restriction_time],
                                    tuning = tuning,
                                    produce_fit = TRUE,
                                    subsample_n = subsample_n,
                                    V = V,
                                    params = params)

  dtest <- data.frame(X_holdout)
  dtrain <- data.frame(X)
  f0_hat <- -stats::predict(boost_results$opt_model, newdata = dtest)[,1]
  f0_hat_train <- -stats::predict(boost_results$opt_model, newdata = dtrain)[,1]

  return(list(f0_hat = f0_hat,
              f0_hat_train = f0_hat_train))

}

