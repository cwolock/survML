#' Estimate AUC VIM
#'
#' @param time \code{n x 1} numeric vector of observed
#' follow-up times If there is censoring, these are the minimum of the
#' event and censoring times.
#' @param event \code{n x 1} numeric vector of status indicators of
#' whether an event was observed. Defaults to a vector of 1s, i.e. no censoring.
#' @param approx_times Numeric vector of length J1 giving times at which to
#' approximate integrals.
#' @param landmark_times Numeric vector of length J2 giving
#' times at which to estimate AUC
#' @param f_hat Full oracle predictions (n x J1 matrix)
#' @param fs_hat Residual oracle predictions (n x J1 matrix)
#' @param S_hat Estimates of conditional event time survival function (n x J2 matrix)
#' @param G_hat Estimate of conditional censoring time survival function (n x J2 matrix)
#' @param folds Numeric vector of length n giving cross-fitting folds
#' @param sample_split Logical indicating whether or not to sample split
#' @param ss_folds Numeric vector of length n giving sample-splitting folds
#' @param scale_est Logical, whether or not to force the VIM estimate to be nonnegative
#' @param alpha The level at which to compute confidence intervals and hypothesis tests. Defaults to 0.05
#' @param robust Logical, whether or not to use the doubly-robust debiasing approach. This option
#' is meant for illustration purposes only --- it should be left as \code{TRUE}.
#'
#' @return data frame giving results
#'
#' @export

vim <- function(type,
                time,
                event,
                X,
                landmark_times = quantile(time[event == 1], probs = c(0.25, 0.5, 0.75)),
                restriction_time = max(time[event == 1]),
                approx_times = NULL,
                large_feature_vector,
                small_feature_vector,
                f_hat = NULL,
                fs_hat = NULL,
                S_hat = NULL,
                G_hat = NULL,
                S_hat_train = NULL,
                G_hat_train = NULL,
                conditional_surv_generator = NULL,
                conditional_surv_generator_control = NULL,
                large_oracle_generator = NULL,
                large_oracle_generator_control = NULL,
                small_oracle_generator = NULL,
                small_oracle_generator_control = NULL,
                cf_folds = NULL,
                cf_fold_num = 5,
                sample_split = TRUE,
                ss_folds = NULL,
                robust = TRUE,
                scale_est = FALSE,
                alpha = 0.05){

  if (is.null(cf_folds)){
    print("Setting up cross-fitting and sample-splitting folds...")
    folds <- generate_folds(n = length(time), V = cf_fold_num, sample_split = sample_split)
    cf_folds <- folds$cf_folds
    ss_folds <- folds$ss_folds
  } else{
    print("Using user-provided folds...")
  }

  if (!is.null(S_hat) & !is.null(G_hat) & !is.null(S_hat_train) & !is.null(G_hat_train) & is.null(approx_times)){
    stop("If using precomputed nuisance estimates, you must provide the grid of approx_times on which they were estimated.")
  }
  if (!is.null(S_hat) & !is.null(G_hat) & !is.null(S_hat_train) & !is.null(G_hat_train) & is.null(cf_folds)){
    stop("If using precomputed nuisance estimates, you must provide cross-fitting fold identifiers.")
  }

  # if (is.null(conditional_surv_generator)){
  #
  # }

  landmark_vims <- c("AUC", "Brier", "accuracy", "R-squared")
  global_vims <- c("C-index", "survival_time_MSE")
  if (!is.null(landmark_times) & type %in% landmark_vims){
    approx_times <- sort(unique(c(time[event == 1 & time <= max(landmark_times)], landmark_times)))
  }
  if (!is.null(restriction_time) & type %in% global_vims){
    approx_times <- sort(unique(c(time[event == 1 & time <= restriction_time], restriction_time)))
  }

  # estimate S and G if any of them are not provided by user
  if ((is.null(S_hat) | is.null(G_hat) | is.null(S_hat_train) | is.null(G_hat_train))){
    print("Estimating conditional survival nuisance functions...")
    if (is.null(conditional_surv_generator)){
      conditional_surv_generator <- generate_nuisance_predictions_stackG
    }
    if (is.null(conditional_surv_generator_control)){
      conditional_surv_generator_control <- list()
    }
    if (is.null(conditional_surv_generator_control$approx_times)){
      conditional_surv_generator_control$approx_times <- approx_times
    }

    generator_args <- formalArgs(conditional_surv_generator)
    conditional_surv_generator_control[which(!(names(conditional_surv_generator_control)%in%generator_args))] <- NULL

    nuisance_preds <- do.call(crossfit_surv_preds,
                              c(list(time = time,
                                     event = event,
                                     X = X,
                                     newtimes = approx_times,
                                     folds = cf_folds,
                                     pred_generator = conditional_surv_generator),
                                conditional_surv_generator_control))

    # nuisance_preds <- crossfit_surv_preds(time = time,
    #                                       event = event,
    #                                       X = X,
    #                                       newtimes = approx_times,
    #                                       folds = cf_folds,
    #                                       pred_generator = conditional_surv_generator,
    #                                       approx_times = approx_times,
    #                                       SL.library = conditional_surv_generator_control$SL.library,
    #                                       V = conditional_surv_generator_control$V,
    #                                       bin_size = conditional_surv_generator_control$bin_size)
    S_hat <- nuisance_preds$S_preds
    G_hat <- nuisance_preds$G_preds
    S_hat_train <- nuisance_preds$S_preds_train
    G_hat_train <- nuisance_preds$G_preds_train
  } else{
    print("Using pre-computed conditional survival function estimates...")
    nuisance_preds <- list(S_preds = S_hat,
                           G_preds = G_hat,
                           S_preds_train = S_hat_train,
                           G_preds_train = G_hat_train)
  }

  # switch <- FALSE
  if (is.null(f_hat)){
    print("Estimating 'big' oracle prediction function...")
    if (is.null(large_oracle_generator_control)){
      large_oracle_generator_control <- list()
    }
    large_oracle_generator_control$indx <- which(!(1:ncol(X) %in% large_feature_vector))

    if (is.null(large_oracle_generator_control$approx_times)){
      large_oracle_generator_control$approx_times <- approx_times
    }
    if (is.null(large_oracle_generator_control$landmark_times)){
      large_oracle_generator_control$landmark_times <- landmark_times
    }
    if (is.null(large_oracle_generator_control$nuisance_preds)){
      large_oracle_generator_control$nuisance_preds <- nuisance_preds
    }

    if (is.null(large_oracle_generator)){
      large_oracle_generator <- generate_oracle_predictions_DR
    }

    generator_args <- formalArgs(large_oracle_generator)
    large_oracle_generator_control[which(!(names(large_oracle_generator_control)%in%generator_args))] <- NULL

    full_preds <- do.call(crossfit_oracle_preds,
                          c(list(time = time,
                                 event = event,
                                 X = X,
                                 folds = cf_folds,
                                 pred_generator = large_oracle_generator),
                            large_oracle_generator_control))
    large_oracle_preds <- full_preds$oracle_preds
    large_oracle_preds_train <- full_preds$oracle_preds_train
  }

  augmented_nuisance_preds <- nuisance_preds
  augmented_nuisance_preds$large_oracle_preds_train <- large_oracle_preds_train

  if (is.null(fs_hat)){
    print("Estimating 'small' oracle prediction function...")

    if (is.null(small_oracle_generator_control)){
      small_oracle_generator_control <- list()
    }

    small_oracle_generator_control$indx <- which(1:ncol(X) %in% large_feature_vector &
                                                   !(1:ncol(X) %in% small_feature_vector))

    if (is.null(small_oracle_generator)){
      small_oracle_generator <- generate_oracle_predictions_SL
    }

    if (is.null(small_oracle_generator_control$approx_times)){
      small_oracle_generator_control$approx_times <- approx_times
    }
    if (is.null(small_oracle_generator_control$landmark_times)){
      small_oracle_generator_control$landmark_times <- landmark_times
    }
    if (is.null(small_oracle_generator_control$nuisance_preds)){
      small_oracle_generator_control$nuisance_preds <- augmented_nuisance_preds
    }
    if (is.null(small_oracle_generator_control$restriction_time)){
      small_oracle_generator_control$restriction_time <- restriction_time
    }



    generator_args <- formalArgs(small_oracle_generator)
    small_oracle_generator_control[which(!(names(small_oracle_generator_control)%in%generator_args))] <- NULL

    reduced_preds <- do.call(crossfit_oracle_preds,
                             c(list(time = time,
                                    event = event,
                                    X = X,
                                    folds = cf_folds,
                                    pred_generator = small_oracle_generator),
                               small_oracle_generator_control))

    small_oracle_preds <- reduced_preds$oracle_preds
    small_oracle_preds_train <- reduced_preds$oracle_preds_train
  }

  if (type == "AUC"){
    # if (switch){
      # large_oracle_preds = lapply(large_oracle_preds, function(x) 1-x)
      # small_oracle_preds = lapply(small_oracle_preds, function(x) 1-x)
    # }
    res <- vim_AUC(time = time,
                   event = event,
                   landmark_times = landmark_times,
                   approx_times = approx_times,
                   f_hat = large_oracle_preds,
                   fs_hat = small_oracle_preds,
                   S_hat = S_hat,
                   G_hat = G_hat,
                   cf_folds = cf_folds,
                   sample_split = sample_split,
                   ss_folds = ss_folds,
                   robust = robust,
                   scale_est = scale_est,
                   alpha = alpha)
  } else if (type == "Brier"){
    res <- vim_brier(time = time,
                     event = event,
                     landmark_times = landmark_times,
                     approx_times = approx_times,
                     f_hat = large_oracle_preds,
                     fs_hat = small_oracle_preds,
                     S_hat = S_hat,
                     G_hat = G_hat,
                     cf_folds = cf_folds,
                     sample_split = sample_split,
                     ss_folds = ss_folds,
                     scale_est = scale_est,
                     alpha = alpha)
  } else if (type == "C_index"){
    res <- vim_cindex(time = time,
                      event = event,
                      restriction_time = restriction_time,
                      approx_times = approx_times,
                      f_hat = large_oracle_preds,
                      fs_hat = small_oracle_preds,
                      S_hat = S_hat,
                      G_hat = G_hat,
                      cf_folds = cf_folds,
                      sample_split = sample_split,
                      ss_folds = ss_folds,
                      scale_est = scale_est,
                      alpha = alpha)
  } else if (type == "survival_time_MSE"){
    res <- vim_rmse_mse(time = time,
                        event = event,
                        restriction_time = restriction_time,
                        approx_times = approx_times,
                        f_hat = large_oracle_preds,
                        fs_hat = small_oracle_preds,
                        S_hat = S_hat,
                        G_hat = G_hat,
                        cf_folds = cf_folds,
                        sample_split = sample_split,
                        ss_folds = ss_folds,
                        scale_est = scale_est,
                        alpha = alpha)
  } else if (type == "R_squared"){
    res <- vim_rsquared(time = time,
                        event = event,
                        landmark_times = landmark_times,
                        approx_times = approx_times,
                        f_hat = large_oracle_preds,
                        fs_hat = small_oracle_preds,
                        S_hat = S_hat,
                        G_hat = G_hat,
                        cf_folds = cf_folds,
                        sample_split = sample_split,
                        ss_folds = ss_folds,
                        scale_est = scale_est,
                        alpha = alpha)
  } else if (type == "accuracy"){
    res <- vim_accuracy(time = time,
                        event = event,
                        landmark_times = landmark_times,
                        approx_times = approx_times,
                        f_hat = large_oracle_preds,
                        fs_hat = small_oracle_preds,
                        S_hat = S_hat,
                        G_hat = G_hat,
                        cf_folds = cf_folds,
                        sample_split = sample_split,
                        ss_folds = ss_folds,
                        scale_est = scale_est,
                        alpha = alpha)
  }

  res$vim <- type
  res$large_feature_vector <- paste0(large_feature_vector, collapse = ",")
  res$small_feature_vector <- paste0(small_feature_vector, collapse = ",")

  return(list(result = res,
              folds = list(cf_folds = cf_folds, ss_folds = ss_folds),
              approx_times = approx_times,
              nuisance_preds = list(S_hat = S_hat,
                                    S_hat_train = S_hat_train,
                                    G_hat = G_hat,
                                    G_hat_train = G_hat_train,
                                    large_oracle_preds = large_oracle_preds,
                                    large_oracle_preds_train = large_oracle_preds_train,
                                    small_oracle_preds = small_oracle_preds,
                                    small_oracle_preds_train = small_oracle_preds_train)))
}
