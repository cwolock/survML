#' Estimate AUC VIM
#'
#' @param type Type of VIM to compute
#' @param time \code{n x 1} numeric vector of observed
#' follow-up times If there is censoring, these are the minimum of the
#' event and censoring times.
#' @param event \code{n x 1} numeric vector of status indicators of
#' whether an event was observed. Defaults to a vector of 1s, i.e. no censoring.
#' @param approx_times Numeric vector of length J1 giving times at which to
#' approximate integrals.
#' @param landmark_times Numeric vector of length J2 giving
#' times at which to estimate AUC
#' @param restriction_time Maximum follow-up time for calculation of C-index or restricted survival time
#' @param conditional_surv_preds
#' @param large_oracle_preds
#' @param small_oracle_preds
#' @param cf_folds Numeric vector of length n giving cross-fitting folds
#' @param cf_fold_num The number of cross-fitting folds, if not providing \code{cf_folds}
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
                conditional_surv_preds = NULL,
                large_oracle_preds = NULL,
                small_oracle_preds = NULL,
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

  if (!is.null(conditional_surv_preds$S_hat) & !is.null(conditional_surv_preds$G_hat) & !is.null(conditional_surv_preds$S_hat_train) & !is.null(conditional_surv_preds$G_hat_train) & is.null(approx_times)){
    stop("If using precomputed nuisance estimates, you must provide the grid of approx_times on which they were estimated.")
  }
  if (!is.null(conditional_surv_preds$S_hat) & !is.null(conditional_surv_preds$G_hat) & !is.null(conditional_surv_preds$S_hat_train) & !is.null(conditional_surv_preds$G_hat_train) & (is.null(cf_folds)) | is.null(ss_folds)){
    stop("If using precomputed nuisance estimates, you must provide cross-fitting fold identifiers.")
  }

  landmark_vims <- c("AUC", "Brier", "accuracy", "R-squared")
  global_vims <- c("C-index", "survival_time_MSE")
  if (!is.null(landmark_times) & type %in% landmark_vims){
    approx_times <- sort(unique(c(time[event == 1 & time <= max(landmark_times)], landmark_times)))
  }
  if (!is.null(restriction_time) & type %in% global_vims){
    approx_times <- sort(unique(c(time[event == 1 & time <= restriction_time], restriction_time)))
  }

  # estimate S and G if any of them are not provided by user
  if ((is.null(conditional_surv_preds$S_hat) | is.null(conditional_surv_preds$G_hat) | is.null(conditional_surv_preds$S_hat_train) | is.null(conditional_surv_preds$G_hat_train))){
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

    conditional_surv_preds <- do.call(crossfit_surv_preds,
                               c(list(time = time,
                                      event = event,
                                      X = X,
                                      newtimes = approx_times,
                                      folds = cf_folds,
                                      pred_generator = conditional_surv_generator),
                                 conditional_surv_generator_control))

  } else{
    print("Using pre-computed conditional survival function estimates...")
  }

  # switch <- FALSE
  if (is.null(large_oracle_preds)){
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
      large_oracle_generator_control$nuisance_preds <- conditional_surv_preds
    }

    if (is.null(large_oracle_generator)){
      large_oracle_generator <- generate_oracle_predictions_DR
    }

    generator_args <- formalArgs(large_oracle_generator)
    large_oracle_generator_control[which(!(names(large_oracle_generator_control)%in%generator_args))] <- NULL

    large_oracle_preds <- do.call(crossfit_oracle_preds,
                          c(list(time = time,
                                 event = event,
                                 X = X,
                                 folds = cf_folds,
                                 pred_generator = large_oracle_generator),
                            large_oracle_generator_control))
  } else{
    print("Using pre-computed 'big' oracle prediction function estimates...")
  }

  augmented_nuisance_preds <- c(conditional_surv_preds, large_oracle_preds)

  if (is.null(small_oracle_preds)){
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

    small_oracle_preds <- do.call(crossfit_oracle_preds,
                             c(list(time = time,
                                    event = event,
                                    X = X,
                                    folds = cf_folds,
                                    pred_generator = small_oracle_generator),
                               small_oracle_generator_control))
  } else{
    print("Using pre-computed 'small' oracle prediction function estimates...")
  }

  print("Estimating variable importance...")
  if (type == "AUC"){
    res <- vim_AUC(time = time,
                   event = event,
                   landmark_times = landmark_times,
                   approx_times = approx_times,
                   f_hat = large_oracle_preds$f_hat,
                   fs_hat = small_oracle_preds$f_hat,
                   S_hat = conditional_surv_preds$S_hat,
                   G_hat = conditional_surv_preds$G_hat,
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
                     f_hat = large_oracle_preds$f_hat,
                     fs_hat = small_oracle_preds$f_hat,
                     S_hat = conditional_surv_preds$S_hat,
                     G_hat = conditional_surv_preds$G_hat,
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
                      f_hat = large_oracle_preds$f_hat,
                      fs_hat = small_oracle_preds$f_hat,
                      S_hat = conditional_surv_preds$S_hat,
                      G_hat = conditional_surv_preds$G_hat,
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
                        f_hat = large_oracle_preds$f_hat,
                        fs_hat = small_oracle_preds$f_hat,
                        S_hat = conditional_surv_preds$S_hat,
                        G_hat = conditional_surv_preds$G_hat,
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
                        f_hat = large_oracle_preds$f_hat,
                        fs_hat = small_oracle_preds$f_hat,
                        S_hat = conditional_surv_preds$S_hat,
                        G_hat = conditional_surv_preds$G_hat,
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
                        f_hat = large_oracle_preds$f_hat,
                        fs_hat = small_oracle_preds$f_hat,
                        S_hat = conditional_surv_preds$S_hat,
                        G_hat = conditional_surv_preds$G_hat,
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
              conditional_surv_preds = conditional_surv_preds,
              large_oracle_preds = large_oracle_preds,
              small_oracle_preds = small_oracle_preds))
}
