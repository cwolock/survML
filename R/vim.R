#' Estimate AUC VIM
#'
#' @param type Type of VIM to compute. Options include \code{"accuracy"}, \code{"AUC"}, \code{"Brier"}, \code{"R-squared"}
#' \code{"C-index"}, and \code{"survival_time_MSE"}.
#' @param time \code{n x 1} numeric vector of observed
#' follow-up times. If there is censoring, these are the minimum of the
#' event and censoring times.
#' @param event \code{n x 1} numeric vector of status indicators of
#' whether an event was observed.
#' @param X \code{n x p} data.frame of observed covariate values
#' @param landmark_times Numeric vector of length J1 giving
#' landmark times at which to estimate VIM (\code{"accuracy"}, \code{"AUC"}, \code{"Brier"}, \code{"R-squared"}).
#' @param restriction_time Maximum follow-up time for calculation of \code{"C-index"} and \code{"survival_time_MSE"}.
#' @param approx_times Numeric vector of length J2 giving times at which to
#' approximate integrals.
#' @param large_feature_vector Numeric vector giving indices of features to include in the 'large' prediction model.
#' @param small_feature_vector Numeric vector giving indices of features to include in the 'small' prediction model. Must be a
#' subset of \code{large_feature_vector}.
#' @param conditional_surv_preds User-provided estimates of the conditional survival functions of the event and censoring
#' variables given the full covariate vector (if not using the \code{vim()} function to compute these nuisance estimates).
#' Must be a named list of lists with elements \code{S_hat}, \code{S_hat_train}, \code{G_hat}, and \code{G_hat_train}. Each of these is itself
#' a list of length \code{K}, where \code{K} is the number of cross-fitting folds. Each element of these lists is a matrix with J2 columns and number of rows
#' equal to either the number of samples in the \code{k}th fold (for \code{S_hat} or \code{G_hat}) or the number of samples used to compute the nuisance estimator
#' for the \code{k}th fold.
#' @param large_oracle_preds User-provided estimates of the oracle prediction function using \code{large_feature_vector}. Must be a named list of lists
#' with elements \code{f_hat} and \code{f_hat_train}. Each of these is iself a list of length \code{K}. Each element of these lists is a matrix with J1 columns
#' (for landmark time VIMs) or 1 column (for \code{"C-index"} and \code{"survival_time_MSE"}).
#' @param small_oracle_preds User-provided estimates of the oracle prediction function using \code{small_feature_vector}. Must be a named list of lists
#' with elements \code{f_hat} and \code{f_hat_train}. Each of these is iself a list of length \code{K}. Each element of these lists is a matrix with J1 columns
#' (for landmark time VIMs) or 1 column (for \code{"C-index"} and \code{"survival_time_MSE"}).
#' @param conditional_surv_generator A user-written function to estimate the conditional survival functions of the event and censoring variables. Must take arguments
#' \code{time}, \code{event}, \code{folds} (cross-fitting fold identifiers), and
#' \code{newtimes} (times at which to generate predictions).
#' @param conditional_surv_generator_control A list of arguments to pass to \code{conditional_surv_generator}.
#' @param large_oracle_generator A user-written function to estimate the oracle prediction function using \code{large_feature_vector}.Must take arguments
#' \code{time}, \code{event}, and \code{folds} (cross-fitting fold identifiers).
#' @param large_oracle_generator_control A list of arguments to pass to \code{large_oracle_generator}.
#' @param small_oracle_generator  A user-written function to estimate the oracle prediction function using \code{small_feature_vector}.Must take arguments
#' \code{time}, \code{event}, and \code{folds} (cross-fitting fold identifiers).
#' @param small_oracle_generator_control A list of arguments to pass to \code{small_oracle_generator}.
#' @param cf_folds Numeric vector of length \code{n} giving cross-fitting folds
#' @param cf_fold_num The number of cross-fitting folds, if not providing \code{cf_folds}
#' @param sample_split Logical indicating whether or not to sample split
#' @param ss_folds Numeric vector of length \code{n} giving sample-splitting folds
#' @param scale_est Logical, whether or not to force the VIM estimate to be nonnegative
#' @param alpha The level at which to compute confidence intervals and hypothesis tests. Defaults to 0.05
#' @param robust Logical, whether or not to use the doubly-robust debiasing approach. This option
#' is meant for illustration purposes only --- it should be left as \code{TRUE}.
#' @param verbose Whether to print progress messages.
#'
#' @return Named list with the following elements:
#' \item{result}{Data frame giving results. See the documentation of the individual \code{vim_*} functions for details.}
#' \item{folds}{A named list giving the cross-fitting fold IDs (\code{cf_folds}) and sample-splitting fold IDs (\code{ss_folds}).}
#' \item{approx_times}{A vector of times used to approximate integrals appearing in the form of the VIM estimator.}
#' \item{conditional_surv_preds}{A named list containing the estimated conditional event and censoring survival functions.}
#' \item{large_oracle_preds}{A named list containing the estimated large oracle prediction function.}
#' \item{small_oracle_preds}{A named list containing the estimated small oracle prediction function.}
#'
#' @seealso [vim_accuracy] [vim_AUC] [vim_brier] [vim_cindex] [vim_rsquared] [vim_survival_time_mse]
#'
#' @export
#'
#' @examples
#' # This is a small simulation example
#' set.seed(123)
#' n <- 200
#' X <- data.frame(X1 = rnorm(n), X2 = rbinom(n, size = 1, prob = 0.5))
#'
#' T <- rexp(n, rate = exp(-2 + X[,1] - X[,2] + .5 *  X[,1] * X[,2]))
#'
#' C <- rexp(n, exp(-2 -.5 * X[,1] - .25 * X[,2] + .5 * X[,1] * X[,2]))
#' C[C > 15] <- 15
#'
#' time <- pmin(T, C)
#' event <- as.numeric(T <= C)
#'
#' # landmark times for AUC
#' landmark_times <- c(1,3,5)
#'
#' output <- vim(type = "AUC",
#'               time = time,
#'               event = event,
#'               X = X,
#'               landmark_times = landmark_times,
#'               large_feature_vector = 1:2,
#'               small_feature_vector = 2,
#'               conditional_surv_generator_control = list(SL.library = c("SL.mean", "SL.glm")),
#'               large_oracle_generator_control = list(SL.library = c("SL.mean", "SL.glm")),
#'               small_oracle_generator_control = list(SL.library = c("SL.mean", "SL.glm")),
#'               cf_fold_num = 2,
#'               sample_split = TRUE,
#'               scale_est = TRUE)
#'
#' print(output$result)

vim <- function(type,
                time,
                event,
                X,
                landmark_times = stats::quantile(time[event == 1], probs = c(0.25, 0.5, 0.75)),
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
                alpha = 0.05,
                verbose = FALSE){

  if (is.null(cf_folds)){
    if (verbose){print("Setting up cross-fitting and sample-splitting folds...")}
    folds <- generate_folds(n = length(time), V = cf_fold_num, sample_split = sample_split)
    cf_folds <- folds$cf_folds
    ss_folds <- folds$ss_folds
  } else{
    if (verbose){print("Using user-provided folds...")}
  }

  if (!is.null(conditional_surv_preds$S_hat) & !is.null(conditional_surv_preds$G_hat) & !is.null(conditional_surv_preds$S_hat_train) & !is.null(conditional_surv_preds$G_hat_train) & is.null(approx_times)){
    stop("If using precomputed nuisance estimates, you must provide the grid of approx_times on which they were estimated.")
  }
  if (!is.null(conditional_surv_preds$S_hat) & !is.null(conditional_surv_preds$G_hat) & !is.null(conditional_surv_preds$S_hat_train) & !is.null(conditional_surv_preds$G_hat_train) & (is.null(cf_folds)) | is.null(ss_folds)){
    stop("If using precomputed nuisance estimates, you must provide cross-fitting fold identifiers.")
  }

  landmark_vims <- c("AUC", "Brier", "accuracy", "R-squared")
  global_vims <- c("C-index", "survival_time_MSE")

  if (type %in% landmark_vims){
    outcome <- "survival_probability"
  } else if (type == "survival_time_MSE"){
    outcome <- "restricted_survival_time"
  } else{
    outcome <- NA
  }

  if (!is.null(landmark_times) & type %in% landmark_vims){
    approx_times <- sort(unique(c(time[event == 1 & time <= max(landmark_times)], landmark_times)))
  }
  if (!is.null(restriction_time) & type %in% global_vims){
    approx_times <- sort(unique(c(time[event == 1 & time <= restriction_time], restriction_time)))
  }

  # estimate S and G if any of them are not provided by user
  if ((is.null(conditional_surv_preds$S_hat) | is.null(conditional_surv_preds$G_hat) | is.null(conditional_surv_preds$S_hat_train) | is.null(conditional_surv_preds$G_hat_train))){
    if (verbose){print("Estimating conditional survival nuisance functions...")}
    if (is.null(conditional_surv_generator)){
      conditional_surv_generator <- generate_nuisance_predictions_stackG
    }
    if (is.null(conditional_surv_generator_control)){
      conditional_surv_generator_control <- list()
    }
    if (is.null(conditional_surv_generator_control$approx_times)){
      conditional_surv_generator_control$approx_times <- approx_times
    }

    generator_args <- methods::formalArgs(conditional_surv_generator)
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
    if (verbose){print("Using pre-computed conditional survival function estimates...")}
  }

  # switch <- FALSE
  if (is.null(large_oracle_preds)){
    if (verbose){print("Estimating 'big' oracle prediction function...")}
    if (is.null(large_oracle_generator_control)){
      large_oracle_generator_control <- list()
    }
    large_oracle_generator_control$indx <- which(!(1:ncol(X) %in% large_feature_vector))
    large_oracle_generator_control$outcome <- outcome

    if (is.null(large_oracle_generator_control$approx_times)){
      large_oracle_generator_control$approx_times <- approx_times
    }
    if (is.null(large_oracle_generator_control$landmark_times)){
      large_oracle_generator_control$landmark_times <- landmark_times
    }
    if (is.null(large_oracle_generator_control$nuisance_preds)){
      large_oracle_generator_control$nuisance_preds <- conditional_surv_preds
    }
    if (is.null(large_oracle_generator_control$restriction_time)){
      large_oracle_generator_control$restriction_time <- restriction_time
    }

    if (is.null(large_oracle_generator) & type != "C-index"){
      large_oracle_generator <- generate_oracle_predictions_DR
    } else if (is.null(large_oracle_generator) & type == "C-index"){
      large_oracle_generator <- generate_oracle_predictions_boost
    }

    generator_args <- methods::formalArgs(large_oracle_generator)
    large_oracle_generator_control[which(!(names(large_oracle_generator_control)%in%generator_args))] <- NULL

    large_oracle_preds <- do.call(crossfit_oracle_preds,
                                  c(list(time = time,
                                         event = event,
                                         X = X,
                                         folds = cf_folds,
                                         pred_generator = large_oracle_generator),
                                    large_oracle_generator_control))
  } else{
    if (verbose){print("Using pre-computed 'big' oracle prediction function estimates...")}
  }

  augmented_nuisance_preds <- c(conditional_surv_preds, large_oracle_preds)

  if (is.null(small_oracle_preds)){
    if (verbose){print("Estimating 'small' oracle prediction function...")}

    if (is.null(small_oracle_generator_control)){
      small_oracle_generator_control <- list()
    }

    small_oracle_generator_control$indx <- which(1:ncol(X) %in% large_feature_vector &
                                                   !(1:ncol(X) %in% small_feature_vector))
    small_oracle_generator_control$outcome <- outcome

    if (is.null(small_oracle_generator) & type != "C-index"){
      small_oracle_generator <- generate_oracle_predictions_SL
    } else if (is.null(small_oracle_generator) & type == "C-index"){
      small_oracle_generator <- generate_oracle_predictions_boost
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

    generator_args <- methods::formalArgs(small_oracle_generator)
    small_oracle_generator_control[which(!(names(small_oracle_generator_control)%in%generator_args))] <- NULL

    small_oracle_preds <- do.call(crossfit_oracle_preds,
                                  c(list(time = time,
                                         event = event,
                                         X = X,
                                         folds = cf_folds,
                                         pred_generator = small_oracle_generator),
                                    small_oracle_generator_control))
  } else{
    if (verbose){print("Using pre-computed 'small' oracle prediction function estimates...")}
  }

  if (verbose){print("Estimating variable importance...")}
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
  } else if (type == "C-index"){
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
    res <- vim_survival_time_mse(time = time,
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
  } else if (type == "R-squared"){
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
