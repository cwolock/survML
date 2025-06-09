#' Estimate variable importance
#'
#' Compute estimates of and confidence intervals for nonparametric variable importance
#' based on the difference predictiveness obtained with and without the feature of interest.
#' Designed for use with time-to-event outcomes subject to right censoring that may be informed
#' by measured covariates.
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
#' @param restriction_time Maximum follow-up time for calculation of \code{"C-index"} and \code{"survival_time_MSE"}. Essentially, this time
#' should be chosen such that the conditional survival function is identified at this time for all covariate values \code{X} present in the data.
#' Choosing the restriction time such that roughly 10% of individuals remain at-risk at that time has been shown to work reasonably well in simulations.
#' @param approx_times Numeric vector of length J2 giving times at which to
#' approximate integrals. Defaults to a grid of 100 timepoints, evenly spaced on the quantile scale of the distribution of observed event times.
#' @param large_feature_vector Numeric vector giving indices of features to include in the 'large' prediction model.
#' @param small_feature_vector Numeric vector giving indices of features to include in the 'small' prediction model. Must be a
#' subset of \code{large_feature_vector}.
#' @param conditional_surv_generator A function to estimate the conditional survival functions of the event and censoring variables. Must take arguments
#' (\code{time}, \code{event}, \code{X}) (for training purposes) and (\code{X_holdout} and \code{newtimes}) (covariate values and times at which to generate predictions). Defaults to [generate_nuisance_predictions_stackG], a pre-built generator function based on the [stackG] function.
#' Alternatively, the user can provide their own function for this argument, or provide pre-computed estimates to
#' \code{conditional_surv_preds} in lieu of this argument.
#' @param conditional_surv_generator_control A list of arguments to pass to \code{conditional_surv_generator}.
#' @param large_oracle_generator A function to estimate the oracle prediction function using \code{large_feature_vector}.
#' Must take arguments \code{time}, \code{event}, \code{X}, \code{X_holdout}, and \code{nuisance_preds}.
#' For all VIM types except for \code{"C-index"}, defaults to [generate_oracle_predictions_DR], a pre-built generator
#' function using doubly-robust pseudo-outcome regression. For \code{"C-index"}, defaults
#' to [generate_oracle_predictions_boost], a pre-built generator function using doubly-robust gradient boosting.
#' Alternatively, the user can provide their own function, or provide pre-computed estimates to \code{large_oracle_preds} in lieu of this argument.
#' @param large_oracle_generator_control A list of arguments to pass to \code{large_oracle_generator}.
#' @param small_oracle_generator  A function to estimate the oracle prediction function using \code{small_feature_vector}. Must take arguments
#' \code{time}, \code{event}, \code{X}, \code{X_holdout}, and \code{nuisance_preds}. For all VIM types except for \code{"C-index"}, defaults to
#' [generate_oracle_predictions_SL], a pre-built generator function based on regression the large oracle predictions on the small feature vector.
#' For \code{"C-index"}, defaults to [generate_oracle_predictions_boost], a pre-built generator function using doubly-robust gradient boosting.
#' Alternatively, the user can provide their own function, or provide pre-computed estimates to \code{small_oracle_preds} in lieu of this argument.
#' @param small_oracle_generator_control A list of arguments to pass to \code{small_oracle_generator}.
#' @param conditional_surv_preds User-provided estimates of the conditional survival functions of the event and censoring
#' variables given the full covariate vector (if not using the \code{conditional_surv_generator} functionality to compute these nuisance estimates).
#' Must be a named list of lists with elements \code{S_hat}, \code{S_hat_train}, \code{G_hat}, and \code{G_hat_train}. If using sample splitting, each of these is
#' itself a list of length \code{2K}, where \code{K} is the number of cross-fitting folds (if not using sample splitting, each is a list of length \code{K}).
#' Each element of these lists is a matrix with J2 columns and number of rows equal to either the number of samples in the \code{k}th fold (for \code{S_hat} and
#' \code{G_hat}) or the number of samples used to compute the nuisance estimates for the \code{k}th fold (for \code{S_hat_train} and \code{G_hat_train}).
#' @param large_oracle_preds User-provided estimates of the oracle prediction function using \code{large_feature_vector} (if not using the \code{large_oracle_generator} functionality to compute these nuisance estimates). Must be a named list of lists
#' with elements \code{f0_hat} and \code{f0_hat_train}. If using sample splitting, each of these is itself a list of length \code{2K} (if not using sample
#' splitting, each is a list of length \code{K}). Each element of these lists is a matrix with J1 columns (for landmark time VIMs) or 1 column
#' (for \code{"C-index"} and \code{"survival_time_MSE"}) and number of rows equal to either the number of samples in the \code{k}th fold (for \code{f0_hat})
#' or the number of samples used to compute the nuisance estimates for the \code{k}th fold (for \code{f0_hat_train}).
#' @param small_oracle_preds User-provided estimates of the oracle prediction function using \code{small_feature_vector} (if not using the \code{small_oracle_generator} functionality to compute these nuisance estimates). Must be a named list of lists
#' with elements \code{f0_hat} and \code{f0_hat_train}. If using sample splitting, each of these is itself a list of length \code{2K} (if not using sample
#' splitting, each is a list of length \code{K}). Each element of these lists is a matrix with J1 columns (for landmark time VIMs) or 1 column
#' (for \code{"C-index"} and \code{"survival_time_MSE"}) and number of rows equal to either the number of samples in the \code{k}th fold (for \code{f0_hat})
#' or the number of samples used to compute the nuisance estimates for the \code{k}th fold (for \code{f0_hat_train}).
#' @param cf_folds Numeric vector of length \code{n} giving cross-fitting folds, if specifying the folds explicitly. This is required if you are providing pre-computed nuisance estimations --- if providing a nuisance generator function, the \code{vim()} will assign folds.
#' @param cf_fold_num The number of cross-fitting folds, if not providing \code{cf_folds}. Note that with samples-splitting, the data will be split into \code{2 x cf_fold_num} folds (i.e., there will be \code{cf_fold_num} folds within each half of the data).
#' @param sample_split Logical indicating whether or not to sample split. Sample-splitting is required for valid hypothesis testing of null importance and is generally recommended. Defaults to \code{TRUE}.
#' @param ss_folds Numeric vector of length \code{n} giving sample-splitting folds, if specifying the folds explicitly. This is required if you are providing pre-computed nuisance estimations --- if providing a nuisance generator function, the \code{vim()} will assign folds.
#' @param scale_est Logical, whether or not to force the VIM estimate to be nonnegative.
#' @param alpha The level at which to compute confidence intervals and hypothesis tests. Defaults to 0.05.
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
#' @details For nuisance estimation, it is generally advisable to use the pre-built nuisance generator functions provided by \code{survML}. See the ''Variable importance in survival analysis'' vignette, or the \href{https://cwolock.github.io/survML/}{package website} for an illustration.
#'
#' @seealso [vim_accuracy] [vim_AUC] [vim_brier] [vim_cindex] [vim_rsquared] [vim_survival_time_mse]
#'
#' @export
#'
#' @examples
#' # This is a small simulation example
#' set.seed(123)
#' n <- 100
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
#' landmark_times <- c(3)
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
#'               sample_split = FALSE,
#'               scale_est = TRUE)
#'
#' print(output$result)
#'
#' @references Wolock C.J., Gilbert P.B., Simon N., and Carone, M. (2025).
#'   "Assessing variable importance in survival analysis using machine learning."

vim <- function(type,
                time,
                event,
                X,
                landmark_times = stats::quantile(time[event == 1], probs = c(0.25, 0.5, 0.75)),
                restriction_time = max(time[event == 1]),
                approx_times = NULL,
                large_feature_vector,
                small_feature_vector,
                conditional_surv_generator = NULL,
                conditional_surv_generator_control = NULL,
                large_oracle_generator = NULL,
                large_oracle_generator_control = NULL,
                small_oracle_generator = NULL,
                small_oracle_generator_control = NULL,
                conditional_surv_preds = NULL,
                large_oracle_preds = NULL,
                small_oracle_preds = NULL,
                cf_folds = NULL,
                cf_fold_num = 5,
                sample_split = TRUE,
                ss_folds = NULL,
                robust = TRUE,
                scale_est = FALSE,
                alpha = 0.05,
                verbose = FALSE){

  precomputed_SG <- FALSE
  precomputed_f <- FALSE
  precomputed_fs <- FALSE
  cf_folds_nuisance <- NULL # cf fold identifiers for nuisance estimation will not always be equal to this for VIM estimation
  # since with sample splitting + no cross-fitting, you use the whole sample to estimate nuisances

  # Deal with folds...
  if (sample_split){ # if sample splitting
    if (is.null(ss_folds)){ # sample splitting folds must be generated.
      if (is.null(cf_folds)){ # cross-fitting folds also blank, generate everything ourselves
        if (!is.numeric(cf_fold_num) | cf_fold_num <= 0 | cf_fold_num%%1 != 0){
          stop("The number of cross-fitting folds must be a positive integer.")
        }
        if (cf_fold_num == 1){ # sample-splitting, no cross-fitting
          cf_folds_nuisance <- rep(1, length(time)) # make cross-fitting folds just a vector of 1s for nuisance purposes
          folds <- generate_folds(n = length(time), V = 1, sample_split = TRUE)
          ss_folds <- folds$ss_folds
          cf_folds <- folds$cf_folds # cross-fitting folds for vim estimation are equal to sample-splitting folds
        } else if (cf_fold_num > 1){ # sample-splitting and cross-fitting, generate both
          if (verbose){print("Setting up cross-fitting and sample-splitting folds...")}
          folds <- generate_folds(n = length(time), V = cf_fold_num, sample_split = TRUE)
          ss_folds <- folds$ss_folds
          cf_folds <- folds$cf_folds
        }
      } else{ # cross-fitting folds were provided --- generate only ss_folds
        if (verbose){print("Setting up sample-splitting folds to accompany user-specified cross-fitting folds...")}
        folds <- generate_folds(n = length(time),
                                V = length(unique(cf_folds))/2,
                                sample_split = TRUE,
                                cf_folds = cf_folds)
        ss_folds <- folds$ss_folds
      }
    } else{ # if sample splitting folds provided
      if (verbose){print("Using user-provided folds...")}
      if (is.null(cf_folds)){ # cf_folds were not provided. This is okay if no cross-fitting intended, otherwise it's an error
        if (cf_fold_num == 1){
          cf_folds_nuisance <- rep(1, length(time)) # all 1s for nuisance estimation
          cf_folds <- ss_folds # equal to ss_folds for vim estimaton
        } else{
          stop("You have provided sample-spliting folds without cross-fitting folds. This is not allowed.")
        }
      } # otherwise cf folds were provided, check they are okay below
    }
  } else{ # otherwise
    if (is.null(cf_folds)){ # if cf_folds aren't provided, generate some (even if all 1s)
      if (!is.numeric(cf_fold_num) | cf_fold_num <= 0 | cf_fold_num%%1 != 0){
        stop("The number of cross-fitting folds must be a positive integer.")
      }
      if (cf_fold_num == 1){
        if (verbose){print("No cross-fitting and no sample-splitting indicated; no folds to set up...")}
        cf_folds = rep(1, length(time))
      } else if (cf_fold_num > 1){
        if (verbose){print("Setting up cross-fitting folds...")}
        folds <- generate_folds(n = length(time), V = cf_fold_num, sample_split = FALSE)
        cf_folds <- folds$cf_folds
      }
    } else{ # if cross-fitting folds are provided, make sure they're valid below
      if (verbose){print("Using user-provided folds...")}
    }
    ss_folds <- rep(1, length(time))
  }

  if (!all(sort(unique(cf_folds)) == 1:max(cf_folds)) | any(cf_folds <= 0)){ # check validity of cf_folds
    stop("Cross-fiting fold identifiers must be consecutive, positive integers.")
  }
  if (!all(sort(unique(ss_folds)) == 1:max(ss_folds)) | any(ss_folds <= 0) | any(ss_folds >= 3)){ # check validity of ss_folds
    stop("Sample-fiting fold identifiers must be either 1 or 2.")
  }

  if (is.null(cf_folds_nuisance)){ # if not setting to all 1s in the no xfit + sample split case, then just set equal to cf_folds
    cf_folds_nuisance <- cf_folds
  }

  if (!is.null(conditional_surv_preds$S_hat) & !is.null(conditional_surv_preds$G_hat) & !is.null(conditional_surv_preds$S_hat_train) & !is.null(conditional_surv_preds$G_hat_train) & is.null(approx_times)){
    stop("If using precomputed nuisance estimates, you must provide the grid of approx_times on which they were estimated.")
  }
  if (!is.null(conditional_surv_preds$S_hat) & !is.null(conditional_surv_preds$G_hat) & !is.null(conditional_surv_preds$S_hat_train) & !is.null(conditional_surv_preds$G_hat_train) & (is.null(cf_folds)) | is.null(ss_folds)){ # this *should* never trigger given that I have manually set up the folds by this point, regardless of the path taken.
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
    if (is.null(approx_times)){
      approx_times <- sort(unique(c(stats::quantile(time[event == 1 & time <= max(landmark_times)],
                                                    probs = seq(0, 1, by = 0.01)),
                                    landmark_times)))
    } else{
      approx_times <- sort(unique(c(approx_times, landmark_times)))
    }
  }
  if (!is.null(restriction_time) & type %in% global_vims){
    if (is.null(approx_times)){
      approx_times <- sort(unique(c(stats::quantile(time[event == 1 & time <= restriction_time],
                                                    probs = seq(0, 1, by = 0.01)),
                                    restriction_time)))
    } else{
      approx_times <- sort(unique(c(approx_times, restriction_time)))
    }
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
                                             folds = cf_folds_nuisance,
                                             pred_generator = conditional_surv_generator),
                                        conditional_surv_generator_control))

  } else{
    precomputed_SG <- TRUE
    if (verbose){print("Using pre-computed conditional survival function estimates...")}
  }

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
                                         folds = cf_folds_nuisance,
                                         pred_generator = large_oracle_generator),
                                    large_oracle_generator_control))
  } else{
    precomputed_f <- TRUE
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
                                         folds = cf_folds_nuisance,
                                         pred_generator = small_oracle_generator),
                                    small_oracle_generator_control))
  } else{
    precomputed_fs <- TRUE
    if (verbose){print("Using pre-computed 'small' oracle prediction function estimates...")}
  }
  conditional_surv_preds_temp <- conditional_surv_preds
  large_oracle_preds_temp <- large_oracle_preds
  small_oracle_preds_temp <- small_oracle_preds
  if (sample_split & length(unique(cf_folds_nuisance)) == 1){ # if sample splitting without cross-fitting, then cf_folds_nuisance will be a vector of 1s for
    # nuisance estimation, but the nuisance estimates need to be split up for vim estimation
    # if user has provided pre-computed estimates, don't do this splitting
    if (!precomputed_SG){
      conditional_surv_preds_temp$S_hat <- list(conditional_surv_preds$S_hat[[1]][cf_folds == 1,],
                                           conditional_surv_preds$S_hat[[1]][cf_folds == 2,])
      conditional_surv_preds_temp$S_hat_train <- list(conditional_surv_preds$S_hat[[1]],
                                                 conditional_surv_preds$S_hat[[1]])
      conditional_surv_preds_temp$G_hat <- list(conditional_surv_preds$G_hat[[1]][cf_folds == 1,],
                                           conditional_surv_preds$G_hat[[1]][cf_folds == 2,])
      conditional_surv_preds_temp$G_hat_train <- list(conditional_surv_preds$G_hat[[1]],
                                                 conditional_surv_preds$G_hat[[1]])
    }
    if (!precomputed_f){
      large_oracle_preds_temp$f_hat <- list(large_oracle_preds$f_hat[[1]][cf_folds == 1,],
                                       large_oracle_preds$f_hat[[1]][cf_folds == 2,])
      large_oracle_preds_temp$f_hat_train <- list(large_oracle_preds$f_hat[[1]],
                                             large_oracle_preds$f_hat[[1]])
    }
    if (!precomputed_fs){
      small_oracle_preds_temp$f_hat <- list(small_oracle_preds$f_hat[[1]][cf_folds == 1,],
                                       small_oracle_preds$f_hat[[1]][cf_folds == 2,])
      small_oracle_preds_temp$f_hat_train <- list(small_oracle_preds$f_hat[[1]],
                                             small_oracle_preds$f_hat[[1]])
    }
    conditional_surv_preds <- conditional_surv_preds_temp
    large_oracle_preds <- large_oracle_preds_temp
    small_oracle_preds <- small_oracle_preds_temp
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
