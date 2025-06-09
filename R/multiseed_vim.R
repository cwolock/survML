#' Estimate variable importance with multiple seeds
#'
#' Repeat the VIM estimation procedure multiple times and aggregate the results,
#' mitigating the additional randomness introduced by sample-splitting and cross-fitting.
#'
#' @inheritParams vim
#' @param n_seed Number of iterations (seeds) to perform the VIM estimation procedure.
#' These will be aggregated into a single result.
#' @param agg_method P-value aggregation method use to combine results from different seeds.
#' Current options are \code{"bonferroni"} (Bonferroni's method),
#' \code{"hommel"} (Hommel's method), \code{"arithmetic"} (arithmetic mean), \code{"geometric"}
#' (geometric mean), \code{"harmonic"} (harmonic mean),
#' \code{"compound_bg"} (compound Bonferroni and geometric mean), and \code{"compound_ba"}
#' (compound Bonferroni and arithmetic mean). These approaches
#' are discussed at length in Vovk and Wang (2020). Defaults to \code{"compound_bg"},
#' which has been shown to work well in many settings.
#' @param ci_grid Grid of VIM values over which to construct a confidence interval by
#' inverting a hypothesis test. The aggregation works by constructing
#' hypothesis tests (at level \code{alpha}) of the null corresponding to each value in
#' \code{ci_grid}, and then inverting these tests to yield a
#' 1 - \code{alpha} confidence interval. For example, for \code{"AUC"} importance, the VIM takes
#' values in (0,1), so a grid of values between 0 and 1
#' would be a reasonable choice.
#'
#' @return Named list with the following elements:
#' \item{agg_result}{Data frame giving results aggregated over seeds.}
#' \item{agg_method}{P-value aggregation method used.}
#' \item{n_seed}{Number of iterations (seeds) used to perform the VIM estimation procedure.}
#' \item{vim_objects}{A list of \code{vim} return objects, each corresponding to a different seed.}
#'
#' @details Using a larger value of \code{n_seed} will result in more stable results, at a greater computational cost.
#'
#' @seealso [vim]
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
#' output <- multiseed_vim(n_seed = 2,
#'               agg_method = "compound_bg",
#'               ci_grid = seq(0, 1, by = 0.01),
#'               type = "AUC",
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
#'
#' @references Vovk V. and Wang R. (2020). "Combining p-values via averaging."
#' @references Wolock C.J., Gilbert P.B., Simon N., and Carone, M. (2025).
#'   "Assessing variable importance in survival analysis using machine learning."
multiseed_vim <- function(n_seed,
                          agg_method = "compound_bg",
                          ci_grid,
                          type,
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
                          cf_fold_num = 5,
                          sample_split = TRUE,
                          scale_est = FALSE,
                          alpha = 0.05,
                          verbose = FALSE){
  result_list <- list()
  vim_list <- list()
  for (i in 1:n_seed){
    vim_i <- vim(type = type,
                 time = time,
                 event = event,
                 X = X,
                 landmark_times = landmark_times,
                 restriction_time = restriction_time,
                 approx_times = approx_times,
                 large_feature_vector = large_feature_vector,
                 small_feature_vector = small_feature_vector,
                 conditional_surv_generator = conditional_surv_generator,
                 conditional_surv_generator_control = conditional_surv_generator_control,
                 large_oracle_generator = large_oracle_generator,
                 large_oracle_generator_control = large_oracle_generator_control,
                 small_oracle_generator = small_oracle_generator,
                 small_oracle_generator_control = small_oracle_generator_control,
                 cf_fold_num = cf_fold_num,
                 sample_split = sample_split,
                 scale_est = scale_est,
                 alpha = alpha,
                 verbose = verbose)
    vim_list[[i]] <- vim_i
    result_list[[i]] <- vim_i$result
  }

  agg_result <- aggregate_vim(result_list = result_list,
                              agg_method = agg_method,
                              ci_grid = ci_grid,
                              n_eff = ifelse(sample_split, length(time)/2, length(time)),
                              alpha = alpha)

  return(list(agg_result = agg_result,
              agg_method = agg_method,
              n_seed = n_seed,
              vim_objects = vim_list))
}

#' Aggregate multiseed VIM results
#'
#' @param result_list List of result data frames return by the \code{vim} function.
#' @param agg_method P-value aggregation method use to combine results from different seeds. Current options are \code{"bonferroni"}
#' (Bonferroni's method), \code{"hommel"} (Hommel's method), \code{"arithmetic"} (arithmetic mean), \code{"geometric"} (geometric mean),
#' \code{"harmonic"} (harmonic mean), \code{"compound_bg"} (compound Bonferroni and geometric mean), and \code{"compound_ba"}
#' (compound Bonferroni and arithmetic mean). These approaches are discussed at length in Vovk and Wang (2020). Defaults to \code{"compound_bg"}, which has been shown to work well in many settings.
#' @param ci_grid Grid of VIM values over which to construct a confidence interval by inverting a hypothesis test. The aggregation works by constructing
#' hypothesis tests (at level \code{alpha}) of the null corresponding to each value in \code{ci_grid}, and then inverting these tests to yield a
#' 1 - \code{alpha} confidence interval. For example, for \code{"AUC"} importance, the VIM takes values in (0,1), so a grid of values between 0 and 1
#' would be a reasonable choice.
#' @param n_eff The effective sample size. Without sample-splitting, this is simply the sample size. With sample-splitting, this is the sample size divided by
#' two (i.e., the size of each of the two halves of the data).
#' @param alpha The level at which to compute confidence intervals and hypothesis tests. Defaults to 0.05.
#'
#' @return Named list with the following elements:
#' \item{agg_result}{Data frame giving results aggregated over seeds.}
#' \item{agg_method}{P-value aggregation method used.}
#' \item{n_seed}{Number of iterations (seeds) used to perform the VIM estimation procedure.}
#' \item{vim_objects}{A list of \code{vim} return objects, each corresponding to a different seed.}
#'
#' @export
aggregate_vim <- function(result_list,
                          agg_method,
                          ci_grid,
                          n_eff,
                          alpha = 0.05){
  res1 <- result_list[[1]]
  static_col_names <- names(res1)[which(!(names(res1) %in% c("est", "var_est", "cil", "ciu",
                                                             "cil_1sided", "p", "large_predictiveness", "small_predictiveness")))]
  avg_col_names <- c("est", "var_est", "large_predictiveness", "small_predictiveness")
  static_cols <- res1[,names(res1) %in% static_col_names]
  result_list <- lapply(result_list, FUN = function(x) x[,!(names(res1) %in% static_col_names)])
  seeds <- 1:length(result_list)

  n1 <- nrow(res1)
  n2 <- length(ci_grid)
  n3 <- length(result_list)

  p_array_twosided <- array(NA, dim = c(n1, n2, n3))
  p_array_onesided <- array(NA, dim = c(n1, n2, n3))

  # calculate p-value for every point in ci_grid and every seed
  for (i in 1:n3){
    dat <- result_list[[i]]
    for (j in 1:n2){
      for (k in 1:n1){
        z <- (dat$est[k] - ci_grid[j]) / sqrt(dat$var_est[k]/n_eff)
        p <- stats::pnorm(abs(z), lower.tail = FALSE)*2
        p_array_twosided[k, j, i] <- p
        p <- stats::pnorm(z, lower.tail = FALSE)
        p_array_onesided[k, j, i] <- p
      }
    }
  }

  p_array_twosided_corrected <- array(NA, dim = c(n1, n2))
  p_array_onesided_corrected <- array(NA, dim = c(n1, n2))

  # aggregate across seeds
  for (j in 1:n2){
    for (k in 1:n1){
      p_two <- aggregate_p(p_array_twosided[k,j,],method=agg_method)
      p_one <- aggregate_p(p_array_onesided[k,j,],method=agg_method)
      p_array_twosided_corrected[k,j] <- p_two
      p_array_onesided_corrected[k,j] <- p_one
    }
  }

  cil <- rep(NA, n1)
  ciu <- rep(NA, n1)
  cil_1sided <- rep(NA, n1)
  p <- rep(NA, n1)

  for (k in 1:n1){
    curr_ps_twosided <- p_array_twosided_corrected[k,]
    curr_ps_onesided <- p_array_onesided_corrected[k,]
    p[k] <- curr_ps_onesided[1]
    if (curr_ps_twosided[1] < alpha){
      lower_index <- min(which(curr_ps_twosided > alpha))
      cil[k] <- ci_grid[lower_index]
      reject_indices <- which(curr_ps_twosided < alpha)
      upper_index <-  min(reject_indices[reject_indices > lower_index])
      ciu[k]<- ci_grid[upper_index]
    } else{
      cil[k] <- ci_grid[1]
      ciu[k]<- ci_grid[min(which(curr_ps_twosided < alpha))]
    }
    if (curr_ps_onesided[1] < alpha){
      cil_1sided[k] <- ci_grid[min(which(curr_ps_onesided > alpha))]
    } else{
      cil_1sided[k] <- ci_grid[1]
    }
  }

  result_list <- lapply(result_list, FUN = function(x) x[,avg_col_names])
  avg_results <- t(apply(matrix(1:nrow(res1)),
                         MARGIN = 1,
                         FUN = function(q) colMeans(dplyr::bind_rows(lapply(result_list, FUN = function(x) x[q,])))))
  agg_results <- cbind(static_cols, avg_results, data.frame(cil = cil, ciu = ciu, cil_1sided = cil_1sided, p = p))
  # cil <- ifelse(cil >= 0, cil, 0)
  agg_results <- agg_results[,names(res1)] # make sure column order is maintained
  return(agg_results)
}

#' Aggregate p-values according to one of several aggregation functions
#' @noRd
aggregate_p <- function(ps, method = "bonferroni"){
  K <- length(ps)
  if (method == "bonferroni"){ # bonferroni
    p_agg <- min(ps) * K
  } else if (method == "hommel"){ # hommel
    constant <- sum(apply(matrix(1:K), MARGIN = 1, FUN = function(x) 1/x))
    sorted_ps <- sort(ps)
    corrected_ps <- apply(matrix(1:K),
                          MARGIN = 1,
                          FUN = function(x) K/x * sorted_ps[x])
    p_agg <- constant * min(corrected_ps)
  } else if (method == "arithmetic"){ # arithmetic mean
    p_agg <- mean(ps) * 2
  } else if (method == "geometric"){ # geometric mean
    p_agg <- exp(mean(log(ps))) * exp(1)
  } else if (method == "harmonic"){ # harmonic mean
    p_agg <- 1/mean(1/ps) * exp(1) * log(K)
  } else if (method == "compound_bg"){ # compound bonferroni-geometric mean
    bonf <- min(ps) * K
    geo <- exp(mean(log(ps))) * exp(1)
    p_agg <- 2 * min(c(bonf, geo))
  } else if (method == "compound_ba"){ # compound bonferroni-arithmetic mean
    bonf <- min(ps) * K
    arith <- mean(ps) * 2
    p_agg <- 2 * min(c(bonf, arith))
  }
  return(min(c(p_agg, 1)))
}

