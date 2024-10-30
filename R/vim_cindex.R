#' Estimate concordance index VIM
#'
#' @param time \code{n x 1} numeric vector of observed
#' follow-up times If there is censoring, these are the minimum of the
#' event and censoring times.
#' @param event \code{n x 1} numeric vector of status indicators of
#' whether an event was observed. Defaults to a vector of 1s, i.e. no censoring.
#' @param approx_times Numeric vector of length J1 giving times at which to
#' approximate integrals.
#' @param restriction_time Restriction time (upper bound for event times to be compared in computing the C-index)
#' @param f_hat Full oracle predictions (n x J1 matrix)
#' @param fs_hat Residual oracle predictions (n x J1 matrix)
#' @param S_hat Estimates of conditional event time survival function (n x J2 matrix)
#' @param G_hat Estimate of conditional censoring time survival function (n x J2 matrix)
#' @param cf_folds Numeric vector of length n giving cross-fitting folds
#' @param sample_split Logical indicating whether or not to sample split
#' @param ss_folds Numeric vector of length n giving sample-splitting folds
#' @param scale_est Logical, whether or not to force the VIM estimate to be nonnegative
#' @param alpha The level at which to compute confidence intervals and hypothesis tests. Defaults to 0.05
#'
#' @return A data frame giving results, with the following columns:
#' \item{restriction_time}{Restriction time (upper bound for event times to be compared in computing the C-index).}
#' \item{est}{VIM point estimate.}
#' \item{var_est}{Estimated variance of the VIM estimate.}
#' \item{cil}{Lower bound of the VIM confidence interval.}
#' \item{ciu}{Upper bound of the VIM confidence interval.}
#' \item{cil_1sided}{Lower bound of a one-sided confidence interval.}
#' \item{p}{p-value corresponding to a hypothesis test of null importance.}
#' \item{large_predictiveness}{Estimated predictiveness of the large oracle prediction function.}
#' \item{small_predictiveness}{Estimated predictiveness of the small oracle prediction function.}
#' \item{vim}{VIM type.}
#' \item{large_feature_vector}{Group of features available for the large oracle prediction function.}
#' \item{small_feature_vector}{Group of features available for the small oracle prediction function.}
#'
#' @seealso [vim] for example usage
#'
#' @export

vim_cindex <- function(time,
                       event,
                       approx_times,
                       restriction_time,
                       f_hat,
                       fs_hat,
                       S_hat,
                       G_hat,
                       cf_folds,
                       sample_split,
                       ss_folds,
                       scale_est = FALSE,
                       alpha = 0.05){

  V <- length(unique(cf_folds))

  # CV_full_plug_ins <- rep(NA, V)
  # CV_reduced_plug_ins <- rep(NA, V)
  CV_var_ests <- rep(NA, V)
  CV_full_numerators <- rep(NA, V)
  CV_reduced_numerators <- rep(NA, V)
  CV_denominators <- rep(NA, V)
  split_numerator_fulls <- rep(NA, V)
  split_denominator_fulls <- rep(NA, V)
  split_numerator_reduceds <- rep(NA, V)
  split_denominator_reduceds <- rep(NA, V)
  # split_plug_in_fulls <- rep(NA, V)
  # split_plug_in_reduceds <- rep(NA, V)
  split_var_est_fulls <- rep(NA, V)
  split_var_est_reduceds <- rep(NA, V)
  for (j in 1:V){
    time_holdout <- time[cf_folds == j]
    event_holdout <- event[cf_folds == j]

    V_0 <- estimate_cindex(time = time_holdout,
                           event = event_holdout,
                           approx_times = approx_times,
                           tau =  restriction_time,
                           preds = f_hat[[j]],
                           S_hat = S_hat[[j]],
                           G_hat = G_hat[[j]])
    V_0s <- estimate_cindex(time = time_holdout,
                            event = event_holdout,
                            approx_times = approx_times,
                            tau =  restriction_time,
                            preds = fs_hat[[j]],
                            S_hat = S_hat[[j]],
                            G_hat = G_hat[[j]])

    # CV_full_plug_ins[j] <- V_0$plug_in
    # CV_reduced_plug_ins[j] <- V_0s$plug_in
    CV_full_numerators[j] <- V_0$numerator
    CV_reduced_numerators[j] <- V_0s$numerator
    CV_denominators[j] <- V_0$denominator
    split_numerator_fulls[j] <- V_0$numerator
    split_denominator_fulls[j] <- V_0$denominator
    split_numerator_reduceds[j] <- V_0s$numerator
    split_denominator_reduceds[j] <- V_0s$denominator
    # split_plug_in_fulls[j] <- V_0$plug_in
    # split_plug_in_reduceds[j] <- V_0s$plug_in
    split_var_est_fulls[j] <- mean(V_0$EIF^2)
    split_var_est_reduceds[j] <- mean(V_0s$EIF^2)
    EIF <- V_0$EIF - V_0s$EIF
    CV_var_ests[j] <- mean(EIF^2)
  }

  if (sample_split){
    folds_0 <- sort(unique(cf_folds[ss_folds == 0]))
    folds_1 <- sort(unique(cf_folds[ss_folds == 1]))
    one_step <- mean(split_numerator_fulls[folds_0])/mean(split_denominator_fulls[folds_0]) -
      mean(split_numerator_reduceds[folds_1])/mean(split_denominator_reduceds[folds_1])
    full_one_step <- mean(split_numerator_fulls[folds_0])/mean(split_denominator_fulls[folds_0])
    # full_plug_in <- mean(split_plug_in_fulls[folds_0])
    reduced_one_step <- mean(split_numerator_reduceds[folds_1])/mean(split_denominator_reduceds[folds_1])
    # reduced_plug_in <- mean(split_plug_in_reduceds[folds_1])
    var_est <- mean(split_var_est_fulls[folds_0]) +
      mean(split_var_est_reduceds[folds_1])
  } else{
    one_step <- mean(CV_full_numerators - CV_reduced_numerators)/mean(CV_denominators)#mean(CV_one_steps)
    var_est <- mean(CV_var_ests)
    full_one_step <- mean(CV_full_numerators)/mean(CV_denominators)#mean(CV_full_one_steps)
    reduced_one_step <- mean(CV_reduced_numerators)/mean(CV_denominators)#mean(CV_reduced_one_steps)
    # full_plug_in <- mean(CV_full_plug_ins)
    # reduced_plug_in <- mean(CV_reduced_plug_ins)
  }

  n_eff <- ifelse(sample_split, length(time)/2, length(time)) # for sample splitting
  cil <- one_step - stats::qnorm(1-alpha/2)*sqrt(var_est/n_eff)
  ciu <- one_step + stats::qnorm(1-alpha/2)*sqrt(var_est/n_eff)
  cil_1sided <- one_step - stats::qnorm(1-alpha)*sqrt(var_est/n_eff)
  p <- ifelse(sample_split,
                 stats::pnorm(one_step/sqrt(var_est/n_eff), lower.tail = FALSE),
                 NA)
  one_step <- ifelse(scale_est, max(c(one_step, 0)), one_step)
  cil <- ifelse(scale_est, max(c(cil, 0)), cil)
  ciu <- ifelse(scale_est, max(c(ciu, 0)), ciu)
  cil_1sided <- ifelse(scale_est, max(c(cil_1sided, 0)), cil_1sided)

  return(data.frame(restriction_time =  restriction_time,
                    est = one_step,
                    var_est = var_est,
                    cil = cil,
                    ciu = ciu,
                    cil_1sided = cil_1sided,
                    p = p,
                    large_predictiveness = full_one_step,
                    small_predictiveness = reduced_one_step))
}
