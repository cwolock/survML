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
#' @param cf_folds Numeric vector of length n giving cross-fitting folds
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

vim_AUC <- function(time,
                    event,
                    approx_times,
                    landmark_times,
                    f_hat,
                    fs_hat,
                    S_hat,
                    G_hat,
                    cf_folds,
                    sample_split,
                    ss_folds,
                    robust = TRUE,
                    scale_est = FALSE,
                    alpha = 0.05){

  J1 <- length(landmark_times)
  V <- length(unique(cf_folds))
  one_step <- rep(NA, J1)
  var_est <- rep(NA, J1)
  full_one_step <- rep(NA, J1)
  reduced_one_step <- rep(NA, J1)
  # full_plug_in <- rep(NA, J1)
  # reduced_plug_in <- rep(NA, J1)
  cil <- rep(NA, J1)
  ciu <- rep(NA, J1)
  cil_1sided <- rep(NA, J1)
  p <- rep(NA, J1)
  for(i in 1:J1) {
    tau <- landmark_times[i]
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

      V_0 <- estimate_AUC(time = time_holdout,
                          event = event_holdout,
                          approx_times = approx_times,
                          tau = tau,
                          preds = f_hat[[j]][,i],
                          S_hat = S_hat[[j]],
                          G_hat = G_hat[[j]],
                          robust = robust)
      V_0s <- estimate_AUC(time = time_holdout,
                           event = event_holdout,
                           approx_times = approx_times,
                           tau = tau,
                           preds = fs_hat[[j]][,i],
                           S_hat = S_hat[[j]],
                           G_hat = G_hat[[j]],
                           robust = robust)

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
      one_step[i] <- mean(split_numerator_fulls[folds_0])/mean(split_denominator_fulls[folds_0]) -
        mean(split_numerator_reduceds[folds_1])/mean(split_denominator_reduceds[folds_1])
      full_one_step[i] <- mean(split_numerator_fulls[folds_0])/mean(split_denominator_fulls[folds_0])
      # full_plug_in[i] <- mean(split_plug_in_fulls[folds_0])
      reduced_one_step[i] <- mean(split_numerator_reduceds[folds_1])/mean(split_denominator_reduceds[folds_1])
      # reduced_plug_in[i] <- mean(split_plug_in_reduceds[folds_1])
      var_est[i] <- mean(split_var_est_fulls[folds_0]) +
        mean(split_var_est_reduceds[folds_1])
    } else{
      one_step[i] <- mean(CV_full_numerators - CV_reduced_numerators)/mean(CV_denominators)#mean(CV_one_steps)
      var_est[i] <- mean(CV_var_ests)
      full_one_step[i] <- mean(CV_full_numerators)/mean(CV_denominators)#mean(CV_full_one_steps)
      reduced_one_step[i] <- mean(CV_reduced_numerators)/mean(CV_denominators)#mean(CV_reduced_one_steps)
      # full_plug_in[i] <- mean(CV_full_plug_ins)
      # reduced_plug_in[i] <- mean(CV_reduced_plug_ins)
    }
    n_eff <- ifelse(sample_split, length(time)/2, length(time)) # for sample splitting
    cil[i] <- one_step[i] - qnorm(1-alpha/2)*sqrt(var_est[i]/n_eff)
    ciu[i] <- one_step[i] + qnorm(1-alpha/2)*sqrt(var_est[i]/n_eff)
    cil_1sided[i] <- one_step[i] - qnorm(1-alpha)*sqrt(var_est[i]/n_eff)
    p[i] <- ifelse(sample_split,
                   stats::pnorm(one_step[i]/sqrt(var_est[i]/n_eff), lower.tail = FALSE),
                   NA)
    one_step[i] <- ifelse(scale_est, max(c(one_step[i], 0)), one_step[i])
    cil[i] <- ifelse(scale_est, max(c(cil[i], 0)), cil[i])
    ciu[i] <- ifelse(scale_est, max(c(ciu[i], 0)), ciu[i])
    cil_1sided[i] <- ifelse(scale_est, max(c(cil_1sided[i], 0)), cil_1sided[i])
  }

  return(data.frame(landmark_time = landmark_times,
                    est = one_step,
                    var_est = var_est,
                    cil = cil,
                    ciu = ciu,
                    cil_1sided = cil_1sided,
                    p = p,
                    large_predictiveness = full_one_step,
                    small_predictiveness = reduced_one_step))
}
