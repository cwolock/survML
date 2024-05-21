#' Estimate restricted prediction time MSE VIM
#'
#' @param time \code{n x 1} numeric vector of observed
#' follow-up times If there is censoring, these are the minimum of the
#' event and censoring times.
#' @param event \code{n x 1} numeric vector of status indicators of
#' whether an event was observed. Defaults to a vector of 1s, i.e. no censoring.
#' @param approx_times Numeric vector of length J1 giving times at which to
#' approximate integrals.
#' @param tau restriction time
#' @param f_hat Full oracle predictions (n x J1 matrix)
#' @param fs_hat Residual oracle predictions (n x J1 matrix)
#' @param S_hat Estimates of conditional event time survival function (n x J2 matrix)
#' @param G_hat Estimate of conditional censoring time survival function (n x J2 matrix)
#' @param folds Numeric vector of length n giving cross-fitting folds
#' @param sample_split Logical indicating whether or not to sample split
#' @param ss_folds Numeric vector of length n giving sample-splitting folds
#' @param scale_est Logical, whether or not to force the VIM estimate to be nonnegative
#'
#' @return data frame giving results
#'
#' @export

vim_rmst_mse <- function(time,
                         event,
                         approx_times,
                         tau,
                         f_hat,
                         fs_hat,
                         S_hat,
                         G_hat,
                         folds,
                         sample_split,
                         ss_folds,
                         scale_est = FALSE){

  V <- length(unique(folds))

  CV_full_plug_ins <- rep(NA, V)
  CV_reduced_plug_ins <- rep(NA, V)
  CV_full_one_steps <- rep(NA, V)
  CV_reduced_one_steps <- rep(NA, V)
  CV_one_steps <- rep(NA, V)
  CV_plug_ins <- rep(NA, V)
  CV_var_ests <- rep(NA, V)
  split_one_step_fulls <- rep(NA, V)
  split_plug_in_fulls <- rep(NA, V)
  split_one_step_reduceds <- rep(NA, V)
  split_plug_in_reduceds <- rep(NA, V)
  split_var_est_fulls <- rep(NA, V)
  split_var_est_reduceds <- rep(NA, V)
  for (j in 1:V){

    time_holdout <- time[folds == j]
    event_holdout <- event[folds == j]

    V_0 <- estimate_rmst_mse(time = time_holdout,
                             event = event_holdout,
                             approx_times = approx_times,
                             tau = tau,
                             preds = f_hat[[j]],
                             S_hat = S_hat[[j]],
                             G_hat = G_hat[[j]])
    V_0s <- estimate_rmst_mse(time = time_holdout,
                              event = event_holdout,
                              approx_times = approx_times,
                              tau = tau,
                              preds = fs_hat[[j]],
                              S_hat = S_hat[[j]],
                              G_hat = G_hat[[j]])

    CV_full_one_steps[j] <- V_0$one_step
    CV_full_plug_ins[j] <- V_0$plug_in
    CV_reduced_one_steps[j] <- V_0s$one_step
    CV_reduced_plug_ins[j] <- V_0s$plug_in
    CV_one_steps[j] <-  V_0$one_step -V_0s$one_step
    CV_plug_ins[j] <-  V_0$plug_in -V_0s$plug_in
    split_one_step_fulls[j] <- V_0$one_step
    split_plug_in_fulls[j] <- V_0$plug_in
    split_one_step_reduceds[j] <- V_0s$one_step
    split_plug_in_reduceds[j] <- V_0s$plug_in
    split_var_est_fulls[j] <- mean(V_0$EIF^2)
    split_var_est_reduceds[j] <- mean(V_0s$EIF^2)
    EIF <- V_0$EIF - V_0s$EIF
    CV_var_ests[j] <- mean(EIF^2)
  }

  if (sample_split){
    folds_0 <- sort(unique(folds[ss_folds == 0]))
    folds_1 <- sort(unique(folds[ss_folds == 1]))
    one_step <- mean(split_one_step_fulls[folds_0]) -
      mean(split_one_step_reduceds[folds_1])
    full_one_step <- mean(split_one_step_fulls[folds_0])
    full_plug_in <- mean(split_plug_in_fulls[folds_0])
    reduced_one_step <- mean(split_one_step_reduceds[folds_1])
    reduced_plug_in <- mean(split_plug_in_reduceds[folds_1])
    var_est <- mean(split_var_est_fulls[folds_0]) +
      mean(split_var_est_reduceds[folds_1])
  } else{
    one_step <- mean(CV_one_steps)
    var_est <- mean(CV_var_ests)
    full_one_step <- mean(CV_full_one_steps)
    reduced_one_step <- mean(CV_reduced_one_steps)
    full_plug_in <- mean(CV_full_plug_ins)
    reduced_plug_in <- mean(CV_reduced_plug_ins)
  }

  n_eff <- ifelse(sample_split, length(time)/2, length(time)) # for sample splitting
  cil <- one_step - 1.96*sqrt(var_est/n_eff)
  ciu <- one_step + 1.96*sqrt(var_est/n_eff)
  cil_1sided <- one_step - 1.645*sqrt(var_est/n_eff)
  p <- ifelse(sample_split,
                 stats::pnorm(one_step/sqrt(var_est/n_eff), lower.tail = FALSE),
                 NA)
  one_step <- ifelse(scale_est, max(c(one_step, 0)), one_step)
  cil <- ifelse(scale_est, max(c(cil, 0)), cil)
  ciu <- ifelse(scale_est, max(c(ciu, 0)), ciu)
  cil_1sided <- ifelse(scale_est, max(c(cil_1sided, 0)), cil_1sided)

  return(data.frame(restriction_time = tau,
                    # full_one_step = full_one_step,
                    # reduced_one_step = reduced_one_step,
                    est = one_step,
                    # full_plug_in = full_plug_in,
                    # reduced_plug_in = reduced_plug_in,
                    var_est = var_est,
                    cil = cil,
                    ciu = ciu,
                    cil_1sided = cil_1sided,
                    p = p))
}
