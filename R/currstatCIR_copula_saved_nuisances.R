#' Estimate a survival function under current status sampling with copula
#' @export
currstatCIR_copula_saved_nuisances <- function(time,
                                               event,
                                               X,
                                               mu_n,
                                               # g_n,
                                               Riemann_grid,
                                               f_sIx_n,
                                               f_s_n,
                                               F_sIx_n,
                                               eval_region,
                                               n_eval_pts = 101,
                                               theta = 0.5,
                                               copula){

  construct_g_n_saved <- function(f_sIx_n, f_s_n){
    fnc <- function(y,w){
      f_sIx_n(y,w) / f_s_n(y)
    }
    return(fnc)
  }
  # construct_F_sIx_n_saved <- function(f_sIx_n, Riemann_grid){
  #   uniq_y <- Riemann_grid[-2] # take out the first time in the haldensify grid (leave 0 in there)
  #   fnc <- function(y, w){
  #     F_sIx_ns <- sapply(uniq_y, function(y){
  #       f_sIx_n(y, as.numeric(w))
  #     })
  #     nbins <- sum(uniq_y <= y)
  #     last_bin <- which.max(uniq_y[uniq_y <= y])
  #     last_cutpoint <- max(uniq_y[uniq_y <= y])
  #     last_dens <- F_sIx_ns[last_bin]
  #     sum(F_sIx_ns[1:nbins] * diff(Riemann_grid)[1:nbins]) + last_dens * (y - last_cutpoint)
  #   }
  # }
  g_n <- construct_g_n_saved(f_sIx_n = f_sIx_n, f_s_n = f_s_n)
  # F_sIx_n <- construct_F_sIx_n_saved(f_sIx_n = f_sIx_n, Riemann_grid = Riemann_grid)

  s <- as.numeric(!is.na(event))

  time[s == 0] <- max(time, na.rm = TRUE)

  dat <- list(delta = event,
              y = time,
              s = s,
              w = X)

  dat$w <- data.frame(stats::model.matrix(stats::as.formula(paste("~",
                                                                  paste(names(dat$w),
                                                                        collapse =  "+"))),
                                          dat$w)[,-1])
  names(dat$w) <- paste("w", 1:ncol(dat$w), sep="")

  F_n <- stats::ecdf(dat$y)
  F_n_inverse <- function(t){
    stats::quantile(dat$y, probs = t, type = 1)
  }

  y_vals <- sort(unique(dat$y))

  if (copula == "clayton"){
    Gamma_n <- construct_Gamma_n_copula_clayton(dat=dat, mu_n=mu_n, g_n=g_n, F_sIx_n = F_sIx_n,
                                                Riemann_grid = Riemann_grid,
                                                theta = theta)
  } else if (copula == "frank"){
    Gamma_n <- construct_Gamma_n_copula_frank(dat=dat, mu_n=mu_n, g_n=g_n, F_sIx_n = F_sIx_n,
                                              Riemann_grid = Riemann_grid,
                                              theta = theta)
  }
  # kappa_n <- construct_kappa_n(dat = dat, mu_n = mu_n, g_n = g_n)

  # only estimate in the evaluation region, which doesn't include the upper bound
  # gcm_x_vals <- sapply(y_vals[y_vals <= eval_region[2]], F_n)
  gcm_x_vals <- sapply(y_vals[y_vals >= eval_region[1] & y_vals <= eval_region[2]], F_n)
  inds_to_keep <- !base::duplicated(gcm_x_vals)
  gcm_x_vals <- gcm_x_vals[inds_to_keep]
  # gcm_y_vals <- sapply(y_vals[y_vals <= eval_region[2]][inds_to_keep], Gamma_n)
  gcm_y_vals <- sapply(y_vals[y_vals >= eval_region[1] & y_vals <= eval_region[2]][inds_to_keep], Gamma_n)
  # check with avi here
  if (!any(gcm_x_vals==0)) {
    gcm_x_vals <- c(0, gcm_x_vals)
    gcm_y_vals <- c(0, gcm_y_vals)
  }
  gcm <- fdrtool::gcmlcm(x=gcm_x_vals, y=gcm_y_vals, type="gcm")
  theta_n <- stats::approxfun(
    x = gcm$x.knots[-length(gcm$x.knots)],
    y = gcm$slope.knots,
    method = "constant",
    rule = 2,
    f = 0
  )

  eval_cdf_upper <- mean(dat$y <= eval_region[2])
  eval_cdf_lower <- mean(dat$y <= eval_region[1])
  # theta_prime <- construct_deriv(r_Mn = theta_n,
  #                                deriv_method = deriv_method,
  #                                # y = seq(0, eval_cdf_upper, length.out = n_eval_pts))
  #                                y = seq(eval_cdf_lower, eval_cdf_upper, length.out = n_eval_pts))

  # Compute estimates
  # ests <- sapply(seq(0,eval_cdf_upper,length.out = n_eval_pts), theta_n)
  # deriv_ests <- sapply(seq(0, eval_cdf_upper, length.out = n_eval_pts), theta_prime)
  ests <- sapply(seq(eval_cdf_lower,eval_cdf_upper,length.out = n_eval_pts), theta_n)
  # deriv_ests <- sapply(seq(eval_cdf_lower, eval_cdf_upper, length.out = n_eval_pts), theta_prime)
  # kappa_ests <- sapply(y_vals, kappa_n)
  # transform kappa to quantile scale
  # kappa_ests_rescaled <- sapply(seq(eval_cdf_lower, eval_cdf_upper, length.out = n_eval_pts), function(x){
  # kappa_ests_rescaled <- sapply(seq(0, eval_cdf_upper, length.out = n_eval_pts), function(x){
  # ind <- which(y_vals == F_n_inverse(x))
  # kappa_ests[ind]
  # })
  # tau_ests <- deriv_ests * kappa_ests_rescaled
  # q <- ChernoffDist::qChern(p = alpha/2)
  # half_intervals <- sapply(1:n_eval_pts, function(x){
  # (4*tau_ests[x]/length(dat$y))^{1/3}*q
  # })
  # cils <- ests - half_intervals
  # cius <- ests + half_intervals

  ests[ests < 0] <- 0
  ests[ests > 1] <- 1
  # cils[cils < 0] <- 0
  # cils[cils > 1] <- 1
  # cius[cius < 0] <- 0
  # cius[cius > 1] <- 1

  ests <- Iso::pava(ests)
  # cils <- Iso::pava(cils)
  # cius <- Iso::pava(cius)

  results <- data.frame(
    t = sapply(seq(eval_cdf_lower, eval_cdf_upper, length.out = n_eval_pts), F_n_inverse),
    # t = sapply(seq(0, eval_cdf_upper, length.out = n_eval_pts), F_n_inverse),
    S_hat_est = 1-ests,
    S_hat_cil = NA,
    S_hat_ciu = NA
  )

  results <- results[results$t >= eval_region[1] & results$t <= eval_region[2],]
  rownames(results) <- NULL

  return(results)

}
