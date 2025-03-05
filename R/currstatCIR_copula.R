#' Estimate a survival function under current status sampling with copula
#'
#' @param time \code{n x 1} numeric vector of observed monitoring times. For individuals that were never
#' monitored, this can be set to any arbitrary value, including \code{NA}, as long as the corresponding
#' \code{event} variable is \code{NA}.
#' @param event \code{n x 1} numeric vector of status indicators of
#' whether an event was observed prior to the monitoring time. This value must be \code{NA} for
#' individuals that were never monitored.
#' @param X \code{n x p} dataframe of observed covariate values.
#' @param SL_control List of \code{SuperLearner} control parameters. This should be a named list; see
#' \code{SuperLearner} documentation for further information.
#' @param HAL_control List of \code{haldensify} control parameters. This should be a named list; see
#' \code{haldensify} documentation for further information.
#' @param deriv_method Method for computing derivative. Options are \code{"m-spline"} (the default,
#' fit a smoothing spline to the estimated function and differentiate the smooth approximation),
#' \code{"linear"} (linearly interpolate the estimated function and use the slope of that line), and
#' \code{"line"} (use the slope of the line connecting the endpoints of the estimated function).
#' @param eval_region Region over which to estimate the survival function.
#' @param n_eval_pts Number of points in grid on which to evaluate survival function.
#' The points will be evenly spaced, on the quantile scale, between the endpoints of \code{eval_region}.
#' @param alpha The level at which to compute confidence intervals and hypothesis tests. Defaults to 0.05
#'
#' @return Data frame giving results, with columns:
#' \item{t}{Time at which survival function is estimated}
#' \item{S_hat_est}{Survival function estimate}
#' \item{S_hat_cil}{Lower bound of confidence interval}
#' \item{S_hat_ciu}{Upper bound of confidence interval}
#'
#' @examples
#' \dontrun{# This is a small simulation example
#' set.seed(123)
#' n <- 300
#' x <- cbind(2*rbinom(n, size = 1, prob = 0.5)-1,
#'            2*rbinom(n, size = 1, prob = 0.5)-1)
#' t <- rweibull(n,
#'               shape = 0.75,
#'               scale = exp(0.4*x[,1] - 0.2*x[,2]))
#' y <- rweibull(n,
#'               shape = 0.75,
#'               scale = exp(0.4*x[,1] - 0.2*x[,2]))
#'
#' # round y to nearest quantile of y, just so there aren't so many unique values
#' quants <- quantile(y, probs = seq(0, 1, by = 0.05), type = 1)
#' for (i in 1:length(y)){
#'   y[i] <- quants[which.min(abs(y[i] - quants))]
#' }
#' delta <- as.numeric(t <= y)
#'
#' dat <- data.frame(y = y, delta = delta, x1 = x[,1], x2 = x[,2])
#'
#' dat$delta[dat$y > 1.8] <- NA
#' dat$y[dat$y > 1.8] <- NA
#' eval_region <- c(0.05, 1.5)
#' res <- survML::currstatCIR_copula(time = dat$y,
#'                            event = dat$delta,
#'                            X = dat[,3:4],
#'                            SL_control = list(SL.library = c("SL.mean", "SL.glm"),
#'                                              V = 3),
#'                            HAL_control = list(n_bins = c(5),
#'                                               grid_type = c("equal_mass"),
#'                                               V = 3),
#'                            eval_region = eval_region)
#'
#' xvals = res$t
#' yvals = res$S_hat_est
#' fn=stepfun(xvals, c(yvals[1], yvals))
#' plot.function(fn, from=min(xvals), to=max(xvals))}
#'
#' @export
currstatCIR_copula <- function(time,
                        event,
                        X,
                        SL_control = list(SL.library = c("SL.mean", "SL.glm"),
                                          V = 3),
                        HAL_control = list(n_bins = c(5),
                                           grid_type = c("equal_mass"),
                                           V = 3),
                        deriv_method = "m-spline",
                        # missing_method = "extended",
                        eval_region,
                        n_eval_pts = 101,
                        alpha = 0.05,
                        theta = 0.5){

  s <- as.numeric(!is.na(event))

  time[s == 0] <- max(time, na.rm = TRUE)

  # if (missing_method == "cc"){
  #   time <- time[s == 1]
  #   event <- event[s == 1]
  #   W <- W[s == 1,]
  #   s <- s[s == 1]
  # }

  dat <- list(delta = event,
              y = time,
              s = s,
              w = X)

  dat$w <- data.frame(stats::model.matrix(stats::as.formula(paste("~",
                                                                  paste(names(dat$w),
                                                                        collapse =  "+"))),
                                          dat$w)[,-1])
  names(dat$w) <- paste("w", 1:ncol(dat$w), sep="")

  # estimate conditional density (only among observed)
  cond_density_fit <- construct_f_sIx_n(dat = dat,
                                        HAL_control = HAL_control)
  f_sIx_n <- cond_density_fit$fnc
  Riemann_grid <- c(0, cond_density_fit$breaks)
  # estimate marginal density (marginalizing the conditional density over whole sample)
  f_s_n <- construct_f_s_n(dat = dat, f_sIx_n = f_sIx_n)
  # estimate density ratio
  g_n <- construct_g_n(f_sIx_n = f_sIx_n, f_s_n = f_s_n)
  # estimate conditional CDF of response times
  F_sIx_n <- construct_F_sIx_n(dat = dat, f_sIx_n = f_sIx_n, Riemann_grid = Riemann_grid)

  # estimate outcome regression (only among observed)
  mu_n <- construct_mu_n(dat = dat, SL_control = SL_control, Riemann_grid = Riemann_grid)
  # mu_n <- construct_mu_n(dat = dat, SL_control = SL_control, Riemann_grid = Riemann_grid)

  y_vals <- sort(unique(dat$y))

  # Use the empirical cdf for F_n
  F_n <- stats::ecdf(dat$y)
  F_n_inverse <- function(t){
    stats::quantile(dat$y, probs = t, type = 1)
  }

  Gamma_n <- construct_Gamma_n_copula(dat=dat, mu_n=mu_n, g_n=g_n, F_sIx_n = F_sIx_n,
                               Riemann_grid = Riemann_grid,
                               theta = theta)
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


#' Estimate primitive
#' @noRd
construct_Gamma_n_copula <- function(dat, mu_n, g_n, F_sIx_n,
                              Riemann_grid, theta) {

  m <- function(u,v){
    (1 + v^(-theta)*(u^(-theta/(theta+1))-1))^(-1/theta)
  }

  n_orig <- length(dat$y)
  dim_w <- length(dat$w)
  mu_ns <- apply(as_df(dat), 1, function(r) {
    y <- r[["y"]]
    w <- as.numeric(r)[1:dim_w]
    mu_n(y=y, w=w)
  })

  g_ns <- apply(as_df(dat), 1, function(r) {
    y <- r[["y"]]
    w <- as.numeric(r)[1:dim_w]
    g_n(y=y, w=w)
  })

  F_ns <- apply(as_df(dat), 1, function(r) {
    y <- r[["y"]]
    w <- as.numeric(r)[1:dim_w]
    F_sIx_n(y = y, w = w)
  })

  # ms <- apply(cbind(mu_ns, F_ns),
  #             MARGIN = 1,
  #             FUN = function(r){
  #               m(r[1], r[2])
  #             })

  m_thetas <- apply(cbind(mu_ns, F_ns),
                    MARGIN = 1,
                    FUN = function(r){
                      m(r[1], r[2])^(theta + 1)
                    })

  # will later be multiplied by indicator
  piece_1 <- (dat$delta - mu_ns) / g_ns / (1 + theta) * m_thetas * (F_ns)^(-theta) * mu_ns^(-(1 + 2*theta)/(1 + theta))
  piece_1[is.na(piece_1)] <- 0
  # piece_2 <-

  print("made it past piece 1")

  # piece 1 maps to (\Delta - \mu_n(Y_i, W_i))/g_n(Y_i, W_i)
  # piece_1 <- (dat$delta-mu_ns) / g_ns
  # piece_1[is.na(piece_1)] <- 0 # there are NAs for missing values, but these get
  # multiplied by 0 later anyway

  # often there aren't so many unique monitoring times, and we can save a lot of
  # time by only computing piece_2 on the unique values
  unique_y <- sort(unique(dat$y))
  uniq_w <- dplyr::distinct(dat$w)
  print(length(unique_y))
  print(nrow(uniq_w))

  unique_mus <- lapply(unique_y, function(y){
    unique_mus <- apply(uniq_w, 1, function(w) { mu_n(y=y, w=as.numeric(w)) })
    unique_mus
  })

  unique_Fs <- lapply(unique_y, function(y){
    unique_Fs <- apply(uniq_w, 1, function(w) { F_sIx_n(y=y, w=as.numeric(w)) })
    unique_Fs
  })

  unique_lambda <- sapply(unique_y, function(y) {
    print(y)
    # uniq_w <- dplyr::distinct(dat$w)
    # unique_mus <- apply(uniq_w, 1, function(w) { mu_n(y=y, w=as.numeric(w)) })
    # unique_Fs <- apply(uniq_w, 1, function(w) { F_sIx_n(y = y, w = as.numeric(w))})
    this_mus <- unique_mus[[which(y == unique_y)]]
    this_Fs <- unique_Fs[[which(y == unique_y)]]
    # mus <- apply(dat$w, 1, function(w) { mu_n(y=y, w=as.numeric(w)) })
    # Fs <- apply(dat$w, 1, function(w) { F_sIx_n(y = y, w = as.numeric(w))})
    mus <- apply(dat$w, 1, function(w) {this_mus[which(colSums(t(uniq_w) == w) == ncol(uniq_w))]})
    Fs <- apply(dat$w, 1, function(w) {this_Fs[which(colSums(t(uniq_w) == w) == ncol(uniq_w))]})
    ms <- apply(cbind(mus, Fs),
                MARGIN = 1,
                FUN = function(r){
                  m(r[1], r[2])
                })
    mean(ms)
  })

  piece_2 <- sapply(dat$y, function(y) {
    unique_lambda[unique_y == y]
  })

  unique_piece_3 <- sapply(unique_y, function(y) {
    this_mus <- unique_mus[[which(y == unique_y)]]
    this_Fs <- unique_Fs[[which(y == unique_y)]]
    mus <- apply(dat$w, 1, function(w) {this_mus[which(colSums(t(uniq_w) == w) == ncol(uniq_w))]})
    Fs <- apply(dat$w, 1, function(w) {this_Fs[which(colSums(t(uniq_w) == w) == ncol(uniq_w))]})
    # mus <- apply(dat$w, 1, function(w) { mu_n(y=y, w=as.numeric(w)) })
    # Fs <- apply(dat$w, 1, function(w) { F_sIx_n(y = y, w = as.numeric(w))})
    indicators <- (dat$y <= y)
    mart <- indicators - Fs
    ms <- apply(cbind(mus, Fs),
                MARGIN = 1,
                FUN = function(r){
                  m(r[1], r[2])^(theta + 1)
                })
    mean(ms * Fs^(-theta - 1)*(mus^(-theta/(theta + 1)) - 1) * mart)
  })

  piece_3 <- sapply(dat$y, function(y) {
    unique_piece_3[unique_y == y]})


  # unique_piece_2 <- sapply(unique_y, function(y) {
  #   mean(apply(dat$w, 1, function(w) { mu_n(y=y, w=as.numeric(w)) }))
  # })
  #
  # # match to pre-computed values
  # # piece 2 maps to \theta_n(Y_i)
  # piece_2 <- sapply(dat$y, function(y) {
  #   unique_piece_2[unique_y == y]
  # })
  # w_distinct <- dplyr::distinct(dat$w)
  # if there actually is missingness to deal with
  # if (sum(dat$s) != length(dat$s)){
  #   piece_4 <- sapply(Riemann_grid, function(y) {
  #     mean(apply(dat$w, 1, function(w) { mu_n(y=y, w=as.numeric(w)) }))
  #   })
  #
  #   # calculate the Riemann integral w.r.t. the conditional density
  #   # for each unique value of w
  #   unique_Riemann_integrals <- apply(w_distinct, MARGIN = 1, function(w){
  #     density_vals <- sapply(Riemann_grid, function(y){
  #       f_sIx_n(y = y, w = as.numeric(w))
  #     })
  #     Riemann_integrals <- sapply(Riemann_grid, function(y){
  #       sum(diff(Riemann_grid[Riemann_grid <= y]) * piece_4[Riemann_grid <= y][-1] * density_vals[Riemann_grid <= y][-1])
  #     })
  #     Riemann_integrals
  #   })
  #
  #   unique_Riemann_integrals <- t(unique_Riemann_integrals)
  # } else{
  #   unique_Riemann_integrals <- matrix(1, nrow = nrow(w_distinct), ncol = length(Riemann_grid))
  # }

  fnc <- function(y) {
    piece_4 <- as.integer(dat$y<=y) * dat$s
    obs_pieces <- mean(piece_4 * piece_1) + mean(piece_4 * piece_2) + mean(piece_4*piece_3)
    # obs_pieces <- (sum(piece_3*piece_1) + mean(piece_3*piece_2))/n_orig
    return(obs_pieces)
  }

  return(fnc)
}

#' Estimate the conditional cdf
#' @noRd
construct_F_sIx_n <- function(dat, f_sIx_n, Riemann_grid){
  uniq_y <- Riemann_grid[-2] # take out the first time in the haldensify grid (leave 0 in there)


  fnc <- function(y, w){
    F_sIx_ns <- sapply(uniq_y, function(y){
      f_sIx_n(y, as.numeric(w))
    })
    nbins <- sum(uniq_y <= y)
    if (nbins == 0){
      print(y)
      print(uniq_y)
    }
    last_bin <- which.max(uniq_y[uniq_y <= y])
    last_cutpoint <- max(uniq_y[uniq_y <= y])
    last_dens <- F_sIx_ns[last_bin]


    sum(F_sIx_ns[1:nbins] * diff(Riemann_grid)[1:nbins]) + last_dens * (y - last_cutpoint)
  }

  # fnc <- function(y, w){
  #   pweibull(y, shape = 0.75, scale = exp(0.4*w[1] - 0.2*w[2] + 0.1*w[3]))
  # }
}

#' Estimate outcome regression
#' @noRd
construct_mu_n_theta <- function(dat, SL_control, Riemann_grid, theta) {

  fnc <- function(y, w){
    weib_scale <- exp(0.4*w[1] - 0.2*w[2] + 0.1*w[3])
    F_Y_of_y <- pweibull(y,
                         shape = 0.75,
                         scale = weib_scale)
    F_T_of_y <- pweibull(y,
                         shape = 0.75,
                         scale = weib_scale)
    cutoff <- -(1/theta)
    inside <- (1/theta)*(F_T_of_y^(-theta) + F_Y_of_y^(-theta) - 2)
    if (theta < 0 & inside >= cutoff){
      ret <- 0
    } else{
      ret <-  ((F_Y_of_y)^(-theta - 1))/((F_T_of_y^(-theta) + F_Y_of_y^(-theta) - 1)^((theta + 1)/theta))
    }
    ret
  }


  return(fnc)

}


