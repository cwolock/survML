#' Estimate a survival function under current status sampling
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
#' @param sample_split Logical indicating whether to perform inference using sample splitting
#' @param m Number of sample-splitting folds, defaults to 5.
#' @param eval_region Region over which to estimate the survival function.
#' @param n_eval_pts Number of points in grid on which to evaluate survival function.
#' The points will be evenly spaced, on the quantile scale, between the endpoints of \code{eval_region}.
#' @param alpha The level at which to compute confidence intervals and hypothesis tests. Defaults to 0.05
#' @param sensitivity_analysis Logical, whether to perform a copula-based sensitivity analysis. Defaults to \code{FALSE}
#' @param copula_control A named list of control parameters for the copula-based sensitivity analysis.
#' This should be a named list.
#'
#' @return List of data frames giving results. If not performing a sensitivity analysis, a single data frame is returned; if
#' performing a sensitivity analysis, a separate data frame will be returned for each value of the copula association parameter.
#' The results data frames have columns:
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
#' res <- survML::currstatCIR(time = dat$y,
#'                            event = dat$delta,
#'                            X = dat[,3:4],
#'                            SL_control = list(SL.library = c("SL.mean", "SL.glm"),
#'                                              V = 3),
#'                            HAL_control = list(n_bins = c(5),
#'                                               grid_type = c("equal_mass"),
#'                                               V = 3),
#'                            sensitivity_analysis = FALSE,
#'                            eval_region = eval_region)$primary_results
#'
#' xvals = res$t
#' yvals = res$S_hat_est
#' fn=stepfun(xvals, c(yvals[1], yvals))
#' plot.function(fn, from=min(xvals), to=max(xvals))}
#'
#' @export
currstatCIR <- function(time,
                        event,
                        X,
                        SL_control = list(SL.library = c("SL.mean", "SL.glm"),
                                          V = 3),
                        HAL_control = list(n_bins = c(5),
                                           grid_type = c("equal_mass"),
                                           V = 3),
                        deriv_method = "m-spline",
                        sample_split = FALSE,
                        m = 5,
                        eval_region,
                        n_eval_pts = 101,
                        alpha = 0.05,
                        sensitivity_analysis = FALSE,
                        copula_control = list(taus = c(-0.1, -0.05, 0.05, 0.1))){

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

  if (sample_split){
    folds <- sample(rep(seq_len(m), length = length(dat$y)))
  } else{
    folds <- rep(1, length(dat$y))
  }
  K <- length(unique(folds))
  y_vals <- sort(unique(dat$y))
  est_matrix <- matrix(NA, nrow = n_eval_pts, ncol = K)

  for (k in 1:K){
    dat_k <- list(delta = dat$delta[folds == k],
                  y = dat$y[folds == k],
                  s = dat$s[folds == k],
                  w = dat$w[folds == k,])
    # estimate conditional density (only among observed)
    cond_density_fit <- construct_f_sIx_n(dat = dat_k,
                                          HAL_control = HAL_control,
                                          SL_control = SL_control)
    f_sIx_n <- cond_density_fit$fnc
    Riemann_grid <- c(0, cond_density_fit$breaks)
    # estimate marginal density (marginalizing the conditional density over whole sample)
    f_s_n <- construct_f_s_n(dat = dat, f_sIx_n = f_sIx_n)
    # estimate density ratio
    g_n <- construct_g_n(f_sIx_n = f_sIx_n, f_s_n = f_s_n)

    # estimate outcome regression (only among observed)
    mu_n <- construct_mu_n(dat = dat_k, SL_control = SL_control, Riemann_grid = Riemann_grid)

    # Use the empirical cdf for F_n
    F_n <- stats::ecdf(dat_k$y)
    F_n_inverse <- function(t){
      stats::quantile(dat_k$y, probs = t, type = 1)
    }

    Gamma_n <- construct_Gamma_n(dat=dat_k, mu_n=mu_n, g_n=g_n,
                                 f_sIx_n = f_sIx_n, Riemann_grid = Riemann_grid)
    kappa_n <- construct_kappa_n(dat = dat_k, mu_n = mu_n, g_n = g_n)

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

    # Compute estimates
    ests <- sapply(seq(eval_cdf_lower,eval_cdf_upper,length.out = n_eval_pts), theta_n)
    est_matrix[,k] <- ests
  }

  if (!sample_split){
    ests <- est_matrix[,1]
    theta_prime <- construct_deriv(r_Mn = theta_n,
                                   deriv_method = deriv_method,
                                   # y = seq(0, eval_cdf_upper, length.out = n_eval_pts))
                                   y = seq(eval_cdf_lower, eval_cdf_upper, length.out = n_eval_pts))

    deriv_ests <- sapply(seq(eval_cdf_lower, eval_cdf_upper, length.out = n_eval_pts), theta_prime)
    kappa_ests <- sapply(y_vals, kappa_n)
    # transform kappa to quantile scale
    kappa_ests_rescaled <- sapply(seq(eval_cdf_lower, eval_cdf_upper, length.out = n_eval_pts), function(x){
      # kappa_ests_rescaled <- sapply(seq(0, eval_cdf_upper, length.out = n_eval_pts), function(x){
      ind <- which(y_vals == F_n_inverse(x))
      kappa_ests[ind]
    })
    tau_ests <- deriv_ests * kappa_ests_rescaled
    q <- ChernoffDist::qChern(p = alpha/2)
    half_intervals <- sapply(1:n_eval_pts, function(x){
      (4*tau_ests[x]/length(dat$y))^{1/3}*q
    })
    cils <- ests - half_intervals
    cius <- ests + half_intervals
  } else{
    avg_ests <- rowMeans(est_matrix)
    sigmas <- sapply(1:n_eval_pts, function(x) sqrt((1/(m-1)) * sum((est_matrix[x,] - avg_ests[x])^2)))
    q <- stats::qt(p = alpha/2, df = m-1)
    half_intervals <- sapply(1:n_eval_pts, function(x){
      sigmas[x]/(sqrt(m))*q
    })
    ests <- avg_ests
    cils <- ests - half_intervals
    cius <- ests + half_intervals
  }

  ests[ests < 0] <- 0
  ests[ests > 1] <- 1
  cils[cils < 0] <- 0
  cils[cils > 1] <- 1
  cius[cius < 0] <- 0
  cius[cius > 1] <- 1

  ests <- Iso::pava(ests)
  cils <- Iso::pava(cils)
  cius <- Iso::pava(cius)

  results <- data.frame(
    t = sapply(seq(eval_cdf_lower, eval_cdf_upper, length.out = n_eval_pts), F_n_inverse),
    # t = sapply(seq(0, eval_cdf_upper, length.out = n_eval_pts), F_n_inverse),
    S_hat_est = 1-ests,
    S_hat_cil = 1-cils,
    S_hat_ciu = 1-cius
  )

  results <- results[results$t >= eval_region[1] & results$t <= eval_region[2],]
  rownames(results) <- NULL

  if (sensitivity_analysis){
    sensitivity_results <- vector(mode = "list", length = length(copula_control$taus))
    names(sensitivity_results) <- as.character(copula_control$taus)
    F_sIx_n <- construct_F_sIx_n(dat = dat, SL_control = SL_control)
    for (k in 1:length(copula_control$taus)){
      tau <- copula_control$taus[k]
      tau_to_theta <- function(tau){
        thetas <- c(seq(-10, -0.01, by = 0.01), seq(0.01, 10, by = 0.01))
        taus <- rep(NA, length(thetas))
        for (i in 1:length(thetas)){
          mycop <- copula::frankCopula(param = thetas[i])
          taus[i] <- copula::tau(mycop)
        }
        return(thetas[which.min(abs(taus - tau))])
      }
      theta <- tau_to_theta(tau)
      Gamma_n <- construct_Gamma_n_copula(dat=dat, mu_n=mu_n, g_n=g_n, F_sIx_n = F_sIx_n,
                                          Riemann_grid = Riemann_grid,
                                          theta = theta)
      gcm_x_vals <- sapply(y_vals[y_vals >= eval_region[1] & y_vals <= eval_region[2]], F_n)
      inds_to_keep <- !base::duplicated(gcm_x_vals)
      gcm_x_vals <- gcm_x_vals[inds_to_keep]
      gcm_y_vals <- sapply(y_vals[y_vals >= eval_region[1] & y_vals <= eval_region[2]][inds_to_keep], Gamma_n)
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
      ests <- sapply(seq(eval_cdf_lower,eval_cdf_upper,length.out = n_eval_pts), theta_n)
      ests[ests < 0] <- 0
      ests[ests > 1] <- 1
      ests <- Iso::pava(ests)
      sensitivity_results[[k]] <- data.frame(
        t = sapply(seq(eval_cdf_lower, eval_cdf_upper, length.out = n_eval_pts), F_n_inverse),
        S_hat_est = 1-ests,
        tau = tau
      )
    }
  }

  if (sensitivity_analysis){
    all_results <- list(primary_results = results, sensitivity_results = sensitivity_results)

  } else{
    all_results <- list(primary_results = results)
  }
  return(all_results)

}

#' Estimate outcome regression
#' @noRd
construct_mu_n <- function(dat, SL_control, Riemann_grid) {
  # Construct newX (all distinct combinations of X and S)
  w_distinct <- dplyr::distinct(dat$w)
  w_distinct <- cbind("w_index"=c(1:nrow(w_distinct)), w_distinct)
  y_distinct <- sort(unique(round(c(dat$y, Riemann_grid), digits = 5)))
  newW <- expand.grid(w_index=w_distinct$w_index, Yprime=y_distinct)
  newW <- dplyr::inner_join(w_distinct, newW, by="w_index")
  newW$w_index <- NULL

  model_sl <- SuperLearner::SuperLearner(
    Y = dat$delta[dat$s == 1],
    X = cbind(dat$w[dat$s == 1,], Yprime=dat$y[dat$s == 1]),
    newX = newW,
    family = "binomial",
    method = "method.NNLS",
    SL.library = SL_control$SL.library,
    cvControl = list(V = SL_control$V)
  )
  pred <- as.numeric(model_sl$SL.predict)
  newW$index <- c(1:nrow(newW))

  fnc <- function(y,w) {
    cond <- paste0("round(Yprime,5)==",round(y,5))
    for (i in c(1:length(w))) {
      cond <- paste0(cond," & round(w",i,",5)==",round(w[i],5))
    }
    index <- (dplyr::filter(newW, eval(parse(text=cond))))$index
    if (length(index)!=1) {
      stop(paste0("y=",y,", w=c(",paste(w, collapse=","),")"))
    }
    return(pred[index])
  }

  return(fnc)

}

#' Estimate conditional density
#' @noRd
construct_f_sIx_n <- function(dat, HAL_control, SL_control){

  # fit hal
  haldensify_fit <- haldensify::haldensify(A = dat$y[dat$s == 1],
                                           W = dat$w[dat$s == 1,],
                                           n_bins = HAL_control$n_bins,
                                           grid_type = HAL_control$grid_type,
                                           cv_folds = HAL_control$V)

  w_distinct <- dplyr::distinct(dat$w)

  if (all(dat$s == 1)){
    binary_pred <- rep(1, nrow(w_distinct))
  } else{
    binary_fit <- SuperLearner::SuperLearner(
      Y = dat$s,
      X = dat$w,
      newX = w_distinct,
      family = "binomial",
      method = "method.NNLS",
      SL.library = SL_control$SL.library,
      cvControl = list(V = SL_control$V)
    )
    binary_pred <- as.numeric(binary_fit$SL.predict)
  }

  w_distinct <- cbind("w_index"=c(1:nrow(w_distinct)), w_distinct)
  # only get predictions at the breakpoints, since estimator is piecewise constant
  y_distinct <- haldensify_fit$breaks
  newW <- expand.grid(w_index=w_distinct$w_index, y=y_distinct)
  newW <- dplyr::inner_join(w_distinct, newW, by="w_index")
  newW$w_index <- NULL

  pred <- stats::predict(haldensify_fit, new_A = newW$y, new_W = newW[,-ncol(newW)])

  newW$index <- c(1:nrow(newW))

  breaks <- haldensify_fit$breaks

  fnc <- function(y,w) {
    left_y <- max(breaks[breaks <= max(y, min(breaks))])
    cond <- paste0("round(y,5)==",round(left_y,5))
    for (i in c(1:length(w))) {
      cond <- paste0(cond," & round(w",i,",5)==",round(w[i],5))
    }
    index <- (dplyr::filter(newW, eval(parse(text=cond))))$index
    dens_pred <- pred[index]
    for (i in c(1:length(w))) {
      if (i == 1){
        cond <- paste0("round(w1,5) == ", round(w[i],5))
      } else{
        cond <- paste0(cond," & round(w",i,",5)==",round(w[i],5))
      }
    }
    index <- (dplyr::filter(w_distinct, eval(parse(text=cond))))$w_index
    binary_pred <- binary_pred[index]
    return(dens_pred * binary_pred)
  }

  return(list(fnc = fnc, breaks = haldensify_fit$breaks))
}

#' Estimate marginal density
#' @noRd
construct_f_s_n <- function(dat, f_sIx_n){
  uniq_y <- sort(unique(dat$y))
  f_sIx_ns <- sapply(uniq_y, function(y){
    mean(apply(dat$w, 1, function(w){
      f_sIx_n(y,as.numeric(w))
    }))})

  fnc <- function(y){
    f_sIx_ns[uniq_y == y]
  }
}

#' Estimate density ratio
#' @noRd
construct_g_n <- function(f_sIx_n, f_s_n){
  function(y,w){
    f_sIx_n(y,w) / f_s_n(y)
  }
}

#' Estimate primitive
#' @noRd
construct_Gamma_n <- function(dat, mu_n, g_n, f_sIx_n, Riemann_grid) {
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

  # piece 1 maps to (\Delta - \mu_n(Y_i, W_i))/g_n(Y_i, W_i)
  piece_1 <- (dat$delta-mu_ns) / g_ns
  piece_1[is.na(piece_1)] <- 0 # there are NAs for missing values, but these get
  # multiplied by 0 later anyway

  # often there aren't so many unique monitoring times, and we can save a lot of
  # time by only computing piece_2 on the unique values
  unique_y <- sort(unique(dat$y))
  uniq_w <- dplyr::distinct(dat$w)

  unique_mus <- lapply(unique_y, function(y){
    unique_mus <- apply(uniq_w, 1, function(w) { mu_n(y=y, w=as.numeric(w)) })
    unique_mus
  })

  unique_piece_2 <- sapply(unique_y, function(y){
    this_mus <- unique_mus[[which(y == unique_y)]]
    mus <- apply(dat$w, 1, function(w) {this_mus[which(colSums(t(uniq_w) == w) == ncol(uniq_w))]})
    mean(mus)
  })

  # match to pre-computed values
  # piece 2 maps to \theta_n(Y_i)
  piece_2 <- sapply(dat$y, function(y) {
    unique_piece_2[unique_y == y]
  })

  fnc <- function(y) {
    piece_3 <- as.integer(dat$y<=y) * dat$s
    obs_pieces <- mean(piece_3*piece_1) + mean(piece_3*piece_2)
    return(obs_pieces)
  }

  return(fnc)
}

#' Estimate primitive
#' @noRd
construct_Gamma_n_copula <- function(dat, mu_n, g_n, F_sIx_n,
                                     Riemann_grid, theta) {

  m <- function(u,v){
    -(1/theta)*log(1 - (u*(1 - exp(-theta)))/(exp(-theta * v) + u*(1 - exp(-theta*v))))
  }

  partial_u <- function(u,v){
    -(1/theta) * (exp(-theta * v)*(exp(-theta) - 1))/((exp(-theta*v) + u*(exp(-theta)-exp(-theta*v)))*(exp(-theta*v) + u*(1 - exp(-theta*v))))
  }

  partial_v <- function(u,v){
    (u*exp(-theta*v)*(1-u)*(1-exp(-theta)))/((exp(-theta*v) + u*(exp(-theta)-exp(-theta*v)))*(exp(-theta*v) + u*(1 - exp(-theta*v))))
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

  partial_us <- apply(cbind(mu_ns, F_ns),
                      MARGIN = 1,
                      FUN = function(r){
                        partial_u(r[1], r[2])
                      })

  # will later be multiplied by indicator
  piece_1 <- (dat$delta - mu_ns) / g_ns * partial_us
  piece_1[is.na(piece_1)] <- 0

  # often there aren't so many unique monitoring times, and we can save a lot of
  # time by only computing piece_2 on the unique values
  unique_y <- sort(unique(dat$y))
  uniq_w <- dplyr::distinct(dat$w)

  unique_mus <- lapply(unique_y, function(y){
    unique_mus <- apply(uniq_w, 1, function(w) { mu_n(y=y, w=as.numeric(w)) })
    unique_mus
  })

  unique_Fs <- lapply(unique_y, function(y){
    unique_Fs <- apply(uniq_w, 1, function(w) { F_sIx_n(y=y, w=as.numeric(w)) })
    unique_Fs
  })

  unique_lambda <- sapply(unique_y, function(y) {
    this_mus <- unique_mus[[which(y == unique_y)]]
    this_Fs <- unique_Fs[[which(y == unique_y)]]
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
    indicators <- (dat$y <= y)
    mart <- indicators - Fs
    partial_vs <- apply(cbind(mus, Fs),
                        MARGIN = 1,
                        FUN = function(r){
                          partial_v(r[1], r[2])
                        })
    mean(partial_vs * mart)
  })

  piece_3 <- sapply(dat$y, function(y) {
    unique_piece_3[unique_y == y]})

  fnc <- function(y) {
    piece_4 <- as.integer(dat$y<=y) * dat$s
    obs_pieces <- mean(piece_4 * piece_1) + mean(piece_4 * piece_2) + mean(piece_4*piece_3)
    return(obs_pieces)
  }

  return(fnc)
}


#' Estimate part of scale factor
#' @noRd
construct_kappa_n <- function(dat, mu_n, g_n){
  fnc <- function(y){
    mean(apply(dat$w, 1, function(w) {
      numer1 <- mu_n(y = y, w = as.numeric(w))
      numer <- numer1 * (1 - numer1)
      denom1 <- g_n(y = y, w = as.numeric(w))
      numer / denom1
    }))
  }
  return(fnc)
}

#' Estimate derivative
#' @noRd
construct_deriv <- function(deriv_method="m-spline", r_Mn, y) {
  # r_Mn should be the theta_n function
  # grid should be something like seq(0, 1, by = 0.01)
  if (deriv_method=="line") {
    fnc <- function(y) { r_Mn(1)-r_Mn(0) }
  } else {
    if (r_Mn(0)==r_Mn(1)) {
      fnc <- function(y) { 0 }
      warning("Estimated function is flat; variance estimation not possible.")
    }
    # Estimate entire function on grid
    r_Mns <- sapply(y, r_Mn)

    # Compute set of midpoints of jump points (plus endpoints)
    points_x <- y[1]
    points_y <- r_Mns[1]
    for (i in 2:length(y)) {
      if (r_Mns[i]-r_Mns[round(i-1)]!=0) {
        points_x <- c(points_x, (y[i]+y[round(i-1)])/2)
        points_y <- c(points_y, mean(c(r_Mns[i],r_Mns[round(i-1)])))
      }
    }
    points_x <- c(points_x, y[length(y)])
    points_y <- c(points_y, r_Mns[length(y)])

    if (deriv_method=="linear") {
      fnc_pre <- stats::approxfun(x=points_x, y=points_y, method="linear", rule=2)
    }

    if (deriv_method=="m-spline") {
      fnc_pre <- stats::splinefun(x=points_x, y=points_y, method="monoH.FC")
    }

    # Construct numerical derivative
    fnc <- function(y) {
      width <- 0.1 # Note: may need to play around with this value
      x1 <- y - width/2
      x2 <- y + width/2
      if (x1<0) { x2<-width; x1<-0; }
      if (x2>1) { x1<-1-width; x2<-1; }
      return(max((fnc_pre(x2)-fnc_pre(x1))/width,0))
    }
  }
  return(fnc)
}

#' Estimate the conditional cdf
#' @noRd
construct_F_sIx_n <- function(dat, SL_control){

  newX <- dplyr::distinct(dat$w)
  newtimes <- sort(unique(dat$y))
  fit <- stackG(time = dat$y,
                X = dat$w,
                newX = newX,
                newtimes = newtimes,
                SL_control = SL_control,
                bin_size = 0.05,
                time_basis = "continuous",
                time_grid_approx = stats::quantile(dat$y, probs = seq(0, 1, by = 0.01)))
  preds <- 1-fit$S_T_preds

  fnc <- function(y, w){
    return(preds[which(colSums(t(newX) == w) == ncol(newX)),which(newtimes == y)])
  }
}

#' Convert list data to df
#' @noRd
as_df <- function(dat) {
  cbind(dat$w, y=dat$y)
}
