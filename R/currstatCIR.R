#' Estimate a survival function under current status sampling
#'
#'
#' @export
currstatCIR <- function(time,
                        event,
                        W,
                        direction = "increasing",
                        SL_control = list(SL.library = c("SL.mean"),
                                          V = 5,
                                          method = "method.NNLS"),
                        HAL_control = list(n_bins = c(5,10),
                                           grid_type = c("equal_range", "equal_mass"),
                                           V = 5),
                        deriv_method = "m-spline",
                        missing_method = "extended",
                        n_eval_pts = 101){

  s <- ifelse(is.na(time) & is.na(event), 0, 1)

  event[s == 0] <- 0
  time[s == 0] <- 0

  if (missing_method == "cc"){
    time <- time[s == 1]
    event <- time[s == 1]
    W <- W[s == 1,]
  }
  dat <- list(delta = event,
              y = time,
              s = s,
              w = W)

  dat$w <- data.frame(model.matrix(as.formula(paste("~",
                                                    paste(names(dat$w),
                                                          collapse =  "+"))),
                                   dat$w)[,-1])
  names(dat$w) <- paste("w", 1:ncol(dat$w), sep="")

  any_missing <- (sum(dat$s) != length(dat$s))

  # estimate conditional density (only among observed)
  cond_density_fit <- construct_f_sIx_n(dat=dat, HAL_control = HAL_control)
  f_sIx_n <- cond_density_fit$fnc
  Riemann_grid <- c(0, cond_density_fit$breaks)
  # estimate marginal density (marginalizing the conditional density over whole sample)
  f_s_n <- construct_f_s_n(dat=dat, f_sIx_n=f_sIx_n)
  # estimate density ratio
  g_n <- construct_g_n(f_sIx_n=f_sIx_n, f_s_n=f_s_n)

  # estimate outcome regression (only among observed)
  mu_n <- construct_mu_n(dat=dat, SL_control = SL_control, Riemann_grid = Riemann_grid)

  y_vals <- sort(unique(dat$y))

  if (any_missing){
    alpha_n <- construct_alpha_n(dat = dat, SL_control = SL_control) # sampling weights
    F_n <- construct_Phi_n(dat = dat,
                           alpha_n = alpha_n,
                           f_sIx_n = f_sIx_n,
                           Riemann_grid = Riemann_grid) # debiased cdf estimate
    F_ns <- sapply(y_vals, FUN = F_n)
    F_n_inverse <- function(t){
      if (any(F_ns >= t)){
        val <- y_vals[min(which(F_ns >= t))]
      } else{
        val <- max(y_vals)
      }
    }
  } else{
    # if no missingness, can use the empirical cdf for F_n and the sampling weights are just 1
    F_n <- ecdf(dat$y)
    F_n_inverse <- function(t){
      quantile(dat$y, probs = t, type = 1)
    }
    alpha_n <- function(w) return(1)
  }

  Gamma_n <- construct_Gamma_n(dat=dat, mu_n=mu_n, g_n=g_n,
                               alpha_n = alpha_n, f_sIx_n = f_sIx_n, Riemann_grid = Riemann_grid)
  kappa_n <- construct_kappa_n(dat = dat, mu_n = mu_n, g_n = g_n, alpha_n = alpha_n)

  gcm_x_vals <- sapply(y_vals, F_n)
  inds_to_keep <- !base::duplicated(gcm_x_vals)
  gcm_x_vals <- gcm_x_vals[inds_to_keep]
  gcm_y_vals <- sapply(y_vals[inds_to_keep], Gamma_n)
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

  theta_prime <- construct_deriv(r_Mn = theta_n,
                                 deriv_method = deriv_method,
                                 dir = direction,
                                 y = seq(0, 1, length.out = n_eval_pts))

  # Compute estimates
  ests <- sapply(seq(0,1,length.out = n_eval_pts), theta_n)
  deriv_ests <- sapply(seq(0, 1, length.out = n_eval_pts), theta_prime)
  kappa_ests <- sapply(y_vals, kappa_n)
  # transform kappa to quantile scale
  kappa_ests_rescaled <- sapply(seq(0, 1, length.out = n_eval_pts), function(x){
    ind <- which(y_vals == F_n_inverse(x)) # this one is more general b/c any form of F_n_inverse can be given
    kappa_ests[ind]
  })
  tau_ests <- deriv_ests * kappa_ests_rescaled
  q <- ChernoffDist::qChern(p = 0.975)
  half_intervals <- sapply(1:n_eval_pts, function(x){
    (4*tau_ests[x]/length(dat$y))^{1/3}*q
  })
  cils <- ests - half_intervals
  cius <- ests + half_intervals

  # Plot estimates vs true values
  results <- data.frame(
    x = sapply(seq(0, 1, length.out = n_eval_pts), F_n_inverse), # this form is more general b/c any form of F_n_inverse can be given
    x_quants = seq(0, 1, length.out = n_eval_pts),
    y = ests,
    y_low = cils,
    y_hi = cius
  )

  return(results)

}

#' Estimate missingness probabilities
#' @noRd
construct_alpha_n <- function(dat, SL_control){
  w_distinct <- dplyr::distinct(dat$w)
  w_distinct <- cbind("w_index"=c(1:nrow(w_distinct)), w_distinct)
  newW <- w_distinct
  newW$w_index <- NULL
  model_sl <- SuperLearner::SuperLearner(
    Y = dat$s,
    X = dat$w,
    newX = newW,
    family = "binomial",
    method = SL_control$method,
    SL.library = SL_control$SL.library,
    cvControl = list(V = SL_control$V),
  )
  pred <- as.numeric(model_sl$SL.predict)

  fnc <- function(w) {
    cond <- paste0("round(w1,5)==",round(w[1],5))
    for (i in c(2:length(w))) {
      cond <- paste0(cond," & round(w",i,",5)==",round(w[i],5))
    }
    filtered <- dplyr::filter(w_distinct, eval(parse(text=cond)))
    index <- filtered$w_index
    if (length(index)!=1) {
      stop(paste0("w=c(",paste(w, collapse=","),")"))
    }
    return(pred[index])
  }

  return(fnc)
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
    method = SL_control$method,
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
construct_f_sIx_n <- function(dat, HAL_control){

  # fit hal
  haldensify_fit <- haldensify::haldensify(A = dat$y[dat$s == 1],
                                           W = dat$w[dat$s == 1,],
                                           n_bins = HAL_control$n_bins,
                                           grid_type = HAL_control$grid_type,
                                           cv_folds = HAL_control$V)

  w_distinct <- dplyr::distinct(dat$w)
  w_distinct <- cbind("w_index"=c(1:nrow(w_distinct)), w_distinct)
  # only get predictions at the breakpoints, since estimator is piecewise constant
  y_distinct <- haldensify_fit$breaks
  newW <- expand.grid(w_index=w_distinct$w_index, y=y_distinct)
  newW <- dplyr::inner_join(w_distinct, newW, by="w_index")
  newW$w_index <- NULL

  pred <- predict(haldensify_fit, new_A = newW$y, new_W = newW[,-ncol(newW)])

  newW$index <- c(1:nrow(newW))

  breaks <- haldensify_fit$breaks

  fnc <- function(y,w) {
    # if (s <= min(haldensify_fit$breaks)){
    # left_s <- min(haldensify_fit$breaks)
    # } else{
    left_y <- max(breaks[breaks <= max(y, min(breaks))])
    # }

    cond <- paste0("round(y,5)==",round(left_y,5))
    for (i in c(1:length(w))) {
      cond <- paste0(cond," & round(w",i,",5)==",round(w[i],5))
    }
    index <- (dplyr::filter(newW, eval(parse(text=cond))))$index
    # if (length(index)!=1) {
    # stop(paste0("s=",s,", x=c(",paste(x, collapse=","),")"))
    # }
    return(pred[index])
  }

  return(list(fnc = fnc, breaks = haldensify_fit$breaks))
}

#' Estimate marginal density
#' @noRd
construct_f_s_n <- function(dat, f_sIx_n) {
  uniq_y <- sort(unique(dat$y))
  f_sIx_ns <- sapply(uniq_y, function(y){
    mean(apply(dat$w, 1, function(w) {
      f_sIx_n(y,as.numeric(w))
    }))})

  fnc <- function(y){
    f_sIx_ns[uniq_y == y]
  }
}

construct_g_n <- function(f_sIx_n, f_s_n) {
  function(y,w) { f_sIx_n(y,w) / f_s_n(y) }
}

#' Estimate primitive
#' @noRd
construct_Gamma_n <- function(dat, mu_n, g_n, alpha_n, f_sIx_n, Riemann_grid) {
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

  alpha_ns <- apply(as_df(dat), 1, function(r){
    w <- as.numeric(r)[1:dim_w]
    alpha_n(w = w)
  })

  # piece 1 maps to (\Delta - \mu_n(Y_i, W_i))/g_n(Y_i, W_i)
  piece_1 <- (dat$delta-mu_ns) / g_ns

  # since there aren't so many unique s values in my application
  # we can save a lot of time by only computing piece_2 on those ~100 unique
  # s values rather than all 1800 non-unique ones
  unique_y <- sort(unique(dat$y))

  unique_piece_2 <- sapply(unique_y, function(y) {
    sum(apply(dat$w, 1, function(w) { mu_n(y=y, w=as.numeric(w)) }))
  })

  # match to pre-computed values
  # piece 2 maps to \theta_n(Y_i)
  piece_2 <- sapply(dat$y, function(y) {
    # sum(apply(dat$x, 1, function(x) { mu_n(s=s, x=as.numeric(x)) }))
    unique_piece_2[unique_y == y]
  })
  w_distinct <- dplyr::distinct(dat$w)
  # if there actually is missingness to deal with
  if (sum(dat$s) != length(dat$s)){
    piece_4 <- sapply(Riemann_grid, function(y) {
      mean(apply(dat$w, 1, function(w) { mu_n(y=y, w=as.numeric(w)) }))
    })

    # calculate the Riemann integral w.r.t. the conditional density
    # for each unique value of w
    unique_Riemann_integrals <- apply(w_distinct, MARGIN = 1, function(w){
      density_vals <- sapply(Riemann_grid, function(y){
        f_sIx_n(y = y, w = as.numeric(w))
      })
      Riemann_integrals <- sapply(Riemann_grid, function(y){
        sum(diff(Riemann_grid[Riemann_grid <= y]) * piece_4[Riemann_grid <= y][-1] * density_vals[Riemann_grid <= y][-1])
      })
      Riemann_integrals
    })

    unique_Riemann_integrals <- t(unique_Riemann_integrals)
  } else{
    unique_Riemann_integrals <- matrix(1, nrow = nrow(w_distinct), ncol = length(Riemann_grid))
  }

  fnc <- function(y) {
    piece_3 <- as.integer(dat$y<=y) * dat$s/alpha_ns
    obs_pieces <- (sum(piece_3*piece_1) + mean(piece_3*piece_2))/n_orig
    piece_4_integrals <- apply(dat$w, MARGIN = 1, FUN = function(w){
      unique_Riemann_integrals[which(colSums(t(w_distinct) == w) == ncol(w_distinct)),max(which(Riemann_grid <= y))]
    })
    unobs_pieces <- sum((1 - dat$s/alpha_ns) * piece_4_integrals)/n_orig
    return(obs_pieces + unobs_pieces)
  }

  return(fnc)
}

#' Estimate marginal CDF
#' @noRd
construct_Phi_n <- function(dat, alpha_n, f_sIx_n, Riemann_grid){

  w_distinct <- dplyr::distinct(dat$w)

  # calculate the Riemann integral w.r.t. the conditional density
  # for each unique value of w
  unique_Riemann_integrals <- apply(w_distinct, MARGIN = 1, function(w){
    density_vals <- sapply(Riemann_grid, function(y){
      f_sIx_n(y = y, w = as.numeric(w))
    })
    Riemann_integrals <- sapply(Riemann_grid, function(y){
      sum(diff(Riemann_grid[Riemann_grid <= y]) * density_vals[Riemann_grid <= y][-1])
    })
    Riemann_integrals
  })

  unique_Riemann_integrals <- t(unique_Riemann_integrals)

  alpha_ns <- apply(dat$w, MARGIN = 1, FUN = alpha_n)

  piece_1 <- dat$s/alpha_ns
  piece_2 <- (1 - dat$s/alpha_ns)

  fnc <- function(y){

    piece_1_ind <- piece_1 * (dat$y <= y)
    piece_2_integral <- apply(dat$w, MARGIN = 1, FUN = function(w){
      unique_Riemann_integrals[which(colSums(t(w_distinct) == w) == ncol(w_distinct)),max(which(Riemann_grid <= y))]
    })
    piece_2_ind <- piece_2 * piece_2_integral
    return(mean(piece_1_ind + piece_2_ind))
  }

  return(fnc)
}

#' Estimate part of scale factor
#' @noRd
construct_kappa_n <- function(dat, mu_n, g_n, alpha_n){

  fnc <- function(y){
    mean(apply(dat$w, 1, function(w) {
      numer1 <- mu_n(y = y, w = as.numeric(w))
      numer <- numer1 * (1 - numer1)
      denom1 <- g_n(y = y, w = as.numeric(w))
      denom2 <- alpha_n(w = as.numeric(w))
      numer / (denom1 * denom2)
    }))
  }

  return(fnc)
}

#' Estimate derivative
#' @noRd
construct_deriv <- function(deriv_method="m-spline", r_Mn, dir, y) {
  # r_Mn should be the theta_n function
  # dir should be ?
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
      if (dir=="decr") {
        return(min((fnc_pre(x2)-fnc_pre(x1))/width,0))
      } else {
        return(max((fnc_pre(x2)-fnc_pre(x1))/width,0))
      }
    }
  }
  return(fnc)
}

#' Convert list data to df
#' @noRd
as_df <- function(dat) { cbind(dat$w, y=dat$y) }
