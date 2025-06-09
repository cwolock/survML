#' Gradient boosting for C-index
#'
#' Using doubly-robust gradient boosting to generate estimates of the prediction function that maximizes the C-index
#'
#' @param time \code{n x 1} numeric vector of observed
#' follow-up times. If there is censoring, these are the minimum of the
#' event and censoring times.
#' @param event \code{n x 1} numeric vector of status indicators of
#' whether an event was observed.
#' @param X \code{n x p} data.frame of observed covariate values
#' @param newX \code{m x p} data.frame of new observed covariate
#' values at which to obtain \code{m} predictions for the estimated algorithm.
#' Must have the same names and structure as \code{X}.
#' @param approx_times Numeric vector of length J2 giving times at which to
#' approximate C-index integral. Note that the last time in \code{approx_times} is taken to be the restriction time (i.e., the maximum follow-up)
#' for comparison of pairs of individuals. Essentially, this time should be chosen such that the conditional survival function is identified at
#' this time for all covariate values \code{X} present in the data. Choosing the restriction time such that roughly 10% of individuals remain at-risk
#' at that time has been shown to work reasonably well in simulations.
#' @param S_hat \code{n x J2} matrix of conditional event time survival function estimates
#' @param G_hat \code{n x J2} matrix of conditional censoring time survival function estimates
#' @param V Number of cross-validation folds for selection of tuning parameters
#' @param tuning Logical, whether or not to use cross-validation to select tuning parameters
#' @param produce_fit Logical, whether to produce a fitted prediction function using the selected optimal parameters.
#' @param subsample_n Number of samples to use for boosting procedure. Using a subsample of the full sample can greatly reduce runtime
#' @param boosting_params Named list of parameter values for the boosting procedure. Elements of this list include \code{mstop} (number of
#' boosting iterations), \code{nu} (learning rate), \code{sigma} (smoothness parameter for sigmoid approximation, with smaller meaning
#' less smoothing), and \code{learner} (base learner, can take values \code{"glm"}, \code{"gam"}, or \code{"tree"})
#'
#' @return Vector of predictions corresponding to \code{newX}
#'
#' @examples
#' # This is a small simulation example
#' set.seed(123)
#' n <- 250
#' X <- data.frame(X1 = rnorm(n), X2 = rbinom(n, size = 1, prob = 0.5))
#'
#' T <- rexp(n, rate = exp(-2 + X[,1] - X[,2] + .5 *  X[,1] * X[,2]))
#' C <- rexp(n, exp(-2 -.5 * X[,1] - .25 * X[,2] + .5 * X[,1] * X[,2]))
#' C[C > 15] <- 15
#'
#' time <- pmin(T, C)
#' event <- as.numeric(T <= C)
#'
#' # Note that this a very small Super Learner library, for computational purposes.
#' SL.library <- c("SL.mean", "SL.glm")
#'
#' # Note that we do not use times beyond the 90th percentile of observed follow-up times
#' approx_times <- c(0, unique(quantile(time, probs = seq(0, 0.9, by = 0.01))))
#'
#' # estimate conditional survival functions at approx_times
#' fit <- stackG(time = time,
#'               event = event,
#'               X = X,
#'               newX = X,
#'               newtimes = approx_times,
#'               direction = "prospective",
#'               bin_size = 0.1,
#'               time_basis = "continuous",
#'               surv_form = "PI",
#'               learner = "SuperLearner",
#'               time_grid_approx = approx_times,
#'               SL_control = list(SL.library = SL.library,
#'                                 V = 3))
#'
#' # use boosting to estimate optimal (according to C-index) prediction function
#' boosted_preds <- boost_c_index(time = time,
#'                                event = event,
#'                                X = X,
#'                                newX = X,
#'                                approx_times = approx_times,
#'                                S_hat = fit$S_T_preds,
#'                                G_hat = fit$S_C_preds,
#'                                V = 3,
#'                                tuning = TRUE,
#'                                produce_fit = TRUE,
#'                                subsample_n = 200,
#'                                boosting_params = list(mstop = c(100, 200),
#'                                                       nu = 0.1,
#'                                                       sigma = 0.1,
#'                                                       learner = "glm"))
#' boosted_preds
#'
#' @export
boost_c_index <- function(time,
                          event,
                          X,
                          newX,
                          S_hat,
                          G_hat,
                          V,
                          approx_times,
                          tuning,
                          produce_fit = TRUE,
                          subsample_n,
                          boosting_params){

  if (subsample_n < length(time)){
    subsample_inds <- sample(1:length(time), size = subsample_n, replace = FALSE)
    time <- time[subsample_inds]
    event <- event[subsample_inds]
    X <- X[subsample_inds,,drop=FALSE]
    S_hat <- S_hat[subsample_inds,,drop=FALSE]
    G_hat <- G_hat[subsample_inds,,drop=FALSE]
  }

  KM_IFs <- calc_KM_IF(time = time,
                       event = event,
                       S_hat = S_hat,
                       G_hat = G_hat,
                       approx_times = approx_times)

  pi <- S_hat + KM_IFs

  reverse_sorted_mstop <- boosting_params$mstop[order(boosting_params$mstop, decreasing = TRUE)]
  tau <- max(approx_times)

  k <- length(approx_times)
  n <- length(time)
  folds <- sample(rep(seq_len(V), length = n))

  param_grid <- expand.grid(mstop = max(boosting_params$mstop),
                            nu = boosting_params$nu,
                            sigma = boosting_params$sigma,
                            learner = boosting_params$learner)
  param_grid_full <- expand.grid(mstop = reverse_sorted_mstop,
                                 nu = boosting_params$nu,
                                 sigma = boosting_params$sigma,
                                 learner = boosting_params$learner)
  param_grid_full$CV_risk <- NA
  mod_list <- vector("list", nrow(param_grid))

  K <- length(boosting_params$mstop)
  for (i in 1:nrow(param_grid)){

    mstop_curr = param_grid[i,1]
    nu_curr <- param_grid[i,2]
    sigma_curr <- param_grid[i,3]
    learner_curr <- param_grid[i,4]

    Sweights_train <- function(i1, i2){
      # prob that i2 exceeds i1
      exceed_prob <- -sum(pi_train[i2,2:k] * diff(pi_train[i1,1:k]))
      return(exceed_prob)
    }

    # the i,j entry of this matrix is P(i fails before j)
    my_Cindex <- function (sigma = 0.1) {
      approxGrad <- function(x) { ## sigmoid function for gradient
        exp(x/sigma) / (sigma * (1 + exp(x/sigma))^2)
      }
      approxLoss <- function(x) { ## sigmoid function for loss
        1 / (1 + exp(x / sigma))
      }
      ngradient = function(y, f, w = 1) { ## negative gradient
        if (!all(w %in% c(0,1)))
          stop(sQuote("weights"), " must be either 0 or 1 for family ",
               sQuote("UnoC"))
        survtime <- y[,1]
        event <- y[,2]
        if (length(w) == 1) w <- rep(1, length(event))
        if (length(f) == 1) {
          f <- rep(f, length(survtime))
        }
        n <- length(survtime)
        etaj <- matrix(f, nrow = n, ncol = n, byrow = TRUE)
        etak <- matrix(f, nrow = n, ncol = n)
        etaMat <- etak - etaj
        rm(etaj); rm(etak);
        weights_out <- wweights
        M1 <- approxGrad(etaMat) * weights_out
        ng <- colSums(M1) - rowSums(M1)
        return(ng)
      }
      risk = function(y, f, w = 1) { ## empirical risk
        survtime <- y[,1]
        event <- y[,2]
        if (length(f) == 1) {
          f <- rep(f, length(y))
        }
        n <- length(survtime)
        etaj <- matrix(f, nrow = n, ncol = n, byrow = TRUE)
        etak <- matrix(f, nrow = n, ncol = n)
        etaMat <- (etak - etaj)
        rm(etaj); rm(etak);
        weights_out <- wweights
        M1 <- approxLoss(etaMat) * weights_out
        return(- sum(M1))
      }
      mboost::Family( ## build the family object
        ngradient = ngradient,
        risk = risk,
        weights = "zeroone",
        offset = function(y, w = 1) {0},
        check_y = function(y) {
          if (!inherits(y, "Surv"))
            stop("response is not an object of class ", sQuote("Surv"),
                 " but ", sQuote("family = UnoC()"))
          y},
        rclass = function(f){},
        name = paste("Efron-type c-index boosting")
      )
    }

    if (tuning){ # cross-validation tuning
      CV_risks <- matrix(NA, nrow = V, ncol = K)
      for (j in 1:V){
        time_train <- time[folds != j]
        event_train <- event[folds != j]
        X_train <- X[folds != j,,drop=FALSE]
        S_hat_train <- S_hat[folds != j,]
        pi_train <- pi[folds != j,]
        X_test <- X[folds == j,,drop=FALSE]
        S_hat_test <- S_hat[folds == j,]
        pi_test <- pi[folds == j]
        time_test <- time[folds == j]
        event_test <- event[folds == j]

        w <- rep(1, length(time_train))
        index_grid <- gtools::combinations(n = length(time_train),
                                           r = 2,
                                           v = 1:length(time_train),
                                           repeats.allowed = TRUE)
        weight_vec <- mapply(FUN = Sweights_train, index_grid[,1], index_grid[,2])
        wweights <- matrix(NA, length(time_train), length(time_train))
        wweights[lower.tri(wweights, diag=TRUE)] <- weight_vec
        wweights <- t(wweights)
        wweights[lower.tri(wweights)] <- 1-t(wweights)[lower.tri(wweights)]
        Wmat <- w %o% w
        wweights <- wweights * Wmat
        diag(wweights) <- 0.5
        wweights <- wweights / sum(wweights)

        dtrain <- data.frame(time = time_train,
                             event = event_train,
                             X_train)
        dtest <- data.frame(X_test)

        feature_names <- names(X_train)

        feature_names <- apply(as.matrix(feature_names),
                               MARGIN = 1,
                               FUN = function(x) paste0("", x, ""))
        feature_form <- paste(feature_names, collapse = "+")
        formula_text <- paste0("Surv(time, event) ~ ", feature_form)
        if (learner_curr == "tree"){
          mod <- mboost::blackboost(survival::Surv(time, event) ~ .,
                                    family = my_Cindex(sigma = sigma_curr),
                                    control = mboost::boost_control(mstop = mstop_curr,
                                                                    trace = FALSE,
                                                                    nu = nu_curr),
                                    data = dtrain)
        } else if (learner_curr == "glm"){
          mod <- mboost::glmboost(survival::Surv(time, event) ~ .^2,
                                  family = my_Cindex(sigma = sigma_curr),
                                  control = mboost::boost_control(mstop = mstop_curr,
                                                                  trace = FALSE,
                                                                  nu = nu_curr),
                                  data = dtrain)
        } else if (learner_curr == "gam"){
          mod <- mboost::gamboost(survival::Surv(time, event) ~ .,
                                  family = my_Cindex(sigma = sigma_curr),
                                  control = mboost::boost_control(mstop = mstop_curr,
                                                                  trace = FALSE,
                                                                  nu = nu_curr),
                                  data = dtrain)
        }

        for (l in 1:K){
          sub_mstop = reverse_sorted_mstop[l]
          mboost::mstop(mod) <- sub_mstop
          preds <- -stats::predict(mod, newdata = dtest)[,1]

          risk_j <- -estimate_cindex(time = time_test,
                                     event = event_test,
                                     approx_times = approx_times,
                                     preds = preds,
                                     S_hat = S_hat_test,
                                     G_hat = matrix(0.5, nrow = nrow(S_hat_test),
                                                    ncol = ncol(S_hat_test)),
                                     tau = tau)$plug_in
          CV_risks[j,l] <- risk_j
        }
      }
      param_grid_full$CV_risk[(((i-1)*K)+1):(i*K)] <- colMeans(CV_risks)

    } else{ # no tuning

      time_train <- time
      event_train <- event
      X_train <- X
      S_hat_train <- S_hat
      pi_train <- pi
      dtrain <- data.frame(time = time_train,
                           event = event_train,
                           X_train)

      #Survival function version
      w <- rep(1, length(time_train))
      index_grid <- gtools::combinations(n = length(time_train),
                                         r = 2,
                                         v = 1:length(time_train),
                                         repeats.allowed = TRUE)
      weight_vec <- mapply(FUN = Sweights_train, index_grid[,1], index_grid[,2])
      wweights <- matrix(NA, length(time_train), length(time_train))
      wweights[lower.tri(wweights, diag=TRUE)] <- weight_vec
      wweights <- t(wweights)
      wweights[lower.tri(wweights)] <- 1-t(wweights)[lower.tri(wweights)]
      Wmat <- w %o% w
      wweights <- wweights * Wmat
      diag(wweights) <- 0.5
      wweights <- wweights / sum(wweights)

      feature_names <- names(X_train)
      feature_names <- apply(as.matrix(feature_names),
                             MARGIN = 1,
                             FUN = function(x) paste0("", x, ""))
      feature_form <- paste(feature_names, collapse = "+")
      formula_text <- paste0("Surv(time, event) ~ ", feature_form)

      if (learner_curr == "tree"){
        mod <- mboost::blackboost(survival::Surv(time, event) ~ .,
                                  family = my_Cindex(sigma = sigma_curr),
                                  control = mboost::boost_control(mstop = mstop_curr,
                                                                  trace = FALSE,
                                                                  nu = nu_curr),
                                  data = dtrain)
      } else if (learner_curr == "glm"){
        mod <- mboost::glmboost(survival::Surv(time, event) ~ .^2,
                                family = my_Cindex(sigma = sigma_curr),
                                control = mboost::boost_control(mstop = mstop_curr,
                                                                trace = FALSE,
                                                                nu = nu_curr),
                                data = dtrain)
      } else if (learner_curr == "gam"){
        mod <- mboost::gamboost(survival::Surv(time, event) ~ .,
                                family = my_Cindex(sigma = sigma_curr),
                                control = mboost::boost_control(mstop = mstop_curr,
                                                                trace = FALSE,
                                                                nu = nu_curr),
                                data = dtrain)
      }

      mod_list[[i]] <- mod
      param_grid$CV_risk[i] <- 999
    }
  }

  if (!tuning){
    param_grid_full <- param_grid
  }

  param_grid_full <- param_grid_full[order(param_grid_full$mstop), ]
  opt_index <- which.min(round(param_grid_full$CV_risk, digits = 3))

  if (tuning & produce_fit){
    mstop_curr = param_grid_full[opt_index,1]
    nu_curr <- param_grid_full[opt_index,2]
    sigma_curr <- param_grid_full[opt_index,3]
    learner_curr <- param_grid_full[opt_index,4]
    time_train <- time
    event_train <- event
    X_train <- X
    S_hat_train <- S_hat
    pi_train <- pi
    dtrain <- data.frame(time = time_train,
                         event = event_train,
                         X_train)

    #Survival function version
    w <- rep(1, length(time_train))
    index_grid <- gtools::combinations(n = length(time_train),
                                       r = 2,
                                       v = 1:length(time_train),
                                       repeats.allowed = TRUE)
    weight_vec <- mapply(FUN = Sweights_train, index_grid[,1], index_grid[,2])
    wweights <- matrix(NA, length(time_train), length(time_train))
    wweights[lower.tri(wweights, diag=TRUE)] <- weight_vec
    wweights <- t(wweights)
    wweights[lower.tri(wweights)] <- 1-t(wweights)[lower.tri(wweights)]
    Wmat <- w %o% w
    wweights <- wweights * Wmat
    diag(wweights) <- 0.5 # diagonal should be 1/2, I think?
    wweights <- wweights / sum(wweights)

    feature_names <- names(X_train)
    feature_names <- apply(as.matrix(feature_names),
                           MARGIN = 1,
                           FUN = function(x) paste0("btree(", x, ")"))
    feature_form <- paste(feature_names, collapse = "+")
    formula_text <- paste0("Surv(time, event) ~ ", feature_form)

    if (learner_curr == "tree"){
      mod <- mboost::blackboost(survival::Surv(time, event) ~ .,
                                family = my_Cindex(sigma = sigma_curr),
                                control = mboost::boost_control(mstop = mstop_curr,
                                                                trace = FALSE,
                                                                nu = nu_curr),
                                data = dtrain)
    } else if (learner_curr == "glm"){
      mod <- mboost::glmboost(survival::Surv(time, event) ~ .^2,
                              family = my_Cindex(sigma = sigma_curr),
                              control = mboost::boost_control(mstop = mstop_curr,
                                                              trace = FALSE,
                                                              nu = nu_curr),
                              data = dtrain)
    } else if (learner_curr == "gam"){
      mod <- mboost::gamboost(survival::Surv(time, event) ~ .,
                              family = my_Cindex(sigma = sigma_curr),
                              control = mboost::boost_control(mstop = mstop_curr,
                                                              trace = FALSE,
                                                              nu = nu_curr),
                              data = dtrain)
    }
    opt_model <- mod
  } else if (!tuning & produce_fit){
    opt_model <- mod_list[[opt_index]]
  } else{
    opt_model <- NULL
  }

  dtest <- data.frame(newX)
  f_hat <- -stats::predict(opt_model, newdata = dtest)[,1]

  return(f_hat)
}
