#' Estimate a conditional survival function via stacking
#'
#' @param time Observed time
#' @param event Indicator of event (vs censoring)
#' @param X Covariate matrix
#' @param newX Values of covariates at which to make a prediction
#' @param newtimes Times at which to make the survival function prediction
#' @param bin_size Quantiles on which to grid times. If NULL, defaults to every observed time
#' @param SL.library SuperLearner library
#'
#' @return An object of class \code{stack}
#'
#' @export
#'
#' @examples
stackSurv <- function(time,
                      event,
                      X,
                      newX,
                      newtimes,
                      bin_size = NULL,
                      test_time = NULL,
                      test_event = NULL,
                      test_X = NULL,
                      SL.library){

  X <- as.matrix(X)
  time <- as.matrix(time)
  event <- as.matrix(event)
  dat <- data.frame(X, time, event)

  # if user gives bin size, set time grid based on quantiles. otherwise, every observed time
  if (!is.null(bin_size)){
    time_grid <- quantile(dat$time, probs = seq(0, 1, by = bin_size))
    time_grid[1] <- 0 # manually set first point to 0, instead of first observed time
  } else{
    time_grid <- sort(unique(time))
    time_grid <- c(0, time_grid)
  }

  trunc_time_grid <- time_grid[-length(time_grid)]

  ncol_stacked <- ncol(X) + length(trunc_time_grid) + 1 # covariates, risk set dummies, binary outcome
  stacked <- matrix(NA, ncol = ncol_stacked, nrow = 1)
  for (i in 1:(length(trunc_time_grid))){
    risk_set <- dat[dat$time > time_grid[i],]
    risk_set_covariates <- risk_set[,1:ncol(X)]
    event_indicators <- matrix(ifelse(risk_set$time <= time_grid[i + 1 ] & risk_set$event == 1, 1, 0))
    dummies <- matrix(0, ncol = length(trunc_time_grid), nrow = nrow(risk_set))
    dummies[,i] <- 1
    newdata <- as.matrix(cbind(dummies, risk_set_covariates, event_indicators))
    stacked <- rbind(stacked, newdata)
  }

  stacked <- stacked[-1,]
  risk_set_names <- paste0("risk_set_", seq(1, (length(trunc_time_grid))))
  colnames(stacked)[1:(length(trunc_time_grid))] <- risk_set_names
  stacked <- data.frame(stacked)
  Y <- stacked$event_indicators
  X <- stacked[,-ncol(stacked)]

  fit <- SuperLearner::SuperLearner(Y = Y,
                                    X = X,
                                    family = binomial(),
                                    SL.library = SL.library,
                                    method = "method.NNLS",
                                    verbose = FALSE)


  get_stacked_pred <- function(t){
    n_bins <- sum(t > trunc_time_grid)
    new_stacked <- matrix(rep(t(newX), n_bins), ncol = ncol(newX) , byrow = TRUE)
    dummy_stacked <- matrix(NA, ncol = length(trunc_time_grid), nrow = 1)
    for (i in 1:n_bins){
      dummies <- matrix(0, ncol = length(trunc_time_grid), nrow = nrow(newX))
      if (t > trunc_time_grid[i]){
        dummies[,i] <- 1
      }
      dummy_stacked <- rbind(dummy_stacked, dummies)
    }
    dummy_stacked <- dummy_stacked[-1,]
    new_stacked <- cbind(dummy_stacked, new_stacked)
    risk_set_names <- paste0("risk_set_", seq(1, (length(trunc_time_grid))))
    colnames(new_stacked) <- c(risk_set_names, colnames(newX))
    new_stacked <- data.frame(new_stacked)
    haz_preds <- predict(fit, newdata=new_stacked)$pred
    haz_preds <- matrix(haz_preds, nrow = nrow(newX))
    surv_preds <- 1 - haz_preds
    total_surv_preds <- apply(surv_preds, prod, MARGIN = 1)
    return(total_surv_preds)
  }

  predictions <- apply(X = matrix(newtimes), FUN = get_stacked_pred, MARGIN = 1)

  res <- list(S_T_preds = predictions,
              fit = fit)
  class(res) <- "stackSurv"
  return(res)
}
