#' Estimate a conditional survival function via stacking
#'
#' @param time Observed time
#' @param event Indicator of event (vs censoring)
#' @param X Covariate matrix
#' @param newX Values of covariates at which to make a prediction
#' @param newtimes Times at which to make the survival function prediction
#' @param bin_size Quantiles on which to grid times. If NULL, defaults to every observed time
#' @param V Number of cross validation folds
#' @param time_basis How to treat time (continuous or dummy)
#' @param entry Study entry variable, if applicable
#' @param direction Prospective or retrospective study
#' @param SL.library SuperLearner library
#'
#' @return An object of class \code{survMLs}
#'
#' @export
#'
#' @examples
survMLs <- function(time,
                    event,
                    X,
                    newX,
                    newtimes,
                    bin_size = NULL,
                    V = 10,
                    time_basis = "continuous",
                    entry = NULL,
                    direction = "prospective",
                    SL.library = NULL,
                    tau = NULL){

  if (direction == "retrospective"){
    if (is.null(tau)){
      tau <- max(entry)
    }
    time <- tau - time
    newtimes <- tau - newtimes
    entry <- tau - entry
    event <- rep(1, length(time))
  }

  X <- as.matrix(X)
  time <- as.matrix(time)
  event <- as.matrix(event)
  dat <- data.frame(X, time, event)

  # if user gives bin size, set time grid based on quantiles. otherwise, every observed time
  if (!is.null(bin_size)){
    time_grid <- quantile(dat$time[dat$event == 1], probs = seq(0, 1, by = bin_size))
    time_grid[1] <- 0
  } else{
    time_grid <- sort(unique(dat$time[dat$event == 1]))
    time_grid <- c(0, time_grid)
  }

  # create stacked dataset
  stacked <- survML:::stack_haz(time = time,
                                event = event,
                                X = X,
                                time_grid = time_grid,
                                entry = entry,
                                time_basis = time_basis)

  .Y <- stacked[,ncol(stacked)]
  .X <- data.frame(stacked[,-ncol(stacked)])
  # fit Super Learner
  fit <- SuperLearner::SuperLearner(Y = .Y,
                                    X = .X,
                                    SL.library = SL.library,
                                    family = binomial(),
                                    method = 'method.NNLS',
                                    verbose = FALSE,
                                    cvControl = list(V = V))

  # create function to get discrete hazard predictions
  if (time_basis == "continuous"){
    get_hazard_preds <- function(t){
      new_stacked <- data.frame(t = t, newX)
      preds <- predict(fit, newdata=new_stacked)$pred
      return(preds)
    }
  } else if (time_basis == "dummy"){
    get_hazard_preds <- function(t){
      dummies <- matrix(0, ncol = length(time_grid), nrow = nrow(newX))
      index <- max(which(time_grid <= t))
      dummies[,index] <- 1
      new_stacked <- cbind(dummies, newX)
      risk_set_names <- paste0("risk_set_", seq(1, (length(time_grid))))
      colnames(new_stacked)[1:length(time_grid)] <- risk_set_names
      new_stacked <- data.frame(new_stacked)
      preds <- predict(fit, newdata=new_stacked)$pred
      return(preds)
    }
  }

  # don't estimate hazard at t =0
  #hazard_preds <- apply(X = matrix(time_grid), FUN = get_hazard_preds, MARGIN = 1)
  hazard_preds <- apply(X = matrix(time_grid[-1]), FUN = get_hazard_preds, MARGIN = 1)

  get_surv_preds <- function(t){
    if (sum(time_grid[-1] <= t) != 0){ # if you don't fall before the first time in the grid
      final_index <- max(which(time_grid[-1] <= t))
    # if (sum(time_grid <= t) != 0){ # if you don't fall before the first time in the grid
    #   final_index <- max(which(time_grid <= t))
      haz <- as.matrix(hazard_preds[,1:final_index])
      anti_haz <- 1 - haz
      surv <- apply(anti_haz, MARGIN = 1, prod)
    } else{
      surv <- rep(1, nrow(hazard_preds))
    }
    return(surv)
  }

  surv_preds <- apply(X = matrix(newtimes), FUN = get_surv_preds, MARGIN = 1)


  if (direction == "retrospective"){
    surv_preds <- 1 - surv_preds
  }

  res <- list(S_T_preds = surv_preds,
              fit = fit)
  class(res) <- "survMLs"
  return(res)
}
