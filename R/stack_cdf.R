#' Stack a dataset
#'
#' @return A stacked dataset
#' @noRd
stack_cdf <- function(time, X, time_grid, time_basis){

  if (time_basis == "continuous"){
    # we will treat time as continuous
    ncol_stacked <- ncol(X) + 2 # covariates, time, binary outcome
    stacked <- matrix(NA, ncol = ncol_stacked, nrow = 1)
    for (i in 1:(length(time_grid))){ # can change this to not do anything in last time bin
      event_indicators <- matrix(ifelse(time <= time_grid[i], 1, 0))
      t <- time_grid[i]
      newdata <- as.matrix(cbind(t, X, event_indicators))
      stacked <- rbind(stacked, newdata)
    }
    stacked <- stacked[-1,]
    colnames(stacked)[ncol(stacked)] <- "event_indicators"
  } else if (time_basis == "dummy"){
    dat <- data.frame(X, time)
    ncol_stacked <- ncol(X) + length(time_grid) + 1 # covariates, risk set dummies, binary outcome
    stacked <- matrix(NA, ncol = ncol_stacked, nrow = 1)
    for (i in 1:(length(time_grid))){
      risk_set <- dat
      risk_set_covariates <- risk_set[,1:ncol(X)]
      event_indicators <- matrix(ifelse(risk_set$time <= time_grid[i], 1, 0))
      dummies <- matrix(0, ncol = length(time_grid), nrow = nrow(risk_set))
      dummies[,i] <- 1
      newdata <- as.matrix(cbind(dummies, risk_set_covariates, event_indicators))
      stacked <- rbind(stacked, newdata)
    }

    stacked <- stacked[-1,]
    risk_set_names <- paste0("risk_set_", seq(1, (length(time_grid))))
    colnames(stacked)[1:(length(time_grid))] <- risk_set_names
  }

  stacked <- data.frame(stacked)
  return(stacked)
}

