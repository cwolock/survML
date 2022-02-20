#' Stack a dataset, using the hazard (to be merged with stack.R)
#'
#' @return A stacked dataset
#' @noRd
stack_haz <- function(time, event, X, time_grid, entry){
  trunc_time_grid <- time_grid[-length(time_grid)] # do I need to truncate if treating time as continuous? look at this later
  dat <- data.frame(X, event = event, time = time)
  # we will treat time as continuous
  ncol_stacked <- ncol(X) + 2 # covariates, time, binary outcome
  stacked <- matrix(NA, ncol = ncol_stacked, nrow = 1)
  for (i in 1:(length(trunc_time_grid))){ # can change this to not do anything in last time bin
    if (is.null(entry)){ # no entry (truncation variable) given
      risk_set <- dat[dat$time > time_grid[i],]
    } else{ # entry given
      risk_set <- dat[dat$time > time_grid[i] & entry <= time_grid[i+1],]
    }
    risk_set_covariates <- risk_set[,1:ncol(X)]
    event_indicators <- matrix(ifelse(risk_set$time <= time_grid[i + 1 ] & risk_set$event == 1, 1, 0))
    t <- rep(time_grid[i + 1], nrow(event_indicators))
    newdata <- as.matrix(cbind(t, risk_set_covariates, event_indicators))
    stacked <- rbind(stacked, newdata)
  }
  stacked <- stacked[-1,]
  colnames(stacked)[ncol(stacked)] <- "event_indicators"
  return(stacked)
}

