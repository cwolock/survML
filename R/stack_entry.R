#' Stack a dataset, using the entry time conditional on event time (to be merged with stack.R)
#'
#' @return A stacked dataset
#' @noRd
stack_entry <- function(time, entry, X, time_grid){
  trunc_time_grid <- time_grid[-length(time_grid)] # do I need to truncate if treating time as continuous? look at this later
  dat <- data.frame(X, time = time, entry = entry)
  # we will treat time as continuous
  ncol_stacked <- ncol(X) + 2 # covariates, time, binary outcome
  stacked <- matrix(NA, ncol = ncol_stacked, nrow = 1)
  for (i in 1:(length(trunc_time_grid))){ # can change this to not do anything in last time bin
    risk_set <- dat[dat$time >= time_grid[i+1],]
    risk_set_covariates <- risk_set[,1:ncol(X)]
    event_indicators <- matrix(ifelse(risk_set$entry <= time_grid[i + 1 ], 1, 0))
    t <- rep(time_grid[i + 1], nrow(risk_set_covariates))
    newdata <- as.matrix(cbind(t, risk_set_covariates, event_indicators))
    stacked <- rbind(stacked, newdata)
  }
  stacked <- stacked[-1,]
  print(stacked)
  colnames(stacked)[ncol(stacked)] <- "event_indicators"
  return(stacked)
}

