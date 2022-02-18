#' Stack a dataset
#'
#' @return A stacked dataset
#' @noRd
stack <- function(time, X, time_grid){
  trunc_time_grid <- time_grid[-length(time_grid)] # do I need to truncate if treating time as continuous? look at this later
  # we will treat time as continuous
  ncol_stacked <- ncol(X) + 2 # covariates, time, binary outcome
  stacked <- matrix(NA, ncol = ncol_stacked, nrow = 1)
  # for (i in 1:(length(trunc_time_grid))){ # can change this to not do anything in last time bin
  #   event_indicators <- matrix(ifelse(time <= time_grid[i + 1], 1, 0))
  #   t <- time_grid[i + 1]
  #   newdata <- as.matrix(cbind(t, X, event_indicators))
  #   stacked <- rbind(stacked, newdata)
  # }
  for (i in 1:(length(time_grid))){ # can change this to not do anything in last time bin
    event_indicators <- matrix(ifelse(time <= time_grid[i], 1, 0))
    t <- time_grid[i]
    newdata <- as.matrix(cbind(t, X, event_indicators))
    stacked <- rbind(stacked, newdata)
  }
  stacked <- stacked[-1,]
  colnames(stacked)[ncol(stacked)] <- "event_indicators"
  return(stacked)
}

