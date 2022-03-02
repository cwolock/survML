#' Stack a dataset, using the entry time conditional on event time (to be merged with stack.R)
#'
#' @return A stacked dataset
#' @noRd
stack_entry <- function(time, entry, X, time_grid, direction){
  trunc_time_grid <- time_grid[-length(time_grid)] # do I need to truncate if treating time as continuous? look at this later
  dat <- data.frame(X, time = time, entry = entry)
  # we will treat time as continuous
  ncol_stacked <- ncol(X) + 2 # covariates, time, binary outcome
  stacked <- matrix(NA, ncol = ncol_stacked, nrow = 1)
  for (i in 1:(length(trunc_time_grid))){
    if (direction == "forward"){
      #for (i in 1:(length(time_grid))){# can change this to not do anything in last time bin
      # if (i > 1){
      #   risk_set <- dat[dat$time > time_grid[i-1],]
      # } else{
      #   risk_set <- dat
      # }
      risk_set <- dat[dat$time > time_grid[i],]# maybe this should be >= i+1? Need to think more carefully about this. obv in the limit, doesn't matter
      risk_set_covariates <- risk_set[,1:ncol(X)]
      #event_indicators <- matrix(ifelse(risk_set$entry <= time_grid[i], 1, 0))
      event_indicators <- matrix(ifelse(risk_set$entry <= time_grid[i + 1 ], 1, 0))
      #t <- rep(time_grid[i], nrow(risk_set_covariates))
      t <- rep(time_grid[i + 1], nrow(risk_set_covariates))
      newdata <- as.matrix(cbind(t, risk_set_covariates, event_indicators))
      stacked <- rbind(stacked, newdata)
    } else if (direction == "reverse"){
      risk_set <- dat[dat$time < time_grid[i+1],]
      risk_set_covariates <- risk_set[,1:ncol(X)]
      #event_indicators <- matrix(ifelse(risk_set$entry <= time_grid[i], 1, 0))
      event_indicators <- matrix(ifelse(risk_set$entry >= time_grid[i ], 1, 0))
      #t <- rep(time_grid[i], nrow(risk_set_covariates))
      t <- rep(time_grid[i], nrow(risk_set_covariates))
      newdata <- as.matrix(cbind(t, risk_set_covariates, event_indicators))
      stacked <- rbind(stacked, newdata)
    }

  }
  stacked <- stacked[-1,]
  colnames(stacked)[ncol(stacked)] <- "event_indicators"
  return(stacked)
}

