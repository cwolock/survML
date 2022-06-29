#' Stack a dataset, using the entry time conditional on event time (to be merged with stack.R)
#'
#' @return A stacked dataset
#' @noRd
stack_entry <- function(time, entry, X, time_grid, time_basis, ids = FALSE){
  trunc_time_grid <- time_grid[-length(time_grid)] # do I need to truncate if treating time as continuous? look at this later
  dat <- data.frame(X, time = time, entry = entry)
  if (ids){
    id_vec <- NA
  }
  if (time_basis == "continuous"){
    ncol_stacked <- ncol(X) + 2 # covariates, time, binary outcome
    stacked <- matrix(NA, ncol = ncol_stacked, nrow = 1)
    for (i in 1:(length(trunc_time_grid))){
      #for (i in 1:(length(time_grid))){# can change this to not do anything in last time bin
      # if (i > 1){
      #   risk_set <- dat[dat$time > time_grid[i-1],]
      # } else{
      #   risk_set <- dat
      # }
      risk_set <- dat[dat$time > time_grid[i],]# maybe this should be >= i+1? Need to think more carefully about this. obv in the limit, doesn't matter
      if (ids){
        id_i <- which(dat$time > time_grid[i])
        id_vec <- c(id_vec, id_i)
      }
      risk_set_covariates <- risk_set[,1:ncol(X)]
      #event_indicators <- matrix(ifelse(risk_set$entry <= time_grid[i], 1, 0))
      event_indicators <- matrix(ifelse(risk_set$entry <= time_grid[i + 1 ], 1, 0))
      #t <- rep(time_grid[i], nrow(risk_set_covariates))
      t <- rep(time_grid[i + 1], nrow(risk_set_covariates))
      newdata <- as.matrix(cbind(t, risk_set_covariates, event_indicators))
      stacked <- rbind(stacked, newdata)
    }
    stacked <- stacked[-1,]
    colnames(stacked)[ncol(stacked)] <- "event_indicators"
  } else if (time_basis ==  "dummy"){
    ncol_stacked <- ncol(X) + length(trunc_time_grid) + 1
    stacked <- matrix(NA, ncol = ncol_stacked, nrow = 1)
    for (i in 1:(length(trunc_time_grid))){
      risk_set <- dat[dat$time > time_grid[i],]# maybe this should be >= i+1? Need to think more carefully about this. obv in the limit, doesn't matter
      if (ids){
        id_i <- which(dat$time > time_grid[i])
        id_vec <- c(id_vec, id_i)
      }
      risk_set_covariates <- risk_set[,1:ncol(X)]
      event_indicators <- matrix(ifelse(risk_set$entry <= time_grid[i + 1 ], 1, 0))
      dummies <- matrix(0, ncol = length(trunc_time_grid), nrow = nrow(risk_set))
      dummies[,i] <- 1
      newdata <- as.matrix(cbind(dummies, risk_set_covariates, event_indicators))
      stacked <- rbind(stacked, newdata)
    }
    stacked <- stacked[-1,]
    risk_set_names <- paste0("risk_set_", seq(1, (length(time_grid))))
    colnames(stacked)[1:(length(time_grid))] <- risk_set_names
    colnames(stacked)[ncol(stacked)] <- "event_indicators"
  }
  if (ids){
    id_vec <- id_vec[-1]
    ids <- id_vec
  }
  stacked <- data.frame(stacked)
  return(list(stacked = stacked, ids = ids))
}

