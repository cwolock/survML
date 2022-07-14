#' Stack a dataset, using the hazard (to be merged with stack.R)
#'
#' @return A stacked dataset
#' @noRd
stack_haz <- function(time, event, X, time_grid, entry, time_basis){
  trunc_time_grid <- time_grid[-length(time_grid)] # do I need to truncate if treating time as continuous? look at this later
  dat <- data.frame(X, event = event, time = time)

  if (time_basis == "continuous"){
    ncol_stacked <- ncol(X) + 2 # covariates, time, binary outcome
  } else if (time_basis == "dummy"){
    ncol_stacked <- ncol(X) + length(trunc_time_grid) + 1 # covariates, risk set dummies, binary outcome
  }
  stacked <- matrix(NA, ncol = ncol_stacked, nrow = 1)
  for (i in 1:(length(trunc_time_grid))){ # can change this to not do anything in last time bin
    if (is.null(entry)){ # no entry (truncation variable) given
      #risk_set <- dat[dat$time > time_grid[i],]
      risk_set <- dat[dat$time >= time_grid[i],]
    } else{ # entry given
      #risk_set <- dat[dat$time > time_grid[i] & entry <= time_grid[i+1],]
      risk_set <- dat[dat$time >= time_grid[i] & entry <= time_grid[i],]
    }
    risk_set_covariates <- risk_set[,1:ncol(X)]
    #event_indicators <- matrix(ifelse(risk_set$time <= time_grid[i + 1 ] & risk_set$event == 1, 1, 0))
    event_indicators <- matrix(ifelse(risk_set$time < time_grid[i + 1 ] & risk_set$event == 1, 1, 0))
    if (time_basis == "continuous"){
      #t <- rep(time_grid[i + 1], nrow(event_indicators))
      t <- rep(time_grid[i], nrow(event_indicators))
    } else if (time_basis == "dummy"){
      t <- matrix(0, ncol = length(trunc_time_grid),
                  nrow = nrow(risk_set))
      t[,i] <- 1
    }
    newdata <- as.matrix(cbind(t, risk_set_covariates, event_indicators))
    stacked <- rbind(stacked, newdata)
    colnames(stacked)[ncol(stacked)] <- "event_indicators"
    if (time_basis == "dummy"){
      risk_set_names <- paste0("risk_set_", seq(1, (length(trunc_time_grid))))
      colnames(stacked)[1:(length(trunc_time_grid))] <- risk_set_names
    }
  }
  stacked <- stacked[-1,]
  stacked <- data.frame(stacked)
  return(stacked)
}

