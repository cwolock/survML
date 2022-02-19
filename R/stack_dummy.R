#' Stack a dataset
#'
#' @return A stacked dataset
#' @noRd
stack_dummy <- function(time, X, time_grid){
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
  stacked <- data.frame(stacked)
  return(stacked)
}

