#' Generate cross-fitting and sample-splitting folds
#'
#' @param n Total sample size
#' @param V Number of cross-fitting folds to use
#' @param sample_split Logical, whether or not sample-splitting is being used
#'
#' @return Named list of cross-fitting and sample-splitting folds
#'
#' @export
generate_folds <- function(n, V, sample_split){
  cf_fold_num <- V
  ss_fold_num <- 2*cf_fold_num
  .V <- ifelse(sample_split, ss_fold_num, cf_fold_num)
  cf_folds <- sample(rep(seq_len(.V), length = n)) # 2V of them
  if (sample_split){
    ss_folds <- c(rep(1, .V/2), rep(2, .V/2))
  } else{
    ss_folds <- rep(1, .V)
  }
  ss_folds <- as.numeric(cf_folds %in% which(ss_folds == 2))
  return(list(cf_folds = cf_folds, ss_folds = ss_folds))
}
