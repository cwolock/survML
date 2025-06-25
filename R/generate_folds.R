#' Generate cross-fitting and sample-splitting folds
#'
#' @param n Total sample size
#' @param V Number of cross-fitting folds to use
#' @param sample_split Logical, whether or not sample-splitting is being used
#' @param cf_folds Cross-fitting folds, if already provided by user.
#' @param cluster_IDs Cluster membership IDs for each of the \code{n} observations. If provided, cross-fitting and
#' sample-splitting folds will respect cluster membership, so each member of a cluster is assigned to the same fold.
#'
#' @return Named list of cross-fitting and sample-splitting folds
#'
#' @export
generate_folds <- function(n, V, sample_split, cf_folds = NULL, cluster_IDs = NULL){
  cf_fold_num <- V
  ss_fold_num <- 2*cf_fold_num
  .V <- ifelse(sample_split, ss_fold_num, cf_fold_num)
  if (is.null(cf_folds)){ # if cf_folds are not provided, generate them
    if (!is.null(cluster_IDs)){ # if cluster IDs are provided, make cf_folds at the cluster level
      # i.e., every observation in a cluster is assigned to the same cf_fold
      unique_cluster_IDs <- unique(cluster_IDs)
      clust_cf_folds <- sample(rep(seq_len(.V), length = length(unique_cluster_IDs)))
      cf_folds <- unlist(lapply(cluster_IDs, FUN = function(x) clust_cf_folds[which(unique_cluster_IDs == x)]))
    } else{ # if no cluster IDs are provided, make cf_folds at the observation level
      cf_folds <- sample(rep(seq_len(.V), length = n)) # 2V of them
    }
  }

  if (sample_split){
    ss_folds <- c(rep(1, .V/2), rep(2, .V/2))
  } else{
    ss_folds <- rep(1, .V)
  }

  # if (sample_split & cf_fold_num == 1){
  # ss_folds <- sample(c(rep(1, ceiling(n/2)), rep(2, floor(n/2))))
  # } else{
  ss_folds <- as.numeric(cf_folds %in% which(ss_folds == 2)) + 1
  # }

  return(list(cf_folds = cf_folds, ss_folds = ss_folds))
}
