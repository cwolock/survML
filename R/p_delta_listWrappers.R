#' List wrappers for p_delta prediction
#'
#' This function lists all prediction algorithms built into the \code{conSurv} package.
#'
#' @usage \code{p_delta_listWrappers()}
#'
#' @return Invisible character vector with the requested functions.
#'
#' @export

p_delta_listWrappers <- function() {
  everything <- sort(unclass(lsf.str(envir = asNamespace("conSurv"), all.names = T)))
  message("All prediction algorithm wrappers in conSurv:\n")
  funs <- everything[grepl(pattern = "^p_delta", everything)]
  funs <- funs[!funs %in% c("p_delta_listWrappers")]
  print(funs)
  invisible(everything)
}
