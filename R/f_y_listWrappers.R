#' List wrappers for f_y prediction
#'
#' This function lists all prediction algorithms built into the \code{conSurv} package.
#'
#' @usage \code{f_y_listWrappers()}
#'
#' @return Invisible character vector with the requested functions.
#'
#' @export

f_y_listWrappers <- function() {
  everything <- sort(unclass(lsf.str(envir = asNamespace("conSurv"), all.names = T)))
  message("All prediction algorithm wrappers in conSurv:\n")
  funs <- everything[grepl(pattern = "^f_y", everything)]
  funs <- funs[!funs %in% c("f_y_listWrappers")]
  print(funs)
  invisible(everything)
}
