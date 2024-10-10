#' Conformal prediction
#'
#' Generic function for performing conformal prediction of real-valued outcomes.
#'
#' @param x training covariates
#' @param y training observations
#' @param x_out new covariates
#' @param y_out new observations
#'
#' @return
#' Conformal prediction object
#'
#' @name cops
NULL

#' @rdname cops
#' @export
cops <- function(x, y, x_out, y_out = NULL, method = c("cidr", "lspm", "locb"), ...) {
  method <- match.arg(method)
  if (method == "cidr") {
    conformal_idr(x, y, x_out, y_out, ...)
  } else if (method == "lspm") {
    conformal_lspm(x, y, x_out, y_out, ...)
  } else if (method == "locb") {
    conformal_bin(x, y, x_out, y_out, ...)
  }
}

