#' Least squares prediction machine
#'
#' Computes conformal prediction set using the least squares prediction algorithm.
#'
#' @param y training observations
#' @param X training covariates
#' @param X_ts test covariates
#'
#' @return
#' LSPM fit based on y and X, applied to X_ts
#'
#' @name lspm
NULL

#' @rdname lspm
#' @export
fit_lspm <- function(y, X, X_ts) {
  fit <- lspm(y, as.matrix(X), as.matrix(X_ts))
  return(fit)
}
