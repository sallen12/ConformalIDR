#' Local binning
#'
#' Computes conformal prediction set using local binning.
#'
#' @param y training observations
#' @param X training covariates
#' @param X_ts test covariates
#' @param k number of bins
#'
#' @return
#' Local binning fit based on y and X, applied to X_ts
#'
#' @name locb
NULL

#' @rdname locb
#' @export
fit_locb <- function(y, X, X_ts, k = 1) {
  out <- kmeans(X, k)
  tr_cl <- out$cluster
  ts_cl <- clue::cl_predict(out, as.matrix(X_ts))
  return(list(tr = data.frame(y = y, cl = tr_cl), ts = as.vector(ts_cl)))
}
