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
conformal_lspm <- function(x, y, x_out) {
  if (is.vector(x)) x <- as.matrix(x)
  if (is.vector(x_out)) x_out <- as.matrix(x_out)

  fit <- lspm(y, x, x_out)

  n <- length(y) + 1
  points <- cbind(-Inf, fit, Inf) |> t()
  cdf <- sample_to_bounds(n)
  cdf_lcnf <- array(cdf$lcnf, dim(points))
  cdf_ucnf <- array(cdf$ucnf, dim(points))
  cdf_oos <- array(cdf$oos, dim(points))

  out <- list(points = points, cdf_oos = cdf_oos, cdf_lcnf = cdf_lcnf,
              cdf_ucnf = cdf_ucnf, x = x, y = y, x_out = x_out)
  structure(out, class = "conformal_fit")
}


sample_to_bounds <- function(n) {
  lcnf <- c((0:(n - 1))/n, (n - 1)/n)
  ucnf <- c(1/n, (1:n)/n)
  oos <- c(0, (1:n)/n)
  return(list(lcnf = lcnf, ucnf = ucnf, oos = oos))
}

