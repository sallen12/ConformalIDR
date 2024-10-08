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
conformal_bin <- function(x, y, x_out, k = 1) {
  if (is.vector(x)) x <- as.matrix(x)
  if (is.vector(x_out)) x_out <- as.matrix(x_out)

  out <- kmeans(x, k)
  tr_cl <- out$cluster
  ts_cl <- clue::cl_predict(out, x_out) |> as.vector()

  n <- length(y) + 1
  points <- sapply(1:k, function(i) y[tr_cl == i] |> sort())

  out <- list(x = x, y = y, x_out = x_out)
  out <- c(out, bins = list(train = tr_cl, test = ts_cl))
  structure(out, class = "conformal_fit")
}


