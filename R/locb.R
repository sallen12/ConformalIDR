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
conformal_bin <- function(x, y, x_out, y_out = NULL, k = 1) {

  out <- kmeans(as.matrix(x), k)
  tr_cl <- out$cluster
  ts_cl <- clue::cl_predict(out, as.matrix(x_out)) |> as.vector()

  points <- sapply(1:k, function(i) y[tr_cl == i] |> sort())
  points <- lapply(points, function(x) c(-Inf, x, Inf))
  cdfs <- lapply(points, function(x) sample_to_bounds(length(x) - 1))

  eval <- !is.null(y_out)
  out_list <- lapply(seq_along(x_out), function(i) {
    bin <- ts_cl[i]
    out <- c(points = points[bin], cdfs[[bin]],
             list(x = x, y = y, x_out = x_out[i]))
    out <- structure(out, class = "conformal_fit")

    if (eval) {
      pcal <- pit(out, y_out[i])
      score <- crps(out, y_out[i])
      thick <- thickness(out)
      out <- c(out, list(y_out = y_out[i], pit = pcal, crps = score, thick = thick))
    }
    out <- c(out, bins = list(x = tr_cl, x_out = ts_cl[i]))
    structure(out, class = "conformal_fit")
  })

}
