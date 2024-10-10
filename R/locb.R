#' Local binning
#'
#' Performs conformal prediction using local binning and Mondrian predictive systems.
#' k-means clustering is used to bin the covariates, and separate predictive
#' distributions are then constructed for observations in each bin.
#'
#' @inheritParams cops
#' @param k number of bins
#'
#' @return
#' An object of class \code{conformal_fit} containing the fit of the local binning.
#'
#'
#' @references
#'
#' \emph{Conformal predictive systems:}
#'
#' Vovk, V., Gammerman, A. and G. Shafer (2022):
#' `Algorithmic learning in a random world',
#' Second Series, Chapter 7.
#' \doi{10.1007/978-3-031-06649-8}
#'
#' \emph{Local binning:}
#'
#' Bostroem, H., Johansson, U. and T. Loefstroem (2021):
#' `Mondrian conformal predictive systems',
#' \emph{Proceedings of Machine Learning Research} 152: 1--15.
#'
#'
#' @seealso \link{cops} \link{cidr} \link{locb}
#'
#' @author Sam Allen
#'
#' @name locb
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
