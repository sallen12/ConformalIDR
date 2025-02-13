#' Conformal binning
#'
#' Derive conformal predictive systems based on binning. k-means clustering is used
#' to bin the covariates, and separate predictive systems are then constructed for each bin.
#'
#'
#' @inheritParams cops
#' @param k number of bins
#'
#' @return
#' \code{conformal_bin()} returns a list of objects of class \code{"cops"} containing
#' the conformal binning fit. See \code{\link{cops}} for details.
#'
#' @details
#' Details to be added here
#'
#'
#' @references
#'
#' \emph{Conformal binning:}
#'
#' Allen, S., Gavrilopolous, G., Henzi, A. and J. Ziegel (2024+):
#' `Conformal isotonic distributional regression'.
#'
#'
#' \emph{Mondrian predictive systems:}
#'
#' Bostroem, H., Johansson, U. and T. Loefstroem (2021):
#' `Mondrian conformal predictive systems',
#' \emph{Proceedings of Machine Learning Research} 152: 1--15.
#'
#'
#' @seealso \code{\link{cops}} \code{\link{cidr}} \code{\link{lspm}}
#'
#' @author Sam Allen
#'
#' @examples
#'
#' n <- 1000
#' x <- rnorm(n)
#' y <- rnorm(n, x, exp(x))
#'
#' N <- 100
#' x_out <- rnorm(N)
#' y_out <- rnorm(N, x_out, exp(x_out))
#' plot(x_out, y_out)
#'
#' fit <- conformal_bin(x, y, x_out, y_out, k = 10)
#' fit1 <- conformal_bin(x, y, x_out, y_out)
#' fit5 <- conformal_bin(x, y, x_out, y_out, k = 5)
#' fit100 <- conformal_bin(x, y, x_out, y_out, k = 100)
#'
#' # compare for different k
#' crps_vec <- c("k = 1" = sapply(fit1, function(x) x$crps) |> mean(),
#'               "k = 5" = sapply(fit5, function(x) x$crps) |> mean(),
#'               "k = 10" = sapply(fit, function(x) x$crps) |> mean(),
#'               "k = 100" = sapply(fit100, function(x) x$crps) |> mean())
#' print(crps_vec)
#'
#' @importFrom stats kmeans
#' @name cbin
#' @export
conformal_bin <- function(x, y, x_out, y_out = NULL, online = FALSE, weights = NULL, k = 1) {

  check_cops_args(x, y, x_out, y_out, "cbin", online, weights, k = k)

  eval <- !is.null(y_out)

  out_list <- lapply(seq_along(x_out), function(i) {

    xn1 <- c(x, x_out[i])
    n1 <- length(xn1)

    out <- kmeans(as.matrix(xn1), k)
    tr_cl <- out$cluster[-n1]
    ts_cl <- out$cluster[n1]

    points <- sapply(1:k, function(cl) y[tr_cl == cl] |> sort())
    points <- lapply(points, function(z) c(-Inf, z, Inf))
    cdfs <- lapply(points, function(z) sample_to_bounds(length(z) - 1))

    out <- c(points = points[ts_cl], cdfs[[ts_cl]],
             list(x = x, y = y, x_out = x_out[i]))
    out <- structure(out, class = "cops")

    if (eval) {
      pcal <- pit(out, y_out[i])
      score <- crps(out, y_out[i])
      thick <- thickness(out)
      out <- c(out, list(y_out = y_out[i], pit = pcal, crps = score, thick = thick))
    }
    out <- c(out, bins = list(x = tr_cl, x_out = ts_cl))
    structure(out, class = "cops")
  })

}
