#' Conformal IDR
#'
#' Perform conformal isotonic distributional regression (IDR) to generate a probabilistic
#' prediction for a real-valued outcome variable. IDR assumes an increasing relationship
#' between the outcome and the covariate.
#'
#' @inheritParams cops
#' @param online logical specifying whether to sequentially update the training
#'  data set. Default is \code{FALSE}. If \code{FALSE}, the same training data is
#'  used for prediction of all the \code{x_out}.
#' @param weights non-negative weights assigned to each covariate and
#'  the new covariate. If omitted, equal weights are used. If \code{online} is \code{FALSE},
#'  \code{weights} is a vector of length \code{length(x) + 1} with the last weight
#'  corresponding to the new covariate. If \code{online} is \code{TRUE}, \code{weights} is a list
#'  of vectors of weights, with lengths \code{length(x) + 1} up to \code{length(x) + length(x_out)}.
#'
#' @return
#' An object of class \code{conformal_fit} containing the conformal IDR fit.
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
#' \emph{Isotonic distributional regression:}
#'
#' Henzi, A., Ziegel, J. and T. Gneiting (2022):
#' `Isotonic distributional regression',
#' \emph{Journal of the Royal Statistical Society: Series B} 83, 963--993.
#' \doi{10.1111/rssb.12450}
#'
#' \emph{Conformal IDR:}
#'
#' Allen, S., Gavrilopolous, G., Henzi, A. and J. Ziegel (2024+):
#' `Conformal isotonic distributional regression'.
#'
#'
#' @seealso \link{cops} \link{lspm} \link{locb}
#'
#' @author Sam Allen
#'
#' @note
#' The sequential computation of ranks requires that g++ version >= 3.3.1 is
#' installed. The code was run on a Ubuntu system and adaptations to lines 5
#' and 6 of online_idr_computation may be necessary on different machines.
#' @name cidr
#' @export
conformal_idr <- function(x, y, x_out, y_out = NULL, online = FALSE, weights = NULL) {
  x_order <- order(x)
  x <- x[x_order]
  y <- y[x_order]
  n <- length(x)
  n_out <- length(x_out)
  if (online) {
    w_names <- paste0("w", seq_len(n_out))
    w_out <- W <- vector("list", n_out)
    if (is.null(weights)) {
      weights <- vector("list", n_out)
      for (j in seq_len(n_out)) {
        weights[[j]] <- rep(1, n + 1)
        W[[j]] <- rep(1, n)
        w_out[[j]] <- rep(1, j)
      }
    } else {
      for (j in seq_len(n_out)) {
        w_out[[j]] <- weights[[j]][(n + 1):(n + j)]
        weights[[j]] <- W[[j]] <-
          weights[[j]][seq_len(n)][x_order]
      }
    }
    names(weights) <- names(W) <- names(w_out) <- w_names
    unique_x <- aggregate(
      lapply(weights, function(x) x[seq_len(n)]),
      by = list(x = x),
      FUN = sum
    )
    w <- as.list(unique_x[names(unique_x) != "x"])
    x_r <- unique_x$x
    rle_x <- rle(x)
    pos_x <- aggregate(
      c(list(ind = rep(x = 0:(length(x_r) - 1), times = rle_x$lengths)), W),
      by = list(y = y),
      FUN = list
    )
    y_unique_r <- pos_x$y
    W <- pos_x[!(names(pos_x) %in% c("y", "ind"))]
    pos_x <- pos_x$ind
    n_thr <- length(unique(c(y, NA)))
    n_x <- length(unique(c(x, x_out)))
    out <- cidr_sequential(
      x_r = x_r,
      x_out = x_out,
      w = w,
      W = W,
      w_out = w_out,
      pos_x = pos_x,
      y_unique_r = y_unique_r,
      y_out = NA,
      n_thr = n_thr,
      n_x = n_x
    )
  } else {
    if (is.null(weights)) {
      weights <- rep(1, n)
      W <- rep(1, n)
      w_out <- 1
    } else {
      w_out <- weights[n + 1]
      weights <- W <- weights[-(n + 1)][x_order]
    }
    unique_x <- aggregate(
      list(w = weights),
      by = list(x = x),
      FUN = sum
    )
    w <- unique_x$w
    x_r <- unique_x$x
    rle_x <- rle(x)
    pos_x <- aggregate(
      list(
        W = W,
        ind = rep(x = 0:(length(x_r) - 1), times = rle_x$lengths)
      ),
      by = list(y = y),
      FUN = list
    )
    W <- pos_x$W
    y_unique_r <- pos_x$y
    pos_x <- pos_x$ind
    n_thr <- length(y_unique_r)
    n_x <- length(x_r)
    out <- cidr_static(
      x_r = x_r,
      x_out = x_out,
      w = w,
      W = W,
      w_out = w_out,
      pos_x = pos_x,
      y_unique_r = y_unique_r,
      y_out = NA,
      n_thr = n_thr,
      n_x = n_x
    )
  }

  out$points <- replicate(length(x_out), out$points)
  out <- c(out, list(x = x, y = y, x_out = x_out))
  out <- structure(out, class = "conformal_fit")

  if (!is.null(y_out)) {
    pcal <- pit(out, y_out)
    score <- crps(out, y_out)
    thick <- thickness(out)
    out <- c(out, list(y_out = y_out, pit = pcal, crps = score, thick = thick))
    out <- structure(out, class = "conformal_fit")
  }

  return(out)
}

