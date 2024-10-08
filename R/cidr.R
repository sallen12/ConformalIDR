#' Online conformal IDR
#'
#' Computes IDR on a training data set and sequentially adds observations and
#' updates the fit.
#'
#' @param x training covariates
#' @param y training observations
#' @param x_out new covariates
#' @param y_out new observations
#' @param online whether to sequentially update the training data set. Default
#'     is \code{TRUE}. If \code{FALSE}, the same training data is used for
#'     predicting with all the \code{x_out}.
#' @param weights if \code{online} is \code{FALSE}, a vector of the same length
#'     as \code{x} plus one containing non-negative weights, with the last
#'     weight for the new covariate. If \code{online} is \code{TRUE}, a list
#'     of vectors of weights, with lengths \code{length(x) + 1} up to
#'     \code{length(x) + length(x_out)}. If omitted, equal weights are used.
#'
#' @return
#' Returns an object of class \code{conformal_idr}. It contains the jump points
#' of the conditional CDFs, which are either a matrix with \code{length(x_out)}
#' columns for \code{online == TRUE}, and a vector otherwise. The values of the
#' conditional CDFs with the corresponding jump points are returned as matrices,
#' where \code{cdf_lwr}, \code{cdf_upr}, and \code{cdf_oos} are the CDFs at the
#' lower bound, upper bound, and interpolation of the two for the given value
#' of \code{x_out[i]}, and \code{cdf_lcnf}, \code{cdf_ucnf} are the lower and
#' upper bounds from the conformal IDR.
#'
#' @note
#' The sequential computation of ranks requires that g++ version >= 3.3.1 is
#' installed. The code was run on a Ubuntu system and adaptations to lines 5
#' and 6 of online_idr_computation may be necessary on different machines.
#' @name cidr
NULL

n <- 1000
n_out <- 10
x <- round(rnorm(n), 2)
y <- round(rnorm(n, x), 2)
x_out <- round(rnorm(n_out), 2)
y_out <- round(rnorm(n_out, x_out), 2)
online <- TRUE
weights <- vector("list", n_out)
for (j in seq_along(weights)) weights[[j]] <- rep(1, n + j)


#' @rdname cidr
#' @export
conformal_idr <- function(x, y, x_out, online = FALSE, weights = NULL) {
  y_out <- NA
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
    n_thr <- length(unique(c(y, y_out)))
    n_x <- length(unique(c(x, x_out)))
    out <- cidr_sequential(
      x_r = x_r,
      x_out = x_out,
      w = w,
      W = W,
      w_out = w_out,
      pos_x = pos_x,
      y_unique_r = y_unique_r,
      y_out = y_out,
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
      y_out = y_out,
      n_thr = n_thr,
      n_x = n_x
    )
  }

  out <- c(out, list(x = x, y = y, x_out = x_out))

  structure(out, class = "conformal_fit")
}
