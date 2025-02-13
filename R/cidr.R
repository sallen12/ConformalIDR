#' Conformal IDR
#'
#' Perform conformal isotonic distributional regression (IDR). Returns a conformal
#' predictive system for a real-valued outcome variable, under the assumption that
#' there is an increasing relationship between the outcome and the covariates.
#'
#'
#' @inheritParams cops
#'
#' @return
#' \code{conformal_cidr()} returns an object of class \code{"cops"} containing the
#' conformal IDR fit. If \code{online = TRUE}, then a list of \code{"cops"} objects are returned.
#' See \code{\link{cops}} for details.
#'
#' @details
#' Details to be added here
#'
#'
#' @references
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
#' @seealso \code{\link{cops}} \code{\link{lspm}} \code{\link{cbin}}
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
#' fit <- conformal_idr(x, y, x_out, y_out)
#' plot(fit)
#'
#' # assign heigher weights to more recent covariates
#' weights <- c((1:n)/n, rep(1, N))
#' fit <- conformal_idr(x, y, x_out, y_out, weights = weights)
#' plot(fit)
#'
#'
#' @note
#' The sequential computation of ranks requires that g++ version >= 3.3.1 is
#' installed.
#'
#' @importFrom stats aggregate
#' @name cidr
#' @export
conformal_idr <- function(x, y, x_out, y_out = NULL, online = FALSE, weights = NULL) {

  check_cops_args(x, y, x_out, y_out, "cidr", online, weights)

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
  out <- structure(out, class = "cops")

  if (!is.null(y_out)) {
    pcal <- pit(out, y_out)
    score <- crps(out, y_out)
    thick <- thickness(out)
    out <- c(out, list(y_out = y_out, pit = pcal, crps = score, thick = thick))
    out <- structure(out, class = "cops")
  }

  return(out)
}

