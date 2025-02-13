#' Evaluation of Conformal Predictive Systems
#'
#' @description
#' Generic functions for evaluating conformal predictive systems. This includes functions
#' to calculate the Continuous Ranked Probability Score (CRPS) and thickness of
#' the predictive systems, obtain Probability Integral Transform (PIT) values, assess
#' threshold calibration, and to plot the predictive systems.
#'
#'
#' @param fit object of class \code{"cops"} to be evaluated
#' @param y labels that are predicted by \code{fit}. Default is to use \code{fit$y_out}, if available
#' @param thresholds thresholds at which to evaluate threshold calibration
#'
#'
#' @returns
#'
#' \code{pit()} returns a vector of PIT values
#' \code{crps()} returns a vector of CRPS values
#' \code{thickness()} returns a vector of thickness values
#' \code{threshcal()} returns a matrix of forecast exceedance probabilities at the
#'  thresholds in \code{thresholds}
#' \code{plot()} returns a base plot object
#'
#'
#' @details
#' Details to be added here
#'
#' \code{plot()} plots a conformal predictive system in \code{fit}, with \code{index}
#'
#' @references
#'
#' Allen, S., Gavrilopolous, G., Henzi, A. and J. Ziegel (2024+):
#' `Conformal isotonic distributional regression'.
#'
#'
#' @seealso
#' \code{\link{cops}}, \code{\link{cidr}}, \code{\link{lspm}}, \code{\link{cbin}}
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
#'
#' fit <- conformal_idr(x, y, x_out, y_out)
#'
#' ## plot conformal predictive system
#' plot(fit)
#' plot(fit, index = 2)
#' plot(fit, main = "Example Conformal Predictive System", ylim = c(-0.1, 1.1))
#'
#'
#' ## calculate CRPS values
#' crps(fit) |> mean()
#'
#' fit2 <- conformal_idr(x, y, x_out) # if y_out not provided, y must be given as an argument in crps()
#' crps(fit2, y_out) |> mean()
#'
#'
#' ## calculate PIT values
#' pit(fit) |> mean()
#' pit(fit2, y_out) |> mean()
#'
#'
#' ## calculate thickness
#' thickness(fit) |> mean()
#' thickness(fit2) |> mean()
#'
#'
#' ## threshold calibration
#' thresholds <- quantile(y, seq(0.1, 0.9, 0.1))
#' F_t <- threshcal(fit, thresholds)
#' F_t2 <- threshcal(fit2, thresholds, y_out)
#'
#' @importFrom graphics lines
#' @importFrom stats stepfun
#' @name cops_eval
NULL


#' @rdname cops_eval
#' @export
pit <- function(fit, y) {
  UseMethod("pit")
}


#' @rdname cops_eval
#' @export
crps <- function(fit, y) {
  UseMethod("crps")
}


#' @rdname cops_eval
#' @export
threshcal <- function(fit, thresholds, y) {
  UseMethod("threshcal")
}


#' @rdname cops_eval
#' @export
thickness <- function(fit) {
  UseMethod("thickness")
}


#' @exportS3Method pit cops
plot.cops <- function(fit, index = 1, ...) {
  points <- fit$points
  cdf_crisp <- fit$cdf_crisp
  cdf_lower <- fit$cdf_lower
  cdf_upper <- fit$cdf_upper
  if (is.vector(points)) {n <- 1} else {n <- ncol(points)}
  if (n > 1) {
    points <- points[, index]
    cdf_crisp <- cdf_crisp[, index]
    cdf_lower <- cdf_lower[, index]
    cdf_upper <- cdf_upper[, index]
  }
  mi <- min(points[is.finite(points)]) - 1
  ma <- max(points[is.finite(points)]) + 1
  points <- c(mi - 1, pmin(pmax(points, mi), ma), ma + 1)
  plot(
    points,
    c(0, cdf_crisp, 1),
    type = "s",
    xlim = c(mi + 1, ma - 1),
    xlab = "Threshold",
    ylab = "CDFs",
    ...
  )
  lines(points, c(0, cdf_lower, 1), type = "s", lty = 5, col = "red")
  lines(points, c(0, cdf_upper, 1), type = "s", lty = 5, col = "red")
}


#' @exportS3Method pit cops
pit.cops <- function(fit, obs = fit$y_out) {
  n <- length(obs)
  pit0 <- function(x, y, z) stepfun(x = x, y = c(0, y))(z)
  if (n == 1) {
    out <- pit0(fit$points, fit$cdf_crisp, obs)
  } else {
    out <- sapply(1:length(obs), function(i) pit0(fit$points[, i], fit$cdf_crisp[, i], obs[i]))
  }
  return(out)
}


#' @exportS3Method crps cops
crps.cops <- function(fit, obs = fit$y_out) {
  n <- length(obs)

  if (n > 1) {
    m <- nrow(fit$points)
    ens <- fit$points[-c(1, m), ] |> t()
    w <- pmax(apply(fit$cdf_crisp, 2, diff)[-(m - 1), ], 0) |> t()
  } else {
    m <- length(fit$points)
    ens <- fit$points[-c(1, m)]
    w <- pmax(diff(fit$cdf_crisp)[-(m - 1)], 0)
  }
  out <- scoringRules::crps_sample(obs, ens, w = w)
  return(out)
}


#' @exportS3Method threshcal cops
threshcal.cops <- function(fit, thresholds, obs = fit$y_out) {
  n <- length(obs)
  pit0 <- function(x, y, z) stepfun(x = x, y = c(0, y))(z)
  if (n > 1) {
    out <- sapply(thresholds, function(t) {
      sapply(1:length(obs), function(i) {
        pit0(fit$points[, i], fit$cdf_crisp[, i], t)
      })
    })
  } else {
    out <- sapply(thresholds, function(t) pit0(fit$points, fit$cdf_crisp, t))
  }

  return(out)
}


#' @exportS3Method thickness cops
thickness.cops <- function(fit) {
  if (is.vector(fit$points)) {
    out <- max(abs(fit$cdf_lower - fit$cdf_upper))
  } else {
    out <- apply(abs(fit$cdf_lower - fit$cdf_upper), 2, max)
  }
  return(out)
}

