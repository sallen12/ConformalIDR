#' Conformal predictive systems
#'
#' @description
#' Generic function for performing conformal prediction of real-valued outcomes.
#'
#' \code{ConformalIDR} provides methods to obtain conformal predictive systems using
#' a local binning procedure, least squares predction machine, or conformal isotonic
#' distributional regression.
#'
#' @param x_tr vector of covariate values for training
#' @param y_tr vector of observations corresponding to `x_tr`
#' @param x_ts vector of covariate values for prediction
#' @param online logical specifying whether to perform online or offline prediction
#' @param ... additional arguments to other methods
#'
#' @return
#' An object of class \code{conformal_fit}. This is a list containing
#' the jump points of the conditional CDFs (\code{points}), which are either a matrix with
#' \code{length(x_out)} columns or a vector of this length. The values of the
#' conditional CDFs with the corresponding jump points are returned as matrices,
#' where \code{cdf_lwr}, \code{cdf_upr}, and \code{cdf_oos} are the CDFs at the
#' lower bound, upper bound, and interpolation of the two for the given value
#' of \code{x_out[i]}, and \code{cdf_lcnf}, \code{cdf_ucnf} are the lower and
#' upper bounds from the conformal IDR.
#'
#' @details
#' Details to be added here
#'
#' @references
#' References to be added here
#'
#' @seealso
#' \code{\link{cidr}}, \code{\link{lspm}}, \code{\link{locb}}
#'
#' @name conformal
NULL


#' @rdname conformal
#' @export
pit <- function(fit, y) {
  UseMethod("pit")
}


#' @rdname conformal
#' @export
crps <- function(fit, y) {
  UseMethod("crps")
}


#' @rdname conformal
#' @export
threshcal <- function(fit, y, thresholds) {
  UseMethod("threshcal")
}


#' @rdname conformal
#' @export
thickness <- function(fit) {
  UseMethod("thickness")
}


#' @export
plot.conformal_fit <- function(fit, index = 1, ...) {
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


#' @exportS3Method pit conformal_fit
pit.conformal_fit <- function(fit, obs) {
  n <- length(obs)
  pit0 <- function(x, y, z) stepfun(x = x, y = c(0, y))(z)
  if (n == 1) {
    out <- pit0(fit$points, fit$cdf_crisp, obs)
  } else {
    out <- sapply(1:length(obs), function(i) pit0(fit$points[, i], fit$cdf_crisp[, i], obs[i]))
  }
  return(out)
}


#' @exportS3Method crps conformal_fit
crps.conformal_fit <- function(fit, obs) {
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


#' @exportS3Method threshcal conformal_fit
threshcal.conformal_fit <- function(fit, thresholds, obs = fit$y_out) {
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


#' @exportS3Method thickness conformal_fit
thickness.conformal_fit <- function(fit) {
  if (is.vector(fit$points)) {
    out <- max(abs(fit$cdf_lower - fit$cdf_upper))
  } else {
    out <- apply(abs(fit$cdf_lower - fit$cdf_upper), 2, max)
  }
  return(out)
}

