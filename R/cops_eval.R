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
#' @param alpha level of prediction intervals to be evaluated
#' @param crisp logical specifying whether prediction intervals are to be obtained from
#' the crips CDF (\code{crisp = TRUE}) or from the bounds of the predictive system (\code{crisp = FALSE});
#' default is \code{crisp = FALSE}.
#'
#'
#' @returns
#'
#' \code{pit()} returns a vector of PIT values
#'
#' \code{crps()} returns a vector of CRPS values
#'
#' \code{thickness()} returns a vector of thickness values
#'
#' \code{threshcal()} returns a matrix of forecast exceedance probabilities at the
#'  thresholds in \code{thresholds}
#'
#' \code{plot()} returns a base plot object
#'
#' \code{get_pint()} returns a matrix containing lower and upper bounds of prediction intervals
#'
#' \code{cov()} returns a vector of coverage indicators
#'
#' \code{is()} returns a vector of Interval Score values
#'
#' \code{width()} returns a vector of prediction interval widths
#'
#'
#' @details
#' Details to be added here
#'
#' \code{plot()} plots a conformal predictive system in \code{fit}, with \code{index}
#'
#' @references
#'
#' Allen, S., Gavrilopolous, G., Henzi, A. and J. Ziegel (2025):
#' `In-sample calibration yields conformal calibration guarantees',
#' \emph{arXiv pre-print} arXiv:2503.03841
#' \doi{10.48550/arXiv.2503.03841}
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


#' @rdname cops_eval
#' @export
get_pint <- function(fit, alpha, crisp) {
  UseMethod("get_pint")
}


#' @rdname cops_eval
#' @export
coverage <- function(fit, y, alpha, average, crisp) {
  UseMethod("coverage")
}


#' @rdname cops_eval
#' @export
int_score <- function(fit, y, alpha, crisp) {
  UseMethod("int_score")
}


#' @rdname cops_eval
#' @export
width <- function(fit, alpha, crisp) {
  UseMethod("width")
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
pit.cops <- function(fit, y = fit$y_out) {
  n <- length(y)
  pit0 <- function(x, y, z) stepfun(x = x, y = c(0, y))(z)
  if (n == 1) {
    out <- pit0(fit$points, fit$cdf_crisp, y)
  } else {
    out <- sapply(1:length(y), function(i) pit0(fit$points[, i], fit$cdf_crisp[, i], y[i]))
  }
  return(out)
}


#' @exportS3Method crps cops
crps.cops <- function(fit, y = fit$y_out) {
  n <- length(y)

  if (n > 1) {
    m <- nrow(fit$points)
    ens <- fit$points[-c(1, m), ] |> t()
    w <- pmax(apply(fit$cdf_crisp, 2, diff)[-(m - 1), ], 0) |> t()
  } else {
    m <- length(fit$points)
    ens <- fit$points[-c(1, m)]
    w <- pmax(diff(fit$cdf_crisp)[-(m - 1)], 0)
  }
  out <- scoringRules::crps_sample(y, ens, w = w)
  return(out)
}


#' @exportS3Method threshcal cops
threshcal.cops <- function(fit, thresholds, y = fit$y_out) {
  n <- length(y)
  pit0 <- function(x, y, z) stepfun(x = x, y = c(0, y))(z)
  if (n > 1) {
    out <- sapply(thresholds, function(t) {
      sapply(1:length(y), function(i) {
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


#' @exportS3Method get_pint cops
get_pint.cops <- function(fit, alpha, crisp = FALSE){

  x <- fit$points
  if (crisp) {
    F_l <- F_u <- fit$cdf_crisp
  } else {
    F_l <- fit$cdf_upper
    F_u <- fit$cdf_lower
  }

  cops_pint_i <- function(x_i, F_l_i, F_u_i, a){
    if (any(diff(x_i) < 0)){
      ord <- order(x_i)
      x_i <- x_i[ord]
      F_l_i <- F_l_i[ord]
      F_u_i <- F_u_i[ord]
    }

    ind_low <- which(F_u_i >= a/2)[1]
    ind_low <- pmax(2, ind_low)
    lower <- if (!is.na(ind_low)) x_i[ind_low] else NA

    ind_upp  <- which(F_l_i >= 1 - a/2)[1]
    ind_upp <- pmin(length(x_i) - 1, ind_upp)
    upper <- if (!is.na(ind_upp)) x_i[ind_upp] else NA

    c(Lower = lower, Upper = upper)
  }

  if (is.matrix(x)) {
    pints <- sapply(1:ncol(x), function(i) cops_pint_i(x[, i], F_l[, i], F_u[, i], alpha))
  } else {
    pints <- cops_pint_i(x, F_l, F_u, alpha)
  }
  pints <- t(pints)
  colnames(pints) <- c("Lower", "Upper")
  return(pints)
}


#' @exportS3Method coverage cops
coverage.cops <- function(fit, y = fit$y_out, alpha, average = TRUE, crisp = FALSE){
  pints <- get_pint(fit, alpha, crisp)
  cov <- as.numeric(y >= pints[, 1] & y <= pints[, 2])
  if (average) {
    return(mean(cov, na.rm = TRUE))
  } else {
    return(cov)
  }
}


#' @exportS3Method int_score cops
int_score.cops <- function(fit, y = fit$y_out, alpha, crisp = FALSE){
  pints <- get_pint(fit, alpha, crisp)
  is <- (pints[, 2] - pints[, 1]) + (2/alpha)*(y < pints[, 1])*(pints[, 1] - y) + (2/alpha)*(y > pints[, 2])*(y - pints[, 2])
  return(is)
}


#' @exportS3Method width cops
width.cops <- function(fit, alpha, crisp = FALSE){
  pints <- get_pint(fit, alpha, crisp)
  width <- pints[, 2] - pints[, 1]
  return(width)
}
