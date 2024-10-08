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
#' An object of type 'conformal_fit'
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
  cdf_oos <- fit$cdf_oos
  cdf_lcnf <- fit$cdf_lcnf
  cdf_ucnf <- fit$cdf_ucnf
  if (!is.null(dim(points))) points <- points[, index]
  mi <- min(points[is.finite(points)]) - 1
  ma <- max(points[is.finite(points)]) + 1
  points <- c(mi - 1, pmin(pmax(points, mi), ma), ma + 1)
  plot(
    points,
    c(0, cdf_oos[, index], 1),
    type = "s",
    xlim = c(mi + 1, ma - 1),
    xlab = "Threshold",
    ylab = "CDFs",
    ...
  )
  lines(points, c(0, cdf_lcnf[, index], 1), type = "s", lty = 5, col = "red")
  lines(points, c(0, cdf_ucnf[, index], 1), type = "s", lty = 5, col = "red")
}


#' @exportS3Method pit conformal_fit
pit.conformal_fit <- function(fit, obs) {
  n <- length(obs)
  if (is.vector(fit$points)) fit$points <- replicate(n, fit$points)
  pit0 <- function(x, y, z) stepfun(x = x, y = c(0, y))(z)
  pit <- sapply(1:length(obs), function(i) pit0(fit$points[, i], fit$cdf_oos[, i], obs[i]))
  return(pit)
}


#' @exportS3Method crps conformal_fit
crps.conformal_fit <- function(fit, obs) {
  n <- length(obs)
  if (is.vector(fit$points)) fit$points <- replicate(n, fit$points)
  ens <- fit$points[-c(1, nrow(fit$points)), ] |> t()
  if (length(obs) > 1) {
    w <- pmax(apply(fit$cdf_oos, 2, diff)[-1, ], 0)  |> t()
  } else {
    w <- pmax(diff(fit$cdf_oos)[-1], 0)
  }
  crps <- scoringRules::crps_sample(obs, ens, w = w)
  return(crps)
}


#' @exportS3Method threshcal conformal_fit
threshcal.conformal_fit <- function(fit, obs, thresholds) {
  n <- length(obs)
  if (is.vector(fit$points)) fit$points <- replicate(n, fit$points)
  pit0 <- function(x, y, z) stepfun(x = x, y = c(0, y))(z)
  F_t <- sapply(thresholds, function(t) {
    sapply(1:length(obs), function(i) {
      pit0(fit$points[, i], fit$cdf_oos[, i], t)
    })
  })
  return(F_t)
}


#' @exportS3Method thickness conformal_fit
thickness.conformal_fit <- function(fit) {
  thicc <- apply(abs(fit$cdf_lcnf - fit$cdf_ucnf), 2, max)
  return(thicc)
}



# 1. aggregate weights over training x (single vector w or list)
# 2. for non-online, do weighted IDR (probably no replacement necessary)
# 3. for online, replace aggregated weights in each loop

plot_data_cidr <- function(fit, indices) {
  points <- fit$points
  cdf_oos <- fit$cdf_oos
  cdf_lwr <- fit$cdf_lwr
  cdf_upr <- fit$cdf_upr
  cdf_lcnf <- fit$cdf_lcnf
  cdf_ucnf <- fit$cdf_ucnf
  n_ind <- length(indices)
  out <- vector("list", length = n_ind)
  if (is.matrix(points)) {
    for (i in seq_len(n_ind)) {
      tmp_points <- points[, indices[i]]
      mi <- min(tmp_points[is.finite(tmp_points)]) - 1
      ma <- max(tmp_points[is.finite(tmp_points)]) + 1
      tmp_points <- c(mi - 1, pmin(pmax(tmp_points, mi), ma), ma + 1)
      out[[i]] <- data.frame(
        x_out = unname(fit$x_out[indices[i]]),
        points = tmp_points,
        cdf_lwr = c(0, cdf_lwr[, indices[i]], 1),
        cdf_upr = c(0, cdf_upr[, indices[i]], 1),
        cdf_oos = c(0, cdf_oos[, indices[i]], 1),
        cdf_lcnf = c(0, cdf_lcnf[, indices[i]], 1),
        cdf_ucnf = c(0, cdf_ucnf[, indices[i]], 1)
      )
    }
  } else {
    mi <- min(points[is.finite(points)]) - 1
    ma <- max(points[is.finite(points)]) + 1
    points <- c(mi - 1, pmin(pmax(points, mi), ma), ma + 1)
    for (i in seq_len(n_ind)) {
      out[[i]] <- data.frame(
        x_out = unname(fit$x_out[indices[i]]),
        points = points,
        cdf_lwr = c(0, cdf_lwr[, indices[i]], 1),
        cdf_upr = c(0, cdf_upr[, indices[i]], 1),
        cdf_oos = c(0, cdf_oos[, indices[i]], 1),
        cdf_lcnf = c(0, cdf_lcnf[, indices[i]], 1),
        cdf_ucnf = c(0, cdf_ucnf[, indices[i]], 1)
      )
    }
  }
  do.call(rbind, out)
}

