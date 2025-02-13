#' Conformal Predictive Systems
#'
#' Generic function for generating conformal predictive systems when forecasting real-valued outcomes
#'
#' @param x training covariates
#' @param y training labels
#' @param x_out new covariates
#' @param y_out new labels
#' @param method method used to perform conformal prediction. One of \code{"cidr"}
#' (default), \code{"lspm"}, and \code{"cbin"}.
#' @param online logical specifying whether to sequentially update the training
#'  data set. Default is \code{FALSE}.
#' @param weights non-negative weights assigned to each training covariate and the
#'  new covariate. If omitted, equal weights are used (the default).
#' @param ... additional arguments passed to the chosen method.
#'
#'
#' @return
#' \code{cops()} returns an object of class \code{"cops"}, or a list of \code{"cops"}
#' objects if \code{online = TRUE} or \code{method = "cbin"}.
#'
#' An object of class \code{"cops"} is a list containing:
#' - \code{points}, the jump points of the upper and lower bands of the conformal predictive system
#' - \code{cdf_lwr}, the lower band at the corresponding point value in \code{points}
#' - \code{cdf_upr}, the upper band at the corresponding point value in \code{points}
#' - \code{cdf_cri}, the crisp predictive distribution function between \code{cdf_lwr} and \code{cdf_upr}
#' - \code{thick}, the thickness of the conformal predictive system
#'
#' If \code{y_out} is specified, then the list will additionally contain metrics that
#' evaluate the performance of the predictive distribution \code{cdf_cri} evaluated using \code{y_out}.
#' This includes:
#' - \code{crps}, Continuous Ranked Probability Score (CRPS) values
#' - \code{pit}, Probability Integral Transform (PIT) values
#'
#' The functions \code{crps} and \code{pit} can be applied to \code{"cops"} objects
#' to obtain the CRPS and PIT values of the conformal predictive system. The function
#' \code{threshcal} additionally allows one to assess whether or not the forecast is
#' threshold calibrated. \code{thickness} returns the thickness of the predictive system,
#' while \code{plot} displays the crisp predictive distribution function and corresponding
#' lower and upper bands.
#'
#'
#' @details
#'
#' ## Conformal Predictive Systems
#'
#' Suppose we are interested in forecasting \eqn{Y \in \mathbb{R}}. Conformal predictive
#' systems provide a means to forecast \eqn{Y} with theoretical calibration guarantees.
#'
#' Conformal predictive systems are defined in terms of lower and upper bands,
#' \eqn{\Pi_{l}} and \eqn{\Pi_{u}}, which are increasing functions between 0 and 1,
#' such that \eqn{\Pi_{l}(z) \leq \Pi_{u}(z)} for all \eqn{z \in \mathbb{R}},
#' \eqn{\Pi_{l}(z) \to 0} as \eqn{z \to -\infty}, and \eqn{\Pi_{u}(z) \to 1} as
#' \eqn{z \to \infty}.
#'
#' The bands \eqn{\Pi_{l}} and \eqn{\Pi_{u}} are constructed such that they contain
#' a calibrated predictive distribution \eqn{F}. That is, \eqn{F} is a cumulative
#' distribution function (cdf) that is calibrated according to standard notions of
#' probabilistic forecast calibration, and for which
#' \eqn{\Pi_{l}(z) \leq F(z) \leq \Pi_{u}(z)} for all \eqn{z \in \mathbb{R}}.
#'
#' This calibration guarantee is only useful if the lower and upper bands are close
#' to one another. The \emph{thickness} of a conformal predictive system is defined
#' as the maximum distance between \eqn{\Pi_{l}} and \eqn{\Pi_{u}}.
#'
#' The lower and upper bands will typically depend on covariates \eqn{X}, and
#' the predictive distribution \eqn{F} can be interpreted as an estimate of the
#' conditional cdf of \eqn{Y} given \eqn{X}. In practice, we typically have a sequence
#' of training covariates \eqn{X_{1}, \dots, X_{n}} and corresponding labels \eqn{Y_{1}, \dots, Y_{n}},
#' and wish to predict a new label \eqn{Y_{n+1}} given \eqn{(X_{1}, Y_{1}), \dots, (X_{n}, Y_{n})},
#' and a new covariate \eqn{X_{n+1}}.
#'
#' We may additionally want to predict additional labels \eqn{Y_{n+2}, Y_{n+3}, \dots} given covariates
#' \eqn{X_{n+2}, X_{n+3}, \dots}. This can either be performed sequentially, i.e. \emph{online},
#' so that the training data updates and the prediction for \eqn{Y_{N}}, for \eqn{N > n},
#' depends on all previous covariate-label pairs, \eqn{(X_{1}, Y_{1}), \dots, (X_{N-1}, Y_{N-1})},
#' and \eqn{X_{N}}; or statically, i.e. \emph{offline}, so that the training data does not update
#' and the prediction for \eqn{Y_{N}} depends only on \eqn{(X_{1}, Y_{1}), \dots, (X_{n}, Y_{n})}
#' and \eqn{X_{N}}, irrespective of \eqn{N}.
#'
#' Different methods are available to generate conformal predictive systems, i.e.
#' to construct the lower and upper bands as a function of the covariates. This
#' package provides the functionality for three different methods: the Least
#' Squares Prediction Machine (LSPM), local binning, and conformal isotonic
#' distributional regression (IDR). See the references below for further details.
#'
#'
#' ## Package Functionality
#'
#' The function \code{cops()} provides a generic wrapper that calls one of
#' \code{\link{lspm}}, \code{\link{cbin}}, and \code{\link{cidr}}, depending on
#' the argument \code{method}. \code{method} must be one of \code{"lspm"}, \code{"cbin"},
#' and \code{"cidr"}, corresponding to the LSPM, local binning, and conformal IDR, respectively.
#'
#' All methods take the following arguments:
#' - \code{x}, a vector or matrix of training covariates \eqn{X_{1}, \dots, X_{n}}
#' - \code{y}, a vector of corresponding labels \eqn{Y_{1}, \dots, Y_{n}}
#' - \code{x_out}, a vector or matrix of new covariate(s) \eqn{X_{n+1}, X_{n+2}, \dots}
#' - \code{y_out}, a vector of new labels \eqn{Y_{n+1}, Y_{n+2}, \dots}
#' - \code{online}, a logical specifying whether prediction is performed online (\code{TRUE})
#'  or offline (\code{FALSE})
#'
#' If \code{x} is a vector, it should have the same length as \code{y}, \eqn{n}. If \code{x} is
#' a matrix, it should have \eqn{n} rows, with each row corresponding to the covariate vector
#' \eqn{X_{i}}, for \eqn{i = 1, \dots, n}.
#'
#' Different weights can be assigned to different training covariates, allowing more
#' recent observations to be emphasised during model training, for example. These
#' weights can be entered using the \code{weights} argument. By default, \code{weights = NULL},
#' in which case the same weight is assigned to all covariates.
#'
#' If \code{online} is \code{FALSE}, \code{weights} is a vector of length \code{length(x) + 1}
#' with the last weight corresponding to the new covariate. If \code{online} is \code{TRUE},
#' \code{weights} is a list of vectors of weights, with lengths \code{length(x) + 1} up
#' to \code{length(x) + length(x_out)}.
#'
#' Additional arguments to the methods can be entered as variable arguments via \code{...}.
#' Local binning, for example, requires an additional hyperparameter specifying the number
#' of bins, while the LSPM can be "studentised" by specifying \code{student = TRUE}.
#' See the individual help pages \code{\link{lspm}}, \code{\link{cbin}}, and \code{\link{cidr}}
#' for further details.
#'
#' \code{cops()} outputs an object of class \code{"cops"}, or a list of such objects
#' if \code{online = TRUE} or \code{method = "cbin"}. A \code{"cops"} object is a list
#' containing the lower and upper bands of the estimated conformal predictive system, as
#' well as a \emph{crisp} cdf that lies between the bands. This crisp cdf is essentially
#' an estimate of the calibrated forecast distribution \eqn{F} described above.
#'
#' Note that when forecasting \eqn{Y_{n+1}}, the value of \eqn{Y_{n+1}} is obviously not used to
#' make the prediction. Hence, the input \code{y_out} is not required to derive the
#' conformal predictive systems. The default is therefore that \code{y_out = NULL}.
#'
#' If the new labels \code{y_out} are provided, then the \code{"cops"} object additionally
#' returns evaluation metrics that quantify the performance of the conformal predictive
#' system when predicting \eqn{Y_{n+1}, Y_{n+2}, \dots}. This includes the Continuous
#' Ranked Probability Score (CRPS) of the crisp predictive distribution, as well as
#' Probability Integral Transform (PIT) values to assess its (unconditional) calibration.
#'
#'
#' @references
#'
#' \emph{Conformal predictive systems and the LSPM:}
#'
#' Vovk, V., Gammerman, A. and G. Shafer (2022):
#' `Algorithmic learning in a random world',
#' Second Series, Chapter 7.
#' \doi{10.1007/978-3-031-06649-8}
#'
#'
#' \emph{Conformal IDR and local binning:}
#'
#' Allen, S., Gavrilopolous, G., Henzi, A. and J. Ziegel (2024+):
#' `Conformal isotonic distributional regression'.
#'
#'
#' @seealso \code{\link{cops}} \code{\link{cidr}} \code{\link{cbin}}
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
#' fit <- cops(x, y, x_out, y_out)
#' plot(fit)
#'
#' fit_lspm <- cops(x, y, x_out, y_out, method = "lspm")
#' fit_cbin <- cops(x, y, x_out, y_out, method = "cbin", k = 10)
#'
#' crps_vec <- c(cidr = mean(fit$crps),
#'               lspm = mean(fit_lspm$crps),
#'               cbin = sapply(fit_cbin, function(x) x$crps) |> mean())
#' print(crps_vec)
#'
#' @md
#'
#' @name cops
NULL

#' @rdname cops
#' @export
cops <- function(x, y, x_out, y_out = NULL, method = c("cidr", "lspm", "cbin"),
                 online = FALSE, weights = NULL, ...) {
  method <- match.arg(method)
  if (method == "cidr") {
    conformal_idr(x, y, x_out, y_out, online, weights)
  } else if (method == "lspm") {
    conformal_lspm(x, y, x_out, y_out, online, weights, ...)
  } else if (method == "cbin") {
    conformal_bin(x, y, x_out, y_out, online, weights, ...)
  }
}


check_cops_args <- function(x, y, x_out, y_out, method, online, weights, ...) {

  # y
  if (!is.vector(y) || !is.numeric(y)) stop("'y' must be a numeric vector")

  # x
  if (!(is.vector(x) || is.matrix(x)) || !is.numeric(x))
    stop("'x' must be a numeric vector or matrix")
  if (is.matrix(x)) {
    if (!identical(nrow(x), length(y)))
      stop("'x' must have the same number of rows as the length of 'y'")
  } else {
    if (!identical(length(x), length(y)))
      stop("'x' and 'y' must have the same length")
  }

  # y_out
  if (!is.null(y_out)) {
    if (!is.vector(y_out) || !is.numeric(y_out)) stop("'y_out' must be a numeric vector")
  }

  # x_out
  if (!(is.vector(x_out) || is.matrix(x_out)) || !is.numeric(x_out))
    stop("'x_out' must be a numeric vector or matrix")
  if (is.matrix(x_out)) {
    if (!is.null(y_out) && !identical(nrow(x_out), length(y_out)))
      stop("'x_out' must have the same number of rows as the length of 'y_out'")
    if (is.vector(x))
      stop("'x' and 'x_out' must be of the same type")
    if (!identical(ncol(x), ncol(x_out)))
      stop("'x' and 'x_out' must have the same number of columns")
  } else {
    if (is.matrix(x))
      stop("'x' and 'x_out' must be of the same type")
    if (!is.null(y_out) && !identical(length(x_out), length(y_out)))
      stop("'x_out' and 'y_out' must have the same length")
  }

  # online
  if (!is.logical(online)) stop("'online' must be either TRUE or FALSE")
  if (length(online) > 1) stop("'online' must be either TRUE or FALSE")

  # weights
  if (!is.null(weights)) {

    if (!online) {
      if (!is.vector(weights) || !is.numeric(weights))
        stop("'weights' must be a numeric vector when 'online = FALSE'")
      if (any(weights < 0))
        stop("'weights' must only contain non-negative values")

      if (is.matrix(x)) {
        if (!identical(length(weights), nrow(x) + nrow(x_out)))
          stop("'weights' must have length equal to nrow(x) + nrow(x_out)")
      } else {
        if (!identical(length(weights), length(x) + length(x_out)))
          stop("'weights' must have length equal to length(x) + length(x_out)")
      }
    } else {
      if (!is.list(weights))
        stop("'weights' must be a list when 'online = TRUE'")
      if (!all(sapply(weights, is.numeric)))
        stop("'weights' must be a list containing numeric vectors")
      if (!all(sapply(weights, is.vector)))
        stop("'weights' must be a list containing numeric vectors")
      if (any(sapply(weights, function(w) any(w < 0))))
        stop("'weights' must only contain non-negative values")

      if (is.matrix(x)) {
        n <- nrow(x)
        n_out <- nrow(x_out)
        if (!identical(length(weights), n_out))
          stop("'weights' must have length equal to nrow(x_out)")
        if (!identical(sapply(weights, length), (n + 1):(n + n_out)))
          stop("the elements of 'weights' must have lengths nrow(x) + 1, nrow(x) + 2, ..., nrow(x) + nrow(x_out)")
      } else {
        n <- length(x)
        n_out <- length(x_out)
        if (!identical(length(weights), n_out))
          stop("'weights' must have length equal to length(x_out)")
        if (!identical(sapply(weights, length), (n + 1):(n + n_out)))
          stop("the elements of 'weights' must have lengths length(x) + 1, length(x) + 2, ..., length(x) + length(x_out)")
      }
    }
  }
}

