#' Least Squares Prediction Machine
#'
#' Derive conformal predictive systems using the Least Squares Prediction Machine (LSPM) algorithm.
#'
#'
#' @inheritParams cops
#' @param student logical specifying whether to studentize the least squares residuals. Default is \code{TRUE}.
#'
#' @return
#' \code{conformal_lspm()} returns an object of class \code{"cops"} containing the
#' LSPM fit. If \code{online = TRUE}, then a list of \code{"cops"} objects are returned.
#' See \code{\link{cops}} for details.
#'
#'
#' @details
#' Details to be added here
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
#' fit <- conformal_lspm(x, y, x_out, y_out)
#' plot(fit)
#'
#'
#' @name lspm
#' @export
conformal_lspm <- function(x, y, x_out, y_out = NULL, online = FALSE, weights = NULL, student = TRUE) {

  check_cops_args(x, y, x_out, y_out, "lspm", online, weights, student = student)

  fit <- lspm(y, as.matrix(x), as.matrix(x_out))

  n <- length(y) + 1
  points <- cbind(-Inf, fit, Inf) |> t()
  cdf <- sample_to_bounds(n)
  cdf_lower <- array(cdf$cdf_lower, dim(points))
  cdf_upper <- array(cdf$cdf_upper, dim(points))
  cdf_crisp <- array(cdf$cdf_crisp, dim(points))

  out <- list(points = points, cdf_crisp = cdf_crisp,
              cdf_lower = cdf_lower, cdf_upper = cdf_upper,
              x = x, y = y, x_out = x_out)
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


sample_to_bounds <- function(n) {
  # lower is the same as crisp except for the last entry
  lower <- c((0:(n - 1))/n, (n - 1)/n)
  upper <- c(1/n, (1:n)/n)
  crisp <- c(0, (1:n)/n)
  return(list(cdf_lower = lower, cdf_upper = upper, cdf_crisp = crisp))
}

