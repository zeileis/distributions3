#' Create a reversed Weibull distribution
#'
#' The reversed (or negated) Weibull distribution is a special case of the
#' `\link{GEV}` distribution, obtained when the GEV shape parameter \eqn{\xi}
#' is negative.  It may be referred to as a type III extreme value
#' distribution.
#'
#' @param location The location (maximum) parameter \eqn{m}.
#'   `location` can be any real number.  Defaults to `0`.
#' @param scale The scale parameter \eqn{s}.
#'   `scale` can be any positive number.  Defaults to `1`.
#' @param shape The scale parameter \eqn{\alpha}.
#'   `shape` can be any positive number.  Defaults to `1`.
#'
#' @return A `RevWeibull` object.
#' @export
#'
#' @family continuous distributions
#'
#' @details
#'
#'   We recommend reading this documentation on
#'   <https://zeileis.github.io/distributions3/>, where the math
#'   will render with additional detail and much greater clarity.
#'
#'   In the following, let \eqn{X} be a reversed Weibull random variable with
#'   location parameter  `location` = \eqn{m}, scale parameter `scale` =
#'   \eqn{s}, and shape parameter `shape` = \eqn{\alpha}.
#'   An RevWeibull(\eqn{m, s, \alpha}) distribution is equivalent to a
#'   `\link{GEV}`(\eqn{m - s, s / \alpha, -1 / \alpha}) distribution.
#'
#'   If \eqn{X} has an RevWeibull(\eqn{m, \lambda, k}) distribution then
#'   \eqn{m - X} has a `\link{Weibull}`(\eqn{k, \lambda}) distribution,
#'   that is, a Weibull distribution with shape parameter \eqn{k} and scale
#'   parameter \eqn{\lambda}.
#'
#'   **Support**: \eqn{(-\infty, m)}.
#'
#'   **Mean**: \eqn{m + s\Gamma(1 + 1/\alpha)}.
#'
#'   **Median**: \eqn{m + s(\ln 2)^{1/\alpha}}{m + s(\ln 2)^(1/\alpha)}.
#'
#'   **Variance**:
#'   \eqn{s^2 [\Gamma(1 + 2 / \alpha) - \Gamma(1 + 1 / \alpha)^2]}.
#'
#'   **Probability density function (p.d.f)**:
#'
#'   \deqn{f(x) = \alpha s ^ {-1} [-(x - m) / s] ^ {\alpha - 1}%
#'         \exp\{-[-(x - m) / s] ^ {\alpha} \}}{%
#'        f(x) = (\alpha / s) [-(x - m) / s] ^ (\alpha - 1)%
#'         exp{-[-(x - m) / s] ^ \alpha}}
#'   for \eqn{x < m}.  The p.d.f. is 0 for \eqn{x \geq m}{x >= m}.
#'
#'   **Cumulative distribution function (c.d.f)**:
#'
#'   \deqn{F(x) = \exp\{-[-(x - m) / s] ^ {\alpha} \}}{%
#'        F(x) = exp{-[-(x - m) / s] ^ \alpha}}
#'   for \eqn{x < m}.  The c.d.f. is 1 for \eqn{x \geq m}{x >= m}.
#'
#' @examples
#'
#' set.seed(27)
#'
#' X <- RevWeibull(1, 2)
#' X
#'
#' random(X, 10)
#'
#' pdf(X, 0.7)
#' log_pdf(X, 0.7)
#'
#' cdf(X, 0.7)
#' quantile(X, 0.7)
#'
#' cdf(X, quantile(X, 0.7))
#' quantile(X, cdf(X, 0.7))
RevWeibull <- function(location = 0, scale = 1, shape = 1) {
  if (any(scale <= 0)) {
    stop("scale must be positive")
  }
  if (any(shape <= 0)) {
    stop("shape must be positive")
  }
  stopifnot(
    "parameter lengths do not match (only scalars are allowed to be recycled)" =
      length(location) == length(scale) & length(location) == length(shape) |
        sum(c(length(location) == 1, length(scale) == 1, length(shape) == 1)) >= 2 |
        length(location) == length(scale) & length(shape) == 1 |
        length(location) == length(shape) & length(scale) == 1 |
        length(scale) == length(shape) & length(location) == 1
  )

  d <- data.frame(location = location, scale = scale, shape = shape)
  class(d) <- c("RevWeibull", "distribution")
  d
}

#' @export
mean.RevWeibull <- function(x, ...) {
  rlang::check_dots_used()
  rval <- x$location + x$scale * gamma(1 + 1 / x$shape)
  setNames(rval, names(x))
}

#' @export
variance.RevWeibull <- function(x, ...) {
  rval <- x$scale^2 * gamma(1 + 2 / x$shape) - gamma(1 + 1 / x$shape)^2
  setNames(rval, names(x))
}


#' Draw a random sample from an RevWeibull distribution
#'
#' @inherit RevWeibull examples
#'
#' @param x A `RevWeibull` object created by a call to [RevWeibull()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @return In case of a single distribution object or `n = 1`, either a numeric
#'   vector of length `n` (if `drop = TRUE`, default) or a `matrix` with `n` columns
#'   (if `drop = FALSE`).
#' @export
#'
random.RevWeibull <- function(x, n = 1L, drop = TRUE, ...) {
  n <- make_positive_integer(n)
  if (n == 0L) {
    return(numeric(0L))
  }

  # Convert to the GEV parameterisation
  FUN <- function(at, d) {
    loc <- d$location - d$scale
    scale <- d$scale / d$shape
    shape <- -1 / d$shape
    revdbayes::rgev(n = at, loc = loc, scale = scale, shape = shape)
  }
  apply_dpqr(d = x, FUN = FUN, at = n, type = "random", drop = drop)
}

#' Evaluate the probability mass function of an RevWeibull distribution
#'
#' @inherit RevWeibull examples
#'
#' @param d A `RevWeibull` object created by a call to [RevWeibull()].
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{d} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{d} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[revdbayes]{dgev}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#' @export
#'
pdf.RevWeibull <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  # Convert to the GEV parameterisation
  FUN <- function(at, d) {
    loc <- d$location - d$scale
    scale <- d$scale / d$shape
    shape <- -1 / d$shape
    revdbayes::dgev(x = at, loc = loc, scale = scale, shape = shape, ...)
  }
  apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop, elementwise = elementwise)
}

#' @rdname pdf.RevWeibull
#' @export
#'
log_pdf.RevWeibull <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  # Convert to the GEV parameterisation
  FUN <- function(at, d) {
    loc <- d$location - d$scale
    scale <- d$scale / d$shape
    shape <- -1 / d$shape
    revdbayes::dgev(x = at, loc = loc, scale = scale, shape = shape, log = TRUE)
  }
  apply_dpqr(d = d, FUN = FUN, at = x, type = "logLik", drop = drop, elementwise = elementwise)
}

#' Evaluate the cumulative distribution function of an RevWeibull distribution
#'
#' @inherit RevWeibull examples
#'
#' @param d A `RevWeibull` object created by a call to [RevWeibull()].
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{d} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{d} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[revdbayes]{pgev}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#' @export
#'
cdf.RevWeibull <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  # Convert to the GEV parameterisation
  FUN <- function(at, d) {
    loc <- d$location - d$scale
    scale <- d$scale / d$shape
    shape <- -1 / d$shape
    revdbayes::pgev(q = at, loc = loc, scale = scale, shape = shape, ...)
  }
  apply_dpqr(d = d, FUN = FUN, at = x, type = "probability", drop = drop, elementwise = elementwise)
}

#' Determine quantiles of a RevWeibull distribution
#'
#' `quantile()` is the inverse of `cdf()`.
#'
#' @inherit RevWeibull examples
#' @inheritParams random.RevWeibull
#'
#' @param probs A vector of probabilities.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{x} be evaluated
#'   at all elements of \code{probs} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{x} and \code{probs} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[revdbayes]{qgev}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(probs)` columns (if `drop = FALSE`). In case of a vectorized
#'   distribution object, a matrix with `length(probs)` columns containing all
#'   possible combinations.
#' @export
#'
quantile.RevWeibull <- function(x, probs, drop = TRUE, elementwise = NULL, ...) {
  # Convert to the GEV parameterisation
  FUN <- function(at, d) {
    loc <- d$location - d$scale
    scale <- d$scale / d$shape
    shape <- -1 / d$shape
    revdbayes::qgev(p = at, loc = loc, scale = scale, shape = shape, ...)
  }
  apply_dpqr(d = x, FUN = FUN, at = probs, type = "quantile", drop = drop, elementwise = elementwise)
}

#' Return the support of the RevWeibull distribution
#'
#' @param d An `RevWeibull` object created by a call to [RevWeibull()].
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Currently not used.
#'
#' @return A vector of length 2 with the minimum and maximum value of the support.
#'
#' @export
support.RevWeibull <- function(d, drop = TRUE, ...) {
  rlang::check_dots_used()
  min <- rep(-Inf, length(d))
  max <- d$location
  make_support(min, max, d, drop = drop)
}

#' @exportS3Method
is_discrete.RevWeibull <- function(d, ...) {
  rlang::check_dots_used()
  setNames(rep.int(FALSE, length(d)), names(d))
}

#' @exportS3Method
is_continuous.RevWeibull <- function(d, ...) {
  rlang::check_dots_used()
  setNames(rep.int(TRUE, length(d)), names(d))
}
