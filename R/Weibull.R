#' Create a Weibull distribution
#'
#' Generalization of the gamma distribution. Often used in survival and
#' time-to-event analyses.
#'
#' @param shape The shape parameter \eqn{k}. Can be any positive real number.
#' @param scale The scale parameter \eqn{\lambda}. Can be any positive real
#'   number.
#'
#' @return A `Weibull` object.
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
#'   In the following, let \eqn{X} be a Weibull random variable with
#'   success probability `p` = \eqn{p}.
#'
#'   **Support**: \eqn{R^+} and zero.
#'
#'   **Mean**: \eqn{\lambda \Gamma(1+1/k)}, where \eqn{\Gamma} is
#'   the gamma function.
#'
#'   **Variance**: \eqn{\lambda^2 [ \Gamma (1 + \frac{2}{k} ) - (\Gamma(1+ \frac{1}{k}))^2 ]}
#'
#'   **Probability density function (p.d.f)**:
#'
#'   \deqn{
#'     f(x) = \frac{k}{\lambda}(\frac{x}{\lambda})^{k-1}e^{-(x/\lambda)^k}, x \ge 0
#'   }
#'
#'   **Cumulative distribution function (c.d.f)**:
#'
#'   \deqn{F(x) = 1 - e^{-(x/\lambda)^k}, x \ge 0}
#'
#'   **Moment generating function (m.g.f)**:
#'
#'   \deqn{\sum_{n=0}^\infty \frac{t^n\lambda^n}{n!} \Gamma(1+n/k), k \ge 1}
#'
#' @examples
#'
#' set.seed(27)
#'
#' X <- Weibull(0.3, 2)
#' X
#'
#' random(X, 10)
#'
#' pdf(X, 2)
#' log_pdf(X, 2)
#'
#' cdf(X, 4)
#' quantile(X, 0.7)
Weibull <- function(shape, scale) {
  stopifnot(
    "parameter lengths do not match (only scalars are allowed to be recycled)" =
      length(shape) == length(scale) | length(shape) == 1 | length(scale) == 1
  )
  d <- data.frame(shape = shape, scale = scale)
  class(d) <- c("Weibull", "distribution")
  d
}

#' @export
mean.Weibull <- function(x, ...) {
  rlang::check_dots_used()
  rval <- x$scale * gamma(1 + 1 / x$shape)
  setNames(rval, names(x))
}

#' @export
variance.Weibull <- function(x, ...) {
  rval <- x$scale^2 * (gamma(1 + 2/x$shape) - gamma(1 + 1/x$shape)^2)
  setNames(rval, names(x))
}

#' @export
skewness.Weibull <- function(x, ...) {
  mu <- mean(x)
  sigma <- sqrt(variance(x))
  r <- mu/sigma
  rval <- gamma(1 + 3/x$shape) * (x$scale / sigma)^3 - 3 * r - r^3
  setNames(rval, names(x))
}

#' @export
kurtosis.Weibull <- function(x, ...) {
  mu <- mean(x)
  sigma <- sqrt(variance(x))
  gamma1 <- skewness(x)
  r <- mu/sigma
  rval <- (x$scale/sigma)^4 * gamma(1 + 4/x$shape) - 4 * gamma1 * r - 6 * r^2 - r^4 - 3
  setNames(rval, names(x))
}

#' Draw a random sample from a Weibull distribution
#'
#' @inherit Weibull examples
#'
#' @param x A `Weibull` object created by a call to [Weibull()].
#' @param n The number of samples to draw. Defaults to `1L`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Unused. Unevaluated arguments will generate a warning to
#'   catch mispellings or other possible errors.
#'
#' @family Weibull distribution
#'
#' @return In case of a single distribution object or `n = 1`, either a numeric
#'   vector of length `n` (if `drop = TRUE`, default) or a `matrix` with `n` columns
#'   (if `drop = FALSE`).
#' @export
#'
random.Weibull <- function(x, n = 1L, drop = TRUE, ...) {
  n <- make_positive_integer(n)
  if (n == 0L) {
    return(numeric(0L))
  }
  FUN <- function(at, d) rweibull(n = at, shape = d$shape, scale = d$scale)
  apply_dpqr(d = x, FUN = FUN, at = n, type = "random", drop = drop)
}

#' Evaluate the probability mass function of a Weibull distribution
#'
#' Please see the documentation of [Weibull()] for some properties
#' of the Weibull distribution, as well as extensive examples
#' showing to how calculate p-values and confidence intervals.
#'
#' @inherit Weibull examples
#'
#' @param d A `Weibull` object created by a call to [Weibull()].
#' @param x A vector of elements whose probabilities you would like to
#'   determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{d} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{d} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[stats]{dweibull}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @family Weibull distribution
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#' @export
#'
pdf.Weibull <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) dweibull(x = at, shape = d$shape, scale = d$scale, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop, elementwise = elementwise)
}

#' @rdname pdf.Weibull
#' @export
log_pdf.Weibull <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) dweibull(x = at, shape = d$shape, scale = d$scale, log = TRUE)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "logLik", drop = drop, elementwise = elementwise)
}

#' Evaluate the cumulative distribution function of a Weibull distribution
#'
#' @inherit Weibull examples
#'
#' @param d A `Weibull` object created by a call to [Weibull()].
#' @param x A vector of elements whose cumulative probabilities you would
#'   like to determine given the distribution `d`.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{d} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{d} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[stats]{pweibull}}.
#'   Unevaluated arguments will generate a warning to catch mispellings or other
#'   possible errors.
#'
#' @family Weibull distribution
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of length `probs` (if `drop = TRUE`, default) or a `matrix` with
#'   `length(x)` columns (if `drop = FALSE`). In case of a vectorized distribution
#'   object, a matrix with `length(x)` columns containing all possible combinations.
#' @export
#'
cdf.Weibull <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) pweibull(q = at, shape = d$shape, scale = d$scale, ...)
  apply_dpqr(d = d, FUN = FUN, at = x, type = "probability", drop = drop, elementwise = elementwise)
}

#' Determine quantiles of a Weibull distribution
#'
#' @inherit Weibull examples
#' @inheritParams random.Weibull
#'
#' @param probs A vector of probabilities.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{x} be evaluated
#'   at all elements of \code{probs} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{x} and \code{probs} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Arguments to be passed to \code{\link[stats]{qweibull}}.
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
#' @family Weibull distribution
#'
quantile.Weibull <- function(x, probs, drop = TRUE, elementwise = NULL, ...) {
  FUN <- function(at, d) qweibull(p = at, shape = d$shape, scale = d$scale, ...)
  apply_dpqr(d = x, FUN = FUN, at = probs, type = "quantile", drop = drop, elementwise = elementwise)
}


#' Return the support of the Weibull distribution
#'
#' @param d An `Weibull` object created by a call to [Weibull()].
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param ... Currently not used.
#'
#' @return A vector of length 2 with the minimum and maximum value of the support.
#'
#' @export
support.Weibull <- function(d, drop = TRUE, ...) {
  rlang::check_dots_used()
  min <- rep(0, length(d))
  max <- rep(Inf, length(d))
  make_support(min, max, d, drop = drop)
}

#' @exportS3Method
is_discrete.Weibull <- function(d, ...) {
  rlang::check_dots_used()
  setNames(rep.int(FALSE, length(d)), names(d))
}

#' @exportS3Method
is_continuous.Weibull <- function(d, ...) {
  rlang::check_dots_used()
  setNames(rep.int(TRUE, length(d)), names(d))
}
